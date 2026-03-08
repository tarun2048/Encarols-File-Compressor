
"""
╔══════════════════════════════════════════════════════════════════════╗
║  ENCAROLS  v5.0  —  Adaptive Dual-Engine Compressor + GUI  By Eaevox ║
╠══════════════════════════════════════════════════════════════════════╣
║  ENGINE A  (text files)                                              ║
║    RLE1 → BWT → MTF → RUNA/RUNB → Multi-table Huffman                ║
║    · RLE1: pre-encode runs of ≥4 identical bytes (like bzip2)        ║
║    · BWT:  reorders bytes so contexts cluster across whole file      ║
║    · MTF:  converts clusters to near-zero values                     ║
║    · RUNA/RUNB: zero-run encoding inside Huffman alphabet            ║
║      (no escape byte overhead; RUNA/RUNB get shortest codes)         ║
║    · Multi-table Huffman: 2-6 tables, 50-sym groups                  ║
║    · Compact encoding: lengths stored as bytes (vs 3B/entry old)     ║
║    · Selectors bit-packed (3 bits each vs 1 byte old)                ║
║                                                                      ║
║  ENGINE B  (binary files)                                            ║
║    LZ77 sliding window (32 KB) + Canonical Huffman                   ║
║                                                                      ║
║  GUI                                                                 ║
║    python encarols.py gui  — opens beautiful browser UI              ║
║    Built-in HTTP server, no external dependencies needed             ║
╚══════════════════════════════════════════════════════════════════════╝
"""

import heapq, struct, os, sys, time, json, math
from collections import Counter, defaultdict
from http.server import HTTPServer, BaseHTTPRequestHandler
import threading, webbrowser, base64, traceback, tempfile

__version__ = "5.0.0"
EXT   = ".el"
MAGIC = b'\x45\x4C\x05\x00'

ALGO_BWT  = 0x01
ALGO_LZ77 = 0x02

# LZ77 constants
WIN_BITS  = 15
WIN_SIZE  = 1 << WIN_BITS
MIN_MATCH = 3
MAX_MATCH = 258
MAX_CHAIN = 64
GOOD_LEN  = 64
SYM_MATCH = 256
SYM_EOB   = 257

# BWT constants
BWT_PREFIX  = 512
MTAB_GROUP  = 50    # symbols per Huffman group

# RUNA/RUNB symbols (represent zero runs in MTF output)
RUNA = 0
RUNB = 1


# ═══════════════════════════════════════════════════════
#  BIT I/O
# ═══════════════════════════════════════════════════════

class BitWriter:
    __slots__ = ('_buf','_cur','_n')
    def __init__(self): self._buf=bytearray(); self._cur=0; self._n=0

    def put(self, v:int, n:int):
        self._cur=(self._cur<<n)|(v&((1<<n)-1)); self._n+=n
        while self._n>=8:
            self._n-=8; self._buf.append((self._cur>>self._n)&0xFF)
            self._cur&=(1<<self._n)-1

    def flush(self)->bytes:
        out=bytearray(self._buf)
        if self._n: out.append((self._cur<<(8-self._n))&0xFF); out.append(self._n)
        else: out.append(0)
        return bytes(out)


class BitReader:
    __slots__ = ('_d','_p','_tot')
    def __init__(self, data:bytes):
        pad=data[-1]; self._d=data; self._p=0
        self._tot=(len(data)-1)*8 if pad==0 else (len(data)-2)*8+pad

    def get(self,n:int)->int:
        v=0
        for _ in range(n):
            v=(v<<1)|((self._d[self._p>>3]>>(7-(self._p&7)))&1); self._p+=1
        return v

    @property
    def done(self)->bool: return self._p>=self._tot


# ═══════════════════════════════════════════════════════
#  CANONICAL HUFFMAN
# ═══════════════════════════════════════════════════════

def _build_lengths(freq:dict)->dict:
    if not freq: return {}
    if len(freq)==1: return {next(iter(freq)):1}
    ctr=0; heap=[]
    for sym,f in freq.items(): heapq.heappush(heap,(f,ctr,('L',sym))); ctr+=1
    while len(heap)>1:
        f1,_,n1=heapq.heappop(heap); f2,_,n2=heapq.heappop(heap)
        heapq.heappush(heap,(f1+f2,ctr,('I',n1,n2))); ctr+=1
    L={}
    def walk(node,d):
        if node[0]=='L': L[node[1]]=max(d,1)
        else: walk(node[1],d+1); walk(node[2],d+1)
    walk(heap[0][2],0); return L


def _canon_codes(lengths:dict)->dict:
    if not lengths: return {}
    result={}; code=0; prev=0
    for sym,l in sorted(lengths.items(),key=lambda x:(x[1],x[0])):
        code<<=l-prev; result[sym]=(code,l); code+=1; prev=l
    return result


def _decode_trie(codes:dict)->dict: return {v:k for k,v in codes.items()}


def _decode_sym(br:BitReader, trie:dict, max_l:int)->int:
    code=0
    for n in range(1,max_l+1):
        code=(code<<1)|br.get(1)
        sym=trie.get((code,n))
        if sym is not None: return sym
    raise ValueError("Bad Huffman code")


# ── Compact table encoding ──────────────────────────────
# Store lengths as one byte per symbol (0=absent, 1-20=length)
# Much cheaper than (symbol_index, length) pairs.

def _ser_table_compact(lengths:dict, max_sym:int)->bytes:
    """Serialize table as 2B max_sym + (max_sym+1) bytes of lengths."""
    arr = bytearray(max_sym+1)
    for sym,l in lengths.items():
        if 0<=sym<=max_sym: arr[sym]=l
    return struct.pack('>H',max_sym) + bytes(arr)


def _deser_table_compact(data:bytes, pos:int):
    max_sym=struct.unpack_from('>H',data,pos)[0]; pos+=2
    arr=data[pos:pos+max_sym+1]; pos+=max_sym+1
    lengths={i:arr[i] for i in range(max_sym+1) if arr[i]>0}
    return lengths, pos


# ═══════════════════════════════════════════════════════
#  ENGINE A — RLE1 + BWT + MTF + RUNA/RUNB + MULTI-TABLE HUFFMAN
# ═══════════════════════════════════════════════════════

def _rle1_encode(data:bytes)->bytes:
    """Pre-BWT RLE: encode runs of ≥4 identical bytes as [b,b,b,b,extra]."""
    out=bytearray(); i=0; n=len(data)
    while i<n:
        b=data[i]; j=i
        while j<n and data[j]==b and j-i<258: j+=1
        run=j-i
        if run>=4: out.extend([b,b,b,b,run-4])
        else: out.extend(data[i:j])
        i=j
    return bytes(out)


def _rle1_decode(data:bytes)->bytes:
    out=bytearray(); i=0; n=len(data); run_len=0; prev=-1
    while i<n:
        b=data[i]; i+=1; out.append(b)
        if b==prev: run_len+=1
        else: run_len=1; prev=b
        if run_len==4:
            extra=data[i]; i+=1; out.extend([b]*extra)
            run_len=0; prev=-1
    return bytes(out)


def _bwt_encode(data:bytes)->tuple:
    n=len(data); doubled=data+data; P=min(n,BWT_PREFIX)
    indices=sorted(range(n),key=lambda i:doubled[i:i+P])
    return bytes(data[(i-1)%n] for i in indices), indices.index(0)


def _bwt_decode(bwt:bytes,primary:int)->bytes:
    n=len(bwt)
    count=[0]*256
    for b in bwt: count[b]+=1
    first=[0]*256; c=0
    for i in range(256): first[i]=c; c+=count[i]
    occ=[0]*256; lf=[0]*n
    for i in range(n):
        b=bwt[i]; lf[i]=first[b]+occ[b]; occ[b]+=1
    out=bytearray(n); idx=primary
    for i in range(n-1,-1,-1): out[i]=bwt[idx]; idx=lf[idx]
    return bytes(out)


def _mtf_encode(data:bytes)->bytes:
    table=list(range(256)); out=bytearray()
    for b in data: r=table.index(b); out.append(r); table.pop(r); table.insert(0,b)
    return bytes(out)


def _mtf_decode(data:bytes)->bytes:
    table=list(range(256)); out=bytearray()
    for r in data: b=table[r]; out.append(b); table.pop(r); table.insert(0,b)
    return bytes(out)


def _runa_encode(mtf:bytes)->list:
    """
    RUNA/RUNB zero-run encoding (bzip2-style).
    Zero runs become sequences of RUNA(0)/RUNB(1) symbols representing
    the run length in binary (LSB-first).  Non-zero MTF symbols are
    shifted up by 1 to make room for RUNA/RUNB.
    This way RUNA/RUNB get the SHORTEST Huffman codes since they
    represent the most frequent symbol (0) in the MTF output.
    """
    out=[]; i=0; n=len(mtf)
    while i<n:
        if mtf[i]==0:
            j=i
            while j<n and mtf[j]==0: j+=1
            v=j-i
            while v>0:
                out.append(RUNB if v&1 else RUNA); v>>=1
            i=j
        else: out.append(mtf[i]+1); i+=1
    return out


def _runa_decode(syms:list, expected_len:int)->bytes:
    """Inverse RUNA/RUNB decoding."""
    out=bytearray(); i=0; n=len(syms)
    while i<n:
        if syms[i] in (RUNA,RUNB):
            # Collect RUNA/RUNB run → decode zero-run length
            bits=[]; j=i
            while j<n and syms[j] in (RUNA,RUNB): bits.append(syms[j]); j+=1
            # Reconstruct run length from LSB-first binary
            run=0
            for k,b in enumerate(bits): run|=(b<<k)
            out.extend(bytes(run)); i=j
        else: out.append(syms[i]-1); i+=1
    return bytes(out)


def _multi_huff_compress(symbols:list)->bytes:
    """
    Multi-table Huffman compression of the RUNA/RUNB symbol stream.

    Algorithm (mirrors bzip2):
    1. Choose N tables (2-6 based on stream size)
    2. Initialize table assignments by dividing stream into N equal chunks
    3. Iterate 4×: for each group assign to lowest-cost table; rebuild tables
    4. Encode:
       - 1B: n_tables
       - For each table: compact length array (2B max_sym + bytes)
       - Selectors bit-packed at ceil(log2(n_tables)) bits each
       - Bit-stream coded with correct table per group
    """
    n=len(symbols)
    if n==0: return struct.pack('>I',0)

    # Unique symbols in stream
    unique=sorted(set(symbols))
    max_sym=max(unique) if unique else 0

    groups=[symbols[i:i+MTAB_GROUP] for i in range(0,n,MTAB_GROUP)]
    G=len(groups)

    # Choose number of tables
    n_tables=6 if G>=48 else (5 if G>=24 else (4 if G>=12 else (3 if G>=6 else 2)))
    n_tables=max(1,min(n_tables,G))

    # Initial assignments: spread evenly
    assignments=[gi*n_tables//G for gi in range(G)]

    def rebuild():
        freqs=[defaultdict(int) for _ in range(n_tables)]
        for gi,grp in enumerate(groups):
            for s in grp: freqs[assignments[gi]][s]+=1
        tables=[]
        for t in range(n_tables):
            freq=dict(freqs[t])
            # Ensure every unique symbol has at least freq=1 (table can encode anything)
            for s in unique:
                if s not in freq: freq[s]=1
            l=_build_lengths(freq); c=_canon_codes(l)
            tables.append((l,c))
        return tables

    # Iterate to optimize assignments
    for _ in range(4):
        tables=rebuild()
        changed=False
        for gi,grp in enumerate(groups):
            best_t=0; best_cost=float('inf')
            for t,(_,codes) in enumerate(tables):
                cost=sum(codes.get(s,(0,32))[1] for s in grp)
                if cost<best_cost: best_cost=cost; best_t=t
            if best_t!=assignments[gi]: assignments[gi]=best_t; changed=True
        if not changed: break

    tables=rebuild()

    # Encode bitstream
    bw=BitWriter()
    for gi,grp in enumerate(groups):
        _,codes=tables[assignments[gi]]
        for s in grp: c,nb=codes[s]; bw.put(c,nb)
    bits=bw.flush()

    # Pack selectors: ceil(log2(n_tables)) bits each
    sel_bits=max(1,math.ceil(math.log2(max(n_tables,2))))
    sbw=BitWriter()
    for a in assignments: sbw.put(a,sel_bits)
    sel_data=sbw.flush()

    # Build header: n_tables, tables (compact), G, sel_bits, sel_data, n_symbols, bitstream
    out=bytearray([n_tables])
    for lengths,_ in tables: out+=_ser_table_compact(lengths,max_sym)
    out+=struct.pack('>IHB',n,G,sel_bits)
    out+=struct.pack('>I',len(sel_data))+sel_data
    out+=bits
    return bytes(out)


def _multi_huff_decompress(data:bytes)->list:
    pos=0
    n_tables=data[pos]; pos+=1

    tables=[]
    for _ in range(n_tables):
        lengths,pos=_deser_table_compact(data,pos)
        codes=_canon_codes(lengths); trie=_decode_trie(codes)
        max_l=max(lengths.values()) if lengths else 1
        tables.append((trie,max_l))

    n_syms,G,sel_bits=struct.unpack_from('>IHB',data,pos); pos+=7
    sel_len=struct.unpack_from('>I',data,pos)[0]; pos+=4
    sel_data=data[pos:pos+sel_len]; pos+=sel_len

    # Decode selectors
    sbr=BitReader(sel_data)
    assignments=[sbr.get(sel_bits) for _ in range(G)]

    br=BitReader(data[pos:])
    out=[]; decoded=0; gi=0

    while decoded<n_syms:
        trie,max_l=tables[assignments[gi]]
        group_sz=min(MTAB_GROUP,n_syms-decoded)
        for _ in range(group_sz):
            sym=_decode_sym(br,trie,max_l)
            out.append(sym); decoded+=1
        gi+=1

    return out


def _compress_bwt(data:bytes)->bytes:
    """Engine A pipeline: RLE1 → BWT → MTF → RUNA/RUNB → Multi-table Huffman"""
    rle1=_rle1_encode(data)
    bwt,pidx=_bwt_encode(rle1)
    mtf=_mtf_encode(bwt)
    runa_syms=_runa_encode(mtf)
    payload=_multi_huff_compress(runa_syms)
    return struct.pack('>III',len(data),len(rle1),pidx)+payload


def _decompress_bwt(data:bytes)->bytes:
    orig_sz,rle1_sz,pidx=struct.unpack_from('>III',data,0)
    runa_syms=_multi_huff_decompress(data[12:])
    mtf=_runa_decode(runa_syms,rle1_sz)
    bwt=_mtf_decode(mtf)
    rle1=_bwt_decode(bwt,pidx)
    orig=_rle1_decode(rle1)
    if len(orig)!=orig_sz:
        raise RuntimeError(f"BWT decode size mismatch: {len(orig)} vs {orig_sz}")
    return orig


# ═══════════════════════════════════════════════════════
#  ENGINE B — LZ77 + HUFFMAN
# ═══════════════════════════════════════════════════════

def _lz77_tokens(data:bytes)->list:
    n=len(data); pos=0; out=[]
    head:dict={}; prev:dict={}

    def ins(p):
        if p+2<n:
            k=(data[p],data[p+1],data[p+2])
            prev[p]=head.get(k,-1); head[k]=p

    def mlen(cur,cap):
        lo,hi=0,cap
        while lo<hi:
            mid=(lo+hi+1)>>1
            if data[cur:cur+mid]==data[pos:pos+mid]: lo=mid
            else: hi=mid-1
        return lo

    while pos<n:
        bl=bd=0
        if pos+MIN_MATCH<=n:
            k=(data[pos],data[pos+1],data[pos+2])
            cur=head.get(k,-1); itr=0
            while cur>=0 and itr<MAX_CHAIN:
                d=pos-cur
                if d>WIN_SIZE: break
                ml=mlen(cur,min(MAX_MATCH,n-pos))
                if ml>bl: bl=ml; bd=d
                if bl>=GOOD_LEN: break
                cur=prev.get(cur,-1); itr+=1
        if bl>=MIN_MATCH:
            out.append(('M',bl,bd))
            for i in range(bl): ins(pos+i)
            pos+=bl
        else: out.append(('L',data[pos])); ins(pos); pos+=1
    return out


def _compress_lz77(data:bytes)->bytes:
    tokens=_lz77_tokens(data)
    lf:dict=defaultdict(int); nf:dict=defaultdict(int)
    lf[SYM_EOB]=1
    for t in tokens:
        if t[0]=='L': lf[t[1]]+=1
        else: lf[SYM_MATCH]+=1; nf[t[1]-MIN_MATCH]+=1
    ll=_build_lengths(dict(lf)); nl=_build_lengths(dict(nf)) if nf else {}
    lc=_canon_codes(ll); nc=_canon_codes(nl)
    bw=BitWriter()
    for t in tokens:
        if t[0]=='L': c,n=lc[t[1]]; bw.put(c,n)
        else:
            _,ln,d=t
            c,n=lc[SYM_MATCH]; bw.put(c,n)
            c,n=nc[ln-MIN_MATCH]; bw.put(c,n)
            bw.put(d-1,WIN_BITS)
    c,n=lc[SYM_EOB]; bw.put(c,n)
    bits=bw.flush()
    max_lsym=max(ll.keys()) if ll else 0
    max_nsym=max(nl.keys()) if nl else 0
    return (struct.pack('>I',len(data))
            +_ser_table_compact(ll,max_lsym)
            +_ser_table_compact(nl,max_nsym)
            +bits)


def _decompress_lz77(data:bytes)->bytes:
    orig=struct.unpack_from('>I',data,0)[0]; pos=4
    ll,pos=_deser_table_compact(data,pos)
    nl,pos=_deser_table_compact(data,pos)
    lc=_canon_codes(ll); nc=_canon_codes(nl)
    lt=_decode_trie(lc); nt=_decode_trie(nc)
    ml=max(ll.values()) if ll else 1; mn=max(nl.values()) if nl else 1
    br=BitReader(data[pos:]); buf=bytearray()
    while not br.done:
        sym=_decode_sym(br,lt,ml)
        if sym==SYM_EOB: break
        elif sym<SYM_MATCH: buf.append(sym)
        else:
            ls=_decode_sym(br,nt,mn); ln=ls+MIN_MATCH
            d=br.get(WIN_BITS)+1; start=len(buf)-d
            for i in range(ln): buf.append(buf[start+i])
    if len(buf)!=orig: raise RuntimeError(f"LZ77 size mismatch: {len(buf)} vs {orig}")
    return bytes(buf)


# ═══════════════════════════════════════════════════════
#  AUTO-DETECTION & PUBLIC API
# ═══════════════════════════════════════════════════════

def _detect(data:bytes)->int:
    sample=data[:4096]
    if not sample: return ALGO_LZ77
    p=sum(1 for b in sample if 32<=b<=126 or b in (9,10,13))
    return ALGO_BWT if p/len(sample)>0.85 else ALGO_LZ77


def compress(data:bytes, algo:int=0)->bytes:
    if algo==0: algo=_detect(data)
    payload=(_compress_bwt(data) if algo==ALGO_BWT else _compress_lz77(data))
    return bytes([algo])+payload


def decompress(data:bytes)->bytes:
    algo=data[0]
    if algo==ALGO_BWT:  return _decompress_bwt(data[1:])
    if algo==ALGO_LZ77: return _decompress_lz77(data[1:])
    raise ValueError(f"Unknown algorithm 0x{algo:02X}")


# ═══════════════════════════════════════════════════════
#  FILE I/O
# ═══════════════════════════════════════════════════════

def encode_file(src:str, dst:str, force:int=0)->dict:
    with open(src,'rb') as f: raw=f.read()
    algo=force if force else _detect(raw)
    t0=time.perf_counter()
    payload=compress(raw,algo)
    dt=time.perf_counter()-t0
    with open(dst,'wb') as f: f.write(MAGIC+payload)
    enc=os.path.getsize(dst); orig=len(raw)
    return dict(orig=orig, enc=enc, ratio=enc/orig*100 if orig else 0,
                saved=orig-enc, time=dt, algo=algo,
                algo_name='BWT+MTF+RUNA+MultiHuffman' if algo==ALGO_BWT else 'LZ77+Huffman')


def decode_file(src:str, dst:str)->dict:
    with open(src,'rb') as f: raw=f.read()
    if raw[:4]!=MAGIC: raise ValueError("Not a valid ENCAROLS v5 file.")
    t0=time.perf_counter()
    data=decompress(raw[4:])
    dt=time.perf_counter()-t0
    with open(dst,'wb') as f: f.write(data)
    return dict(size=len(data), time=dt)


# ═══════════════════════════════════════════════════════
#  TESTS
# ═══════════════════════════════════════════════════════

_TESTS=[
    b"",b"x",b"hello world",
    b"The quick brown fox jumps over the lazy dog.",
    b"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    b"paragraph sentence thesis argument "*400,
    bytes(range(256))*20,
    b"\r\n"*5000,
    b"The cat sat on the mat. "*100,
    b"abcdef"*300,
]

def run_tests()->None:
    _banner()
    print("\n  Round-trip tests\n")
    passed=failed=0
    for i,data in enumerate(_TESTS,1):
        label=repr(data[:45])[2:-1]+('…' if len(data)>45 else '')
        algo=_detect(data)
        try:
            enc=compress(data,algo); dec=decompress(enc)
            ok=dec==data; col='\033[92m' if ok else '\033[91m'; rst='\033[0m'
            rat=len(enc)/max(len(data),1)*100
            eng='A' if algo==ALGO_BWT else 'B'
            print(f"  [{col}{'PASS' if ok else 'FAIL'}{rst}]  "
                  f"Eng={eng}  in={len(data):>8,}B  enc={len(enc):>8,}B  "
                  f"{rat:>5.0f}%  \"{label}\"")
            if not ok: print(f"       ✘ got {len(dec)}B, expected {len(data)}B"); failed+=1
            else: passed+=1
        except Exception as ex:
            print(f"  [\033[91mERROR\033[0m] Test {i}: {ex}"); traceback.print_exc(); failed+=1
    print(f"\n  {passed} passed, {failed} failed\n")


# ═══════════════════════════════════════════════════════
#  BENCHMARK
# ═══════════════════════════════════════════════════════

def bench_file(path:str)->None:
    import subprocess, tempfile
    with open(path,'rb') as f: raw=f.read()
    freq=Counter(raw); total=len(raw)
    entropy=(-sum((c/total)*math.log2(c/total) for c in freq.values()) if total else 0)
    algo=_detect(raw)

    _banner()
    print(f"\n  File     : {path}")
    print(f"  Size     : {total:,} bytes")
    print(f"  Entropy  : {entropy:.2f} bits/byte  "
          f"(floor: {int(total*entropy/8):,}B = {100*entropy/8:.1f}%)")
    print(f"  Engine   : {'A — BWT+MTF+RUNA+MultiHuffman' if algo==ALGO_BWT else 'B — LZ77+Huffman'}")
    print()

    print(f"  Compressing with Encarols v5…", end=' ', flush=True)
    t0=time.perf_counter(); enc=compress(raw,algo); t1=time.perf_counter()
    dec=decompress(enc); t2=time.perf_counter()
    ok=dec==raw; esz=len(enc)+4
    rat=esz/total*100; sav=total-esz

    print("done")
    print(f"  Encarols v5 : {esz:>10,}B  ({rat:.1f}%)  "
          f"enc {t1-t0:.1f}s  dec {t2-t1:.1f}s  lossless={'✔' if ok else '✘ BUG'}")
    print()

    with tempfile.NamedTemporaryFile(suffix='.dat',delete=False) as f:
        f.write(raw); fname=f.name
    print("  Competitors:")
    for lbl,cmd,ext in [('gzip  -9',['gzip','-9','-k','-f',fname],'.gz'),
                        ('bzip2 -9',['bzip2','-9','-k','-f',fname],'.bz2'),
                        ('xz    -9',['xz','-9','-k','-f',fname],'.xz')]:
        try:
            t0=time.perf_counter()
            import subprocess
            subprocess.run(cmd,capture_output=True,timeout=120)
            dt=time.perf_counter()-t0
            out=fname+ext
            if os.path.exists(out):
                sz=os.path.getsize(out); r=sz/total*100
                vs=("✔ Encarols wins!" if esz<sz
                    else f"  bzip2 wins by {sz-esz:,}B" if 'bzip2' in lbl
                    else f"  ({sz-esz:+,}B vs Encarols)")
                print(f"  {lbl}     : {sz:>10,}B  ({r:.1f}%)  {dt:.1f}s  {vs}")
                os.unlink(out)
        except: print(f"  {lbl}     : not available")
    os.unlink(fname)
    print()


# ═══════════════════════════════════════════════════════
#  HTML GUI  (built-in server, opens in browser)
# ═══════════════════════════════════════════════════════

HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ENCAROLS v5 — Compressor</title>
<style>
  :root{--bg:#0d1117;--card:#161b22;--border:#30363d;--accent:#58a6ff;
        --green:#3fb950;--red:#f85149;--yellow:#e3b341;--text:#e6edf3;
        --muted:#8b949e;--hover:#21262d}
  *{box-sizing:border-box;margin:0;padding:0}
  body{background:var(--bg);color:var(--text);font-family:'Segoe UI',system-ui,sans-serif;
       min-height:100vh;display:flex;flex-direction:column;align-items:center;padding:24px 16px}
  h1{font-size:1.6rem;font-weight:700;letter-spacing:-0.5px;margin-bottom:4px}
  .subtitle{color:var(--muted);font-size:.85rem;margin-bottom:28px}
  .card{background:var(--card);border:1px solid var(--border);border-radius:12px;
        padding:24px;width:100%;max-width:680px;margin-bottom:20px}
  .card-title{font-size:.75rem;font-weight:600;text-transform:uppercase;
              letter-spacing:.08em;color:var(--muted);margin-bottom:16px}

  /* Drop zone */
  .drop-zone{border:2px dashed var(--border);border-radius:8px;padding:40px 20px;
             text-align:center;cursor:pointer;transition:all .2s;position:relative}
  .drop-zone:hover,.drop-zone.drag{border-color:var(--accent);background:rgba(88,166,255,.05)}
  .drop-zone input{position:absolute;inset:0;opacity:0;cursor:pointer;width:100%}
  .drop-icon{font-size:2rem;margin-bottom:8px}
  .drop-label{color:var(--muted);font-size:.9rem}
  .drop-label b{color:var(--accent)}
  .file-info{margin-top:12px;padding:10px 14px;background:var(--hover);
             border-radius:6px;font-size:.85rem;display:none}
  .file-info.show{display:flex;align-items:center;gap:10px}
  .file-icon{font-size:1.2rem}
  .file-name{flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .file-size{color:var(--muted);white-space:nowrap}

  /* Mode / engine toggles */
  .row{display:flex;gap:12px;flex-wrap:wrap}
  .toggle-group{display:flex;border:1px solid var(--border);border-radius:8px;overflow:hidden;flex:1}
  .toggle-btn{flex:1;padding:9px 8px;background:transparent;border:none;color:var(--muted);
              cursor:pointer;font-size:.85rem;font-weight:500;transition:all .15s;text-align:center}
  .toggle-btn.active{background:var(--accent);color:#0d1117;font-weight:700}
  .toggle-btn:hover:not(.active){background:var(--hover);color:var(--text)}

  /* Big action button */
  #go-btn{width:100%;padding:14px;border-radius:8px;border:none;font-size:1rem;
          font-weight:700;cursor:pointer;letter-spacing:.02em;
          background:var(--accent);color:#0d1117;transition:all .2s}
  #go-btn:hover:not(:disabled){filter:brightness(1.1);transform:translateY(-1px)}
  #go-btn:disabled{opacity:.4;cursor:not-allowed;transform:none}

  /* Progress */
  .progress-wrap{height:4px;background:var(--border);border-radius:2px;overflow:hidden;margin-top:14px}
  .progress-bar{height:100%;background:var(--accent);width:0;transition:width .3s;border-radius:2px}
  .status-msg{font-size:.8rem;color:var(--muted);margin-top:8px;min-height:1.2em}

  /* Results */
  .result-grid{display:grid;grid-template-columns:1fr 1fr;gap:12px}
  .stat{background:var(--hover);border-radius:8px;padding:14px}
  .stat-label{font-size:.7rem;text-transform:uppercase;letter-spacing:.07em;
              color:var(--muted);margin-bottom:4px}
  .stat-value{font-size:1.3rem;font-weight:700}
  .stat-value.green{color:var(--green)}
  .stat-value.blue{color:var(--accent)}
  .stat-value.yellow{color:var(--yellow)}
  .ratio-bar{margin-top:14px;background:var(--border);border-radius:4px;height:8px;overflow:hidden}
  .ratio-fill{height:100%;border-radius:4px;transition:width 1s ease;background:var(--green)}
  .result-card{display:none}
  .result-card.show{display:block}

  /* History */
  .history-list{list-style:none;max-height:180px;overflow-y:auto}
  .history-list::-webkit-scrollbar{width:4px}
  .history-list::-webkit-scrollbar-thumb{background:var(--border);border-radius:2px}
  .h-item{display:flex;align-items:center;gap:10px;padding:8px 10px;
          border-radius:6px;font-size:.82rem;margin-bottom:4px;background:var(--hover)}
  .h-badge{font-size:.65rem;padding:2px 7px;border-radius:10px;font-weight:700;white-space:nowrap}
  .h-badge.comp{background:rgba(63,185,80,.2);color:var(--green)}
  .h-badge.decomp{background:rgba(88,166,255,.2);color:var(--accent)}
  .h-names{flex:1;overflow:hidden}
  .h-from{overflow:hidden;text-overflow:ellipsis;white-space:nowrap;color:var(--muted)}
  .h-ratio{color:var(--yellow);font-weight:600;white-space:nowrap}
  .empty{color:var(--muted);font-size:.85rem;text-align:center;padding:16px 0}

  /* Download button */
  #dl-btn{display:none;margin-top:14px;width:100%;padding:11px;
          border-radius:8px;border:1px solid var(--green);color:var(--green);
          background:transparent;cursor:pointer;font-weight:600;font-size:.9rem;transition:all .2s}
  #dl-btn:hover{background:rgba(63,185,80,.1)}
  #dl-btn.show{display:block}

  .engines-info{font-size:.78rem;color:var(--muted);margin-top:12px;line-height:1.6}
  .engines-info span{color:var(--text)}
</style>
</head>
<body>

<h1>🗜 ENCAROLS</h1>
<p class="subtitle">v5.0 &mdash; Adaptive Dual-Engine File Compressor</p>

<!-- File selector -->
<div class="card">
  <div class="card-title">Input File</div>
  <div class="drop-zone" id="dz">
    <input type="file" id="file-input" onchange="onFile(this)">
    <div class="drop-icon">📂</div>
    <div class="drop-label">Drop a file here or <b>click to browse</b></div>
  </div>
  <div class="file-info" id="fi">
    <span class="file-icon">📄</span>
    <span class="file-name" id="fn"></span>
    <span class="file-size" id="fs"></span>
  </div>
</div>

<!-- Options -->
<div class="card">
  <div class="card-title">Options</div>
  <div class="row">
    <div>
      <div style="font-size:.75rem;color:var(--muted);margin-bottom:6px">Mode</div>
      <div class="toggle-group">
        <button class="toggle-btn active" id="m-comp" onclick="setMode('compress')">⬇ Compress</button>
        <button class="toggle-btn" id="m-decomp" onclick="setMode('decompress')">⬆ Decompress</button>
      </div>
    </div>
    <div>
      <div style="font-size:.75rem;color:var(--muted);margin-bottom:6px">Engine</div>
      <div class="toggle-group">
        <button class="toggle-btn active" id="e-auto" onclick="setEngine('auto')">🔍 Auto</button>
        <button class="toggle-btn" id="e-text" onclick="setEngine('text')">📝 Text (BWT)</button>
        <button class="toggle-btn" id="e-bin"  onclick="setEngine('binary')">💾 Binary (LZ77)</button>
      </div>
    </div>
  </div>
  <div class="engines-info">
    <span>Engine A — BWT+MTF+RUNA+MultiHuffman:</span> best for text, HTML, source code, logs.<br>
    <span>Engine B — LZ77+Huffman:</span> best for binary, executables, already-mixed data.
  </div>
</div>

<!-- Action -->
<div class="card">
  <button id="go-btn" disabled onclick="go()">Select a file first</button>
  <div class="progress-wrap"><div class="progress-bar" id="pb"></div></div>
  <div class="status-msg" id="sm"></div>
  <button id="dl-btn" onclick="download()">⬇ Download Result</button>
</div>

<!-- Results -->
<div class="card result-card" id="res">
  <div class="card-title">Results</div>
  <div class="result-grid">
    <div class="stat">
      <div class="stat-label">Original</div>
      <div class="stat-value blue" id="r-orig">—</div>
    </div>
    <div class="stat">
      <div class="stat-label">Output</div>
      <div class="stat-value blue" id="r-out">—</div>
    </div>
    <div class="stat">
      <div class="stat-label">Ratio</div>
      <div class="stat-value yellow" id="r-rat">—</div>
    </div>
    <div class="stat">
      <div class="stat-label">Saved / Time</div>
      <div class="stat-value green" id="r-sav">—</div>
    </div>
  </div>
  <div class="ratio-bar"><div class="ratio-fill" id="r-bar" style="width:0"></div></div>
  <div style="font-size:.78rem;color:var(--muted);margin-top:8px" id="r-eng">Engine: —</div>
</div>

<!-- History -->
<div class="card">
  <div class="card-title">History</div>
  <ul class="history-list" id="hist"><li class="empty">No operations yet</li></ul>
</div>

<script>
let fileData=null, fileName='', mode='compress', engine='auto';
let resultData=null, resultName='';

const dz=document.getElementById('dz');
dz.addEventListener('dragover',e=>{e.preventDefault();dz.classList.add('drag')});
dz.addEventListener('dragleave',()=>dz.classList.remove('drag'));
dz.addEventListener('drop',e=>{
  e.preventDefault(); dz.classList.remove('drag');
  if(e.dataTransfer.files[0]) loadFile(e.dataTransfer.files[0]);
});

function onFile(inp){ if(inp.files[0]) loadFile(inp.files[0]); }

function loadFile(f){
  fileName=f.name;
  const r=new FileReader();
  r.onload=e=>{
    fileData=new Uint8Array(e.target.result);
    document.getElementById('fn').textContent=f.name;
    document.getElementById('fs').textContent=fmtSz(fileData.length);
    document.getElementById('fi').classList.add('show');
    const btn=document.getElementById('go-btn');
    btn.disabled=false;
    btn.textContent=mode==='compress'?'⬇ Compress File':'⬆ Decompress File';
  };
  r.readAsArrayBuffer(f);
}

function setMode(m){
  mode=m;
  document.getElementById('m-comp').classList.toggle('active',m==='compress');
  document.getElementById('m-decomp').classList.toggle('active',m==='decompress');
  if(fileData){
    const btn=document.getElementById('go-btn');
    btn.textContent=m==='compress'?'⬇ Compress File':'⬆ Decompress File';
  }
}
function setEngine(e){
  engine=e;
  ['auto','text','binary'].forEach(x=>
    document.getElementById('e-'+x).classList.toggle('active',x===e));
}

function fmtSz(n){
  if(n<1024) return n+'B';
  if(n<1048576) return (n/1024).toFixed(1)+'KB';
  return (n/1048576).toFixed(2)+'MB';
}

function setPb(w,msg){
  document.getElementById('pb').style.width=w+'%';
  document.getElementById('sm').textContent=msg;
}

async function go(){
  if(!fileData) return;
  const btn=document.getElementById('go-btn');
  btn.disabled=true; resultData=null;
  document.getElementById('dl-btn').classList.remove('show');
  document.getElementById('res').classList.remove('show');
  setPb(10,'Uploading…');

  /* Safe chunked base64 — avoids call-stack overflow on large files */
  function toB64(u8){
    let s=''; const chunk=8192;
    for(let i=0;i<u8.length;i+=chunk)
      s+=String.fromCharCode(...u8.subarray(i,i+chunk));
    return btoa(s);
  }
  const payload={file_b64:toB64(fileData), filename:fileName, mode, engine};

  setPb(40,'Processing…');
  try{
    const resp=await fetch('/api/process',{method:'POST',
      headers:{'Content-Type':'application/json'},
      body:JSON.stringify(payload)});
    const data=await resp.json();
    setPb(100,'');

    if(data.error){ setPb(0,'❌ '+data.error); btn.disabled=false; return; }

    /* Safe decode of potentially large base64 result */
    const raw=atob(data.result_b64);
    resultData=new Uint8Array(raw.length);
    for(let i=0;i<raw.length;i++) resultData[i]=raw.charCodeAt(i);
    resultName=data.result_name;

    // Show results
    const res=document.getElementById('res');
    res.classList.add('show');
    document.getElementById('r-orig').textContent=fmtSz(data.orig_size);
    document.getElementById('r-out').textContent=fmtSz(data.out_size);

    if(mode==='compress'){
      const pct=data.out_size/data.orig_size*100;
      const saved=data.orig_size-data.out_size;
      document.getElementById('r-rat').textContent=pct.toFixed(1)+'%';
      document.getElementById('r-sav').textContent=
        (saved>0?'-'+fmtSz(saved)+' / ':'+'+ fmtSz(-saved)+' / ')+data.time.toFixed(2)+'s';
      document.getElementById('r-bar').style.width=Math.min(pct,100)+'%';
      document.getElementById('r-bar').style.background=
        pct<30?'var(--green)':pct<60?'var(--yellow)':'var(--red)';
    } else {
      document.getElementById('r-rat').textContent='Restored';
      document.getElementById('r-sav').textContent=data.time.toFixed(2)+'s';
      document.getElementById('r-bar').style.width='100%';
      document.getElementById('r-bar').style.background='var(--accent)';
    }
    document.getElementById('r-eng').textContent='Engine: '+data.engine_name;
    document.getElementById('dl-btn').textContent='⬇ Download  '+resultName;
    document.getElementById('dl-btn').classList.add('show');

    // Add to history
    const hist=document.getElementById('hist');
    const empty=hist.querySelector('.empty');
    if(empty) empty.remove();
    const li=document.createElement('li'); li.className='h-item';
    const pct2=(data.out_size/data.orig_size*100).toFixed(1);
    li.innerHTML=`
      <span class="h-badge ${mode==='compress'?'comp':'decomp'}">${mode==='compress'?'COMP':'DECOMP'}</span>
      <span class="h-names">
        <div class="h-from">${fileName}</div>
        <div style="font-size:.75rem;color:var(--muted)">${resultName}</div>
      </span>
      <span class="h-ratio">${mode==='compress'?pct2+'%':'✔'}</span>`;
    hist.prepend(li);
    setPb(100,'✔ Done');
  }catch(e){ setPb(0,'❌ '+e.message); }
  btn.disabled=false;
  btn.textContent=mode==='compress'?'⬇ Compress File':'⬆ Decompress File';
}

function download(){
  if(!resultData) return;
  const blob=new Blob([resultData]);
  const a=document.createElement('a'); a.href=URL.createObjectURL(blob);
  a.download=resultName; a.click();
}
</script>
</body></html>"""


class _Handler(BaseHTTPRequestHandler):
    def log_message(self, *_): pass   # suppress access logs

    def do_GET(self):
        self.send_response(200)
        self.send_header('Content-Type','text/html; charset=utf-8')
        self.end_headers()
        self.wfile.write(HTML.encode())

    def do_POST(self):
        length=int(self.headers.get('Content-Length',0))
        # Read in chunks to avoid blocking on large files
        chunks=[]; remaining=length
        while remaining>0:
            chunk=self.rfile.read(min(remaining,65536))
            if not chunk: break
            chunks.append(chunk); remaining-=len(chunk)
        body=b''.join(chunks)
        try:
            req=json.loads(body)
            file_bytes=base64.b64decode(req['file_b64'])
            fname=req['filename']
            mode=req.get('mode','compress')
            eng=req.get('engine','auto')

            if mode=='compress':
                algo={'auto':0,'text':ALGO_BWT,'binary':ALGO_LZ77}.get(eng,0)
                if algo==0: algo=_detect(file_bytes)
                t0=time.perf_counter()
                result=MAGIC+compress(file_bytes,algo)
                dt=time.perf_counter()-t0
                rname=os.path.splitext(fname)[0]+'.el'
                eng_name=('BWT+MTF+RUNA+MultiHuffman' if algo==ALGO_BWT else 'LZ77+Huffman')
                resp=dict(result_b64=base64.b64encode(result).decode(),
                          result_name=rname, orig_size=len(file_bytes),
                          out_size=len(result), time=dt, engine_name=eng_name)
            else:
                if file_bytes[:4]!=MAGIC:
                    raise ValueError("Not a valid Encarols v5 file.")
                t0=time.perf_counter()
                result=decompress(file_bytes[4:])
                dt=time.perf_counter()-t0
                rname=fname.replace('.el','') if fname.endswith('.el') else fname+'_decoded'
                algo=file_bytes[4]
                eng_name=('BWT+MTF+RUNA+MultiHuffman' if algo==ALGO_BWT else 'LZ77+Huffman')
                resp=dict(result_b64=base64.b64encode(result).decode(),
                          result_name=rname, orig_size=len(file_bytes),
                          out_size=len(result), time=dt, engine_name=eng_name)
        except Exception as e:
            resp=dict(error=str(e))

        data=json.dumps(resp).encode()
        self.send_response(200)
        self.send_header('Content-Type','application/json')
        self.send_header('Content-Length',str(len(data)))
        self.end_headers()
        self.wfile.write(data)


def launch_gui(port:int=8765)->None:
    _banner()
    print(f"  Starting GUI server on http://localhost:{port}")
    print(f"  Opening browser… (press Ctrl+C to quit)\n")
    server=HTTPServer(('localhost',port),_Handler)
    t=threading.Thread(target=lambda:webbrowser.open(f'http://localhost:{port}'),daemon=True)
    t.start()
    try: server.serve_forever()
    except KeyboardInterrupt: print("\n  Server stopped.")


# ═══════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════

def _banner()->None:
    print("""
╔════════════════════════════════════════════════════════════════╗
║  ENCAROLS  v5.0  —  Adaptive Dual-Engine Compressor  By Eaevox ║
║  Engine A: BWT+MTF+RUNA+MultiHuffman  (text files)             ║
║  Engine B: LZ77+Huffman               (binary files)           ║
╚════════════════════════════════════════════════════════════════╝""")


def _usage()->None:
    _banner()
    print(f"""
  USAGE:
    python encarols.py  gui             — Launch browser GUI
    python encarols.py  encode  <in>  [out.el]  [--text|--binary]
    python encarols.py  decode  <in.el>  [out]
    python encarols.py  bench   <file>
    python encarols.py  test
""")


def main()->None:
    args=sys.argv[1:]
    if not args: launch_gui(); return
    if args[0] in ('-h','--help','help'): _usage(); return

    cmd=args[0].lower()
    rest=[a for a in args[1:] if not a.startswith('--')]
    flags=[a for a in args[1:] if a.startswith('--')]
    force=(ALGO_BWT if '--text' in flags else
           ALGO_LZ77 if '--binary' in flags else 0)

    if cmd=='gui':
        launch_gui()

    elif cmd in ('encode','enc','e','compress','c'):
        if not rest: print("  Usage: python encarols.py encode <file> [out.el]"); sys.exit(1)
        src=rest[0]; dst=rest[1] if len(rest)>1 else os.path.splitext(src)[0]+EXT
        if not os.path.isfile(src): print(f"  Error: {src} not found"); sys.exit(1)
        _banner()
        print(f"  Compressing  {src}  ({os.path.getsize(src):,} bytes)…", end=' ', flush=True)
        r=encode_file(src,dst,force)
        print("done")
        print(f"\n  ✔  {src}  →  {dst}")
        print(f"     Original   : {r['orig']:>10,} bytes")
        print(f"     Compressed : {r['enc']:>10,} bytes  ({r['ratio']:.1f}%)")
        if r['saved']>0:
            print(f"     Saved      : {r['saved']:>10,} bytes  ({100-r['ratio']:.1f}% smaller)")
        print(f"     Engine     : {r['algo_name']}  ({r['time']:.1f}s)\n")

    elif cmd in ('decode','dec','d','decompress','x'):
        if not rest: print("  Usage: python encarols.py decode <file.el> [out]"); sys.exit(1)
        src=rest[0]; dst=rest[1] if len(rest)>1 else os.path.splitext(src)[0]+'_decoded'
        if not os.path.isfile(src): print(f"  Error: {src} not found"); sys.exit(1)
        _banner()
        print(f"  Decompressing  {src}…", end=' ', flush=True)
        r=decode_file(src,dst)
        print("done")
        print(f"\n  ✔  {src}  →  {dst}")
        print(f"     Restored   : {r['size']:,} bytes  ({r['time']:.1f}s)\n")

    elif cmd in ('bench','b'):
        if not rest: print("  Usage: python encarols.py bench <file>"); sys.exit(1)
        if not os.path.isfile(rest[0]): print(f"  Error: {rest[0]} not found"); sys.exit(1)
        bench_file(rest[0])

    elif cmd=='test':
        run_tests()

    else:
        print(f"  Unknown command: '{cmd}'"); _usage(); sys.exit(1)


if __name__=='__main__':
    main()