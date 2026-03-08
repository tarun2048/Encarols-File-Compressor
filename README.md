#  Encarols v5.0

**Adaptive Dual-Engine File Compressor with Browser GUI**

Encarols is a pure-Python file compression tool that automatically selects the best compression algorithm for your file — a sophisticated BWT-based pipeline for text files, and LZ77+Huffman for binary files. No external dependencies required.

---

## ✨ Features

- **Dual-Engine Compression** — automatically detects whether your file is text or binary and picks the optimal engine
- **Engine A (Text):** `RLE1 → BWT → MTF → RUNA/RUNB → Multi-table Huffman` — inspired by bzip2, great for logs, source code, CSV, and other text files
- **Engine B (Binary):** `LZ77 sliding window (32 KB) + Canonical Huffman` — efficient for executables, images, and other binary formats
- **Browser GUI** — a beautiful, built-in web interface with drag-and-drop, progress bar, and session history
- **CLI** — full command-line interface with encode, decode, benchmark, and test commands
- **Zero Dependencies** — runs on Python's standard library only
- **Custom `.el` format** — compact file header with magic bytes for reliable identification

---

## 🚀 Quick Start

**No installation needed.** Just run:

```bash
python encarols.py
```

This launches the browser GUI automatically at `http://localhost:8765`.

---

## 📋 Requirements

- Python 3.7 or higher
- No third-party packages needed

---

## 🖥️ GUI Usage

```bash
python encarols.py gui
```

- Opens your default browser at `http://localhost:8765`
- Drag and drop any file to compress or decompress
- Choose engine manually (Auto / Text / Binary) or let Encarols decide
- Download the result directly from the browser
- Session history shows all operations with compression ratios

---

## ⌨️ CLI Usage

### Compress a file

```bash
python encarols.py encode <input_file> [output.el] [--text|--binary]
```

| Flag | Description |
|------|-------------|
| `--text` | Force Engine A (BWT pipeline) |
| `--binary` | Force Engine B (LZ77+Huffman) |
| *(none)* | Auto-detect best engine |

**Examples:**
```bash
python encarols.py encode document.txt
python encarols.py encode document.txt archive.el --text
python encarols.py encode program.exe program.el --binary
```

### Decompress a file

```bash
python encarols.py decode <input.el> [output_file]
```

**Example:**
```bash
python encarols.py decode archive.el restored.txt
```

### Benchmark a file

```bash
python encarols.py bench <file>
```

Runs both engines and prints a side-by-side performance and ratio comparison.

### Run self-tests

```bash
python encarols.py test
```

Runs built-in correctness tests to verify encode/decode round-trips.

### Help

```bash
python encarols.py --help
```

---

## 🔧 How It Works

### Engine A — Text Files

Designed for high compression of structured text using a bzip2-inspired pipeline:

| Stage | Description |
|-------|-------------|
| **RLE1** | Pre-encodes runs of ≥4 identical bytes |
| **BWT** | Burrows-Wheeler Transform — reorders bytes so similar contexts cluster |
| **MTF** | Move-to-Front encoding — converts clusters to near-zero values |
| **RUNA/RUNB** | Zero-run encoding inside the Huffman alphabet (RUNA/RUNB get shortest codes) |
| **Multi-table Huffman** | 2–6 dynamically optimised Huffman tables across 50-symbol groups |

### Engine B — Binary Files

Uses a classic LZ77 sliding window approach:

| Stage | Description |
|-------|-------------|
| **LZ77** | 32 KB sliding window with chain hashing, min match 3 bytes, max 258 |
| **Canonical Huffman** | Compact code storage using length arrays |

### Auto-Detection

Encarols analyses your file's byte distribution to automatically choose the best engine before compression begins.

---

## 📦 File Format

Compressed files use the `.el` extension and begin with a 4-byte magic header:

```
Magic: 0x45 0x4C 0x05 0x00
```

The fifth byte identifies the engine used (`0x01` = BWT, `0x02` = LZ77), enabling reliable decompression without user input.

---

## 💡 Tips

- Text files (`.txt`, `.log`, `.csv`, `.py`, `.json`, `.xml`) compress best with Engine A
- Binary files (`.exe`, `.bin`, `.mp3`, `.db`) compress better with Engine B
- Use `bench` to compare both engines on a specific file before committing
- The GUI's Auto mode is recommended for most users

---

## 📁 Project Structure

```
encarols.py      # Single-file application — compressor, decompressor, and GUI
```

---

## 👨‍💻 Developer

**Tarun**  
*Original engine design by Eaevox*

---

## 📄 License

This project is open source. Feel free to use, modify, and distribute.

---

> *Encarols v5.0 — Adaptive Dual-Engine Compressor*
