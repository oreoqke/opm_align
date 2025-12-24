# OPM Align - Protein Structure Matcher

A bioinformatics tool for finding the closest structural and sequence matches for new protein structures in the [OPM (Orientations of Proteins in Membranes)](https://opm.phar.umich.edu/) database.

## Overview

OPM Align combines sequence alignment (BLAST) with structural superposition (USalign) to identify and classify protein matches. It's particularly useful for transmembrane protein analysis, automatically identifying the subunit with the most transmembrane regions for optimal matching.

### Features

- **Sequence alignment** using NCBI BLAST for initial screening
- **Structural superposition** using USalign for 3D comparison
- **Transmembrane detection** via UniProt API integration
- **Match classification** into 5 types based on identity, overlap, and RMSD
- **Cross-platform** support (Linux/macOS)
- **Multi-threaded** processing for faster analysis

## Requirements

- Python 3.8+
- [NCBI BLAST 2.15.0+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
- [USalign](https://zhanggroup.org/US-align/)
- PostgreSQL (for OPM database access)

## Installation

### 1. Clone the repository

```bash
git clone git@github.com:oreoqke/opm_align.git
cd opm_align
```

### 2. Install NCBI BLAST

Download from: https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html

**macOS (Homebrew):**
```bash
brew install blast
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install ncbi-blast+
```

Verify installation:
```bash
blastp -version
```

### 3. Install USalign

Download from: https://zhanggroup.org/US-align/

Either add to your PATH or place the `USalign` executable in the project directory. The tool will automatically search for it.

### 4. Set up Python environment

```bash
# Create virtual environment
python3 -m venv env

# Activate (Linux/macOS)
source env/bin/activate

# Install dependencies
pip install -r requirements.txt
```

### 5. Configure database connection

Edit `get_all_pdbid.py` to specify your OPM database name and credentials if needed.

## Usage

### Basic usage

```bash
python main.py -i new_pdb.txt -o result.csv -d results
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | `new_pdb.txt` | File containing list of PDB IDs to match |
| `-o, --output` | `result.csv` | Output file for sequence alignment results |
| `-d, --output-dir` | `results` | Directory for all output files |

### Input format

Create a text file with one PDB ID per line:
```
8A9Y
7XYZ
6ABC
```

### Output files

| File | Description |
|------|-------------|
| `results.csv` | Complete USalign results |
| `selection.csv` | Classified matches with type annotations |
| `top_matches.csv` | Best match per PDB (lowest RMSD) |
| `to_align.txt` | Pairs selected for structural alignment |
| `bad_matches.txt` | PDB IDs with low sequence similarity (<15%) |

### Match classification

| Type | Description |
|------|-------------|
| 1 | Same protein (high identity, overlap, low RMSD) |
| 2 | High identity but lower structural overlap |
| 3 | Same family (medium identity and overlap) |
| 4 | Family match (lower identity, same overlap) |
| 0 | No match (filtered out) |

## Cleanup utility

Use the cross-platform cleanup script to remove generated files:

```bash
python cleanup.py fasta   # Remove FASTA files
python cleanup.py log     # Remove log files
python cleanup.py db      # Remove BLAST database files
python cleanup.py all     # Remove all results
```

## Project structure

```
opm_align/
├── main.py                  # Main entry point
├── align_structures.py      # USalign integration & classification
├── transmembrane_region.py  # UniProt transmembrane detection
├── get_all_fasta.py         # FASTA sequence downloader
├── get_all_pdbid.py         # OPM database connector
├── cleanup.py               # Cross-platform cleanup utility
└── requirements.txt         # Python dependencies
```

## Authors

- **Alexey Kovalenko** ([@oreoqke](https://github.com/oreoqke))
- **Stanislav Cherepanov**

## License

MIT License - see [LICENSE](LICENSE) for details.

## References

- [NCBI BLAST documentation](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
- [USalign](https://zhanggroup.org/US-align/)
- [OPM Database](https://opm.phar.umich.edu/)
