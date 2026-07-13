# SeqViz 🧬

**Pairwise DNA sequence alignment and heatmap visualisation from the command line.**

Python · Needleman-Wunsch · Smith-Waterman · Biopython · Matplotlib

[![Python](https://img.shields.io/badge/python-3.11%2B-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)]()
[![Development Status](https://img.shields.io/badge/status-alpha-orange)]()

---

## What It Does
SeqViz is capable of taking two seperate DNA sequences, either typed into the terminal or loaded from standard .fasta files. Then aligning these sequences using two well known, and tested, algorithms in bioinformatics. It then renders a heatmap image to visualise the difference between your sequences as a heatmap image and prints the alignment to your terminal for quick and easy visual identification of where the match/mismatch and gaps are and how strong the connection is between DNA sequences. You can compare two genes from the same species, or track evolotionary divergence of two seperate species using cytochrome b with this tool. SeqViz gives a numerical and visual answer to you with a simple command in your terminal and does it locally.

---

## Installation

**Requirements:** Python 3.11 or higher.

# 1. Clone the repository
```bash
git clone https://github.com/BarefootRobber/SeqViz.git
cd SeqViz
```
# 2. Create and activate a virtual environment (recommended)
```bash
python -m venv .venv
source .venv/bin/activate        # macOS / Linux
.venv\Scripts\activate           # Windows
```
# 3. Install SeqViz and all dependencies
```bash
pip install -e ".[dev]"
```
# 4. Confirm the install worked
```bash
seqviz --help
```
**Dependencies installed automatically:**
```
| Package | Purpose |
|---|---|
| `biopython` | FASTA file parsing |
| `numpy` | Scoring matrix construction |
| `matplotlib` + `seaborn` | Heatmap rendering |
| `rich` | Formatted terminal output |
| `click` | CLI argument handling |
```
---

## Usage

### Align two raw sequences

# Needleman-Wunsch global alignment (the classic textbook benchmark) (MATCH = 1, MISMATCH = -1, GAP = -1)
```bash
seqviz --seq1 GATTACA --seq2 GCATGCU --algorithm nw
```
```
──────────────────── Needleman-Wunsch (Global) ─────────────────────
Score: 0
G-ATTACA
GCA-TGCU
```
---

# Smith-Waterman local alignment (MATCH = 2, MISMATCH = -1, GAP = -1)
```bash
seqviz --seq1 ACACACTA --seq2 AGCACACA --algorithm sw
```
```
────────────────────── Smith-Waterman (Local) ──────────────────────
Score: 12
A-CACACTA
AGCACAC-A
```
---

### Align from FASTA files (MATCH = 1, MISMATCH = -1, GAP = -1)

# Global alignment of two real cytochrome b gene sequences
```bash
seqviz --file1 examples/gallus_gallus_cytb.fasta \
       --file2 examples/columba_livia_cytb.fasta \
       --algorithm nw \
       --output results/chicken_vs_pigeon.png
```
```
──────────────────── Needleman-Wunsch (Global) ─────────────────────
Score: 221
AGCATGATGAAA-TTTCGGCTCCCTATTA-GCAGTCTGCCT-CATGACCCAAATCCTCACCGGCCTAC
TACTAGCCATGCAC--TACACAGCAGACA-CATCCCTAGCCTTCTCCTCCGTAGCCCACACTTGCCGG
AACGTACAATACGGCTGACTCATCCGGAATCTCCACGCAAACGGCGCCTCATTCTTCTTCATCTGTAT
CTT-CCTTCACATCGGACGAGGCCTATACTACGGCTCCTACCTCTACAAGGAAACCTGAAACACAGGA
GTAATCCTCCTCCTCACACTCATAGCCACCGCCTTTGTGGGCTATGTTCT-CCCATGAGGACAAATA
AGCATGATGAAACTTT-GGGTCCCTACTAGGCA-TTTGCTTGC-TAACTCAAATCCTAACCGGCTTAC
TACTCGCC--GCACATTACACTGCAGACACCA-CCCTAGCCTTTTCATCCGTTGCACACACATGCCGA
AACGTACAGTACGGCTGGCTAATCCGAAACCTCCATGCAAACGGAGCCTCATTTTTCTTCATCTGTAT
-TTACCTACACATCGGACGAGGACTCTACTACGGATCCTACCTCTACAAAGAGACTTGAAACACAGGA
GTCGTCCTCCTACTAACCCTTATAGCCACTGCATTCGTAGGATATGTCCTACCC-TGAGGACAAATA

  Heatmap saved → /Users/ryanyoung/Documents/GitHub/SeqViz/results/chicken_vs_pigeon.png
```
---

### Run both algorithms side by side
```bash
seqviz --file1 examples/gallus_gallus_cytb.fasta \
       --file2 examples/coturnix_japonica_cytb.fasta \
       --algorithm both \
       --output results/comparison.png
```
---

### Full options reference
```
Options:
  --seq1        TEXT     First sequence as a raw string
  --seq2        TEXT     Second sequence as a raw string
  --file1       PATH     FASTA file for sequence 1
  --file2       PATH     FASTA file for sequence 2
  --algorithm   TEXT     nw · sw · both  [default: nw]
  --match       INT      Match score     [default: 1]
  --mismatch    INT      Mismatch penalty [default: -1]
  --gap         INT      Gap penalty     [default: -1]
  --output      PATH     Save heatmap to this path
  --help                 Show this message and exit
```
---

## Screenshots

### Heatmap — Gallus gallus vs Columba livia (cytochrome b, NW)

![Heatmap output](/README_screenshots/chicken_vs_pigeon.png)

The colour gradient runs from **dark blue** (low / negative score) through
**white** (neutral) to **dark red** (high match score). The diagonal streak
of warm colours shows the path of highest similarity — the optimal global
alignment path traced by the Needleman-Wunsch algorithm.

---

### Terminal output

![Terminal output](/README_screenshots/chicken_vs_pigeon_terminal.png)

---

## Background — Needleman-Wunsch vs Smith-Waterman

These algorithms both function as pairwise sequence alignment tools, with a subtle difference between them. Needleman-Wunsch (nw) is allowed to consider negative scores when filling the matrix and the Smith-Waterman (sw) algorithm floors the lowest score to 0. This is what allows the sw algorithm to find the best local alignment between two sequences, by flooring the score sw essentially ignores the areas of the sequences that do not match, only highlighting the regions of high similarity. In contrast the nw algorithim is forced to consider every pairwise cell and assign a score to it, so it can be traced back to the origin cell of the matrix.

Nw was created in the 1970s to perform global alignment of two sequences. It fills the matrix in an `(m+1) × (n+1)` where the intial row and column of each sequence is filled with cumulative gap penalties, from there every cell takes the maxximum of three possible scores; Match, being the highest scoring option for matching nucleotides at that position, Mismatch, the opposite score of the match normally (these occur when the maximum score is sourced from the top left cell). Gap penalties, usually -1 which represents an insert or deletion from sequence 1 (highest score comes from the cell to the left), or sequence 2 (highest score is sourced from the top cell). Optimal alignment is found by tracing back the matrix from the bottom left hand corner to the top right hand corner following the source of each cell.

Sw was created in 1981 in order to perform local alignment between two sequences. It is a derivetive from the nw allgorithim with an elegant change in logic. It floors all cell scores to 0, meaning no cell can ever be negative, this change means that the conserved region of the sequences will not be missed due to the less conserved regions around them. The traceback in a sw algorithim also starts from the highest scoring cell, back up and left until the nearest 0 cell is reached.

Both algorithms use the same user-configurable scoring scheme, found in scoring.py:

| Event | Default score | Biological meaning |
|---|---|---|
| Match | +1 | Same nucleotide at this position |
| Mismatch | -1 | Different nucleotide — possible point mutation |
| Gap | -1 | Insertion or deletion event (indel) |

---

## Example Data

The `examples/` directory contains real cytochrome b (MT-CYB) gene sequences
downloaded from NCBI for four bird species, chosen to span a range of
evolutionary distances:

| File | Species | NCBI Accession | Order |
|---|---|---|---|
| `gallus_gallus_cytb.fasta` | Chicken | KF964328 | Galliformes |
| `coturnix_japonica_cytb.fasta` | Japanese quail | KF964327 | Galliformes |
| `columba_livia_cytb.fasta` | Rock pigeon | KF964326 | Columbiformes |
| `streptopelia_senegalensis_cytb.fasta` | Laughing dove | KF964325 | Columbiformes |

The chicken–quail pair (same order) shows higher sequence similarity than
the chicken–pigeon pair (different orders, ~90 million years of divergence),
which is visible in the heatmap colour gradient without any additional analysis.

To re-download these sequences from NCBI:
```bash
python scripts/download_examples.py
```
---

## Running Tests

# Full test suite with coverage
```bash
pytest -v
```
# FASTA parsing tests only
```shell
pytest tests/test_fasta.py -v
```
# Alignment algorithm tests only
```bash
pytest tests/test_alignment.py -v
```
# Stop at first failure
```shell
pytest -v -x
```
---

## Project Structure
```
SeqViz/
├── src/seqviz/
│   ├── cli.py              # Typer CLI entry point
│   ├── needleman_wunsch.py # NW matrix fill + traceback
│   ├── smith_waterman.py   # SW matrix fill + traceback
│   ├── fasta_utils.py      # FASTA parsing
│   ├── scoring.py          # Shared scoring function + defaults
│   └── visualisation.py    # Heatmap and terminal rendering
├── tests/
│   ├── test_alignment.py
│   ├── test_fasta.py
│   └── fixtures/           ***# Synthetic FASTA files for testing***
├── examples/               # Real NCBI cytochrome b sequences
├── scripts/
│   └── download_examples.py
└── pyproject.toml
```
---

## Future Plans

- **BLOSUM62 support** — extend to protein sequence alignment using the
  standard substitution matrix used in BLAST, enabling amino acid-level
  comparison in addition to nucleotide alignment
- **Multiple sequence alignment** — align more than two sequences
  simultaneously using progressive alignment (ClustalW-style), a natural
  extension once pairwise alignment is solid
- **Web interface** — a lightweight FastAPI or Streamlit front end so the
  tool can be used without installing anything locally; relevant for
  sharing with wet-lab collaborators
- **Affine gap penalties** — replace the current linear gap penalty with a
  gap-open + gap-extend model, which is biologically more realistic (a
  single long indel is more likely than many short ones)
- **E-value estimation** — add statistical significance scoring so users
  know whether an alignment score is meaningful or could arise by chance,
  bringing SeqViz closer to BLAST-level reporting

---

## Author

**Ryan Young** — Zoologist transitioning to bioinformatics.

[GitHub](https://github.com/BarefootRobber/SeqViz) ·
[Issues](https://github.com/BarefootRobber/SeqViz/issues)

---

## License

MIT — see [LICENSE](LICENSE) for details.
