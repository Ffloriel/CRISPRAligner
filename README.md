# CRISPRAligner

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python tool for identifying, extracting, and aligning CRISPR spacer sequences from bacterial genomes, with a focus on *Erwinia amylovora* CRISPR repeat regions (CRRs).

## Overview

CRISPRAligner processes bacterial genome FASTA files to:

1. Identify CRISPR repeat regions (CRR1, CRR2, and CRR4)
2. Extract spacer sequences between repeats
3. Align spacers across different genomes
4. Generate comparative analyses in multiple formats (FASTA, CSV, and binary phylip)

The tool supports both individual analyses of specific CRISPR cassettes and combined analyses for phylogenetic applications.

## Installation

### Prerequisites

- Python 3.6 or higher
- Biopython
- NumPy

### Setup

```bash
# Clone the repository
git clone https://github.com/Ffloriel/CRISPRAligner.git
cd CRISPRAligner

# Install the package
pip install .
```

## Usage

### Command Line Interface

```bash
# Basic usage with a directory containing FASTA files
crispr_aligner PREFIX [--input_dir /path/to/fasta/files]
```

Where:
- `PREFIX` is used for naming output files
- `--input_dir` (optional) specifies the directory containing input FASTA files (defaults to current directory)

### Example

```bash
# Process FASTA files in the current directory with prefix "Erwinia"
crispr_aligner Erwinia
```

### Output

Results are organized in the following directory structure:

```
Results/
├── Erwinia.Results.csv       # Summary of CRISPR spacers by genome
├── Erwinia.Error.fasta       # Potential assembly errors
├── AllCRR/                   # Combined CRR results
│   ├── Erwinia.AllCRRAligned.csv
│   ├── Erwinia.AllCRRAligned.fasta
│   └── Erwinia.AllCRRBinary.phy
├── CRR1/                     # CRR1-specific results
│   ├── Erwinia.CRR1.fasta
│   ├── Erwinia.CRR1Aligned.csv
│   ├── Erwinia.CRR1Aligned.fasta
│   ├── Erwinia.CRR1Binary.phy
│   └── Erwinia.CRR1Unique.fasta
├── CRR2/                     # CRR2-specific results
│   └── [Similar files for CRR2]
└── CRR4/                     # CRR4-specific results
    └── [Similar files for CRR4]
```

## Programmatic Usage

You can also use the package in your own Python scripts:

```python
from crispr_aligner import find_crispr_spacers, align_spacers, combine_alignments

# Find CRISPR spacers in genomes
find_crispr_spacers("combined.fasta", "Erwinia.")

# Align spacers for a specific CRISPR repeat region
align_spacers(
    "CRR1/Erwinia.CRR1.fasta",
    "CRR1/Erwinia.CRR1Aligned.csv",
    "CRR1/Erwinia.CRR1Aligned.fasta",
    "CRR1/Erwinia.CRR1Unique.fasta",
    "CRR1/Erwinia.CRR1Binary.phy",
    "CRR1"
)

# Combine alignments from different CRISPR regions
combine_alignments(
    "CRR1/Erwinia.CRR1Aligned.fasta",
    "CRR2/Erwinia.CRR2Aligned.fasta",
    "CRR4/Erwinia.CRR4Aligned.fasta",
    "AllCRR/Erwinia.AllCRRAligned.fasta",
    # Additional parameters...
)
```

## Project Structure

```
CRISPRAligner/
├── crispr_aligner/
│   ├── __init__.py                # Package initialization
│   ├── cli.py                     # Command-line interface
│   ├── sequences/                 # FASTA file handling
│   │   ├── __init__.py
│   │   └── fasta_utils.py
│   ├── finder/                    # CRISPR spacer identification
│   │   ├── __init__.py
│   │   └── crr_finder.py
│   ├── aligner/                   # Spacer alignment
│   │   ├── __init__.py
│   │   └── crr_aligner.py
│   └── combiner/                  # Result combination
│       ├── __init__.py
│       └── crr_combiner.py
├── tests/                         # Test suite
│   └── ...
├── README.md
└── pyproject.toml                 # Package configuration
```

## Technical Details

### CRISPR Repeat Regions

CRISPRAligner identifies three different CRISPR repeat regions in *Erwinia amylovora*:

1. **CRR1**: Consensus sequence `GTGTTCCCCGCGTGAGCGGGGATAAACCG`
2. **CRR2**: Consensus sequence `GTGTTCCCCGCGTATGCGGGGATAAACCG`
3. **CRR4**: Consensus sequence `GTTCACTGCCGTACAGGCAGCTTAGAAA`

### Analysis Process

1. **Finding**: Sequences are scanned for matches to CRISPR repeat consensus sequences
2. **Extraction**: Spacers between repeats are extracted
3. **Alignment**: Spacers are aligned across different genomes
4. **Combination**: Results from different CRR types are combined for comparative analysis

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this tool in your research, please reference:

```
Parcey M, Gayder S, Castle AJ, Svircev AM. 2022. Function and Application of the CRISPR-Cas System in the Plant Pathogen Erwinia amylovora. Appl Envir Micro. e02513-21. doi: 10.1128/aem.02513-21
```

## Credits

Originally developed by Michael Parcey at Brock University.  
Reference: Parcey M, Gayder S, Castle AJ, Svircev AM. 2022. Function and Application of the CRISPR-Cas System in the Plant Pathogen Erwinia amylovora. Appl Envir Micro. e02513-21. doi: 10.1128/aem.02513-21