# CRISPRAligner

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Python tool for identifying, extracting, and aligning CRISPR spacer sequences from bacterial genomes, with a focus on CRISPR repeat regions (CRRs) in different bacteria.

## Overview

CRISPRAligner processes bacterial genome FASTA files to:

1. Identify CRISPR repeat regions (like CRR1, CRR2, and when available, CRR4)
2. Group CRISPR arrays by spatial proximity to detect separate arrays with letter suffixes (e.g., CRR1A, CRR1B)
3. Extract spacer sequences between repeats
4. Align spacers across different genomes
5. Generate comparative analyses in multiple formats (FASTA, CSV, and binary phylip)

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
crispr_aligner PREFIX [--input_dir /path/to/fasta/files] [--bacteria SPECIES] [--custom_config FILE]
```

Where:
- `PREFIX` is used for naming output files
- `--input_dir` (optional) specifies the directory containing input FASTA files (defaults to current directory)
- `--bacteria` (optional) species name to load built-in configurations (defaults to "erwinia")
- `--custom_config` (optional) path to a custom JSON configuration file
- `--min_crr1_score`, `--min_crr2_score`, `--min_crr4_score` (optional) minimum match scores
- `--max_frame_size` (optional) maximum size in bp for a CRISPR array frame (default: 20000)

### Examples

```bash
# Process FASTA files in the current directory with prefix "Erwinia"
crispr_aligner Erwinia

# Use built-in Phytobacter configuration (no CRR4)
crispr_aligner Phytobacter --bacteria phytobacter

# Use a custom configuration file
crispr_aligner MyBacteria --custom_config my_species_config.json

# Adjust minimum match scores
crispr_aligner Erwinia --min_crr1_score 20 --min_crr4_score 18

# Change max frame size for spatially separating CRISPR arrays
crispr_aligner Erwinia --max_frame_size 15000
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
├── CRR1A/                    # First spatial frame of CRR1 results
│   ├── Erwinia.CRR1A.fasta
│   ├── Erwinia.CRR1AAligned.csv
│   ├── Erwinia.CRR1AAligned.fasta
│   ├── Erwinia.CRR1ABinary.phy
│   └── Erwinia.CRR1AUnique.fasta
├── CRR1B/                    # Second spatial frame of CRR1 (if found)
│   └── [Similar files for CRR1B]
├── CRR2A/                    # First spatial frame of CRR2
│   └── [Similar files for CRR2A]
└── CRR4A/                    # First spatial frame of CRR4 (if available)
    └── [Similar files for CRR4A]
```

## Configuration System

CRISPRAligner supports different bacterial species through a flexible configuration system.

### Built-in Configurations

The tool includes built-in configurations for several bacteria:
- `erwinia` (default): *Erwinia amylovora* (CRR1, CRR2, CRR4)
- `phytobacter`: *Phytobacter* species (CRR1, CRR2 only, no CRR4)
- And others (use `--bacteria` parameter to specify)

### Custom Configuration Files

For other species, you can create a custom JSON configuration file:

```json
{
  "consensus_sequences": {
    "crr1": "GTGTTCCCCGCGTATGCGGGGATAAACCG",
    "crr2": "GTGCTTCACTGCATGCGTTGATAAACCG",
    "crr4": "GTTCACTGCCGTACAGGCAGCTTAGAAA"
  },
  "short_checks": {
    "crr1_crr2_start": "GTGTTC",
    "crr1_crr2_end": "ATAAACC",
    "crr4_start": "GTTCAC",
    "crr4_middle": "GTACAGG"
  },
  "available_crrs": ["crr1", "crr2", "crr4"]
}
```

The configuration contains:
- `consensus_sequences`: The CRR consensus sequences to search for
- `short_checks`: Short sequences used for initial pattern matching
- `available_crrs`: List of CRR types available for this bacteria (omit unavailable ones)

## Programmatic Usage

You can also use the package in your own Python scripts:

```python
from crispr_aligner import find_crispr_spacers, align_spacers
from crispr_aligner.combiner.crr_combiner import combine_fasta, combine_csv, combine_binary
from crispr_aligner.utils import load_bacteria_config

# Load a configuration for a specific bacteria
config = load_bacteria_config("phytobacter")

# Find CRISPR spacers in genomes with specific configuration
find_crispr_spacers(
    "combined.fasta", 
    "MyPrefix.", 
    results_dir="Results",
    consensus_sequences=config["consensus_sequences"],
    short_checks=config["short_checks"],
    available_crrs=config["available_crrs"],
    min_scores={"crr1": 22, "crr2": 22, "crr4": 21},
    max_frame_size=20000
)

# Get all frame directories created by find_crispr_spacers
frame_dirs = []
for item in Path("Results").iterdir():
    if item.is_dir() and item.name not in ["AllCRR"]:
        frame_dirs.append(item.name)

# Align spacers for each frame
for frame in frame_dirs:
    align_spacers(
        f"Results/{frame}/MyPrefix.{frame}.fasta",
        f"Results/{frame}/MyPrefix.{frame}Aligned.csv",
        f"Results/{frame}/MyPrefix.{frame}Aligned.fasta",
        f"Results/{frame}/MyPrefix.{frame}Unique.fasta",
        f"Results/{frame}/MyPrefix.{frame}Binary.phy",
        frame
    )

# Collect files for combination
fasta_files = [f"Results/{frame}/MyPrefix.{frame}Aligned.fasta" for frame in frame_dirs]
csv_files = [f"Results/{frame}/MyPrefix.{frame}Aligned.csv" for frame in frame_dirs]
phy_files = [f"Results/{frame}/MyPrefix.{frame}Binary.phy" for frame in frame_dirs]

# Combine alignments from different frames
combine_fasta(fasta_files, "Results/AllCRR/MyPrefix.AllCRRAligned.fasta")
combine_csv(csv_files, "Results/AllCRR/MyPrefix.AllCRRAligned.csv")
combine_binary(phy_files, "Results/AllCRR/MyPrefix.AllCRRBinary.phy")
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
│   ├── combiner/                  # Result combination
│   │   ├── __init__.py
│   │   └── crr_combiner.py
│   └── utils/                     # Utility functions
│       ├── __init__.py
│       └── config.py              # Configuration management
├── tests/                         # Test suite
│   └── ...
├── README.md
└── pyproject.toml                 # Package configuration
```

## Technical Details

### CRISPR Repeat Regions

CRISPRAligner identifies different CRISPR repeat regions depending on the bacteria:

For *Erwinia amylovora*:
1. **CRR1**: Consensus sequence `GTGTTCCCCGCGTGAGCGGGGATAAACCG`
2. **CRR2**: Consensus sequence `GTGTTCCCCGCGTATGCGGGGATAAACCG`
3. **CRR4**: Consensus sequence `GTTCACTGCCGTACAGGCAGCTTAGAAA`

For *Phytobacter* species:
1. **CRR1**: Consensus sequence `GTGTTCCCCGCGCGAGCGGGGATAAACCG`
2. **CRR2**: Consensus sequence `GTGTTCCCCGCGCCAGCGGGGATAAACCG`
(CRR4 is not present in Phytobacter)

Other bacteria will have different consensus sequences, which can be specified through the configuration system.

### Frame-Based CRISPR Detection

CRISPRAligner now groups CRISPR arrays into frames based on spatial proximity:

1. CRISPR repeats of the same type (e.g., CRR1) are grouped into a frame if they are within a maximum distance (default: 20,000 bp)
2. Repeats separated by more than the maximum distance are considered separate arrays and assigned different letter suffixes (e.g., CRR1A, CRR1B)
3. Each frame is processed separately for alignment and analysis
4. Results can be combined across frames for comprehensive analysis

This approach preserves the biological reality that spatially separate CRISPR arrays of the same type are functionally distinct units.

### Analysis Process

1. **Finding**: Sequences are scanned for matches to CRISPR repeat consensus sequences
2. **Framing**: Repeats are grouped into spatial frames (e.g., CRR1A, CRR1B)
3. **Extraction**: Spacers between repeats are extracted
4. **Alignment**: Spacers are aligned across different genomes, separately for each frame
5. **Combination**: Results from different frames are combined for comparative analysis

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
