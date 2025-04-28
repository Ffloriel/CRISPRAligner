"""CRISPRAligner - Extract and align CRISPR spacers from bacterial genomes.

This package provides tools for processing FASTA files and identifying CRISPR
arrays in bacterial genomes.
"""

__version__ = "0.1.0"

from crispr_aligner.sequences import combine_fasta_files
from crispr_aligner.finder import find_crispr_spacers
from crispr_aligner.aligner import align_spacers
from crispr_aligner.combiner import combine_alignments
