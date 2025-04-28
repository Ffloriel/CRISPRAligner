"""FASTA file processing utilities for CRISPR analysis."""

from pathlib import Path
from typing import Union
from Bio import SeqIO


def combine_fasta_files(
    input_dir: Union[str, Path], output_filename: str, prefix: str
) -> Path:
    """
    Combine multiple FASTA files into a single output file.

    Args:
        input_dir: Directory containing input FASTA files
        output_filename: Name for the combined output file
        prefix: Prefix for output file naming

    Returns:
        Path to the combined file

    Raises:
        FileNotFoundError: If no FASTA files found in input directory
    """
    input_dir = Path(input_dir)
    output_path = input_dir / "Results" / output_filename

    # Make sure we're starting with a fresh file
    if output_path.exists():
        output_path.unlink()

    # Find all FASTA files in the input directory
    fasta_files = list(input_dir.glob("*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {input_dir}")

    # Process each FASTA file
    combined_count = 0
    with output_path.open("w") as outfile:
        for fasta_file in fasta_files:
            print(f"Processing {fasta_file.name}...")
            try:
                with fasta_file.open("r") as infile:
                    for record in SeqIO.parse(infile, "fasta"):
                        outfile.write(record.format("fasta"))
                        combined_count += 1
            except Exception as e:
                print(f"Error processing {fasta_file.name}: {e}")

    print(
        f"Combined {combined_count} sequences from {len(fasta_files)} files into {output_path}"
    )
    return output_path
