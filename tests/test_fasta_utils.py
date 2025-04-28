"""Unit tests for FASTA file processing utilities."""

import pytest
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crispr_aligner.sequences.fasta_utils import combine_fasta_files


@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory with test FASTA files."""
    # Create Results subdirectory
    results_dir = tmp_path / "Results"
    results_dir.mkdir(exist_ok=True)
    
    # Create test FASTA files
    create_test_fasta(tmp_path / "test1.fasta", [
        ("genome1", "ACGTACGT"),
        ("genome2", "TGCATGCA")
    ])
    
    create_test_fasta(tmp_path / "test2.fasta", [
        ("genome3", "GGGCCCAAA"),
        ("genome4", "TTTAAACCC")
    ])
    
    return tmp_path


def create_test_fasta(path, sequences):
    """Helper function to create test FASTA files."""
    with open(path, 'w') as f:
        for seq_id, seq in sequences:
            record = SeqRecord(Seq(seq), id=seq_id, description="")
            f.write(record.format("fasta"))


def test_combine_fasta_files_success(temp_dir):
    """Test successful combination of FASTA files."""
    # Call the function
    output_path = combine_fasta_files(temp_dir, "combined.fasta", "test")
    
    # Verify the output path
    assert output_path == temp_dir / "Results" / "combined.fasta"
    assert output_path.exists()
    
    # Check the content of the combined file
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 4
    
    # Verify record IDs
    record_ids = [record.id for record in records]
    assert "genome1" in record_ids
    assert "genome2" in record_ids
    assert "genome3" in record_ids
    assert "genome4" in record_ids


def test_combine_fasta_files_no_files(tmp_path):
    """Test error when no FASTA files found."""
    # Create Results directory
    results_dir = tmp_path / "Results"
    results_dir.mkdir(exist_ok=True)
    
    # Should raise FileNotFoundError
    with pytest.raises(FileNotFoundError):
        combine_fasta_files(tmp_path, "combined.fasta", "test")


def test_combine_fasta_files_empty_file(tmp_path):
    """Test handling of empty FASTA files."""
    # Create Results directory
    results_dir = tmp_path / "Results"
    results_dir.mkdir(exist_ok=True)
    
    # Create empty FASTA file
    empty_file = tmp_path / "empty.fasta"
    empty_file.touch()
    
    # Call the function
    output_path = combine_fasta_files(tmp_path, "combined.fasta", "test")
    
    # Verify the output
    assert output_path.exists()
    
    # Should be an empty file (no records)
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 0


def test_combine_fasta_files_existing_output(temp_dir):
    """Test that existing output file is overwritten."""
    # Create an existing output file
    output_path = temp_dir / "Results" / "combined.fasta"
    with output_path.open("w") as f:
        f.write(">existing\nACGT\n")
    
    # Call the function
    combine_fasta_files(temp_dir, "combined.fasta", "test")
    
    # Verify that the file was overwritten
    records = list(SeqIO.parse(output_path, "fasta"))
    assert len(records) == 4
    record_ids = [record.id for record in records]
    assert "existing" not in record_ids