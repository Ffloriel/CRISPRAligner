"""Unit tests for the CRISPR repeat region (CRR) finder module."""

import pytest
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from crispr_aligner.finder.crr_finder import (
    find_crr12_match,
    find_crr4_match,
    count_fasta_records,
    process_sequence,
    find_crispr_spacers
)


@pytest.fixture
def perfect_crr1():
    """Return a perfect match sequence for CRR1."""
    return list("GTGTTCCCCGCGTGAGCGGGGATAAACCG")


@pytest.fixture
def perfect_crr2():
    """Return a perfect match sequence for CRR2."""
    return list("GTGTTCCCCGCGTATGCGGGGATAAACCG")


@pytest.fixture
def perfect_crr4():
    """Return a perfect match sequence for CRR4."""
    return list("GTTCACTGCCGTACAGGCAGCTTAGAAA")


@pytest.fixture
def partial_match_crr1():
    """Return a sequence with partial match to CRR1."""
    return list("GTGTTCCCCGCGTGAGCGGGXXTAAACCG")


@pytest.fixture
def sample_fasta_file(tmp_path):
    """Create a sample FASTA file with known CRISPR repeats for testing."""
    # Create a test file with known CRR patterns
    fasta_path = tmp_path / "test_genome.fasta"
    
    # Create a genome with CRR1, CRR2, and CRR4 patterns
    # Including repeats and spacers for each
    with fasta_path.open("w") as f:
        # Write header
        f.write(">TestGenome1\n")
        
        # Create a sequence with CRR1 repeat, spacer, repeat pattern
        genome_seq = "N" * 100  # Prefix
        
        # Add CRR1 pattern with spacer
        genome_seq += "GTGTTCCCCGCGTGAGCGGGGATAAACCG"  # CRR1 repeat
        genome_seq += "ACTAGCTAGCTAGCTAGCTAGCTAGCT"    # Spacer (27bp)
        genome_seq += "GTGTTCCCCGCGTGAGCGGGGATAAACCG"  # CRR1 repeat
        
        # Add some distance between patterns
        genome_seq += "N" * 100
        
        # Add CRR2 pattern with spacer
        genome_seq += "GTGTTCCCCGCGTATGCGGGGATAAACCG"  # CRR2 repeat
        genome_seq += "GCATGCATGCATGCATGCATGCAT"       # Spacer (24bp)
        genome_seq += "GTGTTCCCCGCGTATGCGGGGATAAACCG"  # CRR2 repeat
        
        # Add some distance between patterns
        genome_seq += "N" * 100
        
        # Add CRR4 pattern with spacer
        genome_seq += "GTTCACTGCCGTACAGGCAGCTTAGAAA"   # CRR4 repeat
        genome_seq += "TGACTGACTGACTGACTGACTGACT"      # Spacer (24bp)
        genome_seq += "GTTCACTGCCGTACAGGCAGCTTAGAAA"   # CRR4 repeat
        
        # Add suffix
        genome_seq += "N" * 100
        
        # Write the sequence with line breaks
        for i in range(0, len(genome_seq), 80):
            f.write(genome_seq[i:i+80] + "\n")
    
    return fasta_path


@pytest.fixture
def results_dir(tmp_path):
    """Create a temporary results directory."""
    results_dir = tmp_path / "Results"
    results_dir.mkdir(exist_ok=True)
    
    # Create subdirectories
    for subdir in ["AllCRR", "CRR1", "CRR2", "CRR4"]:
        (results_dir / subdir).mkdir(exist_ok=True)
    
    return results_dir


def test_find_crr12_match_perfect_crr1(perfect_crr1):
    """Test perfect match for CRR1."""
    result = find_crr12_match(perfect_crr1)
    assert result == 29  # Should be a positive score with all 29 matches


def test_find_crr12_match_perfect_crr2(perfect_crr2):
    """Test perfect match for CRR2."""
    result = find_crr12_match(perfect_crr2)
    assert result == -29  # Should be a negative score with all 29 matches


def test_find_crr12_match_partial(partial_match_crr1):
    """Test partial match."""
    result = find_crr12_match(partial_match_crr1)
    assert 0 < result < 29  # Should be positive but less than 29


def test_find_crr12_match_short_sequence():
    """Test with sequence shorter than consensus."""
    short_seq = list("GTGTTCCC")  # Only 8 characters
    result = find_crr12_match(short_seq)
    assert result == 8  # Should match only the available characters


def test_find_crr4_match_perfect(perfect_crr4):
    """Test perfect match for CRR4."""
    result = find_crr4_match(perfect_crr4)
    assert result == 28  # Should match all 28 characters


def test_find_crr4_match_partial():
    """Test partial match for CRR4."""
    partial = list("GTTCACTGCCGTACAGGCAGCTTXXXXX")
    result = find_crr4_match(partial)
    assert 20 < result < 28  # Should be positive but less than 28


def test_count_fasta_records(sample_fasta_file):
    """Test counting records in a FASTA file."""
    count = count_fasta_records(sample_fasta_file)
    assert count == 1  # Sample file has 1 record


def test_count_fasta_records_empty(tmp_path):
    """Test counting records in an empty FASTA file."""
    empty_file = tmp_path / "empty.fasta"
    empty_file.touch()
    count = count_fasta_records(empty_file)
    assert count == 0


def test_process_sequence(tmp_path):
    """Test processing a sequence to find CRR spacers."""
    # Create sequence with known CRR1 pattern
    seq = list("N" * 50 + "GTGTTCCCCGCGTGAGCGGGGATAAACCG" + 
              "ACTAGCTAGCTAGCTAGCTAGCT" + 
              "GTGTTCCCCGCGTGAGCGGGGATAAACCG" + "N" * 50)
    
    # Create mock output files
    output_files = {}
    for key in ["all", "crr1", "crr2", "crr4", "error", "csv"]:
        file_path = tmp_path / f"{key}.fasta"
        output_files[key] = open(file_path, "w")
    
    try:
        # Process the sequence
        crr1_count, crr2_count, crr4_count = process_sequence(
            "TestGenome", seq, output_files, "test"
        )
        
        # Check counts
        assert crr1_count == 1
        assert crr2_count == 0
        assert crr4_count == 0
        
        # Check file content
        output_files["all"].close()
        with open(tmp_path / "all.fasta", "r") as f:
            content = f.read()
            assert "TestGenome CRR1Spacer1" in content
            assert "ACTAGCTAGCTAGCTAGCTAGCT" in content
    
    finally:
        # Close all files
        for file in output_files.values():
            if not file.closed:
                file.close()


def test_find_crispr_spacers(sample_fasta_file, results_dir):
    """Test finding CRISPR spacers in a FASTA file."""
    # Run the finder
    result = find_crispr_spacers(sample_fasta_file, "test.", results_dir)
    
    # Check result counts
    assert result["crr1"] > 0
    assert result["crr2"] > 0
    assert result["crr4"] > 0
    assert result["total"] == result["crr1"] + result["crr2"] + result["crr4"]
    
    # Check output files
    crr1_file = results_dir / "CRR1" / "test.CRR1.fasta"
    assert crr1_file.exists()
    with crr1_file.open("r") as f:
        content = f.read()
        assert "TestGenome1 CRR1Spacer1" in content
    
    # Check CSV file
    csv_file = results_dir / "test.Results.csv"
    assert csv_file.exists()
    with csv_file.open("r") as f:
        content = f.read()
        assert "TestGenome1" in content
        assert str(result["crr1"]) in content


def test_find_crispr_spacers_file_not_found(tmp_path):
    """Test error when file not found."""
    nonexistent_file = tmp_path / "nonexistent.fasta"
    with pytest.raises(FileNotFoundError):
        find_crispr_spacers(nonexistent_file, "test")


def test_find_crispr_spacers_empty_file(tmp_path, results_dir):
    """Test with empty FASTA file."""
    empty_file = tmp_path / "empty.fasta"
    empty_file.touch()
    
    result = find_crispr_spacers(empty_file, "test.", results_dir)
    assert result["total"] == 0