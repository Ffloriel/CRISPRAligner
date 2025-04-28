# Create a minimal test sample with known CRISPR patterns
def create_test_sample():
    """Create a sample file for testing with known CRISPR patterns."""
    from pathlib import Path
    
    test_dir = Path("tests/test_data")
    test_dir.mkdir(parents=True, exist_ok=True)
    
    # Create a FASTA with known CRISPR repeats and spacers
    with open(test_dir / "sample.fasta", "w") as f:
        f.write(">TestGenome1\n")
        # Add a CRR1 array (2 repeats with a spacer)
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")  # Prefix
        f.write("GTGTTCCCCGCGTGAGCGGGGATAAACCG")    # CRR1 repeat
        f.write("ACTAGCTAGCTAGCTAGCTAGCTAGCTAGCT")  # Spacer (30bp)
        f.write("GTGTTCCCCGCGTGAGCGGGGATAAACCG")    # CRR1 repeat
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")  # Suffix
        
        # Add a CRR2 array
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")  # Prefix
        f.write("GTGTTCCCCGCGTATGCGGGGATAAACCG")    # CRR2 repeat
        f.write("GCATGCATGCATGCATGCATGCATGCATGCA")  # Spacer (30bp)
        f.write("GTGTTCCCCGCGTATGCGGGGATAAACCG")    # CRR2 repeat
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")  # Suffix
        
        # Add a CRR4 array
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")  # Prefix
        f.write("GTTCACTGCCGTACAGGCAGCTTAGAAA")     # CRR4 repeat
        f.write("TGACTGACTGACTGACTGACTGACTGACTGA")  # Spacer (30bp)
        f.write("GTTCACTGCCGTACAGGCAGCTTAGAAA")     # CRR4 repeat
        f.write("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")  # Suffix
        
    print(f"Created test sample at {test_dir / 'sample.fasta'}")

# Run this once to create your test data
# create_test_sample()