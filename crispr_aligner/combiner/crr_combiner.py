"""Utilities for combining aligned CRISPR arrays from different regions."""

from pathlib import Path
from typing import List, Union
from Bio import SeqIO


def combine_aligned_fasta(
    crr1_aligned_path: Union[str, Path],
    crr2_aligned_path: Union[str, Path],
    crr4_aligned_path: Union[str, Path],
    output_path: Union[str, Path]
) -> Path:
    """
    Combine aligned FASTA files from different CRISPR repeat regions.
    
    Args:
        crr1_aligned_path: Path to CRR1 aligned FASTA file
        crr2_aligned_path: Path to CRR2 aligned FASTA file 
        crr4_aligned_path: Path to CRR4 aligned FASTA file
        output_path: Path to write combined output
        
    Returns:
        Path to the output file
    """
    crr1_aligned_path = Path(crr1_aligned_path)
    crr2_aligned_path = Path(crr2_aligned_path)
    crr4_aligned_path = Path(crr4_aligned_path)
    output_path = Path(output_path)
    
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with output_path.open("w") as output_file:
        with crr1_aligned_path.open("r") as crr1_file:
            crr1_records = SeqIO.parse(crr1_file, "fasta")
            
            for record1 in crr1_records:
                output_file.write(f">{record1.id}\n{record1.seq}")
                
                # Find matching record in CRR2
                with crr2_aligned_path.open("r") as crr2_file:
                    crr2_records = SeqIO.parse(crr2_file, "fasta")
                    for record2 in crr2_records:
                        if record1.id == record2.id:
                            output_file.write(str(record2.seq))
                            break
                
                # Find matching record in CRR4
                with crr4_aligned_path.open("r") as crr4_file:
                    crr4_records = SeqIO.parse(crr4_file, "fasta")
                    for record4 in crr4_records:
                        if record1.id == record4.id:
                            output_file.write(str(record4.seq))
                            break
                
                output_file.write("\n")
    
    print(f"Combined FASTA alignments written to {output_path}")
    return output_path


def combine_csv_files(
    crr1_csv_path: Union[str, Path],
    crr2_csv_path: Union[str, Path],
    crr4_csv_path: Union[str, Path],
    output_path: Union[str, Path]
) -> Path:
    """
    Combine CSV files from different CRISPR repeat regions.
    
    Args:
        crr1_csv_path: Path to CRR1 CSV file
        crr2_csv_path: Path to CRR2 CSV file
        crr4_csv_path: Path to CRR4 CSV file
        output_path: Path to write combined output
        
    Returns:
        Path to the output file
    """
    crr1_csv_path = Path(crr1_csv_path)
    crr2_csv_path = Path(crr2_csv_path)
    crr4_csv_path = Path(crr4_csv_path)
    output_path = Path(output_path)
    
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Read all CSV files into memory
    crr1_rows = []
    with crr1_csv_path.open("r") as crr1_file:
        for line in crr1_file:
            line = line.strip()
            crr1_rows.append(line.split(","))
    
    crr2_rows = []
    with crr2_csv_path.open("r") as crr2_file:
        for line in crr2_file:
            line = line.strip()
            crr2_rows.append(line.split(","))
    
    crr4_rows = []
    with crr4_csv_path.open("r") as crr4_file:
        for line in crr4_file:
            line = line.strip()
            crr4_rows.append(line.split(","))
    
    # Write combined CSV
    with output_path.open("w") as out_file:
        for i in range(len(crr1_rows)):
            # Write CRR1 data
            out_file.write(",".join(crr1_rows[i]))
            
            # Add CRR2 data (skip first column)
            if i < len(crr2_rows) and len(crr2_rows[i]) > 1:
                out_file.write(",")
                out_file.write(",".join(crr2_rows[i][1:]))
            
            # Add CRR4 data (skip first column)
            if i < len(crr4_rows) and len(crr4_rows[i]) > 1:
                out_file.write(",")
                out_file.write(",".join(crr4_rows[i][1:]))
            
            out_file.write("\n")
    
    print(f"Combined CSV files written to {output_path}")
    return output_path


def combine_binary_files(
    crr1_binary_path: Union[str, Path],
    crr2_binary_path: Union[str, Path],
    crr4_binary_path: Union[str, Path],
    output_path: Union[str, Path]
) -> Path:
    """
    Combine binary phylip files from different CRISPR repeat regions.
    
    Args:
        crr1_binary_path: Path to CRR1 binary phylip file
        crr2_binary_path: Path to CRR2 binary phylip file
        crr4_binary_path: Path to CRR4 binary phylip file
        output_path: Path to write combined output
        
    Returns:
        Path to the output file
    """
    crr1_binary_path = Path(crr1_binary_path)
    crr2_binary_path = Path(crr2_binary_path)
    crr4_binary_path = Path(crr4_binary_path)
    output_path = Path(output_path)
    
    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    first_line = True
    with output_path.open("w") as out_file:
        with crr1_binary_path.open("r") as crr1_file:
            for line1 in crr1_file:
                line1 = line1.strip()
                
                # Special handling for the first line (header)
                if first_line:
                    first_line = False
                    
                    # Parse header values
                    parts1 = line1.split()
                    taxa_count = parts1[0]
                    
                    # Get character counts from all files
                    with crr2_binary_path.open("r") as crr2_file:
                        parts2 = crr2_file.readline().strip().split()
                        char_count2 = int(parts2[1])
                    
                    with crr4_binary_path.open("r") as crr4_file:
                        parts4 = crr4_file.readline().strip().split()
                        char_count4 = int(parts4[1])
                    
                    # Calculate total character count
                    total_chars = int(parts1[1]) + char_count2 + char_count4
                    
                    # Write combined header
                    out_file.write(f"{taxa_count} {total_chars}\n")
                    continue
                
                # For data lines, append sequences from all files
                out_file.write(line1)
                
                # Get strain ID (first 10 chars)
                strain_id = line1[:min(10, len(line1))]
                
                # Find matching lines in other files
                with crr2_binary_path.open("r") as crr2_file:
                    next(crr2_file)  # Skip header
                    for line2 in crr2_file:
                        line2 = line2.strip()
                        if line2.startswith(strain_id):
                            # Append data part (after ID)
                            out_file.write(line2[min(10, len(line2)):])
                            break
                
                with crr4_binary_path.open("r") as crr4_file:
                    next(crr4_file)  # Skip header
                    for line4 in crr4_file:
                        line4 = line4.strip()
                        if line4.startswith(strain_id):
                            # Append data part (after ID)
                            out_file.write(line4[min(10, len(line4)):])
                            break
                
                out_file.write("\n")
    
    print(f"Combined binary files written to {output_path}")
    return output_path


def combine_alignments(
    crr1_aligned_fasta: Union[str, Path],
    crr2_aligned_fasta: Union[str, Path],
    crr4_aligned_fasta: Union[str, Path],
    all_crr_aligned_fasta: Union[str, Path],
    crr1_aligned_csv: Union[str, Path],
    crr2_aligned_csv: Union[str, Path],
    crr4_aligned_csv: Union[str, Path],
    all_crr_aligned_csv: Union[str, Path],
    crr1_binary_phy: Union[str, Path],
    crr2_binary_phy: Union[str, Path],
    crr4_binary_phy: Union[str, Path],
    all_crr_binary_phy: Union[str, Path]
) -> None:
    """
    Combine all alignment files from different CRISPR repeat regions.
    
    Args:
        crr1_aligned_fasta: Path to CRR1 aligned FASTA file
        crr2_aligned_fasta: Path to CRR2 aligned FASTA file
        crr4_aligned_fasta: Path to CRR4 aligned FASTA file
        all_crr_aligned_fasta: Output path for combined FASTA
        crr1_aligned_csv: Path to CRR1 aligned CSV file
        crr2_aligned_csv: Path to CRR2 aligned CSV file
        crr4_aligned_csv: Path to CRR4 aligned CSV file
        all_crr_aligned_csv: Output path for combined CSV
        crr1_binary_phy: Path to CRR1 binary phylip file
        crr2_binary_phy: Path to CRR2 binary phylip file
        crr4_binary_phy: Path to CRR4 binary phylip file
        all_crr_binary_phy: Output path for combined phylip
    """
    # Combine all file types
    combine_aligned_fasta(crr1_aligned_fasta, crr2_aligned_fasta, crr4_aligned_fasta, all_crr_aligned_fasta)
    combine_csv_files(crr1_aligned_csv, crr2_aligned_csv, crr4_aligned_csv, all_crr_aligned_csv)
    combine_binary_files(crr1_binary_phy, crr2_binary_phy, crr4_binary_phy, all_crr_binary_phy)
    
    print("All alignment files successfully combined")