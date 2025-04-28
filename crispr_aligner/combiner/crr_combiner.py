"""Utilities for combining aligned CRISPR arrays from different regions."""

from pathlib import Path
from typing import Union, List
from Bio import SeqIO


def combine_fasta(
    aligned_fastas: List[Union[str, Path]], output_path: Union[str, Path]
) -> Path:
    """
    Combine multiple aligned FASTA files from different CRISPR repeat regions.

    Args:
        aligned_fastas: List of paths to aligned FASTA files
        output_path: Path to write combined output

    Returns:
        Path to the output file
    """
    # Convert all paths to Path objects
    fastas = [Path(f) for f in aligned_fastas]
    output_path = Path(output_path)

    # Ensure we have at least one file
    if not fastas:
        print("No FASTA files provided for combination")
        return output_path

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get all strain IDs from the first file
    strain_ids = set()
    with fastas[0].open("r") as first_file:
        records = SeqIO.parse(first_file, "fasta")
        for record in records:
            strain_ids.add(record.id)

    # Create a mapping of strain_id -> {file_index: sequence}
    combined_seqs = {strain_id: {} for strain_id in strain_ids}

    # Read all sequences from all files
    for i, fasta_path in enumerate(fastas):
        with fasta_path.open("r") as f:
            records = SeqIO.parse(f, "fasta")
            for record in records:
                if record.id in combined_seqs:
                    combined_seqs[record.id][i] = str(record.seq)

    # Write combined file
    with output_path.open("w") as out_file:
        for strain_id in strain_ids:
            out_file.write(f">{strain_id}\n")
            for i in range(len(fastas)):
                if i in combined_seqs[strain_id]:
                    out_file.write(combined_seqs[strain_id][i])
            out_file.write("\n")

    print(f"Combined {len(fastas)} FASTA alignments written to {output_path}")
    return output_path


def combine_csv(
    csv_files: List[Union[str, Path]], output_path: Union[str, Path]
) -> Path:
    """
    Combine multiple CSV files from different CRISPR repeat regions.

    Args:
        csv_files: List of paths to CSV files
        output_path: Path to write combined output

    Returns:
        Path to the output file
    """
    # Convert all paths to Path objects
    csv_paths = [Path(f) for f in csv_files]
    output_path = Path(output_path)

    # Ensure we have at least one file
    if not csv_paths:
        print("No CSV files provided for combination")
        return output_path

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Read all CSV files into memory
    all_rows = []
    for i, csv_path in enumerate(csv_paths):
        rows = []
        with csv_path.open("r") as f:
            for line in f:
                line = line.strip()
                if line:
                    rows.append(line.split(","))
        all_rows.append(rows)

    # Write combined CSV
    with output_path.open("w") as out_file:
        # Determine max number of rows
        max_rows = max(len(rows) for rows in all_rows) if all_rows else 0

        # Process each row
        for i in range(max_rows):
            row_parts = []

            # For first file, include all columns
            if all_rows and i < len(all_rows[0]):
                row_parts.append(",".join(all_rows[0][i]))

            # For other files, skip first column (strain ID)
            for file_idx in range(1, len(all_rows)):
                if i < len(all_rows[file_idx]) and len(all_rows[file_idx][i]) > 1:
                    row_parts.append(",".join(all_rows[file_idx][i][1:]))

            out_file.write(",".join(row_parts) + "\n")

    print(f"Combined {len(csv_paths)} CSV files written to {output_path}")
    return output_path


def combine_binary(
    binary_files: List[Union[str, Path]], output_path: Union[str, Path]
) -> Path:
    """
    Combine multiple binary phylip files from different CRISPR repeat regions.

    Args:
        binary_files: List of paths to binary phylip files
        output_path: Path to write combined output

    Returns:
        Path to the output file
    """
    # Convert all paths to Path objects
    phy_paths = [Path(f) for f in binary_files]
    output_path = Path(output_path)

    # Ensure we have at least one file
    if not phy_paths:
        print("No binary files provided for combination")
        return output_path

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Read headers to calculate total character count
    taxa_count = 0
    char_counts = []

    for phy_path in phy_paths:
        with phy_path.open("r") as f:
            header = f.readline().strip().split()
            if not taxa_count:
                taxa_count = int(header[0])
            char_counts.append(int(header[1]))

    # Read all data lines
    strain_data = {}

    for file_idx, phy_path in enumerate(phy_paths):
        with phy_path.open("r") as f:
            next(f)  # Skip header
            for line in f:
                line = line.strip()
                if not line:
                    continue

                # Get strain ID (first 10 chars)
                strain_id = line[: min(10, len(line))].strip()
                data = line[min(10, len(line)) :]

                if strain_id not in strain_data:
                    strain_data[strain_id] = {}

                strain_data[strain_id][file_idx] = data

    # Write combined file
    with output_path.open("w") as out_file:
        # Write header
        total_chars = sum(char_counts)
        out_file.write(f"{taxa_count} {total_chars}\n")

        # Write data for each strain
        for strain_id, data_dict in strain_data.items():
            out_file.write(f"{strain_id.ljust(10)}")

            for file_idx in range(len(phy_paths)):
                if file_idx in data_dict:
                    out_file.write(data_dict[file_idx])

            out_file.write("\n")

    print(f"Combined {len(phy_paths)} binary files written to {output_path}")
    return output_path


def combine_alignments(
    aligned_fastas: List[Union[str, Path]],
    output_fasta: Union[str, Path],
    aligned_csvs: List[Union[str, Path]],
    output_csv: Union[str, Path],
    binary_phys: List[Union[str, Path]],
    output_phy: Union[str, Path],
) -> None:
    """
    Combine all alignment files from different CRISPR repeat regions.

    Args:
        aligned_fastas: List of paths to aligned FASTA files
        output_fasta: Output path for combined FASTA
        aligned_csvs: List of paths to aligned CSV files
        output_csv: Output path for combined CSV
        binary_phys: List of paths to binary phylip files
        output_phy: Output path for combined phylip
    """
    combine_fasta(aligned_fastas, output_fasta)
    combine_csv(aligned_csvs, output_csv)
    combine_binary(binary_phys, output_phy)

    print("All alignment files successfully combined")
