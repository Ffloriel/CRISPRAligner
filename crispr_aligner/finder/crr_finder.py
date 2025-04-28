"""CRISPR repeat region (CRR) finder module for identifying CRISPR spacers in bacterial genomes."""

from pathlib import Path
from typing import Tuple, Dict, Union, List, TextIO, Optional
from Bio import SeqIO


def find_crr12_match(query: List[str], crr1_seq: str, crr2_seq: str) -> int:
    """
    Find and score matches to CRR1 and CRR2 consensus sequences.

    Args:
        query: List of characters representing sequence to check
        crr1_seq: CRR1 consensus sequence
        crr2_seq: CRR2 consensus sequence

    Returns:
        Positive score for CRR1 match, negative score for CRR2 match
    """
    crr1 = list(crr1_seq)
    crr2 = list(crr2_seq)

    crr1_match = 0
    crr2_match = 0

    for i in range(len(crr1)):
        if i < len(query):
            if query[i] == crr1[i]:
                crr1_match += 1
            if query[i] == crr2[i]:
                crr2_match += 1

    if crr1_match >= crr2_match:
        return crr1_match
    else:
        return -crr2_match


def find_crr4_match(query: List[str], crr4_seq: str) -> int:
    """
    Find and score matches to CRR4 consensus sequence.

    Args:
        query: List of characters representing sequence to check
        crr4_seq: CRR4 consensus sequence

    Returns:
        Score representing number of matches to CRR4 consensus
    """
    crr4 = list(crr4_seq)
    crr4_match = 0

    for i in range(len(crr4)):
        if i < len(query):
            if query[i] == crr4[i]:
                crr4_match += 1

    return crr4_match


def count_fasta_records(file_path: Union[str, Path]) -> int:
    """
    Count the number of records in a FASTA file.

    Args:
        file_path: Path to FASTA file

    Returns:
        Number of records in the file
    """
    file_path = Path(file_path)

    count = 0
    with file_path.open("r") as seq_file:
        for record in SeqIO.parse(seq_file, "fasta"):
            count += 1

    return count


def process_sequence(
    record_id: str,
    sequence: List[str],
    output_files: Dict[str, TextIO],
    prefix: str,
    consensus_sequences: Dict[str, str],
    short_checks: Dict[str, str],
    min_scores: Dict[str, int],
    available_crrs: List[str] = None,
) -> Tuple[int, int, int]:
    """
    Process a sequence to find CRR spacers.

    Args:
        record_id: Identifier for the sequence record
        sequence: Sequence as a list of characters
        output_files: Dictionary of open file handles for writing results
        prefix: Prefix for output file names
        consensus_sequences: Dictionary of consensus sequences for each CRR type
        short_checks: Dictionary of short check sequences for pattern matching
        min_scores: Dictionary of minimum scores for each CRR type
        available_crrs: List of available CRR types to search for

    Returns:
        Tuple of (CRR1_count, CRR2_count, CRR4_count)
    """
    if available_crrs is None:
        available_crrs = ["crr1", "crr2", "crr4"]

    crr1_count = 0
    crr2_count = 0
    crr4_count = 0

    # Check for CRR1 and CRR2 patterns if they're available
    if "crr1" in available_crrs or "crr2" in available_crrs:
        check1 = list(short_checks["crr1_crr2_start"])
        check2 = list(short_checks["crr1_crr2_end"])

        for m in range(len(sequence) - 29):
            if (
                sequence[m : m + len(check1)] == check1
                or sequence[m + 20 : m + 20 + len(check2)] == check2
            ):
                for n in range(100):
                    crr12_match_start = find_crr12_match(
                        sequence[m : m + 29],
                        consensus_sequences["crr1"],
                        consensus_sequences["crr2"],
                    )
                    crr12_match_end = find_crr12_match(
                        sequence[m + 59 + n : m + 88 + n],
                        consensus_sequences["crr1"],
                        consensus_sequences["crr2"],
                    )

                    # CRR1 match (only if crr1 is available)
                    if (
                        "crr1" in available_crrs
                        and crr12_match_start >= min_scores["crr1"]
                        and crr12_match_end >= min_scores["crr1"]
                    ):
                        spacer = "".join(sequence[m + 29 : m + 59 + n])

                        if len(spacer) > 35:
                            output_files["error"].write(
                                f">{record_id} CRR1Spacer{crr1_count} {m+30}:{m+59+n}\n{spacer}\n"
                            )
                            break
                        else:
                            crr1_count += 1
                            for file_key in ["all", "crr1"]:
                                if file_key in output_files:
                                    output_files[file_key].write(
                                        f">{record_id} CRR1Spacer{crr1_count} {m+30}:{m+59+n}\n{spacer}\n"
                                    )
                            break

                    # CRR2 match (only if crr2 is available)
                    elif (
                        "crr2" in available_crrs
                        and crr12_match_start <= -min_scores["crr2"]
                        and crr12_match_end <= -min_scores["crr2"]
                    ):
                        spacer = "".join(sequence[m + 29 : m + 59 + n])

                        if len(spacer) > 35:
                            output_files["error"].write(
                                f">{record_id} CRR2Spacer{crr2_count} {m+30}:{m+59+n}\n{spacer}\n"
                            )
                            break
                        else:
                            crr2_count += 1
                            for file_key in ["all", "crr2"]:
                                if file_key in output_files:
                                    output_files[file_key].write(
                                        f">{record_id} CRR2Spacer{crr2_count} {m+30}:{m+59+n}\n{spacer}\n"
                                    )
                            break

    # Check for CRR4 patterns if it's available
    if "crr4" in available_crrs and "crr4" in consensus_sequences:
        check3 = list(short_checks["crr4_start"])
        check4 = list(short_checks["crr4_middle"])

        for m in range(len(sequence) - 29):
            if (
                sequence[m : m + len(check3)] == check3
                or sequence[m + 10 : m + 10 + len(check4)] == check4
            ):
                for n in range(100):
                    crr4_match_start = find_crr4_match(
                        sequence[m : m + 28], consensus_sequences["crr4"]
                    )
                    crr4_match_end = find_crr4_match(
                        sequence[m + 58 + n : m + 86 + n], consensus_sequences["crr4"]
                    )

                    if (
                        crr4_match_start >= min_scores["crr4"]
                        and crr4_match_end >= min_scores["crr4"]
                    ):
                        spacer = "".join(sequence[m + 28 : m + 58 + n])

                        if len(spacer) > 35:
                            output_files["error"].write(
                                f">{record_id} CRR4Spacer{crr4_count} {m+29}:{m+58+n}\n{spacer}\n"
                            )
                            break
                        else:
                            crr4_count += 1
                            for file_key in ["all", "crr4"]:
                                if file_key in output_files:
                                    output_files[file_key].write(
                                        f">{record_id} CRR4Spacer{crr4_count} {m+29}:{m+58+n}\n{spacer}\n"
                                    )
                            break

    # Add newlines after each record's spacers
    for file_key in ["all"] + [
        crr for crr in ["crr1", "crr2", "crr4"] if crr in available_crrs
    ]:
        if file_key in output_files:
            output_files[file_key].write("\n")

    return crr1_count, crr2_count, crr4_count


def find_crispr_spacers(
    input_file: Union[str, Path],
    prefix: str,
    results_dir: Union[str, Path],
    consensus_sequences: Dict[str, str],
    min_scores: Dict[str, int],
    short_checks: Dict[str, str],
    available_crrs: List[str],
) -> Dict[str, int]:
    """
    Find CRISPR spacers in a FASTA file.

    Args:
        input_file: Path to input FASTA file
        prefix: Prefix for output file names
        results_dir: Directory to store results
        consensus_sequences: Dictionary with custom consensus sequences
        min_scores: Dictionary with minimum match scores for each CRR type
        short_checks: Dictionary with short sequence patterns for initial checks
        available_crrs: List of available CRR types to search for

    Returns:
        Dictionary with counts for each CRR type and total
    """

    input_file = Path(results_dir, input_file)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Create output directories if they don't exist
    # Always create Results and AllCRR directories
    for dir_name in ["", "AllCRR"]:
        (results_dir / dir_name).mkdir(parents=True, exist_ok=True)

    # Create directories only for available CRRs
    for crr in available_crrs:
        crr_upper = crr.upper()
        (results_dir / crr_upper).mkdir(parents=True, exist_ok=True)

    # Prepare output files dictionary
    output_files = {}

    # Open output files
    output_files["all"] = open(results_dir / "AllCRR" / f"{prefix}AllCRR.fasta", "w")

    # Open files only for available CRRs
    for crr in available_crrs:
        crr_upper = crr.upper()
        output_files[crr] = open(
            results_dir / crr_upper / f"{prefix}{crr_upper}.fasta", "w"
        )

    # Always open error and csv files
    output_files["error"] = open(results_dir / f"{prefix}Error.fasta", "w")
    output_files["csv"] = open(results_dir / f"{prefix}Results.csv", "w")

    try:
        # Write headers
        csv_header = "Strain ID"
        for crr in available_crrs:
            csv_header += f", {crr.upper()} Spacers"
        csv_header += ", Total Spacers\n"
        output_files["csv"].write(csv_header)

        output_files["error"].write(
            "Potential assembly errors found in the following spacers:\n\n"
        )

        # Count total records for progress reporting
        total_records = count_fasta_records(input_file)
        print(f"Processing {total_records} records from {input_file}")

        # Process records
        record_count = 0
        all_counts = {crr: 0 for crr in available_crrs}
        all_counts["total"] = 0

        with input_file.open("r") as seq_file:
            for record in SeqIO.parse(seq_file, "fasta"):
                record_count += 1
                print(f"Processing {record.id} ({record_count}/{total_records})")

                # Determine if sequence needs to be reverse complemented
                seq = record.seq

                # Only attempt reverse complement check if crr1 is available
                if "crr1" in available_crrs:
                    reverse_complement_pattern = (
                        consensus_sequences["crr1"][::-1]
                        .replace("G", "C")
                        .replace("C", "G")
                        .replace("A", "T")
                        .replace("T", "A")
                    )

                    if reverse_complement_pattern in seq:
                        seq = list(seq.reverse_complement())
                    else:
                        seq = list(seq)
                else:
                    seq = list(seq)

                # Find spacers
                crr1_count, crr2_count, crr4_count = process_sequence(
                    record.id,
                    seq,
                    output_files,
                    prefix,
                    consensus_sequences,
                    short_checks,
                    min_scores,
                    available_crrs,
                )

                # Calculate total count based on available CRRs
                crr_counts = {
                    "crr1": crr1_count,
                    "crr2": crr2_count,
                    "crr4": crr4_count,
                }
                total_count = sum(crr_counts[crr] for crr in available_crrs)

                # Update counts
                for crr in available_crrs:
                    all_counts[crr] += crr_counts[crr]
                all_counts["total"] += total_count

                # Write results to CSV
                csv_line = record.id
                for crr in available_crrs:
                    csv_line += f", {crr_counts[crr]}"
                csv_line += f", {total_count}\n"
                output_files["csv"].write(csv_line)

                # Print progress with only available CRRs
                status_msg = f"Spacers found: "
                for crr in available_crrs:
                    crr_upper = crr.upper()
                    status_msg += f"{crr_upper} {crr_counts[crr]}, "
                status_msg += f"total {total_count}"
                print(status_msg)

        print(f"Processing complete. Results written to {results_dir} directory")

        # Print summary with only available CRRs
        summary_msg = f"Total spacers found: "
        for crr in available_crrs:
            crr_upper = crr.upper()
            summary_msg += f"{crr_upper} {all_counts[crr]}, "
        summary_msg += f"total {all_counts['total']}"
        print(summary_msg)

    finally:
        # Close all open files
        for file_handle in output_files.values():
            if not file_handle.closed:
                file_handle.close()

    return all_counts
