"""CRISPR repeat region (CRR) finder module for identifying CRISPR spacers in bacterial genomes."""

from pathlib import Path
from typing import Tuple, Dict, Union, List, TextIO, Optional, Any
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
    max_frame_size: int = 20000,
) -> Dict[str, Dict[str, int]]:
    """
    Process a sequence to find CRR spacers, grouping them into frames.

    Args:
        record_id: Identifier for the sequence record
        sequence: Sequence as a list of characters
        output_files: Dictionary of open file handles for writing results
        prefix: Prefix for output file names
        consensus_sequences: Dictionary of consensus sequences for each CRR type
        short_checks: Dictionary of short check sequences for pattern matching
        min_scores: Dictionary of minimum scores for each CRR type
        available_crrs: List of available CRR types to search for
        max_frame_size: Maximum size in bp for a CRISPR array frame

    Returns:
        Dictionary of counts by CRR type and frame: {crr_type: {frame_letter: count}}
    """
    if available_crrs is None:
        available_crrs = ["crr1", "crr2", "crr4"]

    # Store found repeats with their positions for each CRR type
    found_repeats = {
        "crr1": [],  # List of (position, spacer) tuples
        "crr2": [],
        "crr4": [],
    }

    # Step 1: Find all repeats and their spacers, but don't output yet
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
                            if "error" in output_files:
                                output_files["error"].write(
                                    f">{record_id} CRR1Spacer {m+30}:{m+59+n}\n{spacer}\n"
                                )
                            break
                        else:
                            # Store for later processing
                            found_repeats["crr1"].append((m, spacer))
                            break

                    # CRR2 match (only if crr2 is available)
                    elif (
                        "crr2" in available_crrs
                        and crr12_match_start <= -min_scores["crr2"]
                        and crr12_match_end <= -min_scores["crr2"]
                    ):
                        spacer = "".join(sequence[m + 29 : m + 59 + n])

                        if len(spacer) > 35:
                            if "error" in output_files:
                                output_files["error"].write(
                                    f">{record_id} CRR2Spacer {m+30}:{m+59+n}\n{spacer}\n"
                                )
                            break
                        else:
                            # Store for later processing
                            found_repeats["crr2"].append((m, spacer))
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
                            if "error" in output_files:
                                output_files["error"].write(
                                    f">{record_id} CRR4Spacer {m+29}:{m+58+n}\n{spacer}\n"
                                )
                            break
                        else:
                            # Store for later processing
                            found_repeats["crr4"].append((m, spacer))
                            break

    # Step 2: Group repeats into frames based on max_frame_size
    frame_counts = {
        crr: {} for crr in available_crrs
    }  # {crr_type: {frame_letter: count}}

    # Define letter suffixes for frames
    letters = "abcdefghijklmnopqrstuvwxyz"

    # Process each CRR type
    for crr in available_crrs:
        if not found_repeats[crr]:
            continue

        # Sort repeats by position
        found_repeats[crr].sort(key=lambda x: x[0])

        current_frame = 0
        frame_letter = letters[current_frame]
        frame_start = found_repeats[crr][0][0]

        frame_counts[crr][frame_letter] = 0

        # Group repeats into frames
        for pos, spacer in found_repeats[crr]:
            # If this repeat is too far from frame start, start a new frame
            if pos - frame_start > max_frame_size:
                current_frame += 1
                if current_frame < len(letters):
                    frame_letter = letters[current_frame]
                else:
                    frame_letter = "z"  # Use 'z' for any frames beyond 26
                frame_start = pos
                frame_counts[crr][frame_letter] = 0

            # Increment count for this frame
            frame_counts[crr][frame_letter] += 1

            # Write to appropriate output files
            crr_upper = crr.upper()
            frame_name = f"{crr_upper}{frame_letter.upper()}"

            # Write to all file
            if "all" in output_files:
                output_files["all"].write(
                    f">{record_id} {frame_name}Spacer{frame_counts[crr][frame_letter]} {pos+30}:{pos+30+len(spacer)}\n{spacer}\n"
                )

            # Write to specific frame file
            frame_key = f"{crr}_{frame_letter}"
            if frame_key in output_files:
                output_files[frame_key].write(
                    f">{record_id} {frame_name}Spacer{frame_counts[crr][frame_letter]} {pos+30}:{pos+30+len(spacer)}\n{spacer}\n"
                )

    # Add newlines after each record's spacers
    for file_key in ["all", "error"] + [
        f"{crr}_{letter}"
        for crr in available_crrs
        for letter in frame_counts[crr].keys()
        if f"{crr}_{letter}" in output_files
    ]:
        if file_key in output_files:
            output_files[file_key].write("\n")

    return frame_counts


def find_crispr_spacers(
    input_file: Union[str, Path],
    prefix: str,
    results_dir: Union[str, Path],
    consensus_sequences: Dict[str, str],
    min_scores: Dict[str, int],
    short_checks: Dict[str, str],
    available_crrs: List[str],
    max_frame_size: int = 20000,
) -> Dict[str, Any]:
    """
    Find CRISPR spacers in a FASTA file, organizing them into spatial frames.

    Args:
        input_file: Path to input FASTA file
        prefix: Prefix for output file names
        results_dir: Directory to store results
        consensus_sequences: Dictionary with custom consensus sequences
        min_scores: Dictionary with minimum match scores for each CRR type
        short_checks: Dictionary with short sequence patterns for initial checks
        available_crrs: List of available CRR types to search for
        max_frame_size: Maximum size in bp for a CRISPR array frame

    Returns:
        Dictionary with counts for each CRR frame and total
    """
    input_file = Path(results_dir, input_file)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    # Create output directories if they don't exist
    # Always create Results and AllCRR directories
    for dir_name in ["", "AllCRR"]:
        (results_dir / dir_name).mkdir(parents=True, exist_ok=True)

    # First pass: Process sequences to identify all frames
    print("First pass: Identifying CRISPR array frames...")
    all_frames = {crr: set() for crr in available_crrs}

    with input_file.open("r") as seq_file:
        for record in SeqIO.parse(seq_file, "fasta"):
            seq = record.seq

            # Determine if sequence needs to be reverse complemented
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

            # Process sequence to identify frames (dummy output_files dict)
            frame_counts = process_sequence(
                record.id,
                seq,
                {},  # Empty output_files for first pass
                prefix,
                consensus_sequences,
                short_checks,
                min_scores,
                available_crrs,
                max_frame_size,
            )

            # Collect all frame letters
            for crr, frames in frame_counts.items():
                all_frames[crr].update(frames.keys())

    # Create directories for all discovered frames
    print("Creating directories for discovered CRISPR frames...")
    for crr in available_crrs:
        for frame_letter in sorted(all_frames[crr]):
            crr_upper = crr.upper()
            frame_dir = f"{crr_upper}{frame_letter.upper()}"
            (results_dir / frame_dir).mkdir(parents=True, exist_ok=True)
            print(f"Created directory for {frame_dir}")

    # Prepare output files dictionary
    output_files = {}

    # Open output files
    output_files["all"] = open(results_dir / "AllCRR" / f"{prefix}AllCRR.fasta", "w")
    output_files["error"] = open(results_dir / f"{prefix}Error.fasta", "w")
    output_files["csv"] = open(results_dir / f"{prefix}Results.csv", "w")

    # Open files for each frame
    for crr in available_crrs:
        for frame_letter in sorted(all_frames[crr]):
            crr_upper = crr.upper()
            frame_dir = f"{crr_upper}{frame_letter.upper()}"
            frame_key = f"{crr}_{frame_letter}"
            output_files[frame_key] = open(
                results_dir / frame_dir / f"{prefix}{frame_dir}.fasta", "w"
            )

    # Write CSV header with all frames
    csv_header = "Strain ID"
    for crr in available_crrs:
        for frame_letter in sorted(all_frames[crr]):
            crr_upper = crr.upper()
            frame_name = f"{crr_upper}{frame_letter.upper()}"
            csv_header += f", {frame_name} Spacers"
    csv_header += ", Total Spacers\n"
    output_files["csv"].write(csv_header)

    output_files["error"].write(
        "Potential assembly errors found in the following spacers:\n\n"
    )

    # Second pass: Process records and write output
    print("Second pass: Processing sequences and writing output...")
    total_records = count_fasta_records(input_file)
    print(f"Processing {total_records} records from {input_file}")

    # Initialize counters for all frames
    all_counts = {
        f"{crr}_{letter}": 0 for crr in available_crrs for letter in all_frames[crr]
    }
    all_counts["total"] = 0

    # Process each record
    record_count = 0
    with input_file.open("r") as seq_file:
        for record in SeqIO.parse(seq_file, "fasta"):
            record_count += 1
            print(f"Processing {record.id} ({record_count}/{total_records})")

            # Determine if sequence needs to be reverse complemented
            seq = record.seq

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

            # Process sequence and get frame counts
            frame_counts = process_sequence(
                record.id,
                seq,
                output_files,
                prefix,
                consensus_sequences,
                short_checks,
                min_scores,
                available_crrs,
                max_frame_size,
            )

            # Calculate total count
            record_total = sum(
                sum(counts.values())
                for crr, counts in frame_counts.items()
                if crr in available_crrs
            )

            # Update counts
            for crr in available_crrs:
                for frame_letter, count in frame_counts[crr].items():
                    frame_key = f"{crr}_{frame_letter}"
                    all_counts[frame_key] += count

            all_counts["total"] += record_total

            # Write results to CSV
            csv_line = record.id
            for crr in available_crrs:
                for frame_letter in sorted(all_frames[crr]):
                    count = frame_counts[crr].get(frame_letter, 0)
                    csv_line += f", {count}"
            csv_line += f", {record_total}\n"
            output_files["csv"].write(csv_line)

            # Print progress summary
            status_msg = f"Spacers found: "
            for crr in available_crrs:
                for frame_letter in sorted(all_frames[crr]):
                    crr_upper = crr.upper()
                    frame_name = f"{crr_upper}{frame_letter.upper()}"
                    count = frame_counts[crr].get(frame_letter, 0)
                    if count > 0:
                        status_msg += f"{frame_name} {count}, "
            status_msg += f"total {record_total}"
            print(status_msg)

    print(f"Processing complete. Results written to {results_dir} directory")

    # Print summary of all frames
    summary_msg = f"Total spacers found: "
    for crr in available_crrs:
        for frame_letter in sorted(all_frames[crr]):
            crr_upper = crr.upper()
            frame_name = f"{crr_upper}{frame_letter.upper()}"
            frame_key = f"{crr}_{frame_letter}"
            count = all_counts[frame_key]
            if count > 0:
                summary_msg += f"{frame_name} {count}, "
    summary_msg += f"total {all_counts['total']}"
    print(summary_msg)

    # Close all open files
    try:
        for file_handle in output_files.values():
            if not file_handle.closed:
                file_handle.close()
    except Exception as e:
        print(f"Error closing files: {e}")

    return all_counts
