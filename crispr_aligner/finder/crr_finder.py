"""CRISPR repeat region (CRR) finder module for identifying CRISPR spacers in bacterial genomes."""

from pathlib import Path
from typing import Tuple, Dict, Union, List, TextIO, Optional
from Bio import SeqIO
from Bio.Seq import Seq
import re


def find_crr12_match(query: List[str]) -> int:
    """
    Find and score matches to CRR1 and CRR2 consensus sequences.

    Args:
        query: List of characters representing sequence to check

    Returns:
        Positive score for CRR1 match, negative score for CRR2 match
    """
    crr1 = list("GTGTTCCCCGCGTGAGCGGGGATAAACCG")
    crr2 = list("GTGTTCCCCGCGTATGCGGGGATAAACCG")

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


def find_crr4_match(query: List[str]) -> int:
    """
    Find and score matches to CRR4 consensus sequence.

    Args:
        query: List of characters representing sequence to check

    Returns:
        Score representing number of matches to CRR4 consensus
    """
    crr4 = list("GTTCACTGCCGTACAGGCAGCTTAGAAA")
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
    record_id: str, sequence: List[str], output_files: Dict[str, TextIO], prefix: str
) -> Tuple[int, int, int]:
    """
    Process a sequence to find CRR spacers.

    Args:
        record_id: Identifier for the sequence record
        sequence: Sequence as a list of characters
        output_files: Dictionary of open file handles for writing results
        prefix: Prefix for output file names

    Returns:
        Tuple of (CRR1_count, CRR2_count, CRR4_count)
    """
    crr1_count = 0
    crr2_count = 0
    crr4_count = 0

    check1 = list("GTGTTC")
    check2 = list("GATAAACC")
    check3 = list("GTTCAC")
    check4 = list("GTACGGG")

    for m in range(len(sequence) - 29):
        # Check for CRR1 and CRR2 patterns
        if sequence[m : m + 6] == check1 or sequence[m + 20 : m + 28] == check2:

            for n in range(100):
                crr12_match_start = find_crr12_match(sequence[m : m + 29])
                crr12_match_end = find_crr12_match(sequence[m + 59 + n : m + 88 + n])

                # CRR1 match
                if crr12_match_start >= 22 and crr12_match_end >= 22:

                    spacer = "".join(sequence[m + 29 : m + 59 + n])

                    if len(spacer) > 35:
                        output_files["error"].write(
                            f">{record_id} CRR1Spacer{crr1_count} {m+30}:{m+59+n}\n{spacer}\n"
                        )
                        break
                    else:
                        crr1_count += 1
                        for file_key in ["all", "crr1"]:
                            output_files[file_key].write(
                                f">{record_id} CRR1Spacer{crr1_count} {m+30}:{m+59+n}\n{spacer}\n"
                            )
                        break

                # CRR2 match
                elif crr12_match_start <= -22 and crr12_match_end <= -22:

                    spacer = "".join(sequence[m + 29 : m + 59 + n])

                    if len(spacer) > 35:
                        output_files["error"].write(
                            f">{record_id} CRR2Spacer{crr2_count} {m+30}:{m+59+n}\n{spacer}\n"
                        )
                        break
                    else:
                        crr2_count += 1
                        for file_key in ["all", "crr2"]:
                            output_files[file_key].write(
                                f">{record_id} CRR2Spacer{crr2_count} {m+30}:{m+59+n}\n{spacer}\n"
                            )
                        break

        # Check for CRR4 patterns
        elif sequence[m : m + 6] == check3 or sequence[m + 10 : m + 17] == check4:

            for n in range(100):
                crr4_match_start = find_crr4_match(sequence[m : m + 28])
                crr4_match_end = find_crr4_match(sequence[m + 58 + n : m + 86 + n])

                if crr4_match_start >= 21 and crr4_match_end >= 21:

                    spacer = "".join(sequence[m + 28 : m + 58 + n])

                    if len(spacer) > 35:
                        output_files["error"].write(
                            f">{record_id} CRR4Spacer{crr4_count} {m+29}:{m+58+n}\n{spacer}\n"
                        )
                        break
                    else:
                        crr4_count += 1
                        for file_key in ["all", "crr4"]:
                            output_files[file_key].write(
                                f">{record_id} CRR4Spacer{crr4_count} {m+29}:{m+58+n}\n{spacer}\n"
                            )
                        break

    # Add newlines after each record's spacers
    for file_key in ["all", "crr1", "crr2", "crr4"]:
        output_files[file_key].write("\n")

    return crr1_count, crr2_count, crr4_count


def find_crispr_spacers(
    input_file: Union[str, Path],
    prefix: str,
    results_dir: Union[str, Path],
    consensus_sequences: Optional[Dict[str, str]] = None,
    min_scores: Optional[Dict[str, int]] = None,
) -> Dict[str, int]:
    """
    Find CRISPR spacers in a FASTA file.

    Args:
        input_file: Path to input FASTA file
        prefix: Prefix for output file names
        results_dir: Directory to store results
        consensus_sequences: Dictionary with custom consensus sequences
        min_scores: Dictionary with minimum match scores for each CRR type

    Returns:
        Dictionary with counts for each CRR type and total
    """
    # Use default values if not provided
    if consensus_sequences is None:
        consensus_sequences = {
            "crr1": "GTGTTCCCCGCGTGAGCGGGGATAAACCG",
            "crr2": "GTGTTCCCCGCGTATGCGGGGATAAACCG",
            "crr4": "GTTCACTGCCGTACAGGCAGCTTAGAAA",
        }

    if min_scores is None:
        min_scores = {"crr1": 22, "crr2": 22, "crr4": 21}

    input_file = Path(results_dir, input_file)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    if results_dir is None:
        results_dir = Path("Results")
    else:
        results_dir = Path(results_dir)

    # Create output directories if they don't exist
    for dir_name in ["", "AllCRR", "CRR1", "CRR2", "CRR4"]:
        (results_dir / dir_name).mkdir(parents=True, exist_ok=True)

    # Open all output files
    with open(results_dir / "AllCRR" / f"{prefix}AllCRR.fasta", "w") as all_file, open(
        results_dir / "CRR1" / f"{prefix}CRR1.fasta", "w"
    ) as crr1_file, open(
        results_dir / "CRR2" / f"{prefix}CRR2.fasta", "w"
    ) as crr2_file, open(
        results_dir / "CRR4" / f"{prefix}CRR4.fasta", "w"
    ) as crr4_file, open(
        results_dir / f"{prefix}Results.csv", "w"
    ) as csv_file, open(
        results_dir / f"{prefix}Error.fasta", "w"
    ) as error_file:

        # Prepare output files dictionary
        output_files = {
            "all": all_file,
            "crr1": crr1_file,
            "crr2": crr2_file,
            "crr4": crr4_file,
            "error": error_file,
            "csv": csv_file,
        }

        # Write headers
        csv_file.write(
            "Strain ID, CRR1 Spacers, CRR2 Spacers, CRR4 Spacers, Total Spacers\n"
        )
        error_file.write(
            "Potential assembly errors found in the following spacers:\n\n"
        )

        # Count total records for progress reporting
        total_records = count_fasta_records(input_file)
        print(f"Processing {total_records} records from {input_file}")

        # Process records
        record_count = 0
        all_counts = {"crr1": 0, "crr2": 0, "crr4": 0, "total": 0}

        with input_file.open("r") as seq_file:
            for record in SeqIO.parse(seq_file, "fasta"):
                record_count += 1
                print(f"Processing {record.id} ({record_count}/{total_records})")

                # Determine if sequence needs to be reverse complemented
                seq = record.seq
                if "CGGTTTATCCCCGCTCACGCGGGGAACAC" in seq:
                    seq = list(seq.reverse_complement())
                else:
                    seq = list(seq)

                # Find spacers
                crr1_count, crr2_count, crr4_count = process_sequence(
                    record.id, seq, output_files, prefix
                )

                total_count = crr1_count + crr2_count + crr4_count

                # Update counts
                all_counts["crr1"] += crr1_count
                all_counts["crr2"] += crr2_count
                all_counts["crr4"] += crr4_count
                all_counts["total"] += total_count

                # Write results to CSV
                csv_file.write(
                    f"{record.id}, {crr1_count}, {crr2_count}, {crr4_count}, {total_count}\n"
                )

                print(
                    f"Spacers found: CRR1 {crr1_count}, CRR2 {crr2_count}, "
                    f"CRR4 {crr4_count}, total {total_count}"
                )

    print(f"Processing complete. Results written to {results_dir} directory")
    print(
        f"Total spacers found: CRR1 {all_counts['crr1']}, CRR2 {all_counts['crr2']}, "
        f"CRR4 {all_counts['crr4']}, total {all_counts['total']}"
    )

    return all_counts
