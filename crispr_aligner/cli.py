import argparse
import os
import shutil
from pathlib import Path

from crispr_aligner.sequences.fasta_utils import combine_fasta_files
from crispr_aligner.finder.crr_finder import find_crispr_spacers
from crispr_aligner.aligner.crr_aligner import align_spacers
from crispr_aligner.combiner.crr_combiner import (
    combine_fasta,
    combine_csv,
    combine_binary,
)
from crispr_aligner.utils.config import load_bacteria_config, get_available_bacteria


def create_folder_structure(base_dir, available_crrs):
    """Create the required folder structure for results."""
    results_dir = base_dir / "Results"
    if results_dir.exists():
        shutil.rmtree(results_dir)

    # Always create Results and AllCRR directory
    (base_dir / "Results").mkdir(parents=True, exist_ok=True)
    (base_dir / "Results/AllCRR").mkdir(parents=True, exist_ok=True)

    # Create directories only for available CRRs
    for crr in available_crrs:
        crr_upper = crr.upper()
        (base_dir / f"Results/{crr_upper}").mkdir(parents=True, exist_ok=True)

    return results_dir


def main():
    parser = argparse.ArgumentParser(
        description="CRISPRAligner - Extract and align CRISPR spacers"
    )
    parser.add_argument("prefix", help="Prefix for output files")
    parser.add_argument(
        "--input_dir", default=".", help="Directory containing input FASTA files"
    )

    # Add these arguments:
    parser.add_argument(
        "--bacteria",
        default="erwinia",
        choices=get_available_bacteria(),
        help="Bacteria species to analyze (default: erwinia)",
    )
    parser.add_argument(
        "--custom_config", help="Path to custom consensus sequence config file"
    )
    parser.add_argument(
        "--min_crr1_score",
        type=int,
        default=22,
        help="Minimum score for CRR1 match (default: 22)",
    )
    parser.add_argument(
        "--min_crr2_score",
        type=int,
        default=22,
        help="Minimum score for CRR2 match (default: 22)",
    )
    parser.add_argument(
        "--min_crr4_score",
        type=int,
        default=21,
        help="Minimum score for CRR4 match (default: 21)",
    )

    args = parser.parse_args()

    # Load appropriate configuration
    if args.custom_config:
        config = load_bacteria_config(args.custom_config)
    else:
        config = load_bacteria_config(args.bacteria)

    base_dir = Path(args.input_dir)
    available_crrs = config.get("available_crrs", ["crr1", "crr2", "crr4"])
    results_dir = create_folder_structure(base_dir, available_crrs)

    print("Combining FASTA Files...")
    combined_fasta = f"{args.prefix}.fasta"
    combine_fasta_files(base_dir, combined_fasta, args.prefix)

    print("Finding CRISPR Spacers...")
    find_crispr_spacers(
        combined_fasta,
        f"{args.prefix}.",
        results_dir,
        consensus_sequences=config["consensus_sequences"],
        min_scores={
            "crr1": args.min_crr1_score,
            "crr2": args.min_crr2_score,
            "crr4": args.min_crr4_score,
        },
        short_checks=config["short_checks"],
        available_crrs=config["available_crrs"],
    )

    print("Aligning CRISPR Spacers...")
    for crr in available_crrs:
        crr_upper = crr.upper()
        input_file = results_dir / crr_upper / f"{args.prefix}.{crr_upper}.fasta"

        # Skip if the file doesn't exist or is empty
        if not input_file.exists() or input_file.stat().st_size == 0:
            print(f"Skipping alignment for {crr_upper} (no spacers found)")
            continue

        align_spacers(
            input_file,
            f"{results_dir}/{crr_upper}/{args.prefix}.{crr_upper}Aligned.csv",
            f"{results_dir}/{crr_upper}/{args.prefix}.{crr_upper}Aligned.fasta",
            f"{results_dir}/{crr_upper}/{args.prefix}.{crr_upper}Unique.fasta",
            f"{results_dir}/{crr_upper}/{args.prefix}.{crr_upper}Binary.phy",
            crr_upper,
        )

    print("Combining Aligned CRISPR Cassettes")
    # Prepare the input files for combination
    valid_crrs = []
    for crr in available_crrs:
        crr_upper = crr.upper()
        fasta_file = results_dir / crr_upper / f"{args.prefix}.{crr_upper}Aligned.fasta"
        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            valid_crrs.append(crr)

    # Only proceed if we have at least 2 CRRs with data
    if len(valid_crrs) >= 2:
        # Collect files for each type
        fasta_files = []
        csv_files = []
        phy_files = []

        for crr in valid_crrs:
            crr_upper = crr.upper()
            fasta_file = (
                results_dir / crr_upper / f"{args.prefix}.{crr_upper}Aligned.fasta"
            )
            csv_file = results_dir / crr_upper / f"{args.prefix}.{crr_upper}Aligned.csv"
            phy_file = results_dir / crr_upper / f"{args.prefix}.{crr_upper}Binary.phy"

            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                fasta_files.append(fasta_file)
                csv_files.append(csv_file)
                phy_files.append(phy_file)

        # Output files
        fasta_output = results_dir / "AllCRR" / f"{args.prefix}.AllCRRAligned.fasta"
        csv_output = results_dir / "AllCRR" / f"{args.prefix}.AllCRRAligned.csv"
        phy_output = results_dir / "AllCRR" / f"{args.prefix}.AllCRRBinary.phy"

        # Call the dynamic combine functions
        if fasta_files:
            print(f"Combining aligned FASTA files from {len(fasta_files)} CRR regions")
            combine_fasta(fasta_files, fasta_output)

        if csv_files:
            print(f"Combining CSV files from {len(csv_files)} CRR regions")
            combine_csv(csv_files, csv_output)

        if phy_files:
            print(f"Combining binary phylip files from {len(phy_files)} CRR regions")
            combine_binary(phy_files, phy_output)
    else:
        print(
            f"Need at least 2 CRR regions with data to combine. Found {len(valid_crrs)}."
        )

    print("Analysis complete. Results available in the Results directory.")


if __name__ == "__main__":
    main()
