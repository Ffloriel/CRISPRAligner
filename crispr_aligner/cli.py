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


def create_folder_structure(base_dir):
    """Create the basic folder structure for results."""
    results_dir = base_dir / "Results"
    if results_dir.exists():
        shutil.rmtree(results_dir)

    # Always create Results and AllCRR directory
    (base_dir / "Results").mkdir(parents=True, exist_ok=True)
    (base_dir / "Results/AllCRR").mkdir(parents=True, exist_ok=True)

    # Frame-specific directories will be created by find_crispr_spacers
    # based on the actual frames found in the sequences

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

    parser.add_argument(
        "--max_frame_size",
        type=int,
        default=20000,
        help="Maximum size in bp for a CRISPR array frame (default: 20000)",
    )

    args = parser.parse_args()

    # Load appropriate configuration
    if args.custom_config:
        config = load_bacteria_config(args.custom_config)
    else:
        config = load_bacteria_config(args.bacteria)

    base_dir = Path(args.input_dir)
    available_crrs = config.get("available_crrs", ["crr1", "crr2", "crr4"])
    results_dir = create_folder_structure(base_dir)

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
        available_crrs=available_crrs,
        max_frame_size=args.max_frame_size,  # Pass max_frame_size parameter
    )

    # Discover all frame directories created by find_crispr_spacers
    frame_dirs = []
    for item in results_dir.iterdir():
        if item.is_dir() and item.name not in ["AllCRR"]:
            frame_dirs.append(item.name)

    print(f"Found {len(frame_dirs)} CRISPR array frames: {', '.join(frame_dirs)}")

    print("Aligning CRISPR Spacers...")
    for frame_dir in frame_dirs:
        input_file = results_dir / frame_dir / f"{args.prefix}.{frame_dir}.fasta"

        # Skip if the file doesn't exist or is empty
        if not input_file.exists() or input_file.stat().st_size == 0:
            print(f"Skipping alignment for {frame_dir} (no spacers found)")
            continue

        align_spacers(
            input_file,
            results_dir / frame_dir / f"{args.prefix}.{frame_dir}Aligned.csv",
            results_dir / frame_dir / f"{args.prefix}.{frame_dir}Aligned.fasta",
            results_dir / frame_dir / f"{args.prefix}.{frame_dir}Unique.fasta",
            results_dir / frame_dir / f"{args.prefix}.{frame_dir}Binary.phy",
            frame_dir,
        )

    print("Combining Aligned CRISPR Cassettes")
    # Find frames with valid data for combination
    valid_frames = []
    for frame_dir in frame_dirs:
        fasta_file = results_dir / frame_dir / f"{args.prefix}.{frame_dir}Aligned.fasta"
        if fasta_file.exists() and fasta_file.stat().st_size > 0:
            valid_frames.append(frame_dir)

    # Only proceed if we have at least 2 frame directories with data
    if len(valid_frames) >= 2:
        # Collect files for each type
        fasta_files = []
        csv_files = []
        phy_files = []

        for frame_dir in valid_frames:
            fasta_file = (
                results_dir / frame_dir / f"{args.prefix}.{frame_dir}Aligned.fasta"
            )
            csv_file = results_dir / frame_dir / f"{args.prefix}.{frame_dir}Aligned.csv"
            phy_file = results_dir / frame_dir / f"{args.prefix}.{frame_dir}Binary.phy"

            if fasta_file.exists() and fasta_file.stat().st_size > 0:
                fasta_files.append(fasta_file)
                csv_files.append(csv_file)
                phy_files.append(phy_file)

        # Output files
        fasta_output = results_dir / "AllCRR" / f"{args.prefix}.AllCRRAligned.fasta"
        csv_output = results_dir / "AllCRR" / f"{args.prefix}.AllCRRAligned.csv"
        phy_output = results_dir / "AllCRR" / f"{args.prefix}.AllCRRBinary.phy"

        # Call the combine functions
        if fasta_files:
            print(
                f"Combining aligned FASTA files from {len(fasta_files)} CRISPR frames"
            )
            combine_fasta(fasta_files, fasta_output)

        if csv_files:
            print(f"Combining CSV files from {len(csv_files)} CRISPR frames")
            combine_csv(csv_files, csv_output)

        if phy_files:
            print(f"Combining binary phylip files from {len(phy_files)} CRISPR frames")
            combine_binary(phy_files, phy_output)
    else:
        print(
            f"Need at least 2 CRISPR frames with data to combine. Found {len(valid_frames)}."
        )

    print("Analysis complete. Results available in the Results directory.")


if __name__ == "__main__":
    main()
