import argparse
import os
import shutil
from pathlib import Path

from crispr_aligner.sequences.fasta_utils import combine_fasta_files
from crispr_aligner.finder.crr_finder import find_crispr_spacers
from crispr_aligner.aligner.crr_aligner import align_spacers
from crispr_aligner.combiner.crr_combiner import combine_alignments
from crispr_aligner.utils.config import load_bacteria_config, get_available_bacteria


def create_folder_structure(base_dir):
    """Create the required folder structure for results."""
    results_dir = base_dir / "Results"
    if results_dir.exists():
        shutil.rmtree(results_dir)

    # Create directories
    for folder in [
        "Results",
        "Results/AllCRR",
        "Results/CRR1",
        "Results/CRR2",
        "Results/CRR4",
    ]:
        (base_dir / folder).mkdir(parents=True, exist_ok=True)

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
    )

    print("Aligning CRISPR Spacers...")
    for crr in ["CRR1", "CRR2", "CRR4"]:
        align_spacers(
            f"{results_dir}/{crr}/{args.prefix}.{crr}.fasta",
            f"{results_dir}/{crr}/{args.prefix}.{crr}Aligned.csv",
            f"{results_dir}/{crr}/{args.prefix}.{crr}Aligned.fasta",
            f"{results_dir}/{crr}/{args.prefix}.{crr}Unique.fasta",
            f"{results_dir}/{crr}/{args.prefix}.{crr}Binary.phy",
            crr,
        )

    print("Combining Aligned CRISPR Cassettes")
    combine_alignments(
        f"{results_dir}/CRR1/{args.prefix}.CRR1Aligned.fasta",
        f"{results_dir}/CRR2/{args.prefix}.CRR2Aligned.fasta",
        f"{results_dir}/CRR4/{args.prefix}.CRR4Aligned.fasta",
        f"{results_dir}/AllCRR/{args.prefix}.AllCRRAligned.fasta",
        f"{results_dir}/CRR1/{args.prefix}.CRR1Aligned.csv",
        f"{results_dir}/CRR2/{args.prefix}.CRR2Aligned.csv",
        f"{results_dir}/CRR4/{args.prefix}.CRR4Aligned.csv",
        f"{results_dir}/AllCRR/{args.prefix}.AllCRRAligned.csv",
        f"{results_dir}/CRR1/{args.prefix}.CRR1Binary.phy",
        f"{results_dir}/CRR2/{args.prefix}.CRR2Binary.phy",
        f"{results_dir}/CRR4/{args.prefix}.CRR4Binary.phy",
        f"{results_dir}/AllCRR/{args.prefix}.AllCRRBinary.phy",
    )

    print("Analysis complete. Results available in the Results directory.")


if __name__ == "__main__":
    main()
