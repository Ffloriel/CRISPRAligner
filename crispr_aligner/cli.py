import argparse
import os
import shutil
from pathlib import Path

from crispr_aligner.sequences.fasta_utils import combine_fasta_files
from crispr_aligner.finder.crr_finder import find_crispr_spacers
from crispr_aligner.aligner.crr_aligner import align_spacers
from crispr_aligner.combiner.crr_combiner import combine_alignments

def create_folder_structure(base_dir):
    """Create the required folder structure for results."""
    results_dir = base_dir / "Results"
    if results_dir.exists():
        shutil.rmtree(results_dir)
    
    # Create directories
    for folder in ["Results", "Results/AllCRR", "Results/CRR1", "Results/CRR2", "Results/CRR4"]:
        (base_dir / folder).mkdir(parents=True, exist_ok=True)
    
    return results_dir

def main():
    parser = argparse.ArgumentParser(description="CRISPRAligner - Extract and align CRISPR spacers")
    parser.add_argument("prefix", help="Prefix for output files")
    parser.add_argument("--input_dir", default=".", help="Directory containing input FASTA files")
    
    args = parser.parse_args()
    
    base_dir = Path(args.input_dir)
    results_dir = create_folder_structure(base_dir)
    
    print("Combining FASTA Files...")
    combined_fasta = f"{args.prefix}.fasta"
    combine_fasta_files(base_dir, combined_fasta, args.prefix)
    
    print("Finding CRISPR Spacers...")
    find_crispr_spacers(combined_fasta, f"{args.prefix}.")
    
    print("Aligning CRISPR Spacers...")
    for crr in ["CRR1", "CRR2", "CRR4"]:
        align_spacers(
            f"{crr}/{args.prefix}.{crr}.fasta",
            f"{crr}/{args.prefix}.{crr}Aligned.csv",
            f"{crr}/{args.prefix}.{crr}Aligned.fasta",
            f"{crr}/{args.prefix}.{crr}Unique.fasta",
            f"{crr}/{args.prefix}.{crr}Binary.phy",
            crr
        )
    
    print("Combining Aligned CRISPR Cassettes")
    combine_alignments(
        f"CRR1/{args.prefix}.CRR1Aligned.fasta",
        f"CRR2/{args.prefix}.CRR2Aligned.fasta",
        f"CRR4/{args.prefix}.CRR4Aligned.fasta",
        f"AllCRR/{args.prefix}.AllCRRAligned.fasta",
        f"CRR1/{args.prefix}.CRR1Aligned.csv",
        f"CRR2/{args.prefix}.CRR2Aligned.csv",
        f"CRR4/{args.prefix}.CRR4Aligned.csv",
        f"AllCRR/{args.prefix}.AllCRRAligned.csv",
        f"CRR1/{args.prefix}.CRR1Binary.phy",
        f"CRR2/{args.prefix}.CRR2Binary.phy",
        f"CRR4/{args.prefix}.CRR4Binary.phy",
        f"AllCRR/{args.prefix}.AllCRRBinary.phy"
    )
    
    print("Analysis complete. Results available in the Results directory.")

if __name__ == "__main__":
    main()