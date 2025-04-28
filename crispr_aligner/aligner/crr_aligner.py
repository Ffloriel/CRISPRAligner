"""CRISPR repeat region (CRR) aligner module for processing CRISPR spacers."""

from pathlib import Path
from typing import List, Tuple, Dict, Union, Optional
from Bio import SeqIO


def match_crr(query: str, reference: str) -> Tuple[int, int]:
    """
    Calculate forward and reverse matches between query and reference sequences.
    
    Args:
        query: Query sequence to compare
        reference: Reference sequence to compare against
        
    Returns:
        Tuple containing (forward_match_count, reverse_match_count)
    """
    crr_match = 0
    crr_rev_match = 0
    
    # Calculate forward matches
    for i in range(min(len(query), len(reference))):
        if query[i] == reference[i]:
            crr_match += 1
    
    # Calculate reverse matches
    for i in range(1, min(len(query), len(reference)) + 1):
        if query[-i] == reference[-i]:
            crr_rev_match += 1
            
    return (crr_match, crr_rev_match)


def parse_fasta_by_cassette(input_file: Union[str, Path]) -> List[List]:
    """
    Parse FASTA file and organize sequences by cassette.
    
    Args:
        input_file: Path to input FASTA file
        
    Returns:
        List of samples, each containing ID and sequences
    """
    input_file = Path(input_file)
    
    sample_list = []
    start = True
    first = True
    sample = []
    
    with input_file.open("r") as seq_file:
        seq_obj = SeqIO.parse(seq_file, "fasta")
        
        for record in seq_obj:
            if start:
                start = False
                pre_name = record.id
                
            if pre_name in record.id:
                sample.append(str(record.seq))
            else:
                if first:
                    sample.insert(0, pre_name)
                    sample_list.append(sample)
                    first = False
                    pre_name = record.id
                    sample = []
                    sample.append(str(record.seq))
                else:
                    inserted = False
                    sample.insert(0, pre_name)
                    
                    for a in range(len(sample_list)):
                        if len(sample) > len(sample_list[a]):
                            inserted = True
                            sample_list.insert(a, sample)
                            break
                            
                    if not inserted:
                        sample_list.append(sample)
                        
                    pre_name = record.id
                    sample = []
                    sample.append(str(record.seq))
    
    # Process the last sample
    if sample:
        sample.insert(0, pre_name)
        inserted = False
        
        for b in range(len(sample_list)):
            if len(sample) > len(sample_list[b]):
                inserted = True
                sample_list.insert(b, sample)
                break
                
        if not inserted:
            sample_list.append(sample)
    
    return sample_list


def generate_consensus(sample_list: List[List]) -> List[str]:
    """
    Generate consensus sequences from the sample list.
    
    Args:
        sample_list: List of samples organized by cassette
        
    Returns:
        List of consensus sequences with 'Consensus' as first element
    """
    first = True
    consensus = []
    upstream = []
    
    for cassette in sample_list:
        last = False
        
        if first:
            first = False
            for c in range(1, len(cassette)):
                consensus.append(cassette[c])
        else:
            no_crr_match = True
            
            for d in range(1, len(cassette)):
                match = False
                
                for e in range(len(consensus)):
                    consensus_match = match_crr(cassette[d], consensus[e])
                    
                    if (consensus_match[0] >= (len(consensus[e]) - 2) or 
                        consensus_match[1] >= (len(consensus[e]) - 2)):
                        
                        if len(cassette[d]) > len(consensus[e]):
                            consensus[e] = cassette[d]
                            
                        no_crr_match = False
                        match = True
                        
                        if len(upstream) >= 1:
                            upstream = upstream[::-1]
                            for f in range(len(upstream)):
                                consensus.insert(e, upstream[f])
                                
                        upstream = []
                        last = False
                        break
                        
                if not match:
                    upstream.append(cassette[d])
                    last = True
                    
            if no_crr_match:
                consensus.extend(upstream)
                upstream = []
                
        if last and upstream:
            consensus.extend(upstream)
            upstream = []
    
    # Add the 'Consensus' label as the first element
    consensus.insert(0, "Consensus")
    
    return consensus


def write_aligned_files(
    consensus: List[str], 
    sample_list: List[List], 
    csv_path: Union[str, Path],
    align_path: Union[str, Path],
    unique_path: Union[str, Path],
    binary_path: Union[str, Path],
    crr_type: str
) -> None:
    """
    Write aligned files in various formats based on consensus and samples.
    
    Args:
        consensus: List of consensus sequences with 'Consensus' as first element
        sample_list: List of samples organized by cassette
        csv_path: Path to output CSV file
        align_path: Path to output aligned FASTA file
        unique_path: Path to output unique FASTA file
        binary_path: Path to output binary phylip file
        crr_type: CRR type identifier (e.g., 'CRR1', 'CRR2', 'CRR4')
    """
    csv_path = Path(csv_path)
    align_path = Path(align_path)
    unique_path = Path(unique_path)
    binary_path = Path(binary_path)
    
    # Create parent directories if they don't exist
    for path in [csv_path, align_path, unique_path, binary_path]:
        path.parent.mkdir(parents=True, exist_ok=True)
    
    with csv_path.open("w") as csv_file, \
         align_path.open("w") as align_file, \
         unique_path.open("w") as unique_file, \
         binary_path.open("w") as binary_file:
        
        # Write consensus to CSV
        csv_file.write(f"{consensus[0]},")
        for h in range(1, len(consensus)):
            csv_file.write(f"{consensus[h]},")
            unique_file.write(f">Consensus{h}\n{consensus[h]}\n")
        csv_file.write("\n")
        
        # Write consensus to aligned FASTA
        align_file.write(f">{consensus[0]}\n{''.join(consensus[1:])}\n")
        
        # Write header for binary file
        binary_file.write(f"{len(sample_list)} {len(consensus)-1}\n")
        
        # Process each sample
        for i in range(len(sample_list)):
            csv_file.write(f"{sample_list[i][0]},")
            align_file.write(f">{sample_list[i][0]}\n")
            
            # Write first part of sample ID to binary file
            binary_count = 0
            for p in range(3, len(sample_list[i][0])):
                if binary_count == 10:
                    break
                binary_file.write(f"{sample_list[i][0][p]}")
                binary_count += 1
                
            # Pad ID with spaces if needed
            if len(sample_list[i][0]) < 13:
                binary_file.write(" " * (13 - len(sample_list[i][0])))
            
            # Process each consensus sequence
            for k in range(1, len(consensus)):
                hsp = False
                
                # Find matching sequence in sample
                for j in range(1, len(sample_list[i])):
                    consensus_match = match_crr(sample_list[i][j], consensus[k])
                    
                    if (consensus_match[0] >= len(consensus[k]) - 2 or 
                        consensus_match[1] >= len(consensus[k]) - 2):
                        
                        hsp = True
                        csv_file.write(f"{crr_type}S{j},")
                        binary_file.write("1")
                        
                        # Handle alignment gaps based on match orientation
                        if (len(sample_list[i][j]) != len(consensus[k]) and 
                            consensus_match[0] < consensus_match[1]):
                            align_file.write("-" * (len(consensus[k]) - len(sample_list[i][j])))
                        
                        align_file.write(f"{sample_list[i][j]}")
                        
                        if (len(sample_list[i][j]) != len(consensus[k]) and 
                            consensus_match[0] > consensus_match[1]):
                            align_file.write("-" * (len(consensus[k]) - len(sample_list[i][j])))
                        
                        sample_list[i][j] = "Taken"
                        break
                
                # If no match found, write gaps
                if not hsp:
                    csv_file.write(",")
                    binary_file.write("0")
                    align_file.write("-" * len(consensus[k]))
            
            csv_file.write("\n")
            align_file.write("\n")
            binary_file.write("\n")
    
    print(f"Alignment files written to:")
    print(f"  CSV:    {csv_path}")
    print(f"  FASTA:  {align_path}")
    print(f"  Unique: {unique_path}")
    print(f"  Binary: {binary_path}")


def align_spacers(
    input_file: Union[str, Path],
    csv_output: Union[str, Path],
    aligned_output: Union[str, Path],
    unique_output: Union[str, Path],
    binary_output: Union[str, Path],
    crr_type: str
) -> None:
    """
    Align CRISPR spacers and generate output in multiple formats.
    
    Args:
        input_file: Path to input FASTA file with CRISPR spacers
        csv_output: Path to output CSV file
        aligned_output: Path to output aligned FASTA file
        unique_output: Path to output unique FASTA file
        binary_output: Path to output binary phylip file
        crr_type: CRR type identifier (e.g., 'CRR1', 'CRR2', 'CRR4')
    """
    input_file = Path(input_file)
    
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    print(f"Processing {crr_type} spacers from {input_file}")
    
    # Parse FASTA and organize by cassette
    sample_list = parse_fasta_by_cassette(input_file)
    
    if not sample_list:
        print(f"Warning: No cassettes found in {input_file}")
        return
    
    # Generate consensus sequence
    consensus = generate_consensus(sample_list)
    
    # Write aligned files
    write_aligned_files(
        consensus,
        sample_list,
        csv_output,
        aligned_output,
        unique_output,
        binary_output,
        crr_type
    )
    
    print(f"Successfully aligned {crr_type} spacers")