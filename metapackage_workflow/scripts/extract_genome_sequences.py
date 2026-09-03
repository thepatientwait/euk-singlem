#!/usr/bin/env python3

"""
Extract sequences for HMM matches from individual genome.
Replaces the bulk mfqe_all.py processing with per-genome processing.
"""

import sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# Import the logging configuration
sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger

# Configure logger
logger = configure_logger(__file__)


def extract_sequences_for_genome(matches_file, genome_file, output_file):
    """
    Extract protein sequences that have HMM matches.
    
    Args:
        matches_file: Path to matches file (from hmmsearch)
        genome_file: Path to genome protein sequences
        output_file: Path to output extracted sequences
    """
    logger.info(f"Processing genome: {genome_file}")
    logger.info(f"Using matches from: {matches_file}")
    
    # Read matches - format: HMM_name sequence_name E-value
    matches = set()
    if Path(matches_file).exists() and Path(matches_file).stat().st_size > 0:
        with open(matches_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        sequence_name = parts[1]  # sequence name is second column
                        matches.add(sequence_name)
    
    logger.info(f"Found {len(matches)} sequence matches to extract")
    
    # Extract matched sequences
    extracted_count = 0
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as out_f:
        if Path(genome_file).exists():
            for record in SeqIO.parse(genome_file, "fasta"):
                # Check if this sequence has any matches
                if record.id in matches or record.id.split('|')[0] in matches:
                    # Write sequence with genome prefix for uniqueness
                    taxon_id = Path(genome_file).stem  # Get taxon_id from filename
                    unique_id = f"{taxon_id}|{record.id}"
                    out_f.write(f">{unique_id}\n{str(record.seq)}\n")
                    extracted_count += 1
    
    logger.info(f"Extracted {extracted_count} sequences to {output_file}")
    return extracted_count


def main():
    """Main function using snakemake inputs."""
    # Get inputs from snakemake
    matches_file = snakemake.input.matches
    genome_file = snakemake.input.genome
    output_file = snakemake.output.sequences
    
    try:
        extracted_count = extract_sequences_for_genome(matches_file, genome_file, output_file)
        logger.info(f"Successfully extracted {extracted_count} sequences")
        
        if extracted_count == 0:
            logger.warning(f"No sequences extracted for {genome_file}")
            # Create empty output file to satisfy snakemake
            Path(output_file).touch()
            
    except Exception as e:
        logger.error(f"Error processing {genome_file}: {e}")
        raise


if __name__ == "__main__":
    main()