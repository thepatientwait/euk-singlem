#!/usr/bin/env python3

"""
Extract transcript nucleotide sequences from genome using BUSCO coordinates.
BUSCO outputs protein sequences but we need nucleotide transcripts for SingleM.
"""

import sys
import gzip
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from typing import Tuple
import polars as pl

# Import the logging configuration
sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger

# Configure logger
logger = configure_logger(__file__)


def parse_busco_file(busco_file: Path) -> Tuple[str, str, pl.DataFrame]:
    """
    Parse a BUSCO results file and extract group, genome_acc, and complete genes.
    
    Args:
        busco_file: Path to BUSCO full_table.tsv file
        
    Returns:
        Tuple of (group, genome_acc, complete_genes_dataframe)
    """
    # Extract group and genome accession from path.
    # Path format: outputs/busco/{group}/{genome_acc}/run_eukaryota_odb12/full_table.tsv
    path_parts = busco_file.parts
    group = path_parts[-4]
    genome_acc = path_parts[-3]
    
    logger.info(f"Parsing BUSCO file for {group}/{genome_acc}")
    
    # Read BUSCO full table, skipping comment lines that start with #
    busco_df = pl.read_csv(busco_file, separator='\t', skip_rows=2).rename({'# Busco id': 'Busco id'})
    
    # Filter for complete single-copy genes using the correct column name
    complete_genes = busco_df.filter(pl.col('Status') == 'Complete')
    
    logger.info(f"Found {len(complete_genes)} complete genes for {group}/{genome_acc}")
    return group, genome_acc, complete_genes


def load_genome_sequences(genome_file):
    """
    Load genome sequences from FASTA file.
    
    Args:
        genome_file: Path to genome FASTA file (may be gzipped)
        
    Returns:
        dict: Sequence dictionary {seq_id: seq_record}
    """
    sequences = {}
    
    try:
        if genome_file.endswith('.gz'):
            with gzip.open(genome_file, 'rt') as f:
                for record in SeqIO.parse(f, 'fasta'):
                    sequences[record.id] = record
        else:
            for record in SeqIO.parse(genome_file, 'fasta'):
                sequences[record.id] = record
                
        logger.info(f"Loaded {len(sequences)} sequences from {genome_file}")
        
    except Exception as e:
        logger.error(f"Error loading genome file {genome_file}: {e}")
        raise
        
    return sequences


def extract_transcript_sequence(genome_sequences, seqid, start, end, strand, gene_id):
    """
    Extract transcript sequence from genome.
    
    Args:
        genome_sequences: Dictionary of genome sequences
        seqid: Sequence/contig identifier
        start: Start coordinate (1-based)
        end: End coordinate (1-based)
        strand: Strand ('+' or '-')
        gene_id: Gene identifier
        
    Returns:
        SeqRecord: Transcript sequence record
    """
    if seqid not in genome_sequences:
        logger.warning(f"Sequence {seqid} not found in genome for gene {gene_id}")
        return None
        
    genome_seq = genome_sequences[seqid]
    genome_len = len(genome_seq.seq)

    # BUSCO reports 1-based inclusive coordinates.
    left = min(start, end)
    right = max(start, end)

    if left < 1 or right > genome_len:
        logger.warning(
            f"Coordinates out of bounds for gene {gene_id} on {seqid}: "
            f"{start}-{end} (normalized {left}-{right}, contig length {genome_len})"
        )
        return None

    # Convert 1-based inclusive -> Python 0-based end-exclusive.
    transcript_seq = genome_seq.seq[left - 1:right]

    # Return transcript in biological orientation.
    if strand == '-':
        transcript_seq = transcript_seq.reverse_complement()
    
    # Create transcript record
    transcript_record = SeqRecord(
        transcript_seq,
        id=gene_id,
        description=(
            f"BUSCO transcript {gene_id} from {seqid}:{start}-{end} strand {strand} "
            f"(normalized interval {left}-{right})"
        )
    )
    
    return transcript_record


def extract_busco_transcripts(genome_file, full_table_file, output_file):
    """
    Main function to extract BUSCO transcript sequences.
    
    Args:
        genome_file: Path to genome FASTA file
        full_table_file: Path to BUSCO full_table.tsv file
        output_file: Path to output transcript FASTA file
    """
    logger.info(f"Extracting transcripts from {genome_file}")
    logger.info(f"Using BUSCO full table: {full_table_file}")
    
    # Load genome sequences
    genome_sequences = load_genome_sequences(genome_file)
    
    # Parse BUSCO full table using polars
    group, genome_acc, complete_genes = parse_busco_file(Path(full_table_file))
    logger.info(f"Processing {group}/{genome_acc} with {len(complete_genes)} complete BUSCO genes")
    
    if len(complete_genes) == 0:
        logger.warning("No complete BUSCO genes found in full table")
        # Create empty output file
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        Path(output_file).touch()
        return 0
    
    transcript_records = []
    genes_processed = 0
    genes_found = 0
    
    # Extract transcript sequences for each gene
    for row in complete_genes.iter_rows(named=True):
        genes_processed += 1
        
        busco_id = str(row['Busco id'])
        sequence = str(row['Sequence'])
        start = int(row['Gene Start'])
        end = int(row['Gene End'])
        strand = str(row['Strand'])
        
        transcript_record = extract_transcript_sequence(
            genome_sequences, sequence, start, end, strand, busco_id
        )
        
        if transcript_record:
            transcript_records.append(transcript_record)
            genes_found += 1
        
    logger.info(f"Processed {genes_processed} genes, extracted {genes_found} transcripts")
    
    # Write transcript sequences
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        SeqIO.write(transcript_records, f, 'fasta')
    
    logger.info(f"Wrote {len(transcript_records)} transcript sequences to {output_file}")
    
    return len(transcript_records)


def main():
    """Main function using snakemake inputs."""
    genome_file = snakemake.input.genome
    full_table_file = snakemake.input.full_table
    output_file = snakemake.output.transcripts
    
    try:
        transcript_count = extract_busco_transcripts(genome_file, full_table_file, output_file)
        
        if transcript_count == 0:
            logger.warning("No transcripts extracted")
            # Create empty file to satisfy snakemake
            Path(output_file).touch()
        
        logger.info(f"Successfully extracted {transcript_count} transcripts")
        
    except Exception as e:
        logger.error(f"Error extracting transcripts: {e}")
        raise


if __name__ == "__main__":
    main()