#!/usr/bin/env python3

"""
Aggregate sequences across all genomes for specific BUSCO gene.
Replaces the transpose_hmms_with_sequences_all.py functionality.
"""

import sys
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import polars as pl

# Import the logging configuration
sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger

# Configure logger
logger = configure_logger(__file__)


def parse_busco_gene_from_sequence(seq_id, busco_gene):
    """
    Check if sequence ID corresponds to the target BUSCO gene.
    
    Args:
        seq_id: Sequence identifier
        busco_gene: Target BUSCO gene ID
        
    Returns:
        True if sequence belongs to this gene
    """
    # BUSCO sequences often have the gene ID in the sequence name
    return busco_gene in seq_id


def extract_taxonomy_for_taxon(taxon_id, taxonomy_df):
    """
    Get taxonomy string for a taxon ID.
    
    Args:
        taxon_id: Taxon identifier
        taxonomy_df: DataFrame with taxonomy mapping
        
    Returns:
        Taxonomy string in GTDB format
    """
    try:
        # Look up taxonomy for the taxon
        tax_row = taxonomy_df.filter(pl.col('taxon_id') == str(taxon_id))
        if len(tax_row) > 0:
            taxonomy = tax_row['parsed_taxon_string'][0]
            if taxonomy and ';' in taxonomy:
                return taxonomy
            logger.warning(f"Taxonomy too shallow for {taxon_id}: {taxonomy}")
            return None
        else:
            logger.warning(f"No taxonomy found for {taxon_id}")
            return None
    except Exception as e:
        logger.warning(f"Error getting taxonomy for {taxon_id}: {e}")
        return None


def aggregate_gene_sequences(genome_sequences_list, busco_gene, taxonomy_file, 
                           output_sequences, output_taxonomy):
    """
    Aggregate sequences for a specific BUSCO gene across all genomes.
    
    Args:
        genome_sequences_list: List of paths to genome sequence files
        busco_gene: Target BUSCO gene ID
        taxonomy_file: Path to taxonomy mapping file
        output_sequences: Path to output sequences file
        output_taxonomy: Path to output taxonomy file
    """
    logger.info(f"Aggregating sequences for BUSCO gene: {busco_gene}")
    logger.info(f"Processing {len(genome_sequences_list)} genome files")
    
    # Load taxonomy mapping
    try:
        taxonomy_df = pl.read_csv(taxonomy_file, separator='\t')
        logger.info(f"Loaded taxonomy mapping with {len(taxonomy_df)} entries")
    except Exception as e:
        logger.error(f"Failed to load taxonomy file {taxonomy_file}: {e}")
        raise
    
    # Collect sequences for this gene
    gene_sequences = []
    taxonomy_entries = []
    skipped_for_taxonomy = 0
    
    for genome_file in genome_sequences_list:
        genome_path = Path(genome_file)
        if not genome_path.exists():
            logger.debug(f"Genome file does not exist: {genome_file}")
            continue
            
        # Extract taxon_id from filename (e.g., "12345.faa" -> "12345")
        taxon_id = genome_path.stem
        
        # Read sequences and filter for this gene
        try:
            for record in SeqIO.parse(genome_file, "fasta"):
                if parse_busco_gene_from_sequence(record.id, busco_gene):
                    # Extract original taxon_id from sequence ID if prefixed
                    if '|' in record.id:
                        original_taxon_id = record.id.split('|')[0]
                    else:
                        original_taxon_id = taxon_id
                    
                    # Get taxonomy for this sequence
                    taxonomy = extract_taxonomy_for_taxon(original_taxon_id, taxonomy_df)
                    if taxonomy is None:
                        skipped_for_taxonomy += 1
                        continue

                    gene_sequences.append(record)
                    taxonomy_entries.append((record.id, taxonomy))
                    
        except Exception as e:
            logger.warning(f"Error reading {genome_file}: {e}")
            continue
    
    logger.info(f"Found {len(gene_sequences)} sequences for {busco_gene}")
    if skipped_for_taxonomy > 0:
        logger.info(f"Skipped {skipped_for_taxonomy} sequences for {busco_gene} due to unresolved/shallow taxonomy")
    
    # Write sequences
    output_sequences_path = Path(output_sequences)
    output_sequences_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_sequences, 'w') as f:
        for record in gene_sequences:
            f.write(f">{record.id}\n{str(record.seq)}\n")
    
    # Write taxonomy
    output_taxonomy_path = Path(output_taxonomy)
    output_taxonomy_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_taxonomy, 'w') as f:
        f.write("sequence_id\ttaxonomy\n")
        for seq_id, taxonomy in taxonomy_entries:
            f.write(f"{seq_id}\t{taxonomy}\n")
    
    logger.info(f"Wrote {len(gene_sequences)} sequences to {output_sequences}")
    logger.info(f"Wrote {len(taxonomy_entries)} taxonomy entries to {output_taxonomy}")
    
    return len(gene_sequences)


def main():
    """Main function using snakemake inputs."""
    # Get inputs from snakemake
    genome_sequences_list = snakemake.input.genome_sequences
    busco_gene = snakemake.wildcards.busco_gene
    taxonomy_file = snakemake.input.taxonomy
    output_sequences = snakemake.output.sequences
    output_taxonomy = snakemake.output.taxonomy
    
    try:
        sequence_count = aggregate_gene_sequences(
            genome_sequences_list, busco_gene, taxonomy_file,
            output_sequences, output_taxonomy
        )
        
        if sequence_count == 0:
            logger.warning(f"No sequences found for gene {busco_gene}")
            # Create empty files to satisfy snakemake
            Path(output_sequences).touch()
            Path(output_taxonomy).touch()
        
        logger.info(f"Successfully aggregated {sequence_count} sequences for {busco_gene}")
        
    except Exception as e:
        logger.error(f"Error aggregating gene {busco_gene}: {e}")
        raise


if __name__ == "__main__":
    main()
