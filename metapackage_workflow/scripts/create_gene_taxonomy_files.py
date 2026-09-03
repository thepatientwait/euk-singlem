# init
# ----
import polars as pl
from Bio import SeqIO
import os

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
FAA_FILE = snakemake.input['faa']
TAXONOMY_FILE = snakemake.input['taxonomy']

TAX_OUTPUT = snakemake.output['tax']


# functions
# ---------
def extract_sequence_ids(faa_file):
    """Extract sequence IDs and taxon IDs from FASTA file"""
    seq_ids = []
    taxon_ids = []
    
    for seq_record in SeqIO.parse(faa_file, "fasta"):
        seq_id = seq_record.id
        taxon_id = seq_id.split('~')[0]  # Extract taxon_id from "taxon_id~gene_id"
        seq_ids.append(seq_id)
        taxon_ids.append(taxon_id)
    
    return seq_ids, taxon_ids


def create_taxonomy_mapping(seq_ids, taxon_ids, taxonomy_file):
    """Join sequence IDs with taxonomy strings"""
    # Load taxonomy mapping
    taxonomy = pl.read_csv(taxonomy_file, separator='\t')
    
    # Create mapping dataframe
    seq_df = pl.DataFrame({'seq_id': seq_ids, 'taxon_id': taxon_ids})
    
    # Join with taxonomy
    taxonomy_map = seq_df.join(
        taxonomy.select(['taxon_id', 'parsed_taxon_string']),
        on='taxon_id',
        how='left'
    )
    
    return taxonomy_map


# main
# ----
def main():
    logger.info(f"Creating taxonomy file for {FAA_FILE}")
    
    # Extract sequence IDs from FASTA
    seq_ids, taxon_ids = extract_sequence_ids(FAA_FILE)
    logger.info(f"Extracted {len(seq_ids)} sequence IDs")
    
    # Create taxonomy mapping
    taxonomy_map = create_taxonomy_mapping(seq_ids, taxon_ids, TAXONOMY_FILE)
    
    # Write taxonomy file
    with open(TAX_OUTPUT, 'w') as f:
        for row in taxonomy_map.iter_rows(named=True):
            f.write(f"{row['seq_id']}\t{row['parsed_taxon_string']}\n")
    
    logger.info(f"Wrote taxonomy mapping to {TAX_OUTPUT}")


if __name__ == '__main__':
    main()
