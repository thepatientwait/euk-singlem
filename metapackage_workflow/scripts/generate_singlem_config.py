# init
# ----
import yaml
import os
import polars as pl
from pathlib import Path

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
SEQUENCES_DIR = snakemake.input['sequences_dir']
TAXONOMY_FILE = snakemake.input['taxonomy']

CONFIG_OUTPUT = snakemake.output['config']

TARGET_DOMAINS = snakemake.params['target_domains']
SINGLEM_OUTPUT_DIR = snakemake.params['output_directory']
METAPACKAGE_PATH = snakemake.params['metapackage_path']


# functions
# ---------
def get_available_genes(sequences_dir):
    """Get list of genes that have both fasta and taxonomy files"""
    sequences_path = Path(sequences_dir)
    available_genes = []
    
    # Scan directory for .faa files
    candidate_genes = [f.stem for f in sequences_path.glob("*.faa")]
    
    for gene_id in candidate_genes:
        faa_file = sequences_path / f"{gene_id}.faa"
        tax_file = sequences_path / f"{gene_id}_taxonomy.tsv"
        
        if faa_file.exists() and tax_file.exists():
            available_genes.append(gene_id)
    
    return available_genes


def generate_config_dict(genes, target_domains, output_dir, metapackage_path, sequences_dir):
    """Generate SingleM config dictionary"""
    config_dict = {
        'singlem_packages': {},
        'singlem_new_metapackage_working_directory': os.path.abspath(output_dir),
        'metapackage_path': os.path.abspath(metapackage_path)
    }
    
    for gene_id in genes:
        config_dict['singlem_packages'][gene_id] = {
            'aa_sequences_file': os.path.abspath(os.path.join(sequences_dir, f"{gene_id}.faa")),
            'taxonomy_file': os.path.abspath(os.path.join(sequences_dir, f"{gene_id}_taxonomy.tsv")),
            'gene_description': gene_id,  # Using gene ID as description for now
            'target_domains': target_domains
        }
    
    return config_dict


# main
# ----
def main():
    logger.info(f"Generating SingleM config from sequences in {SEQUENCES_DIR}")
    
    # Get available genes
    genes = get_available_genes(SEQUENCES_DIR)
    logger.info(f"Found {len(genes)} genes with both sequence and taxonomy files")
    
    if len(genes) == 0:
        logger.error("No genes found with both .faa and _taxonomy.tsv files")
        raise RuntimeError("No valid genes for SingleM package creation")
    
    # Generate config
    config_dict = generate_config_dict(genes, TARGET_DOMAINS, SINGLEM_OUTPUT_DIR, METAPACKAGE_PATH, SEQUENCES_DIR)
    
    # Write config
    with open(CONFIG_OUTPUT, 'w') as f:
        yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Generated config for {len(genes)} genes at {CONFIG_OUTPUT}")


if __name__ == '__main__':
    main()
