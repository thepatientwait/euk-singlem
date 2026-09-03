# init
# ----
import polars as pl
from pathlib import Path
import os

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
FULL_TABLES = snakemake.input['full_tables']

SUMMARY_OUTPUT = snakemake.output['summary']
COVERAGE_OUTPUT = snakemake.output['coverage']


# functions
# ---------
def parse_path_info(full_table_path):
    """Extract group and taxon_id from BUSCO output path"""
    parts = Path(full_table_path).parts
    group = parts[parts.index('busco') + 1]
    taxon_id = parts[parts.index('busco') + 2]
    return group, taxon_id


def read_busco_table(full_table_path, group, taxon_id):
    """Read a single BUSCO full_table.tsv file"""
    df = pl.read_csv(
        full_table_path,
        separator='\t',
        comment_prefix='#',
        has_header=False,
        new_columns=[
            'busco_id', 'status', 'sequence', 'start', 'end', 
            'strand', 'score', 'length', 'OrthoDB_url', 'description'
        ]
    )
    
    # Add genome info
    df = df.with_columns([
        pl.lit(group).alias('group'),
        pl.lit(taxon_id).alias('taxon_id')
    ])
    
    return df


def calculate_coverage_stats(combined_df):
    """Calculate gene coverage statistics across all genomes"""
    total_genomes = combined_df['taxon_id'].n_unique()
    
    coverage_stats = (
        combined_df
        .filter(pl.col('status') == 'Complete')
        .group_by('busco_id')
        .agg([
            pl.count().alias('n_genomes'),
            (pl.count() / total_genomes * 100).alias('coverage_pct')
        ])
        .sort('coverage_pct', descending=True)
    )
    
    return coverage_stats, total_genomes


# main
# ----
def main():
    logger.info(f"Aggregating {len(FULL_TABLES)} BUSCO full tables...")
    
    all_data = []
    
    for full_table_path in FULL_TABLES:
        group, taxon_id = parse_path_info(full_table_path)
        logger.info(f"Processing {group}/{taxon_id}")
        
        df = read_busco_table(full_table_path, group, taxon_id)
        all_data.append(df)
    
    # Combine all tables
    logger.info("Combining all BUSCO tables...")
    combined = pl.concat(all_data)
    combined.write_csv(SUMMARY_OUTPUT, separator='\t')
    logger.info(f"Wrote summary to {SUMMARY_OUTPUT}")
    
    # Calculate coverage statistics
    logger.info("Calculating coverage statistics...")
    coverage_stats, total_genomes = calculate_coverage_stats(combined)
    coverage_stats.write_csv(COVERAGE_OUTPUT, separator='\t')
    logger.info(f"Wrote coverage stats to {COVERAGE_OUTPUT}")
    
    # Print summary
    logger.info(f"Total genomes: {total_genomes}")
    logger.info(f"Total genes: {coverage_stats.height}")
    logger.info(f"Genes with ≥90% coverage: {coverage_stats.filter(pl.col('coverage_pct') >= 90).height}")
    logger.info(f"Genes with ≥80% coverage: {coverage_stats.filter(pl.col('coverage_pct') >= 80).height}")
    logger.info(f"Genes with ≥70% coverage: {coverage_stats.filter(pl.col('coverage_pct') >= 70).height}")


if __name__ == '__main__':
    main()
