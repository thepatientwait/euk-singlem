#!/usr/bin/env python3
"""
Assign eukaryotic taxonomy to SingleM OTU tables from transcript processing.
Adapted from assign_viral_taxonomy.py for eukaryotic metapackage creation.
"""

# init
# ----
import os
import polars as pl

# logging
# -------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)


# snakefile parameters
# --------------------
OTU_TABLE           = snakemake.input['otu_table']
TAXONOMY_MAPPING    = snakemake.input['taxonomy_mapping']

GENOME_ACC          = snakemake.wildcards['genome_acc']

TAXONOMY_OTU        = snakemake.output['taxonomy_otu']

# functions
# ---------
UNKNOWN_TAXONOMY = (
    "d__Eukaryota; p__unclassified; c__unclassified; "
    "o__unclassified; f__unclassified; g__unclassified; s__unclassified"
)


def _normalize_accession(genome_acc: str) -> str:
    # E.g. GCA_000258005.2 -> GCA_000258005
    return genome_acc.split('.', 1)[0]


def get_taxon_string(file_path, genome_acc):
    mapping_df = pl.read_csv(file_path, separator='\t')

    # 1) Exact accession match first.
    exact = mapping_df.filter(pl.col('acc') == genome_acc).select('parsed_taxon_string').to_series().to_list()
    if len(exact) == 1:
        return exact[0]
    if len(exact) > 1:
        logger.warning(
            f"Multiple taxonomy rows found for exact accession '{genome_acc}'. Using first match."
        )
        return exact[0]

    # 2) Fallback to accession without version suffix.
    normalized = _normalize_accession(genome_acc)
    fallback = (
        mapping_df
        .with_columns(pl.col('acc').map_elements(_normalize_accession, return_dtype=pl.Utf8).alias('acc_normalized'))
        .filter(pl.col('acc_normalized') == normalized)
        .select('parsed_taxon_string')
        .to_series()
        .to_list()
    )

    if len(fallback) == 1:
        logger.warning(
            f"No exact taxonomy row for '{genome_acc}'. Using version-insensitive match '{normalized}'."
        )
        return fallback[0]
    if len(fallback) > 1:
        logger.warning(
            f"Multiple version-insensitive taxonomy rows found for '{normalized}'. Using first match."
        )
        return fallback[0]

    # 3) Final non-fatal fallback to keep pipeline moving.
    logger.warning(
        f"No taxonomy mapping found for accession '{genome_acc}'. Using placeholder taxonomy."
    )
    return UNKNOWN_TAXONOMY

def assign_taxonomy_to_otu_table(otu_table_path, taxon_string, output_path):
    otu_df = pl.read_csv(otu_table_path, separator='\t')

    # Add taxonomy column
    otu_df = otu_df.with_columns(
        pl.lit(taxon_string).alias('taxonomy')
    )

    # Save updated OTU table
    otu_df.write_csv(output_path, separator='\t')
    logger.info(f"Saved taxonomy-assigned OTU table to {output_path}")

def main():
    logger.info("Starting eukaryotic taxonomy assignment")
    logger.info(f"OTU table: {OTU_TABLE}")
    logger.info(f"Taxonomy mapping file: {TAXONOMY_MAPPING}")
    logger.info(f"Output taxonomy OTU table: {TAXONOMY_OTU}")

    taxon_string = get_taxon_string(TAXONOMY_MAPPING, GENOME_ACC)
    logger.debug(f"Retrieved taxon string: {taxon_string}")

    assign_taxonomy_to_otu_table(OTU_TABLE, taxon_string, TAXONOMY_OTU)
    logger.info("Completed eukaryotic taxonomy assignment")

if __name__ == "__main__":
    main()