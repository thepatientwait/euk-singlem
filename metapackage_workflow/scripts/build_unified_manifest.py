#!/usr/bin/env python3
"""
Build unified download manifest from NCBI, P10K, and JGI metadata.

Consolidates metadata from three sources, resolves missing tax IDs, and outputs a
single manifest for the unified download pipeline.

Output columns: genome_acc, tax_id, download_source, group, original_name
"""
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import polars as pl
import os

from log_config import configure_logger
from ncbi_taxid_resolver import resolve_accessions_to_taxids, lookup_taxid

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# Snakemake parameters
JGI_METADATA_FILES = list(snakemake.input.jgi_csvs)
JGI_JSON_FILES = list(snakemake.input.jgi_jsons)
P10K_METADATA_FILE = snakemake.input['p10k_csv']
NCBI_METADATA_FILE = snakemake.input['ncbi_tsv']
OUTPUT_MANIFEST = snakemake.output['manifest']
NCBI_EMAIL = snakemake.params['ncbi_email']


def load_jgi_metadata(metadata_files: List[str], json_files: List[str]) -> pl.DataFrame:
    """Load JGI taxonomy metadata for accessions with downloadable JSON records."""
    records = []

    json_accessions_by_group = {}
    for json_file in json_files:
        group = Path(json_file).stem
        with open(json_file, encoding='utf-8') as handle:
            data = json.load(handle)
        json_accessions_by_group[group] = {
            file_record['portal_detail_id']
            for organism in data.get('organisms', [])
            for file_record in organism.get('files', [])
            if file_record.get('file_size', 0) == max(
                (item.get('file_size', 0) for item in organism.get('files', [])),
                default=0,
            )
            if file_record.get('portal_detail_id')
        }

    for meta_file in metadata_files:
        # Extract group from filename: e.g., "data/fungi_meta.csv" -> "fungi"
        group = Path(meta_file).stem.split('_')[0]

        df = pl.read_csv(meta_file, encoding='latin-1', ignore_errors=True)

        # Expected columns: portal, NCBI Taxon, name
        if 'portal' in df.columns and 'NCBI Taxon' in df.columns:
            json_accessions = json_accessions_by_group.get(group, set())
            csv_accessions = set(df['portal'].drop_nulls().to_list())
            csv_only = csv_accessions - json_accessions
            json_only = json_accessions - csv_accessions
            logger.info(
                "JGI %s metadata: CSV=%d JSON=%d shared=%d CSV-only=%d JSON-only=%d",
                group,
                len(csv_accessions),
                len(json_accessions),
                len(csv_accessions & json_accessions),
                len(csv_only),
                len(json_only),
            )
            if csv_only:
                logger.warning(
                    "JGI %s: excluding %d CSV accessions without JSON download records; sample=%s",
                    group,
                    len(csv_only),
                    sorted(csv_only)[:10],
                )
            if json_only:
                logger.warning(
                    "JGI %s: %d JSON accessions have no CSV taxonomy metadata; sample=%s",
                    group,
                    len(json_only),
                    sorted(json_only)[:10],
                )
            df = df.filter(pl.col('portal').is_in(json_accessions))
            df = df.select([
                pl.col('portal').alias('genome_acc'),
                pl.col('NCBI Taxon').cast(pl.String).alias('tax_id'),
                pl.col('name').alias('original_name'),
                pl.lit('JGI').alias('download_source'),
                pl.lit(group).alias('group')
            ])
            records.append(df)
        else:
            logger.warning(f"Unexpected columns in {meta_file}. Expected: portal, NCBI Taxon, name")
    
    if records:
        return pl.concat(records).unique()
    else:
        return pl.DataFrame(schema={
            'genome_acc': pl.Utf8,
            'tax_id': pl.String,
            'original_name': pl.Utf8,
            'download_source': pl.Utf8,
            'group': pl.Utf8
        })


def load_p10k_metadata(metadata_file: str) -> pl.DataFrame:
    """Load and parse P10K metadata CSV."""
    df = pl.read_csv(metadata_file)
    
    # Expected columns: Accession, Scientific name
    if 'Accession' in df.columns and 'Scientific name' in df.columns:
        df = df.select([
            pl.col('Accession').alias('genome_acc'),
            pl.lit(None, dtype=pl.String).alias('tax_id'),  # Will be filled via Entrez lookup
            pl.col('Scientific name').alias('original_name'),
            pl.col('Source').fill_null('NCBI').str.to_uppercase().alias('download_source'),
            pl.lit('p10k').alias('group')
        ])
        return df.unique()
    else:
        logger.error(f"Unexpected columns in {metadata_file}. Expected: Accession, Scientific name")
        return pl.DataFrame(schema={
            'genome_acc': pl.Utf8,
            'tax_id': pl.String,
            'original_name': pl.Utf8,
            'download_source': pl.Utf8,
            'group': pl.Utf8
        })


def load_ncbi_metadata(metadata_file: str) -> pl.DataFrame:
    """Load the accession list produced by discover_ncbi_genomes."""
    return pl.read_csv(metadata_file, separator='\t', schema_overrides={
        'genome_acc': pl.String,
        'tax_id': pl.String,
        'original_name': pl.String,
        'download_source': pl.String,
        'group': pl.String,
    }).select([
        'genome_acc', 'tax_id', 'original_name', 'download_source', 'group'
    ]).unique()


def resolve_p10k_tax_ids(manifest: pl.DataFrame) -> pl.DataFrame:
    """
    For NCBI-sourced records with null tax_id, resolve tax IDs from accessions
    via the Datasets API, with Entrez fallback for replaced/suppressed assemblies.
    """
    p10k_records = manifest.filter(pl.col('download_source') == 'NCBI')
    
    if p10k_records.height == 0:
        return manifest
    
    p10k_with_null = p10k_records.filter(pl.col('tax_id').is_null())
    
    if p10k_with_null.height == 0:
        logger.info("All P10K records already have tax IDs")
        return manifest
    
    logger.info(f"Looking up {p10k_with_null.height} P10K tax IDs via NCBI Datasets API...")
    
    accessions = p10k_with_null['genome_acc'].to_list()
    acc_to_taxid = resolve_accessions_to_taxids(accessions, NCBI_EMAIL)

    def _tax_id_for(acc):
        hit = lookup_taxid(acc_to_taxid, acc)
        return hit[0] if hit else None

    # Update manifest with looked-up tax IDs
    manifest = manifest.with_columns(
        pl.col('genome_acc').map_elements(
            _tax_id_for,
            return_dtype=pl.String
        ).fill_null(pl.col('tax_id')).alias('tax_id_resolved')
    )
    
    manifest = manifest.with_columns(
        pl.col('tax_id_resolved').alias('tax_id')
    ).drop('tax_id_resolved')
    
    still_null = manifest.filter(
        (pl.col('download_source') == 'NCBI') & pl.col('tax_id').is_null()
    ).height
    logger.info(f"NCBI records still missing tax IDs after lookup: {still_null}")
    
    return manifest


def main():
    logger.info("Building unified download manifest...")
    
    # Load metadata from all sources
    logger.info("Loading JGI metadata...")
    jgi_manifest = load_jgi_metadata(JGI_METADATA_FILES, JGI_JSON_FILES)
    logger.info(f"Loaded {jgi_manifest.height} JGI records")
    
    logger.info("Loading P10K metadata...")
    p10k_manifest = load_p10k_metadata(P10K_METADATA_FILE)
    logger.info(f"Loaded {p10k_manifest.height} P10K records")

    logger.info("Loading NCBI assembly metadata...")
    ncbi_manifest = load_ncbi_metadata(NCBI_METADATA_FILE)
    logger.info(f"Loaded {ncbi_manifest.height} NCBI records")

    p10k_accessions = p10k_manifest.get_column('genome_acc').unique().to_list()
    ncbi_before = ncbi_manifest.height
    ncbi_manifest = ncbi_manifest.filter(
        ~pl.col('genome_acc').is_in(p10k_accessions)
    )
    logger.info(f"Removed {ncbi_before - ncbi_manifest.height} overlapping NCBI records")
    
    # Combine manifests
    unified = pl.concat([jgi_manifest, p10k_manifest, ncbi_manifest]).unique()
    logger.info(f"Combined manifest: {unified.height} records")
    
    # Resolve P10K tax IDs (if needed)
    unified = resolve_p10k_tax_ids(unified)

    # Genomes without a tax ID cannot receive a taxonomy downstream; drop them
    # here so they never enter the download/BUSCO pipeline.
    before = unified.height
    unresolved_mask = (
        pl.col('tax_id').is_null()
        | (pl.col('tax_id').cast(pl.String).str.strip_chars() == '')
    )
    unresolved_sample = (
        unified.filter(unresolved_mask)['genome_acc'].to_list()[:10]
    )
    unified = unified.filter(~unresolved_mask)
    dropped = before - unified.height
    if dropped:
        logger.warning(
            "Dropped %d genomes with unresolved tax IDs from the unified manifest "
            "(examples: %s); they will not be downloaded or processed further.",
            dropped,
            ", ".join(unresolved_sample),
        )

    # Write output
    unified.write_csv(OUTPUT_MANIFEST, separator='\t')
    logger.info(f"Unified manifest written: {unified.height} records to {OUTPUT_MANIFEST}")


if __name__ == '__main__':
    main()
