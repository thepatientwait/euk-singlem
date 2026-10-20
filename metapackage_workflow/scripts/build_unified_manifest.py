#!/usr/bin/env python3
"""
Build unified download manifest from NCBI, P10K, and JGI metadata.

Consolidates metadata from three sources, queries NCBI for lineage info,
filters out animals (Metazoa) and plants (Viridiplantae), and outputs a
single manifest for the unified download pipeline.

Output columns: genome_acc, tax_id, download_source, group, original_name
"""
import json
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import polars as pl
from Bio import Entrez
import os
import re

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


def query_ncbi_datasets_lineage(tax_ids: List[int]) -> Dict[int, str]:
    """
    Query NCBI Datasets CLI for lineage info and return dict of {taxid: lineage_string}.
    Lineage format: d__domain;p__phylum;c__class;o__order;f__family;g__genus;s__species
    """
    lineage_map = {}
    
    if not tax_ids:
        return lineage_map
    
    # Remove duplicates
    tax_ids = list(set(tax_ids))
    
    try:
        # Query NCBI Datasets for organism taxonomy info
        # We'll use datasets summary taxonomy taxon command
        tax_id_str = ','.join(str(tid) for tid in tax_ids)
        
        logger.info(f"Querying NCBI Datasets for {len(tax_ids)} tax IDs...")
        
        cmd = ['datasets', 'summary', 'taxonomy', 'taxon', tax_id_str, '--as-json-lines']
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            try:
                record = json.loads(line)
                tax_record = record.get('taxonomy', {})
                tax_id = int(tax_record.get('tax_id', 0))
                
                if not tax_id:
                    continue
                
                # Datasets emits rank names in taxonomy.classification.
                classification = tax_record.get('classification', {})
                lineage_parts = {}
                
                for rank, item in classification.items():
                    rank = rank.lower()
                    name = item.get('name', '')
                    
                    if rank in ['superkingdom', 'domain']:
                        lineage_parts['domain'] = name
                    elif rank == 'kingdom':
                        lineage_parts['kingdom'] = name
                    elif rank == 'phylum':
                        lineage_parts['phylum'] = name
                    elif rank == 'class':
                        lineage_parts['class'] = name
                    elif rank == 'order':
                        lineage_parts['order'] = name
                    elif rank == 'family':
                        lineage_parts['family'] = name
                    elif rank == 'genus':
                        lineage_parts['genus'] = name
                    elif rank == 'species':
                        lineage_parts['species'] = name
                
                # Format as GTDB-like lineage
                gtdb_lineage = ';'.join([
                    f"d__{lineage_parts.get('domain', '')}",
                    f"k__{lineage_parts.get('kingdom', '')}",
                    f"p__{lineage_parts.get('phylum', '')}",
                    f"c__{lineage_parts.get('class', '')}",
                    f"o__{lineage_parts.get('order', '')}",
                    f"f__{lineage_parts.get('family', '')}",
                    f"g__{lineage_parts.get('genus', '')}",
                    f"s__{lineage_parts.get('species', '')}"
                ])
                
                lineage_map[tax_id] = gtdb_lineage
            except (json.JSONDecodeError, ValueError, KeyError) as e:
                logger.debug(f"Error parsing taxonomy record: {e}")
                continue
        
        logger.info(f"Successfully queried lineage for {len(lineage_map)} tax IDs")
    
    except subprocess.CalledProcessError as e:
        logger.error(f"NCBI Datasets query failed: {e.stderr}")
        # Fallback: try Entrez API
        logger.info("Falling back to NCBI Entrez API...")
        lineage_map = query_ncbi_entrez_lineage(tax_ids, NCBI_EMAIL)
    
    return lineage_map


def query_ncbi_entrez_lineage(tax_ids: List[int], email: str) -> Dict[int, str]:
    """Fallback: Query NCBI Entrez API for lineage info."""
    lineage_map = {}
    
    Entrez.email = email
    
    try:
        # Batch query
        tax_id_str = ','.join(str(tid) for tid in tax_ids)
        handle = Entrez.efetch(db="Taxonomy", id=tax_id_str, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        for record in records:
            tax_id = int(record.get('TaxId', 0))
            if not tax_id:
                continue
            
            lineage_parts = {}
            lineage_ex = record.get('LineageEx', [])
            
            for item in lineage_ex:
                rank = item.get('Rank', '').lower()
                name = item.get('ScientificName', '')
                
                if rank in ['superkingdom', 'domain']:
                    lineage_parts['domain'] = name
                elif rank == 'phylum':
                    lineage_parts['phylum'] = name
                elif rank == 'class':
                    lineage_parts['class'] = name
                elif rank == 'order':
                    lineage_parts['order'] = name
                elif rank == 'family':
                    lineage_parts['family'] = name
                elif rank == 'genus':
                    lineage_parts['genus'] = name
            
            species = record.get('ScientificName', '')
            
            gtdb_lineage = ';'.join([
                f"d__{lineage_parts.get('domain', '')}",
                f"p__{lineage_parts.get('phylum', '')}",
                f"c__{lineage_parts.get('class', '')}",
                f"o__{lineage_parts.get('order', '')}",
                f"f__{lineage_parts.get('family', '')}",
                f"g__{lineage_parts.get('genus', '')}",
                f"s__{species}"
            ])
            
            lineage_map[tax_id] = gtdb_lineage
        
        logger.info(f"Successfully queried Entrez lineage for {len(lineage_map)} tax IDs")
    
    except Exception as e:
        logger.error(f"NCBI Entrez query failed: {e}")
    
    return lineage_map


def filter_excluded_clades(manifest: pl.DataFrame, lineage_map: Dict[int, str]) -> pl.DataFrame:
    """
    Filter out records where lineage contains Metazoa (animals) or Viridiplantae (plants).
    Also drops records with missing lineage info.
    """
    def is_excluded(tax_id):
        lineage = lineage_map.get(tax_id, '')
        if not lineage:
            return True  # Exclude records without lineage info
        
        # Check for excluded clades
        excluded_patterns = [
            r'k__Metazoa',
            r'k__Viridiplantae',
            r'p__Chordata',
            r'p__Streptophyta',
        ]
        
        for pattern in excluded_patterns:
            if re.search(pattern, lineage, re.IGNORECASE):
                return True
        
        return False
    
    # Filter out excluded records
    before_count = manifest.height
    
    manifest = manifest.with_columns(
        pl.col('tax_id').map_elements(is_excluded, return_dtype=pl.Boolean).alias('_is_excluded')
    ).filter(~pl.col('_is_excluded')).drop('_is_excluded')
    
    after_count = manifest.height
    removed_count = before_count - after_count
    
    logger.info(f"Filtered manifest: {before_count} → {after_count} records ({removed_count} removed)")
    
    return manifest


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
