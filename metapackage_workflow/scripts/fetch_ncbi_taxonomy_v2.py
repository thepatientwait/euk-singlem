#!/usr/bin/env python3
"""
Fetch NCBI taxonomy for all genomes in the unified manifest.

Reads unified_download_manifest.tsv, extracts unique tax IDs, queries NCBI Datasets CLI
for lineage information, and outputs taxonomy_mapping.tsv with GTDB-format lineage.
"""
import json
import subprocess
from pathlib import Path
from typing import Dict, List
import polars as pl
from Bio import Entrez
import os

from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# Snakemake parameters
UNIFIED_MANIFEST = snakemake.input['unified_manifest']
TAXONOMY_OUTPUT = snakemake.output['taxonomy']
NCBI_EMAIL = snakemake.params['email']
DATASETS_BATCH_SIZE = 200


def load_unified_manifest(manifest_file: str) -> pl.DataFrame:
    """Load unified manifest and extract tax IDs."""
    df = pl.read_csv(manifest_file, separator='\t')
    
    # Expected columns: genome_acc, tax_id, download_source, group, original_name
    required_cols = {'genome_acc', 'tax_id'}
    if not required_cols.issubset(set(df.columns)):
        raise ValueError(f"Manifest missing required columns. Expected: {required_cols}, Got: {set(df.columns)}")
    
    logger.info(f"Loaded manifest with {df.height} records")
    return df


def query_ncbi_datasets_taxonomy(tax_ids: List[int]) -> Dict[int, Dict]:
    """
    Query NCBI Datasets CLI for taxonomy records.
    
    Returns dict mapping tax_id -> {lineage_dict}
    """
    taxonomy_map = {}
    
    if not tax_ids:
        return taxonomy_map
    
    # Remove duplicates and None values
    tax_ids = [tid for tid in set(tax_ids) if tid]
    
    if not tax_ids:
        logger.warning("No valid tax IDs to query")
        return taxonomy_map
    
    try:
        logger.info(
            "Querying NCBI Datasets for %d tax IDs in batches of %d...",
            len(tax_ids),
            DATASETS_BATCH_SIZE,
        )
        lineage_count = 0
        batches = [
            tax_ids[start:start + DATASETS_BATCH_SIZE]
            for start in range(0, len(tax_ids), DATASETS_BATCH_SIZE)
        ]

        for batch_number, tax_id_batch in enumerate(batches, start=1):
            tax_id_str = ','.join(str(tid) for tid in tax_id_batch)
            logger.debug(
                "NCBI Datasets batch %d/%d: %d IDs, argument length=%d, first IDs=%s, last IDs=%s",
                batch_number,
                len(batches),
                len(tax_id_batch),
                len(tax_id_str),
                tax_id_batch[:10],
                tax_id_batch[-10:],
            )

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

                    lineage_list = tax_record.get('classification', {}).items()
                    lineage_dict = {}

                    for rank, item in lineage_list:
                        rank = rank.lower()
                        name = item.get('name', '')

                        rank_map = {
                            'superkingdom': 'domain',
                            'domain': 'domain',
                            'phylum': 'phylum',
                            'class': 'class',
                            'order': 'order',
                            'family': 'family',
                            'genus': 'genus',
                            'species': 'species'
                        }

                        if rank in rank_map:
                            lineage_dict[rank_map[rank]] = name

                    taxonomy_map[tax_id] = lineage_dict
                    secondary_tax_ids = tax_record.get('secondary_tax_ids', [])
                    for secondary_tax_id in secondary_tax_ids:
                        taxonomy_map[int(secondary_tax_id)] = lineage_dict
                    lineage_count += 1

                except (json.JSONDecodeError, ValueError, KeyError) as e:
                    logger.debug(f"Error parsing taxonomy record: {e}")
                    continue

        logger.info(f"Successfully queried {lineage_count} lineages from NCBI Datasets")

    except subprocess.CalledProcessError as e:
        logger.error(
            "NCBI Datasets query failed (return code %s). stderr=%r stdout=%r",
            e.returncode,
            (e.stderr or '')[:2000],
            (e.stdout or '')[:2000],
        )
        logger.info("Falling back to NCBI Entrez API...")
        taxonomy_map = query_ncbi_entrez_taxonomy(tax_ids, NCBI_EMAIL)
    
    return taxonomy_map


def query_ncbi_entrez_taxonomy(tax_ids: List[int], email: str) -> Dict[int, Dict]:
    """
    Fallback: Query NCBI Entrez API for taxonomy records.
    """
    taxonomy_map = {}
    Entrez.email = email
    
    try:
        tax_id_str = ','.join(str(tid) for tid in tax_ids)
        logger.info(f"Querying NCBI Entrez for {len(tax_ids)} tax IDs...")
        
        handle = Entrez.efetch(db="Taxonomy", id=tax_id_str, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        returned_tax_ids = {
            int(record.get('TaxId', 0))
            for record in records
            if record.get('TaxId')
        }
        missing_tax_ids = sorted(set(tax_ids) - returned_tax_ids)
        logger.info(
            "NCBI Entrez returned %d records for %d requested tax IDs",
            len(records),
            len(set(tax_ids)),
        )
        if missing_tax_ids:
            logger.warning(
                "NCBI Entrez returned no records for %d tax IDs; sample=%s",
                len(missing_tax_ids),
                missing_tax_ids[:25],
            )
        
        lineage_count = 0
        for record in records:
            tax_id = int(record.get('TaxId', 0))
            if not tax_id:
                continue
            
            lineage_dict = {}
            lineage_ex = record.get('LineageEx', [])
            
            for item in lineage_ex:
                rank = item.get('Rank', '').lower()
                name = item.get('ScientificName', '')
                
                rank_map = {
                    'superkingdom': 'domain',
                    'domain': 'domain',
                    'phylum': 'phylum',
                    'class': 'class',
                    'order': 'order',
                    'family': 'family',
                    'genus': 'genus',
                    'species': 'species'
                }
                
                if rank in rank_map:
                    lineage_dict[rank_map[rank]] = name
            
            # Add current record's species name if not already present
            if 'species' not in lineage_dict:
                lineage_dict['species'] = record.get('ScientificName', '')
            
            taxonomy_map[tax_id] = lineage_dict
            for secondary_tax_id in record.get('AkaTaxIds', []):
                taxonomy_map[int(secondary_tax_id)] = lineage_dict
            lineage_count += 1
        
        logger.info(f"Successfully queried {lineage_count} lineages from NCBI Entrez")
    
    except Exception as e:
        logger.error(f"NCBI Entrez query failed: {e}")
    
    return taxonomy_map


def lineage_dict_to_gtdb_string(lineage_dict: Dict) -> str:
    """Convert lineage dictionary to GTDB-format string."""
    ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    
    parts = []
    for rank, prefix in zip(ranks, prefixes):
        value = lineage_dict.get(rank, '')
        parts.append(f"{prefix}{value}")
    
    return ';'.join(parts)


def build_taxonomy_output(manifest: pl.DataFrame, taxonomy_map: Dict[int, Dict]) -> pl.DataFrame:
    """
    Join manifest with taxonomy data and output in expected format:
    acc, ncbi_id, parsed_taxon_string
    """
    output_records = []
    
    for row in manifest.iter_rows(named=True):
        genome_acc = row['genome_acc']
        tax_id = row['tax_id']
        
        if not tax_id:
            logger.debug(f"Skipping {genome_acc}: no tax ID")
            continue
        
        if tax_id not in taxonomy_map:
            logger.warning(f"No taxonomy found for tax_id {tax_id} (genome: {genome_acc})")
            continue
        
        lineage_dict = taxonomy_map[tax_id]
        parsed_taxon_string = lineage_dict_to_gtdb_string(lineage_dict)
        
        output_records.append({
            'acc': genome_acc,
            'ncbi_id': tax_id,
            'parsed_taxon_string': parsed_taxon_string
        })
    
    return pl.DataFrame(output_records)


def main():
    logger.info("Fetching NCBI taxonomy for unified manifest...")
    
    # Load unified manifest
    manifest = load_unified_manifest(UNIFIED_MANIFEST)
    
    # Extract unique tax IDs
    tax_ids = manifest['tax_id'].unique().drop_nulls().to_list()
    logger.info(f"Found {len(tax_ids)} unique tax IDs")
    
    # Query NCBI for taxonomy
    taxonomy_map = query_ncbi_datasets_taxonomy(tax_ids)
    
    if not taxonomy_map:
        raise ValueError("No taxonomy records retrieved; cannot proceed")
    
    # Build output
    taxonomy_df = build_taxonomy_output(manifest, taxonomy_map)
    
    if taxonomy_df.height == 0:
        raise ValueError("No taxonomy records matched manifest; cannot proceed")
    
    # Write output
    taxonomy_df.write_csv(TAXONOMY_OUTPUT, separator='\t')
    logger.info(f"Taxonomy output: {taxonomy_df.height} records written to {TAXONOMY_OUTPUT}")


if __name__ == '__main__':
    main()
