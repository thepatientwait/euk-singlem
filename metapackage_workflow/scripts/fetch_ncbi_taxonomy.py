# init
# ----
import polars as pl
from Bio import Entrez
import os

import regex as re

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
JGI_METADATA_FILES      = snakemake.input['jgi_metadata']
P10K_METADATA_FILE      = snakemake.input['p10k_metadata']
TAXONOMY_OUTPUT         = snakemake.output['taxonomy']

NCBI_EMAIL          = snakemake.params['email']


# functions
# ---------
def get_jgi_taxon_id_map(metadata_files):
    """Load and extract unique NCBI taxon IDs from metadata files"""
    taxon_id_map_list = []

    for meta_file in metadata_files:
        taxon_id_map_list.append(
            pl.read_csv(meta_file, encoding='latin-1', columns=['portal', 'NCBI Taxon', 'name'])
            .select([
                pl.col('portal').alias('acc'),
                pl.col('NCBI Taxon').alias('ncbi_id').cast(pl.String),
                pl.col('name').alias('strain_string')
            ])
        )
    taxon_id_map = pl.concat(taxon_id_map_list).unique()

    logger.info(f"Found {taxon_id_map.height} unique JGI taxon IDs.")
    
    return taxon_id_map


def ncbi_accs_to_taxids(accessions, chunk_size=100):
    """Batch convert assembly accession -> TaxID."""
    acc_to_taxid = {}

    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i:i + chunk_size]

        # Build a single assembly search query for this chunk
        term = " OR ".join([f'"{acc}"[Assembly Accession]' for acc in chunk])

        search_h = Entrez.esearch(db="assembly", term=term, retmax=len(chunk) * 5)
        search_r = Entrez.read(search_h)
        search_h.close()

        uids = search_r.get("IdList", [])
        if not uids:
            continue

        summary_h = Entrez.esummary(db="assembly", id=",".join(uids))
        summary_r = Entrez.read(summary_h)
        summary_h.close()

        docs = summary_r.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        for doc in docs:
            acc = str(doc.get("AssemblyAccession", "")).strip()
            taxid = str(doc.get("Taxid", "")).strip()
            if acc and taxid:
                acc_to_taxid[acc] = taxid

    return acc_to_taxid


def get_p10k_taxon_id_map(metadata_file):
    """Load and extract unique NCBI taxon IDs from P10K metadata file"""
    taxon_id_map = (
        pl.read_csv(metadata_file, columns=['Accession', 'Scientific name'])
        .select([
            pl.col('Accession').alias('acc').cast(pl.String),
            pl.col('Scientific name').alias('strain_string')
        ])
        .unique()
    )

    acc_list = taxon_id_map['acc'].to_list()
    acc_to_taxid = ncbi_accs_to_taxids(acc_list)

    taxon_id_map = (
        taxon_id_map
        .with_columns(
            pl.col('acc').replace(acc_to_taxid).cast(pl.String).alias('ncbi_id'))
        .select(['acc', 'ncbi_id', 'strain_string'])
        .unique()
    )

    logger.info(f"Found {taxon_id_map.height} unique P10K taxon IDs.")
    
    return taxon_id_map


def fetch_ncbi_taxonomy(taxon_id_map):
    """Fetch taxonomy for a batch of NCBI IDs"""
    ncbi_ids = taxon_id_map['ncbi_id'].unique().to_list()
    logger.info(f"Querying Entrez for {len(ncbi_ids)} unique taxon IDs...")
    
    try:
        handle = Entrez.efetch(db="Taxonomy", id=','.join(ncbi_ids), retmode="xml")
        ncbi_records = Entrez.read(handle)
        handle.close()
    except Exception as e:
        logger.error(f"Error fetching batch: {e}")
    
    return ncbi_records


def parse_taxonomy(record, jgi_strain_string):
    """Parse a single NCBI taxonomy record into GTDB-like format"""
    # Build taxonomy dict
    taxonomy = {
        taxon['Rank']: taxon['ScientificName'] 
        for taxon in record['LineageEx'] 
        if taxon['Rank'] in {'domain', 'phylum', 'class', 'order', 'family', 'genus'}
    }

    # remove version suffix from species name
    pattern = r'\sv[0-9.]+$'

    if jgi_strain_string:
        taxonomy['species'] = re.sub(pattern, '', jgi_strain_string)
        taxonomy['genus'] = jgi_strain_string.split(' ')[0]
    else:
        taxonomy['species'] = record['ScientificName']
    
    # Create parsed string in GTDB-like format
    ranks_order = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    parsed_taxon_string = '; '.join([
        f"{rank[0]}__{next((taxonomy.get(ranks_order[j], '') for j in range(i, -1, -1) if taxonomy.get(ranks_order[j])), '')}" 
        for i, rank in enumerate(ranks_order)
    ])
    
    return parsed_taxon_string


def extract_candidate_taxids(record):
    """Return all tax IDs that can map this record back to input metadata."""
    candidate_ids = set()

    primary_taxid = str(record.get('TaxId', '')).strip()
    if primary_taxid:
        candidate_ids.add(primary_taxid)

    for field in ('AkaTaxIds', 'MergedTaxIds'):
        values = record.get(field)
        if not values:
            continue

        if isinstance(values, list):
            for value in values:
                taxid = str(value).strip()
                if taxid:
                    candidate_ids.add(taxid)
        else:
            taxid = str(values).strip()
            if taxid:
                candidate_ids.add(taxid)

    return candidate_ids


def parse_ncbi_records(ncbi_records, taxon_id_map):
    entries_by_ncbi_id = {}
    for entry in taxon_id_map.iter_rows(named=True):
        key = str(entry['ncbi_id'])
        entries_by_ncbi_id.setdefault(key, []).append(entry)

    for record in ncbi_records:
        candidate_ids = extract_candidate_taxids(record)
        matched_entries = []
        seen_accessions = set()

        for candidate_id in candidate_ids:
            for entry in entries_by_ncbi_id.get(candidate_id, []):
                acc = entry['acc']
                if acc in seen_accessions:
                    continue
                seen_accessions.add(acc)
                matched_entries.append(entry)

        if not matched_entries:
            logger.debug(
                "No matching input entries for returned TaxId %s (aliases: %s)",
                record.get('TaxId'),
                sorted(candidate_ids),
            )
            continue

        for jgi_entry in matched_entries:
            acc                 = jgi_entry['acc']
            strain_string       = jgi_entry['strain_string']
            requested_ncbi_id   = str(jgi_entry['ncbi_id'])
            
            parsed_taxon_string = parse_taxonomy(record, strain_string)
            
            yield {
                'acc': acc,
                'ncbi_id': requested_ncbi_id,
                'parsed_taxon_string': parsed_taxon_string
            }


# main
# ----
def main():
    Entrez.email = NCBI_EMAIL
    
    logger.info("Loading NCBI taxon IDs from JGI and P10K metadata files...")
    taxon_id_map = get_jgi_taxon_id_map(JGI_METADATA_FILES)
    taxon_id_map = taxon_id_map.vstack(get_p10k_taxon_id_map(P10K_METADATA_FILE))

    logger.info("Fetching NCBI taxonomy records...")
    ncbi_records = fetch_ncbi_taxonomy(taxon_id_map)

    logger.info("Parsing NCBI taxonomy records...")
    all_parsed_records = list(parse_ncbi_records(ncbi_records, taxon_id_map))
    
    # Write output
    taxonomy_df = pl.DataFrame(all_parsed_records).unique()
    taxonomy_df.write_csv(TAXONOMY_OUTPUT, separator='\t')
    logger.info(f"Taxonomy mapping complete: {taxonomy_df.height} records written to {TAXONOMY_OUTPUT}.")


if __name__ == '__main__':
    main()
