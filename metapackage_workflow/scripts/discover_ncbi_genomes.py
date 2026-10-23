#!/usr/bin/env python3
"""Discover current NCBI assemblies below the configured taxonomy ID."""
import os
from typing import Iterable

import polars as pl
from Bio import Entrez

from log_config import configure_logger

logger = configure_logger(os.path.splitext(os.path.basename(__file__))[0])
P10K_METADATA = str(snakemake.input.p10k_metadata)
FUNGI_METADATA = str(snakemake.input.fungi_metadata)
ALGAE_METADATA = str(snakemake.input.algae_metadata)
TAXON_ID = str(snakemake.params.tax_id)
EMAIL = snakemake.params.email
REFERENCE_ONLY = bool(snakemake.params.reference_only)
INCLUDE_METAZOA_VIRIDIPLANTAE_REFS = bool(
    getattr(snakemake.params, "include_metazoa_viridiplantae_refs", False))
OUTPUT = snakemake.output.metadata
BATCH_SIZE = 200
SUBCLADE_EXCLUSIONS = " NOT txid33208[Organism:exp] NOT txid3193[Organism:exp]"


def chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start:start + size]


def search(term: str) -> list[dict]:
    logger.info("Searching NCBI assembly database: %s", term)
    with Entrez.esearch(db="assembly", term=term, retmax=100000, usehistory="y") as handle:
        search_result = Entrez.read(handle)

    ids = [str(value) for value in search_result.get("IdList", [])]
    logger.info("Found %d current assemblies", len(ids))
    rows = []
    for batch in chunks(ids, BATCH_SIZE):
        with Entrez.esummary(db="assembly", id=",".join(batch), report="full") as handle:
            summary = Entrez.read(handle)
        documents = summary.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        for document in documents:
            accession = str(document.get("AssemblyAccession", "")).strip()
            tax_id = str(document.get("Taxid", "")).strip()
            if accession and tax_id:
                rows.append({
                    "genome_acc": accession,
                    "tax_id": tax_id,
                    "original_name": str(document.get("SpeciesName", "")).strip(),
                    "download_source": "NCBI",
                    "group": "ncbi",
                })
    return rows


def main() -> None:
    Entrez.email = EMAIL
    base = f"txid{TAXON_ID}[Organism:exp] AND latest[filter]"

    rows = []
    ref_term = base
    if not INCLUDE_METAZOA_VIRIDIPLANTAE_REFS:
        ref_term += SUBCLADE_EXCLUSIONS
    ref_term += " AND reference[Properties]"
    rows.extend(search(ref_term))

    if not REFERENCE_ONLY:
        rows.extend(search(base + SUBCLADE_EXCLUSIONS + " NOT reference[Properties]"))

    if not rows:
        raise RuntimeError("NCBI discovery returned no assembly records")

    ncbi_df = pl.DataFrame(rows).unique("genome_acc").sort("genome_acc")

    p10k_df = pl.read_csv(P10K_METADATA).select(pl.col('Accession').alias('genome_acc'))
    
    fungi_df = (
        pl
        .read_csv(FUNGI_METADATA, encoding="latin-1", ignore_errors=True)
        .select(
            pl.col('NCBI Taxon')
            .cast(pl.String)
            .alias('tax_id')))

    algae_df = (
        pl
        .read_csv(ALGAE_METADATA, encoding="latin-1", ignore_errors=True)
        .select(
            pl.col('NCBI Taxon')
            .cast(pl.String)
            .alias('tax_id')))

    ncbi_df = (
        ncbi_df
        .join(p10k_df, on="genome_acc", how="anti")
        .join(fungi_df, on="tax_id", how="anti")
        .join(algae_df, on="tax_id", how="anti"))

    ncbi_df.write_csv(OUTPUT, separator="\t")

    logger.info("Wrote %d NCBI assemblies to %s", len(ncbi_df), OUTPUT)


if __name__ == "__main__":
    main()
