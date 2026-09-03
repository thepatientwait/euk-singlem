#!/usr/bin/env python3

from pathlib import Path
from typing import Iterable

import polars as pl
from Bio import Entrez


def parse_taxonomy(record: dict) -> str:
    taxonomy = {
        taxon["Rank"]: taxon["ScientificName"]
        for taxon in record.get("LineageEx", [])
        if taxon.get("Rank") in {"superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"}
    }

    species = record.get("ScientificName", "")
    taxonomy["species"] = species

    # Keep GTDB-style prefixes used by the existing pipeline.
    ranks_order = [
        ("domain", "superkingdom"),
        ("kingdom", "kingdom"),
        ("phylum", "phylum"),
        ("class", "class"),
        ("order", "order"),
        ("family", "family"),
        ("genus", "genus"),
        ("species", "species"),
    ]

    values = []
    for short, rank in ranks_order:
        value = taxonomy.get(rank, "")
        prefix = short[0]
        values.append(f"{prefix}__{value}")

    return "; ".join(values)


def chunked(values: Iterable[str], chunk_size: int) -> Iterable[list[str]]:
    chunk = []
    for value in values:
        chunk.append(value)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def fetch_taxonomy_strings(ncbi_ids: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}

    for chunk in chunked(ncbi_ids, 200):
        handle = Entrez.efetch(db="Taxonomy", id=",".join(chunk), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        for record in records:
            tax_id = str(record.get("TaxId", "")).strip()
            if not tax_id:
                continue
            mapping[tax_id] = parse_taxonomy(record)

    return mapping


def main() -> None:
    base_taxonomy_file = Path(str(snakemake.input.base_taxonomy))
    external_manifest_table = Path(str(snakemake.input.external_manifest_table))
    external_taxid_map_file = Path(str(snakemake.input.external_taxid_map))
    eukcc_summary_file = Path(str(snakemake.input.eukcc_summary))

    output_file = Path(str(snakemake.output.taxonomy))
    output_file.parent.mkdir(parents=True, exist_ok=True)

    strict_requery = bool(snakemake.params.strict_requery)

    Entrez.email = str(snakemake.params.email)

    base_df = pl.read_csv(base_taxonomy_file, separator="\t")
    external_manifest_df = pl.read_csv(external_manifest_table, separator="\t")

    if external_taxid_map_file.is_file() and external_taxid_map_file.stat().st_size > 0:
        external_taxids_df = pl.read_csv(external_taxid_map_file, separator="\t")
    else:
        external_taxids_df = pl.DataFrame(schema={"acc": pl.Utf8, "ncbi_id": pl.Utf8})

    if eukcc_summary_file.exists() and eukcc_summary_file.stat().st_size > 0:
        eukcc_df = pl.read_csv(eukcc_summary_file, separator="\t")
    else:
        eukcc_df = pl.DataFrame(schema={"acc": pl.Utf8, "eukcc_taxonomy": pl.Utf8})

    base_cols = [col for col in ["acc", "ncbi_id", "parsed_taxon_string"] if col in base_df.columns]
    base_df = base_df.select(base_cols)

    all_accessions = pl.concat(
        [
            base_df.select("acc"),
            external_manifest_df.select("genome_acc").rename({"genome_acc": "acc"}),
        ]
    ).unique()

    id_map_df = pl.concat(
        [
            base_df.select(["acc", "ncbi_id"]).filter(pl.col("ncbi_id").is_not_null()),
            external_taxids_df.select(["acc", "ncbi_id"]).filter(pl.col("ncbi_id").is_not_null()),
        ],
        how="vertical_relaxed",
    ).unique(subset=["acc"], keep="first")

    acc_df = all_accessions.join(id_map_df, on="acc", how="left")

    ncbi_ids = sorted({str(x) for x in acc_df.get_column("ncbi_id").drop_nulls().to_list() if str(x).strip()})
    requery_map = fetch_taxonomy_strings(ncbi_ids)

    taxonomy_rows = []
    missing_requery = []

    eukcc_tax_map = {}
    if "acc" in eukcc_df.columns and "eukcc_taxonomy" in eukcc_df.columns:
        eukcc_tax_map = dict(zip(eukcc_df["acc"].to_list(), eukcc_df["eukcc_taxonomy"].to_list()))

    base_tax_map = {}
    if "acc" in base_df.columns and "parsed_taxon_string" in base_df.columns:
        base_tax_map = dict(zip(base_df["acc"].to_list(), base_df["parsed_taxon_string"].to_list()))

    for row in acc_df.iter_rows(named=True):
        acc = row["acc"]
        ncbi_id = row.get("ncbi_id")
        parsed_taxon_string = None
        source = "missing"

        if ncbi_id is not None:
            parsed_taxon_string = requery_map.get(str(ncbi_id))
            if parsed_taxon_string:
                source = "ncbi_requery"
            else:
                missing_requery.append((acc, str(ncbi_id)))

        if not parsed_taxon_string:
            parsed_taxon_string = base_tax_map.get(acc)
            if parsed_taxon_string:
                source = "fallback_existing"

        if not parsed_taxon_string:
            parsed_taxon_string = eukcc_tax_map.get(acc)
            if parsed_taxon_string:
                source = "fallback_eukcc"

        taxonomy_rows.append(
            {
                "acc": acc,
                "ncbi_id": ncbi_id,
                "parsed_taxon_string": parsed_taxon_string,
                "taxonomy_source": source,
            }
        )

    if strict_requery and missing_requery:
        sample = ", ".join([f"{acc}:{taxid}" for acc, taxid in missing_requery[:10]])
        raise ValueError(
            f"Strict NCBI requery failed for {len(missing_requery)} accessions. Examples: {sample}"
        )

    pl.DataFrame(taxonomy_rows).write_csv(output_file, separator="\t")


if __name__ == "__main__":
    main()
