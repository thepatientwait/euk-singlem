#!/usr/bin/env python3

from pathlib import Path

import polars as pl

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = {
    "domain": "d__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}


def normalize_taxonomy(value: str | None) -> str:
    if value is None:
        return ""

    parsed = {}
    for chunk in str(value).split(";"):
        chunk = chunk.strip()
        for rank, prefix in PREFIXES.items():
            if chunk.startswith(prefix):
                name = chunk[len(prefix):].strip()
                if name:
                    parsed[rank] = name
                break

    if "domain" not in parsed:
        parsed["domain"] = "Eukaryota"

    return "; ".join(
        f"{PREFIXES[rank]}{parsed[rank]}"
        for rank in RANKS
        if rank in parsed
    )


def main() -> None:
    input_path = Path(str(snakemake.input.taxonomy_mapping))
    taxonomy_output = Path(str(snakemake.output.taxonomy))
    trusted_output = Path(str(snakemake.output.trusted_taxa))

    taxonomy_df = pl.read_csv(input_path, separator="\t")
    required = {"acc", "parsed_taxon_string"}
    missing = required.difference(taxonomy_df.columns)
    if missing:
        raise ValueError(f"Taxonomy mapping is missing columns: {', '.join(sorted(missing))}")

    normalized = (
        taxonomy_df.select(
            pl.col("acc").cast(pl.Utf8).str.strip_chars().alias("genome_id"),
            pl.col("parsed_taxon_string").map_elements(
                normalize_taxonomy,
                return_dtype=pl.Utf8,
            ).alias("taxonomy"),
        )
        .filter((pl.col("genome_id") != "") & (pl.col("taxonomy") != ""))
        .unique(subset=["genome_id"], keep="first")
    )

    taxonomy_output.parent.mkdir(parents=True, exist_ok=True)
    trusted_output.parent.mkdir(parents=True, exist_ok=True)
    normalized.select(["genome_id", "taxonomy"]).write_csv(
        taxonomy_output,
        separator="\t",
        include_header=False,
    )

    trusted_taxa = sorted(
        {
            chunk[len(prefix):].strip()
            for taxonomy in normalized["taxonomy"]
            for chunk in taxonomy.split(";")
            for prefix in PREFIXES.values()
            if chunk.strip().startswith(prefix) and chunk[len(prefix):].strip() != "Eukaryota"
        }
    )
    trusted_output.write_text("\n".join(trusted_taxa) + ("\n" if trusted_taxa else ""), encoding="utf-8")


if __name__ == "__main__":
    main()
