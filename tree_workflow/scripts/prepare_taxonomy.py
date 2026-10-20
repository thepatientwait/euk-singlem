#!/usr/bin/env python3

from pathlib import Path
import re
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


def parse_ranks(value: str | None) -> dict:
    if not value:
        return {}

    parsed = {}
    for chunk in str(value).split(";"):
        chunk = chunk.strip()
        for rank, prefix in PREFIXES.items():
            if chunk.startswith(prefix):
                name = chunk[len(prefix):].strip()
                name = re.sub(r"\s*\([^)]*\)", "", name)
                name = name.replace("sp. ", "").replace("sp.", "species")
                name = name.replace(":", "_").strip()
                
                # Filter out naive NCBI placeholders and uncultured tags
                lower_name = name.lower()
                if "eukaryota" in lower_name or "uncultured" in lower_name:
                    name = ""
                    
                if name:
                    parsed[rank] = name
                break

    # Force domain to Eukarya
    parsed["domain"] = "Eukarya"
    return parsed


def format_aggressive(parsed: dict) -> str:
    # Reject if ANY of the 7 required ranks are missing or empty
    if not all(parsed.get(rank) for rank in RANKS):
        return ""
    
    return ";".join(f"{PREFIXES[rank]}{parsed[rank]}" for rank in RANKS)


def format_conservative(parsed: dict) -> str:
    if not parsed:
        return ""
        
    # Join all ranks; missing/filtered ranks are left blank (e.g., "p__")
    return ";".join(f"{PREFIXES[rank]}{parsed.get(rank, '')}" for rank in RANKS)


def main() -> None:
    input_path = Path(str(snakemake.input.taxonomy_mapping))
    taxonomy_output = Path(str(snakemake.output.taxonomy))
    trusted_output = Path(str(snakemake.output.trusted_taxa))
    conservative_output = Path(str(snakemake.output.conservative_taxonomy))

    taxonomy_df = pl.read_csv(input_path, separator="\t")
    required = {"acc", "parsed_taxon_string"}
    missing = required.difference(taxonomy_df.columns)
    if missing:
        raise ValueError(f"Taxonomy mapping is missing columns: {', '.join(sorted(missing))}")

    # 1. Parse all taxonomy strings into dictionary objects
    parsed_df = (
        taxonomy_df.select(
            pl.col("acc").cast(pl.Utf8).fill_null("").str.strip_chars().alias("genome_id"),
            pl.col("parsed_taxon_string").map_elements(
                parse_ranks, 
                return_dtype=pl.Object
            ).alias("parsed_dict")
        )
        .filter(pl.col("genome_id") != "")
        .unique(subset=["genome_id"], keep="first")
    )

    # 2. Map the dictionaries to the aggressive and conservative string formats
    formatted_df = parsed_df.with_columns(
        pl.col("parsed_dict").map_elements(format_aggressive, return_dtype=pl.Utf8).alias("aggressive_tax"),
        pl.col("parsed_dict").map_elements(format_conservative, return_dtype=pl.Utf8).alias("conservative_tax")
    )

    # 3. Filter and select the final outputs
    aggressive_df = formatted_df.filter(pl.col("aggressive_tax") != "").select(
        ["genome_id", pl.col("aggressive_tax").alias("taxonomy")]
    )
    
    conservative_df = formatted_df.filter(pl.col("conservative_tax") != "").select(
        ["genome_id", pl.col("conservative_tax").alias("taxonomy")]
    )

    # Ensure output directories exist
    taxonomy_output.parent.mkdir(parents=True, exist_ok=True)
    trusted_output.parent.mkdir(parents=True, exist_ok=True)
    conservative_output.parent.mkdir(parents=True, exist_ok=True)

    # Write TSVs
    aggressive_df.write_csv(taxonomy_output, separator="\t", include_header=False)
    conservative_df.write_csv(conservative_output, separator="\t", include_header=False)

    # Write trusted taxa list based solely on the perfect aggressive output
    trusted_taxa = list(
        dict.fromkeys(
            chunk.strip()
            for taxonomy in aggressive_df["taxonomy"]
            for chunk in taxonomy.split(";")
        )
    )
    trusted_output.write_text("\n".join(trusted_taxa) + ("\n" if trusted_taxa else ""), encoding="utf-8")


if __name__ == "__main__":
    main()