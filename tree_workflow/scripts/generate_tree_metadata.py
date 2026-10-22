#!/usr/bin/env python3

from pathlib import Path

import polars as pl

from log_redirect import redirect_to_snakemake_log

# send this rule's stdout/stderr to its log: file
redirect_to_snakemake_log()

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


def expand_taxonomy(frame: pl.DataFrame) -> pl.DataFrame:
    expressions = []
    for rank, prefix in PREFIXES.items():
        expressions.append(
            pl.col("taxonomy")
            .str.extract(rf"(?:^|;\s*)({prefix}[^;]*)", 1)
            .fill_null(f"{prefix}")
            .alias(rank)
        )
    return frame.with_columns(expressions)


def determine_novelty(frame: pl.DataFrame) -> pl.DataFrame:
    # 1. Identify placeholders, supporting both un-prefixed (g1) and prefixed (Family_g1) syntax
    placeholders = {
        "species": pl.col("species").str.contains(r" sp_\d+$"),
        "genus": pl.col("genus").str.contains(r"(?:^|_)[g]\d+$"),
        "family": pl.col("family").str.contains(r"(?:^|_)[f]\d+$"),
        "order": pl.col("order").str.contains(r"(?:^|_)[o]\d+$"),
        "class": pl.col("class").str.contains(r"(?:^|_)[c]\d+$"),
        "phylum": pl.col("phylum").str.contains(r"(?:^|_)[p]\d+$"),
    }
    
    # 2. Build a bottom-up unbroken chain. 
    # A rank is only "novel" if it AND all ranks below it are placeholders.
    s_novel = placeholders["species"]
    g_novel = s_novel & placeholders["genus"]
    f_novel = g_novel & placeholders["family"]
    o_novel = f_novel & placeholders["order"]
    c_novel = o_novel & placeholders["class"]
    p_novel = c_novel & placeholders["phylum"]
    
    # 3. Assign just the rank name or "known"
    novelty_expr = (
        pl.when(p_novel).then(pl.lit("phylum"))
        .when(c_novel).then(pl.lit("class"))
        .when(o_novel).then(pl.lit("order"))
        .when(f_novel).then(pl.lit("family"))
        .when(g_novel).then(pl.lit("genus"))
        .when(s_novel).then(pl.lit("species"))
        .otherwise(pl.lit("known"))
    )
    
    return frame.with_columns(novelty_expr.alias("novelty"))


def read_taxonomy(path: Path, source: str, novelty: str | None = None) -> pl.DataFrame:
    frame = pl.read_csv(path, separator="\t")
    expected = {"genome_id", "taxonomy"}
    missing = expected.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")

    result = frame.select("genome_id", "taxonomy")
    if novelty is None or "novelty" not in frame.columns:
        result = result.with_columns(pl.lit(novelty or "known").alias("novelty"))
    else:
        result = result.with_columns(frame.get_column("novelty").fill_null(novelty).alias("novelty"))
    result = result.with_columns(pl.lit(source).alias("source"))
    return result


def main() -> None:
    reference = read_taxonomy(
        Path(str(snakemake.input.reference_taxonomy)),
        source="reference",
    )
    target = read_taxonomy(
        Path(str(snakemake.input.target_taxonomy)),
        source="new",
        novelty="novel",
    )

    metadata = (
        pl.concat([reference, target], how="vertical_relaxed")
        .unique(subset=["genome_id"], keep="last")
        .pipe(expand_taxonomy)
        .pipe(determine_novelty)
        .sort("genome_id")
    )
    
    output_path = Path(str(snakemake.output.metadata))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_csv(output_path, separator=",")


if __name__ == "__main__":
    main()