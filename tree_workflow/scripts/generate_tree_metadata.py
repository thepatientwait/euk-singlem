#!/usr/bin/env python3

from pathlib import Path

import polars as pl

RANKS = ["domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = {
    "domain": "d__",
    "kingdom": "k__",
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
            .str.extract(rf"(?:^|; )({prefix}[^;]*)", 1)
            .fill_null(f"{prefix}")
            .alias(rank)
        )
    return frame.with_columns(expressions)


def read_taxonomy(path: Path, source: str, novelty: str | None = None) -> pl.DataFrame:
    frame = pl.read_csv(path, separator="\t")
    expected = {"genome_id", "taxonomy"}
    missing = expected.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")

    result = frame.select("genome_id", "taxonomy")
    if novelty is None or "novelty" not in frame.columns:
        result = result.with_columns(pl.lit(novelty or "known species").alias("novelty"))
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
        .sort("genome_id")
    )
    output_path = Path(str(snakemake.output.metadata))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_csv(output_path, separator=",")


if __name__ == "__main__":
    main()
