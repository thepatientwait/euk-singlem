#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl


def build_paths(cluster_table: Path, base_dir: str) -> list[str]:
    if not cluster_table.exists() or cluster_table.stat().st_size == 0:
        return []

    df = pl.read_csv(cluster_table, separator="\t", has_header=False)
    first_col = df.columns[0]

    result = (
        df.select(pl.col(first_col).cast(pl.Utf8).str.strip_chars().alias("raw"))
        .filter(pl.col("raw").is_not_null() & (pl.col("raw") != ""))
        .with_columns(
            pl.when(pl.col("raw").str.starts_with("."))
            .then(pl.lit(base_dir.rstrip("/")) + pl.col("raw").str.slice(1))
            .otherwise(pl.col("raw"))
            .alias("genome_path")
        )
        .select("genome_path")
        .unique()
        .sort("genome_path")
    )

    return result["genome_path"].to_list()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert first column of a clusters TSV into normalized genome paths."
    )
    parser.add_argument("--cluster-table", required=True, help="Path to clusters TSV (e.g. genome_clusters_95.tsv)")
    parser.add_argument(
        "--base-dir",
        required=True,
        help="Base directory that replaces a leading '.' in column 1 values",
    )
    parser.add_argument("--output", required=True, help="Output manifest path")

    args = parser.parse_args()

    cluster_table = Path(args.cluster_table)
    output_path = Path(args.output)

    paths = build_paths(cluster_table=cluster_table, base_dir=args.base_dir)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{path}\n")


if __name__ == "__main__":
    main()
