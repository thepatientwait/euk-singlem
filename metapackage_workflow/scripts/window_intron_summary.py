#!/usr/bin/env python3
"""Collect the per-gene window decisions of rule ``intron_window`` into one table.

Logs the genes whose corrected window still crosses an intron in more than ``max_span_fraction``
of the genomes it covers; with ``params.strict`` the run stops there, so a metapackage cannot be
built from such windows unnoticed.
"""
import polars as pl

from log_config import configure_logger

logger = configure_logger(__file__)


def main():
    frames = [pl.read_csv(path, separator="\t") for path in snakemake.input.choices]
    summary = pl.concat(frames).sort("busco_id")
    summary.write_csv(snakemake.output.summary, separator="\t")

    limit = summary["max_span_fraction"][0]
    naive_over = summary.filter(pl.col("naive_span_fraction") > limit)
    moved = summary.filter(pl.col("corrected_start") != pl.col("naive_start"))
    violations = summary.filter(pl.col("violation"))
    logger.info("%d genes: the naive (best-coverage) window crosses an intron in more than %.0f%% of genomes for %d; "
                "%d corrected windows differ from the naive one; %d genes have no corrected window under the limit",
                len(summary), 100 * limit, len(naive_over), len(moved), len(violations))
    for row in violations.iter_rows(named=True):
        logger.warning("%s: corrected window %d spans an intron in %.1f%% of genomes (limit %.1f%%)", row["busco_id"],
                       row["corrected_start"], 100 * row["corrected_span_fraction"], 100 * row["max_span_fraction"])
    if violations.height and snakemake.params.strict:
        raise RuntimeError(f"{violations.height} genes have no window under the intron limit: "
                           f"{violations['busco_id'].to_list()}")


if __name__ == "__main__":
    main()
