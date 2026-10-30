#!/usr/bin/env python3

from pathlib import Path

import polars as pl

from log_redirect import redirect_to_snakemake_log
from source_groups import group_of, read_source_groups

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
    """One column per rank (``p__Ascomycota`` ...): the lineage is ``;``-separated and every item carries its prefix."""
    expressions = []
    for rank, prefix in PREFIXES.items():
        expressions.append(
            pl.col("taxonomy").str.split(";")
            .list.eval(pl.element().str.strip_chars().filter(pl.element().str.starts_with(prefix)))
            .list.first()
            .fill_null(prefix)
            .alias(rank)
        )
    return frame.with_columns(expressions)


def read_taxonomy(path: Path, source: str, novelty: str | None = None) -> pl.DataFrame:
    # assign_taxonomy_v5 writes headerless two-column TSVs (genome_id, taxonomy).
    # Older versions wrote a header (and optionally a "novelty" column); accept both.
    with open(path) as handle:
        first_line = handle.readline().rstrip("\n")
    has_header = first_line.split("\t")[0] == "genome_id"

    if has_header:
        frame = pl.read_csv(path, separator="\t")
    else:
        frame = pl.read_csv(path, separator="\t", has_header=False)
        if frame.width < 2:
            raise ValueError(f"{path} must have at least 2 tab-separated columns (genome_id, taxonomy)")
        names = ["genome_id", "taxonomy"] + [f"extra_{i}" for i in range(frame.width - 2)]
        frame = frame.rename(dict(zip(frame.columns, names)))

    expected = {"genome_id", "taxonomy"}
    missing = expected.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")

    result = frame.select("genome_id", "taxonomy")
    if "novelty" in frame.columns:
        # written by assign_taxonomy: the highest rank from which everything below is a new taxon, else "known"
        result = result.with_columns(frame.get_column("novelty").fill_null(novelty or "known").alias("novelty"))
    else:
        # taxonomy files from before assign_taxonomy wrote novelty: targets are new, references known
        result = result.with_columns(pl.lit(novelty or "known").alias("novelty"))
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

    # download source of each genome (fungi, algae, p10k, ncbi, tsa, external, ...); "unknown" when a genome is in
    # no BUSCO QC manifest. Lets the notebooks judge transcriptome (tsa) tips separately.
    groups = read_source_groups(list(getattr(snakemake.input, "busco_manifests", []) or []))

    metadata = (
        pl.concat([reference, target], how="vertical_relaxed")
        .unique(subset=["genome_id"], keep="last")
        .pipe(expand_taxonomy)
        .with_columns(pl.col("genome_id").map_elements(lambda g: group_of(g, groups), return_dtype=pl.String).alias("source_group"))
        .sort("genome_id")
    )
    
    output_path = Path(str(snakemake.output.metadata))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata.write_csv(output_path, separator=",")


if __name__ == "__main__":
    main()