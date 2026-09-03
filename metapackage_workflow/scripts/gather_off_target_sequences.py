#!/usr/bin/env python3

import polars as pl
from pathlib import Path
from typing import Any

from log_config import configure_logger


logger = configure_logger(__file__)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake: Any = None


def read_genes(ids_path: Path) -> list[str]:
    df = pl.read_csv(ids_path, separator="\t")
    return sorted(set(df["gene"].to_list()))


def load_per_genome_rows(paths: list[Path]) -> dict[str, list[tuple[str, str]]]:
    """Load sequences from per-genome TSV files using polars."""
    all_dfs = []
    
    for path in sorted(paths):
        if not path.exists() or path.stat().st_size == 0:
            continue
        
        try:
            df = pl.read_csv(path, separator="\t")
            # Keep only non-null rows for required columns
            df = df.filter(
                pl.col("gene").is_not_null()
                & pl.col("sequence_id").is_not_null()
                & pl.col("sequence").is_not_null()
            )
            all_dfs.append(df)
        except Exception as e:
            logger.warning(f"Error reading {path}: {e}")
            continue
    
    if not all_dfs:
        return {}
    
    # Concatenate all dataframes
    combined = pl.concat(all_dfs)
    
    # Group by gene and collect (sequence_id, sequence) tuples
    sequences_by_gene: dict[str, list[tuple[str, str]]] = {}
    for gene, group_df in combined.group_by("gene"):
        sequences_by_gene[gene[0]] = [
            (row["sequence_id"], row["sequence"])
            for row in group_df.select(["sequence_id", "sequence"]).iter_rows(named=True)
        ]
    
    return sequences_by_gene


def write_outputs(
    genes: list[str],
    sequences_by_gene: dict[str, list[tuple[str, str]]],
    aggregate_faa_pattern: str,
    aggregate_tax_pattern: str,
    default_taxonomy: str,
) -> None:
    for gene in genes:
        faa_path = Path(aggregate_faa_pattern.format(busco_gene=gene))
        tax_path = Path(aggregate_tax_pattern.format(busco_gene=gene))
        faa_path.parent.mkdir(parents=True, exist_ok=True)
        tax_path.parent.mkdir(parents=True, exist_ok=True)

        entries = sorted(sequences_by_gene.get(gene, []), key=lambda x: x[0])

        with open(faa_path, "w", encoding="utf-8") as faa_out:
            for sequence_id, sequence in entries:
                faa_out.write(f">{sequence_id}\n{sequence}\n")

        with open(tax_path, "w", encoding="utf-8") as tax_out:
            tax_out.write("sequence_id\ttaxonomy\n")
            for sequence_id, _ in entries:
                tax_out.write(f"{sequence_id}\t{default_taxonomy}\n")


def main() -> None:
    hmm_ids = Path(str(snakemake.input.hmm_ids))
    per_genome_results = [Path(str(path)) for path in snakemake.input.per_genome_results]

    aggregate_faa_pattern = str(snakemake.params.aggregate_faa_pattern)
    aggregate_tax_pattern = str(snakemake.params.aggregate_tax_pattern)
    default_taxonomy = str(snakemake.params.default_taxonomy)

    aggregate_dir = Path(str(snakemake.output.aggregate_dir))
    done_path = Path(str(snakemake.output.done))

    aggregate_dir.mkdir(parents=True, exist_ok=True)
    done_path.parent.mkdir(parents=True, exist_ok=True)

    genes = read_genes(hmm_ids)
    logger.info(
        "Gathering off-target sequences for %d genes from %d per-genome files",
        len(genes),
        len(per_genome_results),
    )

    sequences_by_gene = load_per_genome_rows(per_genome_results)

    write_outputs(
        genes=genes,
        sequences_by_gene=sequences_by_gene,
        aggregate_faa_pattern=aggregate_faa_pattern,
        aggregate_tax_pattern=aggregate_tax_pattern,
        default_taxonomy=default_taxonomy,
    )

    done_path.touch()
    logger.info("Completed off-target gather")


if __name__ == "__main__":
    main()
