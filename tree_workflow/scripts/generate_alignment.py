#!/usr/bin/env python3
"""Concatenate the per-gene BMGE alignments into one supermatrix (one row per genome).

* A genome that has several records for a gene (should not happen: the aggregation keeps one per
  BUSCO) is represented by its lowest-numbered copy (``genome|busco|1`` = the top BUSCO hit),
  chosen explicitly rather than by file order.
* Genes a genome lacks are filled with gaps.
* Genomes with fewer than ``min_genes_per_genome`` genes are dropped: with few genes they carry
  little signal and mostly add gaps (and unstable placement) to the tree.
"""
import re
import sys
from pathlib import Path

from Bio import SeqIO

try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:  # imported from tests
    redirect_to_snakemake_log = None


def copy_rank(record_id):
    """Copy number from ``genome|busco|copy`` (1 = best hit); unparsable ids rank last."""
    match = re.search(r"\|(\d+)$", record_id)
    return int(match.group(1)) if match else sys.maxsize


def read_gene_alignment(path):
    """{genome: aligned sequence} for one gene, best copy per genome. Also returns the number of
    genomes that had more than one record."""
    best = {}
    repeated = set()
    for record in SeqIO.parse(str(path), "fasta"):
        genome = record.id.split("|")[0]
        rank = copy_rank(record.id)
        if genome in best:
            repeated.add(genome)
            if rank >= best[genome][0]:
                continue
        best[genome] = (rank, str(record.seq))
    return {genome: seq for genome, (_, seq) in best.items()}, len(repeated)


def concatenate(alignment_files, exclude=frozenset(), min_genes_per_genome=1):
    """Returns ({genome: concatenated sequence}, stats dict)."""
    per_gene = []
    repeated_total = 0
    for path in alignment_files:
        genomes, repeated = read_gene_alignment(path)
        if not genomes:
            continue
        lengths = {len(seq) for seq in genomes.values()}
        if len(lengths) != 1:
            raise ValueError(f"{path}: aligned sequences have different lengths {sorted(lengths)[:5]}")
        per_gene.append((Path(path).name.split("_bmge")[0], lengths.pop(), genomes))
        repeated_total += repeated

    all_genomes = sorted({g for _, _, genomes in per_gene for g in genomes} - set(exclude))
    n_genes = {g: sum(g in genomes for _, _, genomes in per_gene) for g in all_genomes}
    kept = [g for g in all_genomes if n_genes[g] >= min_genes_per_genome]

    sequences = {
        genome: "".join(genomes.get(genome) or "-" * length for _, length, genomes in per_gene)
        for genome in kept
    }
    stats = {
        "genes": len(per_gene), "genomes_seen": len(all_genomes), "genomes_kept": len(kept),
        "genomes_dropped_few_genes": len(all_genomes) - len(kept),
        "genes_with_repeated_genomes": repeated_total,
        "columns": sum(length for _, length, _ in per_gene),
    }
    return sequences, stats


def main():
    redirect_to_snakemake_log()
    inputs = snakemake.input.alignment_files
    inputs = [Path(inputs)] if isinstance(inputs, (str, Path)) else [Path(p) for p in inputs]
    valid = [p for p in inputs if p.is_file()]
    if len(valid) != len(inputs):
        missing = ", ".join(str(p) for p in inputs if not p.is_file())
        if not valid:
            raise FileNotFoundError(f"No valid BMGE alignment files were found. Invalid inputs: {missing}")
        print(f"Skipping invalid BMGE alignment inputs: {missing}")

    sequences, stats = concatenate(
        valid,
        exclude=set(snakemake.params.exclude),
        min_genes_per_genome=int(snakemake.params.min_genes_per_genome),
    )
    print("generate_alignment:", stats)
    if not sequences:
        raise RuntimeError("No genome has enough genes for the concatenated alignment")
    with open(snakemake.output.concatinated_alignment, "w") as handle:
        for genome, sequence in sequences.items():
            handle.write(f">{genome}\n{sequence}\n")


if __name__ == "__main__":
    main()
