#!/usr/bin/env python3
"""Validate the taxonomy assignment by hiding database labels and checking they come back.

A random fraction of the reference genomes loses its label; ``assign_taxonomy`` then names those genomes like it
names any new genome, from the tree and the remaining labels. Comparing the result with the hidden truth measures,
per rank, how often the method gives the right database name, a placeholder (novel), or a *wrong* database name.
The wrong ones are the cost of the majority-name rule for MIXED nodes and of inheriting names from neighbours.

    label_masking.py  (snakemake: red_tree, red_cutoffs, taxonomy -> table; params fraction, replicates, seed)

A hidden genome whose true name has no other visible member at that rank (a species seen once) cannot be recovered
from labels; the right answer there is "novel". Those are counted as ``unseen`` and kept out of ``tested``, so
accuracy is only measured where recovery is possible.

Output rows: replicate, rank, tested (true name still visible elsewhere), correct, novel, wrong, unseen.
Several replicates with different masks give a spread, so the numbers come with an idea of their noise.
"""
import random

import polars as pl

import assign_taxonomy_v5 as at

RANKS = at.RANKS[1:]


def one_replicate(newick_text, known, cutoffs, fraction, seed):
    """Mask ``fraction`` of the labelled genomes in the tree and return per-rank counts."""
    tree = at.read_red_tree(newick_text)
    in_tree = sorted(leaf.taxon.label.split("|RED_")[0].replace("'", "") for leaf in tree.leaf_node_iter())
    labelled = [g for g in in_tree if at.parse_tax(known.get(g, "")).keys() - {"domain"}]
    rng = random.Random(seed)
    masked = set(rng.sample(labelled, max(1, int(len(labelled) * fraction))))
    visible = {g: ("" if g in masked else known[g]) for g in in_tree if g in known}
    reference, target, _ = at.assign_taxonomy(tree, visible, cutoffs)
    predicted = {r["genome_id"]: at.parse_tax(r["taxonomy"]) for r in reference + target}
    database_names = {rank: {at.parse_tax(known[g]).get(rank) for g in labelled if g not in masked} - {None} for rank in RANKS}
    counts = {rank: {"tested": 0, "correct": 0, "novel": 0, "wrong": 0, "unseen": 0} for rank in RANKS}
    for genome in masked:
        truth = at.parse_tax(known[genome])
        got = predicted[genome]
        for rank in RANKS:
            if rank not in truth:
                continue
            if truth[rank] not in database_names[rank]:
                counts[rank]["unseen"] += 1
                continue
            counts[rank]["tested"] += 1
            name = got.get(rank, "")
            if name == truth[rank]:
                counts[rank]["correct"] += 1
            elif name in database_names[rank]:
                counts[rank]["wrong"] += 1          # another existing database name
            else:
                counts[rank]["novel"] += 1          # a placeholder: the method says "new taxon"
    return len(masked), counts


def run(newick_text, known, cutoffs, fraction=0.15, replicates=3, seed=1):
    rows = []
    for rep in range(replicates):
        n_masked, counts = one_replicate(newick_text, known, cutoffs, fraction, seed + rep)
        for rank, c in counts.items():
            rows.append({"replicate": rep, "masked": n_masked, "rank": rank, **c})
    return pl.DataFrame(rows)


def main():
    from log_redirect import redirect_to_snakemake_log
    redirect_to_snakemake_log()
    cutoffs = at.load_red_cutoffs(str(snakemake.input.red_cutoffs))
    known = at.load_known_taxonomy(str(snakemake.input.taxonomy))
    text = open(str(snakemake.input.red_tree)).read()
    table = run(text, known, cutoffs, float(snakemake.params.fraction), int(snakemake.params.replicates),
                int(snakemake.params.seed))
    table.write_csv(str(snakemake.output.table), separator="\t")
    summary = (table.group_by("rank", maintain_order=True)
               .agg(pl.col("tested").sum(), pl.col("correct").sum(), pl.col("novel").sum(), pl.col("wrong").sum(), pl.col("unseen").sum()))
    print(summary)


if __name__ == "__main__" and "snakemake" in globals():
    main()
