#!/usr/bin/env python3
"""How fragmented are the database taxa on a tree?

For each rank and name, a *piece* is a maximal clade whose labelled tips all carry that name
(tips without a label for that rank are ignored). A monophyletic name has one piece; the
``extra`` pieces are the fragmentation. Used to compare trees (e.g. FastTree ``-nt`` vs
``-lg -gamma``) and by the diagnostics notebook.

    python tree_fragmentation.py conservative_taxonomy.tsv TREE.nwk [TREE2.nwk ...]
"""
import sys
from collections import Counter, defaultdict

from Bio import Phylo

RANKS = ["phylum", "class", "order", "family", "genus"]
PREFIX = {"phylum": "p__", "class": "c__", "order": "o__", "family": "f__", "genus": "g__"}


def read_taxonomy(path):
    """{tip: {rank: name}} from a two-column file: accession <tab> d__..;p__..;c__..;..."""
    out = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0] == "acc":
                continue
            ranks = {}
            for part in fields[1].split(";"):
                for rank, prefix in PREFIX.items():
                    if part.startswith(prefix) and len(part) > len(prefix):
                        ranks[rank] = part
            out[fields[0]] = ranks
    return out


def count_pieces(tree, taxonomy):
    """{rank: {name: (n_tips, n_pieces)}}; iterative post-order, labels summarised per node."""
    result = {}
    root = tree.root
    order = list(tree.find_clades(order="postorder"))
    for rank in RANKS:
        tips = Counter()
        pieces = Counter()
        state = {}  # clade id -> label (str), "" = nothing labelled, None = mixed
        for clade in order:
            if clade.is_terminal():
                name = taxonomy.get(clade.name, {}).get(rank, "")
                if name:
                    tips[name] += 1
                state[id(clade)] = name
                continue
            child_states = [state[id(c)] for c in clade.clades]
            labels = {s for s in child_states if s}
            if None in child_states or len(labels) > 1:
                # this node is mixed: every pure child is a maximal piece
                for s in child_states:
                    if s:
                        pieces[s] += 1
                state[id(clade)] = None
            else:
                state[id(clade)] = next(iter(labels), "")
        if state[id(root)]:
            pieces[state[id(root)]] += 1
        result[rank] = {name: (tips[name], pieces[name]) for name in tips}
    return result


def summarise(result):
    rows = []
    for rank in RANKS:
        names = result[rank]
        fragmented = [n for n, (_, p) in names.items() if p > 1]
        rows.append({
            "rank": rank, "names": len(names), "fragmented_names": len(fragmented),
            "extra_pieces": sum(p - 1 for _, p in names.values()),
            "tips_in_fragmented": sum(t for n, (t, p) in names.items() if p > 1),
        })
    return rows


def worst(result, rank, n=10):
    return sorted(((p, t, name) for name, (t, p) in result[rank].items() if p > 1), reverse=True)[:n]


def main(argv):
    taxonomy = read_taxonomy(argv[1])
    for path in argv[2:]:
        tree = Phylo.read(path, "newick")
        result = count_pieces(tree, taxonomy)
        print(f"\n== {path}  ({tree.count_terminals()} tips)")
        print(f"{'rank':8}{'names':>7}{'fragmented':>12}{'extra pieces':>14}{'tips in frag.':>15}")
        for r in summarise(result):
            print(f"{r['rank']:8}{r['names']:>7}{r['fragmented_names']:>12}{r['extra_pieces']:>14}{r['tips_in_fragmented']:>15}")
        for rank in ("phylum", "class"):
            print(f"  worst {rank}: " + ", ".join(f"{name}({t} tips/{p} pieces)" for p, t, name in worst(result, rank, 6)))


if __name__ == "__main__":
    main(sys.argv)
