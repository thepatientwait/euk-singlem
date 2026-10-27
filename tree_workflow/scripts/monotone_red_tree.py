#!/usr/bin/env python3
"""Make PhyloRank's RED trees monotone: no negative branch lengths, RED never below the parent's.

PhyloRank (``outliers``) sets each node's RED to the median over many rootings (one per phylum) and
rewrites branch lengths as ``RED(child) - RED(parent)``. A node only receives a value from the
rootings in which it lies on the ingroup side, so a parent and its child are medians over different
sets of rootings and the child can end up below the parent. That gives a few negative branch lengths
(near the root, in very large nodes).

This script fixes that in place, for ``*.red_decorated.tree`` (``|RED=`` labels) and
``*.scaled.tree`` (no labels):

* RED(node) is rebuilt from the branch lengths (root = 0, RED = sum of lengths from the root),
* each node is raised to at least its parent's RED, root to tip,
* branch length = RED(node) - RED(parent) (so always >= 0),
* ``|RED=`` labels are rewritten from the same values (3 decimals, like PhyloRank).

The untouched PhyloRank file is kept next to the result as ``<name>.phylorank_raw<suffix>`` (created
once; a later call starts from it again, so the script is idempotent) for diagnostics.

    monotone_red_tree.py TREE [TREE ...]
"""
import re
import shutil
import sys
from pathlib import Path

import dendropy

_RED_RE = re.compile(r"\|RED_(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")
EPS = 1e-12


def raw_path(path):
    path = Path(path)
    return path.with_name(f"{path.stem}.phylorank_raw{path.suffix}")


def read_tree(text):
    # '=' inside labels is awkward for newick readers
    return dendropy.Tree.get(data=text.replace("|RED=", "|RED_"), schema="newick",
                             preserve_underscores=True, suppress_internal_node_taxa=True)


def node_label(node):
    return node.taxon.label if node.is_leaf() and node.taxon is not None else node.label


def set_node_label(node, value):
    if node.is_leaf() and node.taxon is not None:
        node.taxon.label = value
    else:
        node.label = value


def make_monotone(tree):
    """Fix ``tree`` in place. Returns a stats dict.

    ``raised``: nodes whose RED had to be raised to their parent's; ``negative_before``: negative
    branches before; ``max_raise``: largest increase of a node's RED; ``negative_after``: should be 0.
    """
    stats = {"nodes": 0, "negative_before": 0, "most_negative_before": 0.0, "raised": 0,
             "max_raise": 0.0, "negative_after": 0}
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        if parent is None:
            node.red_fixed = 0.0
            node.red_raw = 0.0
            continue
        length = node.edge_length if node.edge_length is not None else 0.0
        node.red_raw = parent.red_raw + length
        if length < 0:
            stats["negative_before"] += 1
            stats["most_negative_before"] = min(stats["most_negative_before"], length)
        node.red_fixed = max(node.red_raw, parent.red_fixed)
        if node.red_fixed > node.red_raw + EPS:
            stats["raised"] += 1
            stats["max_raise"] = max(stats["max_raise"], node.red_fixed - node.red_raw)
    for node in tree.preorder_node_iter():
        stats["nodes"] += 1
        parent = node.parent_node
        if parent is not None:
            node.edge_length = max(node.red_fixed - parent.red_fixed, 0.0)
            if node.edge_length < 0:
                stats["negative_after"] += 1
        label = node_label(node)
        if label and _RED_RE.search(label):
            set_node_label(node, _RED_RE.sub(f"|RED_{node.red_fixed:.3f}", label))
    return stats


def write_tree(tree, path):
    tree.write_to_path(str(path), schema="newick", suppress_rooting=True, unquoted_underscores=True)
    text = Path(path).read_text()
    Path(path).write_text(text.replace("|RED_", "|RED="))


def fix_file(path):
    """Make one tree file monotone (starting from its saved raw copy if there is one)."""
    path = Path(path)
    raw = raw_path(path)
    if not raw.exists():
        shutil.copyfile(path, raw)
    tree = read_tree(raw.read_text())
    stats = make_monotone(tree)
    write_tree(tree, path)
    return stats


def main(paths):
    for path in paths:
        stats = fix_file(path)
        print(f"{path}: {stats['negative_before']} negative branches "
              f"(most negative {stats['most_negative_before']:.4f}) before; "
              f"{stats['raised']} nodes raised to their parent's RED (largest raise {stats['max_raise']:.4f}); "
              f"{stats['negative_after']} negative after. Raw PhyloRank tree: {raw_path(path)}", flush=True)
        if stats["negative_after"]:
            raise SystemExit(f"{path}: negative branches remain")


if __name__ == "__main__":
    main(sys.argv[1:])
