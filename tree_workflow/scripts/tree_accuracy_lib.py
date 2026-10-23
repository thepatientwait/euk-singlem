"""Tree-accuracy statistics for notebooks/tree_accuracy.ipynb and notebooks/tree_comparison.ipynb.

* ``fragmentation``: rooting-independent monophyly. A *piece* of a clade name is a maximal part of
  the unrooted tree (one side of an edge) whose labelled tips all carry that name; a monophyletic
  name has one piece. The count does not depend on where the tree is rooted: for a name X the
  tree is rooted on any labelled tip that is not X, and the pure-X sides never contain that tip.
* ``misplaced_tips``: tips whose labelled neighbourhood agrees (>= ``frac``) on another name that
  is *not* a sister (a different parent clade), so sister-taxon swaps are not flagged.
* ``rf_distance``: normalised Robinson-Foulds distance on the shared tips.
* RED helpers: spread of PhyloRank's per-taxon RED per rank and monotonicity of node RED.

Support values are the SH-like values FastTree writes on internal nodes (``clade.confidence``).
Re-rooting with Biopython can misplace the support of the few nodes on the path between the old
and the new root; the effect on the summaries is negligible and is not corrected.
"""
import copy
import json
import re
from collections import Counter, defaultdict

import numpy as np

from tree_fragmentation import PREFIX, RANKS, read_taxonomy  # noqa: F401  (re-exported)

PARENT_RANK = {"class": "phylum", "order": "class", "family": "order", "genus": "family"}
_RED_RE = re.compile(r"\|RED=([0-9.eE+-]+)")


# ----------------------------------------------------------------------------- pieces
def _pieces_rooted(tree, taxonomy, rank):
    """[(name, clade, labelled tips)] of the maximal pure clades of a rooted tree."""
    state, tips, pieces = {}, {}, []
    for clade in tree.find_clades(order="postorder"):
        key = id(clade)
        if clade.is_terminal():
            name = taxonomy.get(clade.name, {}).get(rank, "")
            state[key] = name
            tips[key] = [clade.name] if name else []
            continue
        child_states = [state[id(c)] for c in clade.clades]
        labels = {s for s in child_states if s}
        if None in child_states or len(labels) > 1:
            for child, s in zip(clade.clades, child_states):
                if s:
                    pieces.append((s, child, tips[id(child)]))
            state[key], tips[key] = None, []
        else:
            state[key] = next(iter(labels), "")
            tips[key] = [t for c in clade.clades for t in tips[id(c)]]
    root = state[id(tree.root)]
    if root:
        pieces.append((root, tree.root, tips[id(tree.root)]))
    return pieces


def _reroot_on(tree, tip_name):
    rerooted = copy.deepcopy(tree)
    target = next(t for t in rerooted.get_terminals() if t.name == tip_name)
    rerooted.root_with_outgroup(target)
    return rerooted


def pieces_unrooted(tree, taxonomy, rank):
    """{name: [(clade, tips) ...]} sorted largest first; independent of the tree's root."""
    labelled = [t.name for t in tree.get_terminals() if taxonomy.get(t.name, {}).get(rank)]
    names = {taxonomy[t][rank] for t in labelled}
    if len(names) < 2:
        return {n: [(tree.root, labelled)] for n in names}
    a = labelled[0]
    label_a = taxonomy[a][rank]
    b = next(t for t in labelled if taxonomy[t][rank] != label_a)
    out = defaultdict(list)
    for root_tip, valid in ((a, lambda n: n != label_a), (b, lambda n: n == label_a)):
        rooted = _reroot_on(tree, root_tip)
        for name, clade, tips in _pieces_rooted(rooted, taxonomy, rank):
            if valid(name):
                out[name].append((clade, tips))
    for name in out:
        out[name].sort(key=lambda p: -len(p[1]))
    return dict(out)


def fragmentation(tree, taxonomy, ranks=RANKS, pieces_by_rank=None):
    """Per-rank summary rows plus per-name rows (name, rank, tips, pieces, extra piece sizes)."""
    summary, names = [], []
    for rank in ranks:
        pieces = (pieces_by_rank or {}).get(rank) or pieces_unrooted(tree, taxonomy, rank)
        extra_pieces = 0
        fragmented = 0
        tips_in_fragmented = 0
        for name, ps in pieces.items():
            total = sum(len(t) for _, t in ps)
            names.append({"rank": rank, "name": name, "tips": total, "pieces": len(ps),
                          "extra_sizes": [len(t) for _, t in ps[1:]]})
            if len(ps) > 1:
                fragmented += 1
                extra_pieces += len(ps) - 1
                tips_in_fragmented += total
        summary.append({"rank": rank, "names": len(pieces), "fragmented_names": fragmented,
                        "extra_pieces": extra_pieces, "tips_in_fragmented": tips_in_fragmented})
    return summary, names


def extra_piece_supports(pieces):
    """Support of the edge above each non-largest piece: [(name, tips, support)]."""
    rows = []
    for name, ps in pieces.items():
        for clade, tips in ps[1:]:
            rows.append((name, len(tips), clade.confidence))
    return rows


# --------------------------------------------------------------------------- misplaced
def misplaced_tips(tree, taxonomy, rank, min_neighbours=8, frac=0.9, max_levels=60):
    """Tips labelled A whose nearest >= ``min_neighbours`` labelled neighbours are >= ``frac`` one
    other name B, where B is not a sister of A (different parent clade). Returns dict rows."""
    parent = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parent[child] = clade
    parent_rank = PARENT_RANK.get(rank)
    leaf_label = {t: taxonomy.get(t.name, {}).get(rank, "") for t in tree.get_terminals()}
    # labelled-tip counts per clade, bottom up
    counts = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            label = leaf_label[clade]
            counts[clade] = Counter({label: 1}) if label else Counter()
        else:
            total = Counter()
            for child in clade.clades:
                total.update(counts[child])
            counts[clade] = total
    rows = []
    for tip, label in leaf_label.items():
        if not label:
            continue
        node, levels = tip, 0
        while node in parent and levels < max_levels:
            node = parent[node]
            levels += 1
            if sum(counts[node].values()) - 1 >= min_neighbours:
                break
        neighbours = counts[node].copy()
        neighbours[label] -= 1
        neighbours = +neighbours
        total = sum(neighbours.values())
        if total < min_neighbours:
            continue
        other, n_other = neighbours.most_common(1)[0]
        if other == label or n_other / total < frac:
            continue
        same_parent = (parent_rank is not None
                       and taxonomy.get(tip.name, {}).get(parent_rank) not in (None, "")
                       and _name_parent(taxonomy, rank, other, parent_rank)
                       == taxonomy[tip.name][parent_rank])
        rows.append({"tip": tip.name, "label": label, "neighbours_say": other,
                     "agreement": n_other / total, "n_neighbours": total, "sister": bool(same_parent)})
    return rows


_PARENT_CACHE = {}


def _name_parent(taxonomy, rank, name, parent_rank):
    key = (id(taxonomy), rank, name)
    if key not in _PARENT_CACHE:
        value = None
        for lineage in taxonomy.values():
            if lineage.get(rank) == name:
                value = lineage.get(parent_rank)
                break
        _PARENT_CACHE[key] = value
    return _PARENT_CACHE[key]


# --------------------------------------------------------------------------------- RF
def _bipartitions(tree, index, common_mask):
    masks = {}
    splits = set()
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            masks[id(clade)] = (1 << index[clade.name]) if clade.name in index else 0
        else:
            m = 0
            for child in clade.clades:
                m |= masks[id(child)]
            masks[id(clade)] = m
        m = masks[id(clade)] & common_mask
        size = bin(m).count("1")
        if 2 <= size <= bin(common_mask).count("1") - 2:
            lowest = common_mask & -common_mask
            splits.add(m if not (m & lowest) else common_mask ^ m)
    return splits


def rf_distance(tree_a, tree_b):
    """(RF distance, normalised RF, shared tips) on the tips both trees have."""
    names_a = {t.name for t in tree_a.get_terminals()}
    names_b = {t.name for t in tree_b.get_terminals()}
    common = sorted(names_a & names_b)
    if len(common) < 4:
        return 0, float("nan"), len(common)
    index = {n: i for i, n in enumerate(common)}
    mask = (1 << len(common)) - 1
    sa, sb = _bipartitions(tree_a, index, mask), _bipartitions(tree_b, index, mask)
    rf = len(sa ^ sb)
    return rf, rf / max(1, len(sa) + len(sb)), len(common)


# -------------------------------------------------------------------------------- RED
def red_by_rank(tsv_path, cutoffs_path):
    """{rank: (cutoff, [median RED of every named taxon of that rank])} from PhyloRank's table."""
    with open(cutoffs_path, "r", encoding="utf-8") as handle:
        cutoffs = json.load(handle)
    values = defaultdict(list)
    with open(tsv_path, "r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            for rank, prefix in PREFIX.items():
                if fields[0].startswith(prefix):
                    try:
                        values[rank].append(float(fields[2]))
                    except ValueError:
                        pass
                    break
    return {rank: (cutoffs.get(rank), values.get(rank, [])) for rank in cutoffs}


def red_monotonicity(tree):
    """(internal nodes with RED, nodes below their parent, negative branches, most negative length)
    for a PhyloRank RED tree read with Bio.Phylo (internal names look like ``0.93|RED=0.52``)."""
    red = {}
    for clade in tree.find_clades(order="preorder"):
        match = _RED_RE.search(clade.name or "")
        if match and not clade.is_terminal():
            red[clade] = float(match.group(1))
    parent = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parent[child] = clade
    lower = 0
    for clade, value in red.items():
        anc = parent.get(clade)
        while anc is not None and anc not in red:
            anc = parent.get(anc)
        if anc is not None and value < red[anc] - 1e-12:
            lower += 1
    lengths = [c.branch_length for c in tree.find_clades() if c.branch_length is not None]
    negative = [x for x in lengths if x < 0]
    return len(red), lower, len(negative), (min(negative) if negative else 0.0)


def rank_overlap(red_values):
    """Fraction of taxa at each rank whose RED falls on the wrong side of the neighbouring rank's
    median: {rank: (share below the median of the rank above, share above the median of the rank below)}."""
    medians = {r: float(np.median(v)) for r, (_, v) in red_values.items() if v}
    ranks = [r for r in ("phylum", "class", "order", "family", "genus", "species") if r in medians]
    out = {}
    for i, rank in enumerate(ranks):
        v = np.asarray(red_values[rank][1])
        below = float((v < medians[ranks[i - 1]]).mean()) if i > 0 else float("nan")
        above = float((v > medians[ranks[i + 1]]).mean()) if i + 1 < len(ranks) else float("nan")
        out[rank] = (below, above)
    return out


# ------------------------------------------------------------------------------ BUSCO
def busco_status_counts(full_table):
    """Counter of BUSCO statuses (Complete, Duplicated, Fragmented, Missing) of one genome."""
    counts = Counter()
    with open(full_table, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t", 2)
            if len(fields) > 1:
                counts[fields[1]] += 1
    return counts
