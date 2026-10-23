#!/usr/bin/env python3
"""Assign taxonomy to every leaf of a RED-decorated tree.

Inputs (via snakemake):
  red_tree        phylorank RED-decorated newick (labels contain ``|RED=<value>``)
  red_cutoffs     JSON mapping rank -> median RED (phylorank ``.dict`` file)
  taxonomy        headerless TSV ``genome_id<TAB>taxonomy`` (database labels)

Outputs (TSV with header ``genome_id<TAB>taxonomy``):
  reference_taxonomy   leaves that had a database taxonomy
  target_taxonomy      leaves without one (new genomes), with placeholder names

Method
------
1. Leaves: parse database labels. A leaf with no rank below domain is a *target*.
2. Harmonise higher ranks: every name at rank R maps to a single parent chain
   (the majority among the reference leaves carrying that name; ties break
   alphabetically). Applied species -> class so the result is hierarchy-consistent.
3. Bottom-up: each node records, per rank, the single value shared by all
   reference leaves below it, or MIXED if they disagree. MIXED is sticky: it
   propagates to every ancestor and stops lower ranks from being assigned.
4. Top-down: nodes inherit from their parent, apply their own single-valued
   ranks (clearing lower inherited ranks when a rank changes), drop inherited
   ranks that are deeper than the node's RED allows, and mint placeholder names
   (``Parent_f3``, ``Genus sp_7``) for empty ranks at leaves and at internal
   nodes whose RED exceeds that rank's cutoff.

Known limitation / TODO
-----------------------
Polyphyletic database names are NOT split into ``Name_1``, ``Name_2`` ... any
more. The old splitting logic (exact-purity clades) fragmented clades because of
single mislabelled genomes and gave traversal-order-dependent suffixes. It was
removed pending a redesign; polyphyletic names are currently left as they are
(their node is MIXED, so targets inside them are named by RED placeholders).
"""

import json
import re
from collections import Counter, defaultdict

import dendropy
import polars as pl

from log_redirect import redirect_to_snakemake_log

# send this rule's stdout/stderr to its log: file
redirect_to_snakemake_log()

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = {rank: f"{rank[0]}__" for rank in RANKS}
DEFAULT_DOMAIN = "Eukarya"

MIXED = object()  # sentinel: descendants disagree at this rank

_RED_RE = re.compile(r"\|RED_(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)")


# --------------------------------------------------------------------------
# Taxonomy string helpers
# --------------------------------------------------------------------------
def parse_tax(tax_str):
    """'d__X;p__Y;...' -> {'domain': 'X', 'phylum': 'Y', ...}; blank ranks omitted."""
    parsed = {}
    for chunk in (tax_str or "").split(";"):
        chunk = chunk.strip()
        for rank, prefix in PREFIXES.items():
            if chunk.startswith(prefix):
                value = chunk[len(prefix):].strip()
                if value:
                    parsed[rank] = value
                break
    return parsed


def format_tax(tax):
    return ";".join(f"{PREFIXES[rank]}{tax.get(rank, '')}" for rank in RANKS)


def nearest_ancestor(tax, rank):
    """Closest populated rank above ``rank`` -> (rank, value); (None, 'unk') if none."""
    for i in range(RANKS.index(rank) - 1, -1, -1):
        value = tax.get(RANKS[i])
        if value:
            return RANKS[i], value
    return None, "unk"


# --------------------------------------------------------------------------
# Tree preparation
# --------------------------------------------------------------------------
def parse_node_label(node):
    """Return (name, red) for a node from its label / taxon label."""
    raw = node.label
    if raw is None and node.taxon is not None:
        raw = node.taxon.label
    raw = (raw or "").replace("'", "")
    red_match = _RED_RE.search(raw)
    red = float(red_match.group(1)) if red_match else None
    return raw.split("|RED_")[0], red


def read_red_tree(newick_text):
    # '=' inside labels is awkward for newick readers, so normalise it first.
    return dendropy.Tree.get(
        data=newick_text.replace("|RED=", "|RED_"),
        schema="newick",
        preserve_underscores=True,
    )


def report_negative_branches(tree):
    """Number and extent of negative branch lengths in the RED tree (printed to the log).

    PhyloRank rewrites branch lengths as differences between median RED values over several
    rootings, which can be negative where RED is not monotone along a path. The lengths are
    not divergences; assignment uses the RED value stored on each node, which is unaffected.
    Returns (count, most_negative).
    """
    negative = [e.length for e in tree.preorder_edge_iter() if e.length is not None and e.length < 0]
    if negative:
        print(f"NOTE: {len(negative)} branches have negative length in the PhyloRank RED tree "
              f"(most negative {min(negative):.4f}, total {sum(negative):.4f}). These are artefacts of "
              "PhyloRank's RED rescaling; branch lengths are ignored and node RED values are used (made monotone, see below).")
    else:
        print("No negative branch lengths in the RED tree.")
    return len(negative), (min(negative) if negative else 0.0)


def annotate_leaves(tree, known_tax):
    """Set genome_id / is_target / strict (rank -> value) / red on every leaf."""
    leaves = []
    seen = set()
    for leaf in tree.leaf_node_iter():
        genome_id, _ = parse_node_label(leaf)
        if genome_id in seen:
            raise ValueError(f"Duplicate leaf name in tree: {genome_id}")
        seen.add(genome_id)

        strict = parse_tax(known_tax.get(genome_id, ""))
        strict.pop("domain", None)  # domain is not informative here
        leaf.genome_id = genome_id
        leaf.is_target = not strict
        leaf.strict = strict
        leaf.red = 1.0
        leaves.append(leaf)
    return leaves


# --------------------------------------------------------------------------
# Step 2: hierarchy-consistent higher ranks
# --------------------------------------------------------------------------
def harmonise_higher_ranks(leaves):
    """Give every name at rank R one parent chain (majority among reference leaves).

    Goes species -> class so that after each pass all groups at that rank are
    internally consistent and stay so in later (higher-rank) passes. Blank higher
    ranks of grouped leaves are filled with the majority value.
    Returns {rank: number of leaf labels changed or filled}.
    """
    references = [leaf for leaf in leaves if not leaf.is_target]
    changes = Counter()

    for idx in range(len(RANKS) - 1, 1, -1):  # species .. class (phylum needs a parent: domain)
        rank = RANKS[idx]
        groups = defaultdict(list)
        for leaf in references:
            name = leaf.strict.get(rank)
            if name:
                groups[name].append(leaf)

        for name in sorted(groups):
            members = groups[name]
            for high_rank in RANKS[1:idx]:
                values = Counter(m.strict[high_rank] for m in members if m.strict.get(high_rank))
                if not values:
                    continue
                winner = min(values, key=lambda v: (-values[v], v))
                for member in members:
                    if member.strict.get(high_rank) != winner:
                        member.strict[high_rank] = winner
                        changes[high_rank] += 1
    return dict(changes)


# --------------------------------------------------------------------------
# Step 3: bottom-up strict state
# --------------------------------------------------------------------------
def annotate_internal_nodes(tree):
    """Set node.state (rank -> value | MIXED) and node.red for every node.

    A rank is absent from ``state`` when no reference leaf below has a value
    (no information, distinct from MIXED). MIXED is sticky and, as soon as it
    occurs, lower ranks are left unset.
    """
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.state = dict(node.strict)
            continue

        _, node.red = parse_node_label(node)

        state = {}
        for rank in RANKS:
            merged = None
            for child in node.child_node_iter():
                value = child.state.get(rank)
                if value is None:
                    continue
                if merged is None:
                    merged = value
                elif merged != value:
                    merged = MIXED
                    break
            if merged is None:
                continue
            state[rank] = merged
            if merged is MIXED:
                break
        node.state = state


def enforce_monotone_red(tree):
    """Make RED non-decreasing from the root to the tips; returns the number of nodes raised.

    PhyloRank's RED is the median over several rootings, which need not rise from a parent to
    its child (that is also what produces negative branch lengths in the RED tree). A child
    that is less diverged than its parent is meaningless for rank assignment, so each node's
    RED is raised to at least its nearest ancestor's. Only the node RED values are used here,
    never PhyloRank's rewritten branch lengths.
    """
    changed = 0
    for node in tree.preorder_node_iter():
        if node.is_leaf() or node.red is None:
            continue
        ancestor = node.parent_node
        while ancestor is not None and getattr(ancestor, "red", None) is None:
            ancestor = ancestor.parent_node
        if ancestor is not None and node.red < ancestor.red:
            node.red = ancestor.red
            changed += 1
    return changed


# --------------------------------------------------------------------------
# Step 4: top-down assignment
# --------------------------------------------------------------------------
def apply_state(assigned, state):
    """Apply a node's single-valued ranks; a changed rank invalidates all lower ranks."""
    for idx, rank in enumerate(RANKS):
        value = state.get(rank)
        if value is MIXED:
            break
        if value is None:
            continue
        if assigned.get(rank) != value:
            assigned[rank] = value
            for lower in RANKS[idx + 1:]:
                assigned.pop(lower, None)


def truncate_by_red(assigned, red, red_cutoffs):
    """Drop inherited ranks deeper than the node's RED supports.

    The deepest rank whose cutoff is below ``red`` is kept; ranks below it are
    removed. If no cutoff is below ``red`` (node older than the phylum cutoff)
    nothing is removed, i.e. database labels are kept.
    """
    if red is None:
        return
    allowed = [i for i, rank in enumerate(RANKS)
               if red_cutoffs.get(rank) is not None and red > red_cutoffs[rank]]
    if not allowed:
        return
    for rank in RANKS[allowed[-1] + 1:]:
        assigned.pop(rank, None)


class Placeholder:
    """A not-yet-numbered name for a new taxon (``Parent_f3`` / ``Genus sp_7``).

    Placeholders are created wherever a node needs a new taxon but are only
    numbered when a leaf that actually uses them is written out, so the numbers
    in the output have no gaps. ``parent`` is the nearest populated rank above
    (a str, or another Placeholder that is resolved first).
    """

    __slots__ = ("rank", "ancestor_rank", "parent", "name")

    def __init__(self, rank, ancestor_rank, parent):
        self.rank = rank
        self.ancestor_rank = ancestor_rank
        self.parent = parent
        self.name = None


def resolve_name(value, counters, used_names):
    """str -> itself; Placeholder -> its (lazily allocated) unique name."""
    if not isinstance(value, Placeholder):
        return value
    if value.name is None:
        prefix = resolve_name(value.parent, counters, used_names)
        key = (value.ancestor_rank, prefix, value.rank)
        while True:
            counters[key] += 1
            n = counters[key]
            if value.rank == "species":
                name = f"{prefix} sp_{n}"
            else:
                code = f"{value.rank[0]}{n}"
                name = f"{prefix}_{code}" if prefix != "unk" else code
            if name not in used_names[value.rank]:  # never reuse a database name
                used_names[value.rank].add(name)
                value.name = name
                break
    return value.name


def assign_top_down(tree, red_cutoffs, used_names):
    counters = defaultdict(int)
    reference, target = [], []

    for node in tree.preorder_node_iter():
        parent = node.parent_node
        assigned = dict(parent.assigned) if parent is not None else {}
        apply_state(assigned, node.state)
        if parent is None:
            assigned.setdefault("domain", DEFAULT_DOMAIN)
        truncate_by_red(assigned, node.red, red_cutoffs)

        is_leaf = node.is_leaf()
        for rank in RANKS:
            if assigned.get(rank):
                continue
            cutoff = red_cutoffs.get(rank)
            if is_leaf or (node.red is not None and cutoff is not None and node.red > cutoff):
                ancestor_rank, parent = nearest_ancestor(assigned, rank)
                assigned[rank] = Placeholder(rank, ancestor_rank, parent)

        node.assigned = assigned
        if is_leaf:
            names = {rank: resolve_name(assigned[rank], counters, used_names) for rank in RANKS if rank in assigned}
            row = {"genome_id": node.genome_id, "taxonomy": format_tax(names)}
            (target if node.is_target else reference).append(row)

    return reference, target


# --------------------------------------------------------------------------
# Orchestration / IO
# --------------------------------------------------------------------------
def assign_taxonomy(tree, known_tax, red_cutoffs):
    """Run the full method on a dendropy tree. Returns (reference_rows, target_rows, stats)."""
    leaves = annotate_leaves(tree, known_tax)
    changes = harmonise_higher_ranks(leaves)
    annotate_internal_nodes(tree)
    red_raised = enforce_monotone_red(tree)

    used_names = {rank: set() for rank in RANKS}
    for leaf in leaves:
        for rank, value in leaf.strict.items():
            used_names[rank].add(value)

    reference, target = assign_top_down(tree, red_cutoffs, used_names)
    stats = {
        "leaves": len(leaves),
        "reference": len(reference),
        "target": len(target),
        "harmonised_changes": changes,
        "red_raised_to_parent": red_raised,
        "internal_nodes_without_red": sum(
            1 for n in tree.internal_nodes() if n.red is None
        ),
    }
    return reference, target, stats


def load_red_cutoffs(path):
    with open(path) as handle:
        cutoffs = json.load(handle)
    unknown = set(cutoffs) - set(RANKS)
    if unknown:
        raise ValueError(f"Unknown ranks in RED cutoffs {path}: {sorted(unknown)}")
    return {rank: float(value) for rank, value in cutoffs.items() if value is not None}


def load_known_taxonomy(path):
    frame = pl.read_csv(
        path,
        separator="\t",
        has_header=False,
        new_columns=["genome_id", "taxonomy"],
        schema_overrides={"genome_id": pl.Utf8, "taxonomy": pl.Utf8},
        quote_char=None,
    ).with_columns(pl.col("taxonomy").fill_null(""))
    known = {}
    for genome_id, taxonomy in zip(frame["genome_id"], frame["taxonomy"]):
        if genome_id in known:
            raise ValueError(f"Duplicate genome in taxonomy file {path}: {genome_id}")
        known[genome_id] = taxonomy
    return known


def write_rows(rows, path):
    schema = {"genome_id": pl.Utf8, "taxonomy": pl.Utf8}
    pl.DataFrame(rows, schema=schema).write_csv(path, separator="\t")


def main():
    red_cutoffs = load_red_cutoffs(str(snakemake.input.red_cutoffs))
    known_tax = load_known_taxonomy(str(snakemake.input.taxonomy))
    with open(str(snakemake.input.red_tree)) as handle:
        tree = read_red_tree(handle.read())
    report_negative_branches(tree)

    reference, target, stats = assign_taxonomy(tree, known_tax, red_cutoffs)

    tree_ids = {r["genome_id"] for r in reference} | {r["genome_id"] for r in target}
    missing_from_tree = len(set(known_tax) - tree_ids)
    if missing_from_tree:
        print(f"NOTE: {missing_from_tree} genomes in the taxonomy file are not in the tree.")
    if stats["harmonised_changes"]:
        print(f"Harmonised higher-rank labels (leaf labels changed, by rank): {stats['harmonised_changes']}")
    if stats["red_raised_to_parent"]:
        print(f"NOTE: RED raised to the parent's value on {stats['red_raised_to_parent']} nodes (non-monotone RED).")
    if stats["internal_nodes_without_red"]:
        print(f"WARNING: {stats['internal_nodes_without_red']} internal nodes have no RED value.")

    write_rows(target, str(snakemake.output.target_taxonomy))
    print(f"Wrote {len(target)} target genomes.")
    write_rows(reference, str(snakemake.output.reference_taxonomy))
    print(f"Wrote {len(reference)} reference genomes.")


if __name__ == "__main__":
    main()
