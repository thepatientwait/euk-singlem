"""
Assign placeholder taxonomy + novelty rank to genomes with missing/bad taxonomy,
based on where they attach to the reference tree and phylorank's calibrated RED cutoffs.

Requires: dendropy, polars

USAGE NOTES -- update the CONFIG block before running:
  - RED_TREE_PATH / RED_CUTOFFS_PATH must come from the SAME phylorank run
    (midpoint-rooted, trusted-taxa filtered, --skip_rd_refine).
  - TAXONOMY_PATH is your taxonomy.tsv, read-only. This script NEVER writes
    to this path -- see the assertion in main() that enforces this.

======================================================================
PART 1 -- REFERENCE SUPER-CLADE SUBDIVISION (new)
======================================================================
Problem this solves: RED-labeled nodes don't always sit exactly at their
rank's calibrated threshold. A clade like Ascomycota can have a labeled
root whose OWN red is still below the phylum cutoff, with several
independent branches diverging *before* any of them individually cross
that threshold. If left as-is, any target genome anchoring against the
single "Ascomycota" label could get compared against the wrong reference
point, and worse, if this same effect happens elsewhere it can fragment
one true monophyletic group into several unrelated "novel phylum" calls.

Fix: for EVERY original rank-label in the tree (at every rank, not just
phylum), walk its descendants looking for the branches where RED first
crosses that rank's cutoff (find_crossing_branches). If there's exactly
one crossing branch, the anchor just gets relocated to the true crossing
point, same name. If there are multiple (e.g. 4, as in the Ascomycota
case), each gets a lettered suffix: Ascomycota_A, Ascomycota_B, etc.

This is run tree-wide, at every rank, for every named clade -- so it can
also subdivide classes, orders, families, genera wherever the same
premature-branching pattern occurs, not just phyla.

MIRROR OPERATION -- COMBINING OVER-SPLIT CLADES (new):
The reverse problem also happens: a named clade's node can already sit
PAST its rank's cutoff (over-resolved), while climbing back up from it,
the parent (and the parent's parent, etc.) is ALSO still past that same
cutoff. If that climb reaches a point shared with another separately
-named clade at the same rank whose own ancestry is ALSO entirely above
the cutoff up to that shared point, then by RED these are really one
entity that ended up with two names (most likely from inconsistent
reference curation) -- they get merged back into one.

Mechanism: for each post-split label at rank R, climb upward from its
node as long as the parent's RED is still > cutoff[R], stopping at the
highest ancestor that still qualifies (find_merge_root). Group all rank-R
labels by this shared ancestor; any group of more than one gets merged
into a single name at that ancestor, and the deeper, now-redundant labels
are removed. Naming: mode of the constituent names; ties broken by
whichever tied candidate belongs to the deepest (highest-RED) node, for
reproducibility (documented here in case you'd rather it be random).

Combine runs AFTER split, per rank, using split's output as its input.
Both passes are single-pass (not iterated against each other).

The result is stored as `node.revised_ranks` (separate from the original
parsed labels) and used for everything downstream: both the revised
REFERENCE taxonomy (written to a new file, see below) and the anchor point
for TARGET genome novelty assignment (Part 2).

OUTPUTS from Part 1: SPLIT_MAP_OUTPUT and COMBINE_MAP_OUTPUT record every
subdivision and merge made, for traceability back to the original names.

NOT implemented (deferred on purpose): overriding a RED-based split/merge
when all of a parent node's siblings actually belong to one already-named
clade. Flagged for a future pass, not built here.

======================================================================
PART 2 -- TARGET GENOME NOVELTY + PLACEHOLDER TAXONOMY (as before)
======================================================================
For each maximal target-only clade, find the nearest ancestor with a
REVISED label (not the raw phylorank label), classify novelty from its
RED, then recurse into the clade's own internal topology to find further
resolvable splits (genus within a novel family, species within a novel
genus, etc.), using the same branch-crossing style detection that this
recursion already had (current_rank state tracked through the descent).

Per earlier decisions: intermediate ranks with no real split evidence are
left BLANK rather than given an invented code. Placeholder numbers are
scoped to their parent placeholder. Species is terminal -- an entire
species-level sub-branch is finalized as one unit.

OUTPUTS:
  - REVISED_REFERENCE_TAXONOMY_OUTPUT: every REFERENCE genome's taxonomy,
    rebuilt from the revised/subdivided labels. Where a genome's lineage
    was already clean (no subdivision needed anywhere), this matches its
    original taxonomy.tsv entry exactly.
  - TARGET_NOVELTY_OUTPUT: genome_id, placeholder taxonomy, novelty rank
    for every target genome.
"""

import json
import re
from collections import Counter, defaultdict

import dendropy
import polars as pl

# ---------------- CONFIG ----------------
RED_TREE_PATH = "red_midpoint_trusted_norefine/tree.midpoint.trusted.norefine.decorated.red_decorated.tree"
RED_CUTOFFS_PATH = "red_midpoint_trusted_norefine/tree.midpoint.trusted.norefine.decorated.dict"
TAXONOMY_PATH = "taxonomy.tsv"  # READ-ONLY, never written to
REVISED_REFERENCE_TAXONOMY_OUTPUT = "taxonomy.revised_reference.tsv"
TARGET_NOVELTY_OUTPUT = "target_novelty_taxonomy.tsv"
SPLIT_MAP_OUTPUT = "clade_split_map.tsv"
COMBINE_MAP_OUTPUT = "clade_combine_map.tsv"
PARTIAL_TAXONOMY_MAP_OUTPUT = "partial_taxonomy_recovery_map.tsv"
# -----------------------------------------

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX = {
    "domain": "d__", "phylum": "p__", "class": "c__", "order": "o__",
    "family": "f__", "genus": "g__", "species": "s__",
}
RANK_LETTER = {"phylum": "p", "class": "c", "order": "o", "family": "f", "genus": "g"}
CUTOFF_RANKS = ["phylum", "class", "order", "family", "genus", "species"]

RED_RE = re.compile(r"RED_([0-9.]+)")
TAX_RE = re.compile(r":?((?:d__|p__|c__|o__|f__|g__|s__)[^|]*)")
SUPPORT_RE = re.compile(r"^'?([0-9.]+)")

RED_CUTOFFS = {}  # populated in main()


def load_red_cutoffs(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def classify_rank(red_value: float, red_cutoffs: dict):
    for rank in CUTOFF_RANKS:
        if rank in red_cutoffs and red_value <= red_cutoffs[rank]:
            return rank
    return "species"


def rank_above(rank: str) -> str:
    idx = RANKS.index(rank)
    return RANKS[idx - 1] if idx > 0 else "domain"


def is_deeper(rank_a: str, rank_b: str) -> bool:
    return RANKS.index(rank_a) > RANKS.index(rank_b)


def parse_node_label(raw_label):
    """Extract RED value, the raw (possibly multi-rank) taxonomy string, and support value."""
    if raw_label is None:
        return None, None, None
    red_match = RED_RE.search(raw_label)
    red = float(red_match.group(1)) if red_match else None
    tax_match = TAX_RE.search(raw_label)
    taxon = tax_match.group(1).strip() if tax_match else None
    support_match = SUPPORT_RE.match(raw_label)
    support = float(support_match.group(1)) if support_match else None
    return red, taxon, support


def extract_rank_chunks(raw_taxon_string):
    """Split a possibly multi-rank bundled label ('o__X; f__Y; g__Z') into
    individual (rank, name) pairs."""
    if not raw_taxon_string:
        return []
    chunks = []
    for chunk in raw_taxon_string.split(";"):
        chunk = chunk.strip()
        for rank, prefix in RANK_PREFIX.items():
            if chunk.startswith(prefix):
                name = chunk[len(prefix):].strip()
                if name:
                    chunks.append((rank, name))
                break
    return chunks


def genome_id_from_leaf_label(raw_label):
    if raw_label and "|RED_" in raw_label:
        return raw_label.split("|RED_")[0].replace("'", "")
    return raw_label.replace("'", "") if raw_label else ""


def suffix_letter(n: int) -> str:
    """0-indexed -> 'A', 'B', ..., 'Z', 'AA', 'AB', ..."""
    n += 1
    letters = ""
    while n > 0:
        n, rem = divmod(n - 1, 26)
        letters = chr(65 + rem) + letters
    return letters


# ---------------- PART 1: reference super-clade subdivision ----------------

def find_crossing_branches(node, cutoff_red):
    """Return the nodes where RED first exceeds `cutoff_red` along every
    descendant path from `node` (inclusive). One result = clean single
    crossing (possibly just relocated deeper than `node` itself). Multiple
    results = premature branching before the threshold -- needs subdivision."""
    if node.red is not None and node.red > cutoff_red:
        return [node]
    results = []
    for child in node.child_nodes():
        results.extend(find_crossing_branches(child, cutoff_red))
    return results


def build_revised_labels(tree, red_cutoffs):
    """Populate node.revised_ranks for every node, by re-deriving each
    original label's true crossing point(s) at its own rank's cutoff."""
    for node in tree.preorder_node_iter():
        node.revised_ranks = {}

    original_labels = []
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        for rank, name in extract_rank_chunks(getattr(node, "node_taxon", None)):
            original_labels.append((node, rank, name))

    subdivision_log = []
    for orig_node, rank, name in original_labels:
        cutoff = red_cutoffs.get(rank)
        if cutoff is None:
            continue
        crossing_nodes = find_crossing_branches(orig_node, cutoff)
        if len(crossing_nodes) == 1:
            crossing_nodes[0].revised_ranks[rank] = name
        else:
            new_names = []
            for i, cn in enumerate(crossing_nodes):
                new_name = f"{name}_{suffix_letter(i)}"
                cn.revised_ranks[rank] = new_name
                new_names.append(new_name)
            subdivision_log.append({
                "rank": rank,
                "original_name": name,
                "new_names": ";".join(new_names),
                "num_new_clades": len(new_names),
            })

    return subdivision_log


def find_merge_root(node, cutoff_red):
    """Climb upward from `node` as long as the PARENT's RED is still past
    `cutoff_red`, stopping at the highest ancestor that still qualifies.
    This is the mirror of find_crossing_branches: instead of finding where
    a lineage first crosses INTO past-threshold territory going down, this
    finds the shallowest point from which the whole subtree is guaranteed
    to already be past the threshold."""
    current = node
    while (current.parent_node is not None
           and current.parent_node.red is not None
           and current.parent_node.red > cutoff_red):
        current = current.parent_node
    return current


def apply_combine_pass(tree, red_cutoffs):
    """For each rank, group all currently-revised labels by their shared
    merge root (see find_merge_root). Groups of >1 get merged into a single
    name at that root; the deeper, now-redundant labels are removed."""
    combine_log = []

    for rank in CUTOFF_RANKS:
        cutoff = red_cutoffs.get(rank)
        if cutoff is None:
            continue

        labeled_nodes = [
            node for node in tree.preorder_node_iter()
            if not node.is_leaf() and rank in getattr(node, "revised_ranks", {})
        ]

        groups = defaultdict(list)
        for node in labeled_nodes:
            root = find_merge_root(node, cutoff)
            groups[root].append(node)

        for root, members in groups.items():
            if len(members) <= 1:
                continue

            original_names = [m.revised_ranks[rank] for m in members]
            name_counts = Counter(original_names)
            max_count = max(name_counts.values())
            candidates = [name for name, cnt in name_counts.items() if cnt == max_count]

            if len(candidates) == 1:
                chosen_name = candidates[0]
            else:
                # tie-break: deepest (highest RED) node among tied candidates -- deterministic,
                # not random, so reruns are reproducible. Swap for random.choice if preferred.
                tied_members = [m for m in members if m.revised_ranks[rank] in candidates]
                deepest_member = max(tied_members, key=lambda m: m.red if m.red is not None else -1)
                chosen_name = deepest_member.revised_ranks[rank]

            for m in members:
                if rank in m.revised_ranks:
                    del m.revised_ranks[rank]
            root.revised_ranks[rank] = chosen_name

            genome_ids = sorted({l.genome_id for l in root.leaf_iter()})
            combine_log.append({
                "rank": rank,
                "merged_name": chosen_name,
                "original_names": ";".join(sorted(set(original_names))),
                "num_merged_clades": len(members),
                "genome_count": len(genome_ids),
                "genome_ids": ";".join(genome_ids),
                "example_genomes": ",".join(genome_ids[:3]),
            })

    return combine_log


def find_nearest_revised_labeled_ancestor(node):
    current = node
    while current is not None and not getattr(current, "revised_ranks", {}):
        current = current.parent_node
    return current


def find_nearest_labeled_ancestor(node):
    current = node
    while current is not None and not getattr(current, "node_taxon", None):
        current = current.parent_node
    return current


def collect_revised_lineage(node):
    lineage = {}
    current = node
    while current is not None:
        raw_taxon = getattr(current, "node_taxon", None)
        if raw_taxon:
            for chunk in raw_taxon.split(";"):
                chunk = chunk.strip()
                for rank, prefix in RANK_PREFIX.items():
                    if chunk.startswith(prefix) and rank not in lineage:
                        lineage[rank] = chunk[len(prefix):].strip()
        for rank, name in getattr(current, "revised_ranks", {}).items():
            lineage[rank] = name
        current = current.parent_node
    return lineage


def collect_raw_lineage(node):
    lineage = {}
    current = node
    while current is not None:
        taxon = getattr(current, "node_taxon", None)
        if taxon:
            for chunk in taxon.split(";"):
                chunk = chunk.strip()
                for rank, prefix in RANK_PREFIX.items():
                    if chunk.startswith(prefix) and rank not in lineage:
                        lineage[rank] = chunk[len(prefix):].strip()
        current = current.parent_node
    return lineage


def parse_taxonomy_dict(tax_string: str):
    if not tax_string:
        return {}
    parsed = {}
    for chunk in tax_string.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        for rank, prefix in RANK_PREFIX.items():
            if chunk.startswith(prefix):
                value = chunk[len(prefix):].strip()
                if value:
                    parsed[rank] = value
                break
    return parsed


def clade_group_key(tax_dict: dict, missing_rank: str | None = None):
    ranks = [rank for rank in RANKS if rank in tax_dict and str(tax_dict[rank]).strip()]
    if missing_rank is not None:
        idx = RANKS.index(missing_rank)
        deeper = [rank for rank in ranks if RANKS.index(rank) > idx]
        shallower = [rank for rank in ranks if RANKS.index(rank) < idx]
        # The nearest known downstream rank identifies the fallback clade.
        # Including every deeper rank would split one clade by species name.
        ranks = shallower + deeper[:1]
    return tuple((rank, str(tax_dict[rank]).strip()) for rank in ranks)


def partial_clade_signature(tax_dict: dict):
    return tuple(
        (rank, str(value).strip())
        for rank in RANKS
        if rank in tax_dict and str(tax_dict[rank]).strip()
    )


def synthesize_placeholder_for_rank(rank: str, parent_name: str, counters: defaultdict, group_key=None):
    base = parent_name or "unk"
    if group_key is not None:
        stable_key = ("partial_group", rank, tuple(group_key))
        if stable_key not in counters:
            index_key = ("partial_group_index", rank, base)
            counters[index_key] += 1
            counters[stable_key] = counters[index_key]
        return f"{base}_{RANK_LETTER.get(rank, rank)}{counters[stable_key]}"

    counters[(rank, base)] += 1
    return f"{base}_{RANK_LETTER.get(rank, rank)}{counters[(rank, base)]}"


def find_matching_reference_lineage(tree, tax_dict: dict):
    """Return the unique revised lineage matching every known rank."""
    known = {
        rank: value for rank, value in tax_dict.items()
        if rank in RANKS and str(value).strip()
    }
    if not known:
        return None

    matches = []
    for leaf in tree.leaf_node_iter():
        lineage = collect_revised_lineage(leaf)
        if all(lineage.get(rank) == value for rank, value in known.items()):
            matches.append(lineage)

    if not matches:
        return None

    # Several leaves can represent one clade. Keep the match only when all
    # matching leaves agree on the missing-rank assignments.
    signatures = {
        tuple((rank, lineage.get(rank)) for rank in RANKS if rank not in known)
        for lineage in matches
    }
    if len(signatures) != 1:
        return None
    return matches[0]


def partial_taxonomy_matches_reference_clade(tree, tax_dict: dict):
    """Return whether all known ranks match an existing revised lineage."""
    return find_matching_reference_lineage(tree, tax_dict) is not None


def repair_partial_taxonomy_string(
    tax_string: str,
    counters: defaultdict | None = None,
    tree=None,
    reference_lineage=None,
):
    if not tax_string:
        return tax_string
    if counters is None:
        counters = defaultdict(int)

    parsed = parse_taxonomy_dict(tax_string)
    if not parsed:
        return tax_string

    fixed = {}
    last_real_name = ""
    for rank in RANKS:
        if rank in parsed:
            fixed[rank] = parsed[rank]
            last_real_name = parsed[rank]
            continue

        deeper_present = any(other_rank in parsed for other_rank in RANKS[RANKS.index(rank) + 1:])
        if deeper_present:
            subgroup = clade_group_key(parsed, missing_rank=rank)
            matched_lineage = reference_lineage or (
                find_matching_reference_lineage(tree, parsed)
                if tree is not None else None
            )
            if (
                matched_lineage is not None
                and rank in matched_lineage
                and all(
                    matched_lineage.get(known_rank) == value
                    for known_rank, value in parsed.items()
                    if matched_lineage.get(known_rank) is not None
                )
            ):
                fixed[rank] = matched_lineage[rank]
            else:
                fixed[rank] = synthesize_placeholder_for_rank(
                    rank,
                    last_real_name or "unk",
                    counters,
                    group_key=subgroup,
                )
            last_real_name = fixed[rank]

    ordered = []
    for rank in RANKS:
        if rank in fixed:
            ordered.append(f"{RANK_PREFIX[rank]}{fixed[rank]}")

    return ';'.join(ordered)


def fill_missing_taxonomy_ranks(tax_dict: dict, counters: defaultdict | None = None):
    if counters is None:
        counters = defaultdict(int)

    repaired = dict(tax_dict)
    last_real_name = ""
    for rank in RANKS:
        if rank in repaired:
            last_real_name = str(repaired[rank]).strip()
            continue
        deeper_present = any(other_rank in repaired for other_rank in RANKS[RANKS.index(rank) + 1:])
        if deeper_present:
            subgroup = clade_group_key(repaired, missing_rank=rank)
            repaired[rank] = synthesize_placeholder_for_rank(
                rank,
                last_real_name or "unk",
                counters,
                group_key=subgroup,
            )
            last_real_name = str(repaired[rank]).strip()
    return repaired


def has_partial_taxonomy(tax_string: str):
    if not tax_string:
        return False
    rank_values = parse_taxonomy_dict(tax_string)
    if not rank_values:
        return False
    chunks = [c.strip() for c in tax_string.split(";") if c.strip()]
    return any(
        any(chunk.startswith(prefix) and chunk[len(prefix):].strip() == "" for prefix in RANK_PREFIX.values())
        for chunk in chunks
    )


def recover_partial_taxonomy_map(known_taxonomy: dict, tree=None):
    rows = []
    counters = defaultdict(int)
    tree_lineages = {
        leaf.genome_id: collect_revised_lineage(leaf)
        for leaf in tree.leaf_node_iter()
        if getattr(leaf, "genome_id", None)
    } if tree is not None else {}
    for gid, tax_string in known_taxonomy.items():
        if not tax_string or not has_partial_taxonomy(tax_string):
            continue
        rows.append({
            "genome_id": gid,
            "original_taxonomy": tax_string,
            "recovered_taxonomy": repair_partial_taxonomy_string(
                tax_string,
                counters,
                tree=tree,
                reference_lineage=tree_lineages.get(gid),
            ),
        })
    return rows


def build_reference_taxonomy_rows(tree, known_taxonomy, partial_counters):
    ref_rows = []
    for leaf in tree.leaf_node_iter():
        if leaf.is_target:
            continue
        lineage = collect_revised_lineage(leaf)
        original_tax = known_taxonomy.get(leaf.genome_id, "")
        if original_tax:
            partial = has_partial_taxonomy(original_tax)
            if partial:
                recovered_tax = repair_partial_taxonomy_string(
                    original_tax,
                    partial_counters,
                    tree=tree,
                    reference_lineage=lineage,
                )
                ref_rows.append({"genome_id": leaf.genome_id, "taxonomy": recovered_tax})
                continue
            for chunk in original_tax.split(";"):
                chunk = chunk.strip()
                for rank, prefix in RANK_PREFIX.items():
                    if chunk.startswith(prefix) and rank not in lineage:
                        lineage[rank] = chunk[len(prefix):].strip()
        repaired_lineage = fill_missing_taxonomy_ranks(lineage, partial_counters)
        parts = [
            f"{RANK_PREFIX[rank]}{repaired_lineage[rank]}"
            for rank in RANKS
            if rank in repaired_lineage and str(repaired_lineage.get(rank, '')).strip()
        ]
        ref_rows.append({"genome_id": leaf.genome_id, "taxonomy": ";".join(parts)})
    return ref_rows


# ---------------- PART 2: target genome novelty + placeholder taxonomy ----------------

def dict_to_tax_string(tax_dict, counters=None):
    repaired = fill_missing_taxonomy_ranks(tax_dict, counters)
    return ';'.join(
        f"{RANK_PREFIX[rank]}{repaired[rank]}"
        for rank in RANKS
        if rank in repaired and str(repaired[rank]).strip()
    )


def get_deepest_rank(tax_dict):
    for rank in reversed(RANKS):
        if rank in tax_dict:
            return rank
    return "domain"


def extend_taxonomy_with_placeholders(base_tax, inherited_rank, target_rank, counters):
    if RANKS.index(target_rank) <= RANKS.index(inherited_rank):
        return base_tax

    base_name = base_tax.get(inherited_rank, "unk")
    start_idx = RANKS.index(inherited_rank) + 1
    end_idx = RANKS.index(target_rank)

    prefix_chain = base_name
    for i in range(start_idx, end_idx + 1):
        rank = RANKS[i]
        # Each new branch should get a fresh placeholder chain; reusing the same
        # counter state for earlier ranks collapses distinct clades onto the same
        # phylum/class placeholders.
        key = (inherited_rank, prefix_chain, rank)
        counters[key] += 1

        if rank == "species":
            base_tax[rank] = f"{base_name} sp_{counters[key]}"
            prefix_chain = base_tax[rank]
        else:
            code = f"{RANK_LETTER[rank]}{counters[key]}"
            if i == start_idx:
                prefix_chain = f"{base_name}_{code}"
            else:
                prefix_chain = f"{prefix_chain}_{code}"
            base_tax[rank] = prefix_chain

    return base_tax


def assign_clade_taxonomy(clade_node, base_tax, novelty_rank, counters, results):
    for node in clade_node.preorder_iter():
        node_base = base_tax.copy() if node is clade_node else node.parent_node.assigned_tax.copy()
        inherited_rank = get_deepest_rank(node_base)
        target_rank = classify_rank(node.red, RED_CUTOFFS)

        node.assigned_tax = extend_taxonomy_with_placeholders(node_base, inherited_rank, target_rank, counters)

        if node.is_leaf():
            final_tax = node.assigned_tax.copy()
            deepest = get_deepest_rank(final_tax)

            if deepest != "species":
                base_name = final_tax.get(deepest, "unk")
                key = (deepest, base_name, "species")
                counters[key] += 1
                if base_name == "unk":
                    final_tax["species"] = f"sp_{counters[key]}"
                else:
                    final_tax["species"] = f"{base_name} sp_{counters[key]}"

            results.append({
                "genome_id": node.genome_id,
                "taxonomy": dict_to_tax_string(final_tax),
                "novelty": f"novel {novelty_rank}",
            })


def assign_target_clade(clade_node, entry_rank, anchor_rank, anchor_name, lineage,
                        counters, results):
    base_tax = {}
    for rank in RANKS:
        if rank in lineage:
            base_tax[rank] = lineage[rank]
        if rank == anchor_rank:
            break

    if entry_rank == "species":
        assign_clade_taxonomy(clade_node, base_tax, entry_rank, counters, results)
        return

    cutoff = RED_CUTOFFS.get(entry_rank)
    if cutoff is not None:
        crossing_nodes = []
        for child in clade_node.child_nodes():
            crossing_nodes.extend(find_crossing_branches(child, cutoff))
        if len(crossing_nodes) > 1:
            for branch in set(crossing_nodes):
                assign_clade_taxonomy(branch, base_tax.copy(), entry_rank, counters, results)
            return

    assign_clade_taxonomy(clade_node, base_tax, entry_rank, counters, results)


def main():
    global RED_CUTOFFS

    if "snakemake" in globals():
        RED_TREE_PATH = str(snakemake.input.red_tree)
        RED_CUTOFFS_PATH = str(snakemake.input.red_cutoffs)
        TAXONOMY_PATH = str(snakemake.input.taxonomy)
        REVISED_REFERENCE_TAXONOMY_OUTPUT = str(snakemake.output.reference_taxonomy)
        TARGET_NOVELTY_OUTPUT = str(snakemake.output.target_taxonomy)

    RED_CUTOFFS = load_red_cutoffs(RED_CUTOFFS_PATH)
    print(f"Loaded RED cutoffs: {RED_CUTOFFS}")

    assert REVISED_REFERENCE_TAXONOMY_OUTPUT != TAXONOMY_PATH, \
        "Refusing to write revised taxonomy over the original input file"

    taxonomy_df = pl.read_csv(
        TAXONOMY_PATH, separator="\t", has_header=False, new_columns=["genome_id", "taxonomy"]
    )
    known_taxonomy = dict(zip(taxonomy_df["genome_id"], taxonomy_df["taxonomy"]))

    def is_target_genome(gid: str) -> bool:
        tax_string = known_taxonomy.get(gid)
        return tax_string is None or tax_string.strip() == ""

    with open(RED_TREE_PATH) as fh:
        tree_str = fh.read().replace("|RED=", "|RED_")

    tree = dendropy.Tree.get(
        data=tree_str,
        schema="newick",
        preserve_underscores=True,
        suppress_internal_node_taxa=True,
        suppress_leaf_node_taxa=True,
    )

    for leaf in tree.leaf_node_iter():
        raw_label = leaf.label or (leaf.taxon.label if leaf.taxon else None)
        gid = genome_id_from_leaf_label(raw_label)
        leaf.genome_id = gid
        leaf.is_target = is_target_genome(gid)
        leaf.red = 1.000
        leaf.node_taxon = None
        leaf.support = None

    for node in tree.postorder_node_iter():
        if node.is_leaf():
            continue
        red, taxon, support = parse_node_label(node.label)
        node.red = red
        node.node_taxon = taxon
        node.support = support
        node.is_target = all(getattr(c, "is_target", False) for c in node.child_nodes())

    # ---- PART 1: subdivide over-large reference clades according to RED ----
    subdivision_log = build_revised_labels(tree, RED_CUTOFFS)
    print(f"\n{len(subdivision_log)} named clade(s) required subdivision:")
    for entry in subdivision_log[:20]:
        print(f"  [{entry['rank']}] {entry['original_name']} -> {entry['new_names']}")
    if len(subdivision_log) > 20:
        print(f"  ... and {len(subdivision_log) - 20} more")
    if subdivision_log:
        print(f"Recorded {len(subdivision_log)} subdivision event(s)")

    # ---- PART 2: repair/anchor partial taxon strings in the revised clade structure ----
    partial_counters = defaultdict(int)
    partial_map = recover_partial_taxonomy_map(known_taxonomy, tree=tree)
    if partial_map:
        print(f"Recorded {len(partial_map)} partial-taxonomy recovery record(s)")

    # ---- PART 3: consolidate clades that are too small according to RED ----
    combine_log = apply_combine_pass(tree, RED_CUTOFFS)
    print(f"\n{len(combine_log)} clade group(s) required combining:")
    for entry in combine_log[:20]:
        print(f"  [{entry['rank']}] {entry['original_names']} -> {entry['merged_name']} "
              f"({entry['num_merged_clades']} clades, {entry['genome_count']} genomes)")
    if len(combine_log) > 20:
        print(f"  ... and {len(combine_log) - 20} more")
    if combine_log:
        print(f"Recorded {len(combine_log)} combine event(s)")

    # ---- write revised reference taxonomy (new file, original untouched) ----
    ref_rows = build_reference_taxonomy_rows(tree, known_taxonomy, partial_counters)
    ref_df = pl.DataFrame(ref_rows)
    ref_df.write_csv(REVISED_REFERENCE_TAXONOMY_OUTPUT, separator="\t")
    print(f"\nWrote revised taxonomy for {len(ref_rows)} reference genome(s) to {REVISED_REFERENCE_TAXONOMY_OUTPUT}")

    # ---- PART 4: annotate leaves with no taxon strings using the final revised structure ----
    maximal_clades = []
    for node in tree.preorder_node_iter():
        if getattr(node, "is_target", False):
            parent = node.parent_node
            if parent is None or not getattr(parent, "is_target", False):
                maximal_clades.append(node)

    print(f"\nFound {len(maximal_clades)} maximal target clade(s)")

    counters = defaultdict(int)
    results = []

    for clade_node in maximal_clades:
        genome_ids_preview = [l.genome_id for l in clade_node.leaf_iter()][:3]

        stem_start = clade_node.parent_node
        if stem_start is None:
            print(f"WARNING: clade {genome_ids_preview}... has no parent (tree root) -- skipping")
            continue

        labeled_node = find_nearest_revised_labeled_ancestor(stem_start)
        if labeled_node is None:
            labeled_node = find_nearest_labeled_ancestor(stem_start)
            if labeled_node is None:
                entry_rank = "phylum"
                anchor_rank = "domain"
                lineage = {}
                anchor_name = ""
                print(f"WARNING: no revised or raw labeled ancestor found for clade {genome_ids_preview}... using root-level fallback")
            else:
                entry_red = getattr(labeled_node, "red", None)
                if entry_red is None:
                    entry_rank = "phylum"
                    anchor_rank = "domain"
                    lineage = {}
                    anchor_name = ""
                    print(f"WARNING: nearest raw labeled ancestor has no RED for clade {genome_ids_preview}... using root-level fallback")
                else:
                    entry_rank = classify_rank(entry_red, RED_CUTOFFS)
                    anchor_rank = rank_above(entry_rank)
                    lineage = collect_raw_lineage(labeled_node)
                    anchor_name = lineage.get(anchor_rank) or lineage.get("domain") or ""
        else:
            entry_red = getattr(labeled_node, "red", None)
            if entry_red is None:
                print(f"WARNING: nearest revised labeled ancestor has no RED for clade {genome_ids_preview}... -- skipping")
                continue

            entry_rank = classify_rank(entry_red, RED_CUTOFFS)
            anchor_rank = rank_above(entry_rank)
            lineage = collect_revised_lineage(labeled_node)
            anchor_name = lineage.get(anchor_rank) or "Unclassified"

        base_tax = {}
        for rank in RANKS:
            if rank in lineage:
                base_tax[rank] = lineage[rank]
            if rank == anchor_rank:
                break

        novelty_rank = entry_rank if entry_rank != "species" else "species"
        assign_target_clade(clade_node, entry_rank, anchor_rank, anchor_name, lineage,
                            counters, results)

    result_df = pl.DataFrame(results)
    result_df.write_csv(TARGET_NOVELTY_OUTPUT, separator="\t")
    print(f"Wrote {len(results)} genome(s) across {len(maximal_clades)} top-level clade(s) to {TARGET_NOVELTY_OUTPUT}")


if __name__ == "__main__":
    main()
