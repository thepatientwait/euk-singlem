#!/usr/bin/env python3

import json
import re
from collections import defaultdict, Counter
import dendropy
import polars as pl

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = {r: f"{r[0]}__" for r in RANKS}


def parse_tax(tax_str):
    if not tax_str: 
        return {}
    res = {}
    for chunk in [c.strip() for c in tax_str.split(";")]:
        for r, p in PREFIXES.items():
            if chunk.startswith(p):
                res[r] = chunk[len(p):].strip()
    return res


def format_tax(tax_dict):
    return ";".join(f"{PREFIXES[r]}{tax_dict.get(r, '')}" for r in RANKS)


def get_nearest_ancestor(tax_dict, current_rank):
    """Finds the closest populated taxonomic rank ABOVE the current one."""
    idx = RANKS.index(current_rank)
    for i in range(idx - 1, -1, -1):
        anc = RANKS[i]
        if tax_dict.get(anc):
            return anc, tax_dict[anc]
    return None, "unk"


def main():
    red_tree_path = str(snakemake.input.red_tree)
    red_cutoffs_path = str(snakemake.input.red_cutoffs)
    taxonomy_path = str(snakemake.input.taxonomy)
    revised_reference_output = str(snakemake.output.reference_taxonomy)
    target_novelty_output = str(snakemake.output.target_taxonomy)

    with open(red_cutoffs_path) as fh:
        red_cutoffs = json.load(fh)

    tax_df = pl.read_csv(
        taxonomy_path, separator="\t", has_header=False, new_columns=["genome_id", "taxonomy"]
    ).with_columns(pl.col("taxonomy").fill_null(""))
    known_tax = dict(zip(tax_df["genome_id"], tax_df["taxonomy"]))

    with open(red_tree_path) as fh:
        tree = dendropy.Tree.get(data=fh.read().replace("|RED=", "|RED_"), schema="newick", preserve_underscores=True)

    # Initialize leaves and pre-collect leaf lists per rank for instant access
    leaves_by_rank_val = {r: defaultdict(list) for r in RANKS}
    for node in tree.leaf_node_iter():
        raw_label = node.label or (node.taxon.label if node.taxon else "")
        node.genome_id = raw_label.split("|RED_")[0].replace("'", "")
        node.node_taxon = known_tax.get(node.genome_id, "")
        node.is_target = not node.node_taxon.strip()
        node.red = 1.000
        node.strict_tax = parse_tax(node.node_taxon) if not node.is_target else {}
        
        if not node.is_target:
            for r, val in node.strict_tax.items():
                if val:
                    leaves_by_rank_val[r][val].append(node)

    # Helper to find all leaves in a node's subtree
    def get_subtree_leaves(node):
        if node.is_leaf():
            return [node]
        res = []
        for child in node.child_nodes():
            res.extend(get_subtree_leaves(child))
        return res

    # --- OPTIMIZED PRE-PROCESSING PASS 1: Detect Polyphyly and Split (Edge Case 2) ---
    for rank in RANKS:
        for T, leaves in list(leaves_by_rank_val[rank].items()):
            if len(leaves) <= 1:
                continue
            
            pure_clades = []
            def check_purity(n):
                if n.is_leaf():
                    if not n.is_target:
                        v = n.strict_tax.get(rank, "")
                        return (v == T), (v != "" and v != T)
                    return True, False
                
                all_match = True
                has_conflict = False
                for child in n.child_nodes():
                    m, c = check_purity(child)
                    if not m:
                        all_match = False
                    if c:
                        has_conflict = True
                
                sub_leaves = [l for l in get_subtree_leaves(n) if not l.is_target]
                if sub_leaves:
                    vals = {l.strict_tax.get(rank, "") for l in sub_leaves if l.strict_tax.get(rank, "")}
                    if vals == {T}:
                        pure_clades.append(n)
                        return True, False
                return all_match and not has_conflict, has_conflict

            check_purity(tree.seed_node)
            pure_set = {id(n) for n in pure_clades}
            maximal_clades = [
                n for n in pure_clades 
                if n.parent_node is None or id(n.parent_node) not in pure_set
            ]
            
            if len(maximal_clades) > 1:
                for idx, clade_node in enumerate(maximal_clades, start=1):
                    new_name = f"{T}_{idx}"
                    clade_leaves = [l for l in get_subtree_leaves(clade_node) if not l.is_target and l.strict_tax.get(rank) == T]
                    for leaf in clade_leaves:
                        leaf.strict_tax[rank] = new_name
                        leaves_by_rank_val[rank][new_name].append(leaf)
                del leaves_by_rank_val[rank][T]

    # --- OPTIMIZED PRE-PROCESSING PASS 2: Universal Higher-Rank Consistency (Edge Case 1) ---
    for i, low_rank in enumerate(RANKS[1:], start=1):
        for low_val in list(leaves_by_rank_val[low_rank].keys()):
            low_leaves = leaves_by_rank_val[low_rank][low_val]
            if not low_leaves:
                continue
            
            pure_clades = []
            def check_low_purity(n):
                sub_leaves = [l for l in get_subtree_leaves(n) if not l.is_target]
                if sub_leaves:
                    vals = {l.strict_tax.get(low_rank, "") for l in sub_leaves if l.strict_tax.get(low_rank, "")}
                    if vals == {low_val}:
                        pure_clades.append(n)
                        return True
                for child in n.child_nodes():
                    check_low_purity(child)
                return False

            check_low_purity(tree.seed_node)
            pure_set = {id(n) for n in pure_clades}
            maximal_clades = [
                n for n in pure_clades 
                if n.parent_node is None or id(n.parent_node) not in pure_set
            ]
            
            for clade_node in maximal_clades:
                clade_leaves = [l for l in get_subtree_leaves(clade_node) if not l.is_target]
                for high_rank_idx in range(i):
                    high_rank = RANKS[high_rank_idx]
                    high_vals = [leaf.strict_tax.get(high_rank, "") for leaf in clade_leaves if leaf.strict_tax.get(high_rank, "")]
                    if high_vals:
                        majority_val, _ = Counter(high_vals).most_common(1)[0]
                        for leaf in clade_leaves:
                            leaf.strict_tax[high_rank] = majority_val

    # 1. Bottom-Up: Annotate nodes based STRICTLY on Monophyly
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            continue
            
        red_match = re.search(r"RED_([0-9.]+)", node.label or "")
        node.red = float(red_match.group(1)) if red_match else None
        node.is_target = all(getattr(c, "is_target", False) for c in node.child_nodes())
        
        child_taxes = [c.strict_tax for c in node.child_nodes() if getattr(c, "strict_tax", {})]
        
        node.strict_tax = {}
        if child_taxes:
            for rank in RANKS:
                vals = set()
                for ct in child_taxes:
                    if ct.get(rank):
                        vals.add(ct[rank])
                
                if len(vals) == 1:
                    node.strict_tax[rank] = vals.pop()
                elif len(vals) > 1:
                    break
                else:
                    pass

    # 2. Top-Down: Assign taxonomy with RED Truncation & Placeholder Generation
    counters = defaultdict(int)
    target_results = []
    ref_results = []

    tree.seed_node.assigned_tax = getattr(tree.seed_node, "strict_tax", {})
    if "domain" not in tree.seed_node.assigned_tax or not tree.seed_node.assigned_tax["domain"]:
        tree.seed_node.assigned_tax["domain"] = "Eukarya"

    for node in tree.preorder_node_iter():
        if node != tree.seed_node:
            parent_tax = getattr(node.parent_node, "assigned_tax", {})
            node.assigned_tax = parent_tax.copy()

        if hasattr(node, "strict_tax"):
            for r in RANKS:
                if r in node.strict_tax and node.strict_tax[r]:
                    node.assigned_tax[r] = node.strict_tax[r]

        is_leaf = node.is_leaf()
        node_red = getattr(node, "red", None)

        # RED Cutoff Check: truncate inherited ranks below the node's evolutionary divergence level
        if node_red is not None and red_cutoffs:
            allowed_ranks = []
            for rank in RANKS:
                cutoff = red_cutoffs.get(rank)
                if cutoff is not None and node_red > cutoff:
                    allowed_ranks.append(rank)
            
            if allowed_ranks:
                max_allowed_rank = allowed_ranks[-1]
                max_idx = RANKS.index(max_allowed_rank)
                for r in RANKS[max_idx + 1:]:
                    node.assigned_tax[r] = ""

        for rank in RANKS:
            if rank in node.assigned_tax and node.assigned_tax[rank]:
                continue
            
            cutoff = red_cutoffs.get(rank)
            
            if is_leaf or (node_red is not None and cutoff is not None and node_red > cutoff):
                ancestor_rank, prefix = get_nearest_ancestor(node.assigned_tax, rank)
                
                key = (ancestor_rank, prefix, rank)
                counters[key] += 1
                
                if rank == "species":
                    val = f"{prefix} sp_{counters[key]}"
                else:
                    code = f"{rank[0]}{counters[key]}"
                    val = f"{prefix}_{code}" if prefix != "unk" else code
                    
                node.assigned_tax[rank] = val

        if is_leaf:
            formatted = format_tax(node.assigned_tax)
            if getattr(node, "is_target", False):
                target_results.append({
                    "genome_id": node.genome_id,
                    "taxonomy": formatted
                })
            else:
                ref_results.append({
                    "genome_id": node.genome_id,
                    "taxonomy": formatted
                })

    pl.DataFrame(target_results).write_csv(target_novelty_output, separator="\t", include_header=False)
    print(f"Wrote {len(target_results)} target genome placeholders.")

    pl.DataFrame(ref_results).write_csv(revised_reference_output, separator="\t", include_header=False)
    print(f"Wrote {len(ref_results)} reference genomes.")


if __name__ == "__main__":
    main()