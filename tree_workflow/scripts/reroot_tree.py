#!/usr/bin/env python3
"""Reroot the FastTree tree: on the outgroup's MRCA edge if outgroup accessions are given and
found, else at the midpoint. Needs dendropy, so it runs as a script: rule in the taxonomy env
(a run: block would execute in the driver's Python, where dendropy is not installed)."""
try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:
    redirect_to_snakemake_log = None


def reroot(path, outgroup=()):
    import dendropy
    tree = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True,
                             rooting="force-rooted")
    labels = set(tree.taxon_namespace.labels())
    valid = [acc for acc in outgroup if acc in labels]
    how = "midpoint"
    if len(valid) > 1:
        mrca = tree.mrca(taxon_labels=valid)
        if mrca is not None and mrca.edge is not None:
            length = mrca.edge.length
            half = 0.5 * length if length is not None else None
            tree.reroot_at_edge(mrca.edge, length1=half, length2=half, update_bipartitions=False)
            how = f"outgroup MRCA ({len(valid)} taxa)"
    if how == "midpoint":
        tree.reroot_at_midpoint(update_bipartitions=False)
    return tree, how


def main():
    if redirect_to_snakemake_log:
        redirect_to_snakemake_log()
    tree, how = reroot(str(snakemake.input.tree), list(snakemake.params.outgroup))
    tree.write(path=str(snakemake.output.rerooted_tree), schema="newick", unquoted_underscores=True)
    print(f"reroot_tree: rooted at {how}; {len(tree.leaf_nodes())} tips")


if __name__ == "__main__":
    main()
