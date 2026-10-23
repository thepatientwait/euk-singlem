"""Computations behind the tree diagnostics notebook (notebooks/tree_diagnostics.ipynb).

Questions answered, all from the concatenated alignment, its partitions and the taxonomy:

* occupancy: how many genes does each genome have, how gappy is its row?
* gene overlap: do genomes of the same clade share the same genes (more than genomes of
  different clades)? If not, the clade is held together by too little common data.
* divergence: are same-clade sequences closer than different-clade ones, per rank and per
  gene (``separation`` = probability that a within-clade pair is closer than a between-clade pair)?
* multicopy: how good is the "top-scoring copy" heuristic for Duplicated BUSCO genes?
"""
import random
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from tree_fragmentation import RANKS, read_taxonomy  # noqa: F401  (re-exported)

GAP = ord("-")


# ----------------------------------------------------------------------------- alignment
def load_alignment(path):
    """(ids, uint8 matrix [genomes x columns]) of an aligned FASTA."""
    ids, rows = [], []
    with open(path, "r", encoding="utf-8") as handle:
        name, parts = None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    ids.append(name)
                    rows.append("".join(parts))
                name, parts = line[1:].split()[0], []
            else:
                parts.append(line.strip())
        if name is not None:
            ids.append(name)
            rows.append("".join(parts))
    if not rows:
        return [], np.zeros((0, 0), dtype=np.uint8)
    width = {len(r) for r in rows}
    if len(width) != 1:
        raise ValueError(f"{path}: rows have different lengths")
    return ids, np.frombuffer("".join(rows).encode("ascii", "replace"), dtype=np.uint8).reshape(len(rows), -1)


def read_partitions(path):
    """[(gene, start0, end0_exclusive)] from partitions.tsv (1-based inclusive in the file)."""
    out = []
    with open(path, "r", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            gene, start, end, *_ = line.rstrip("\n").split("\t")
            out.append((gene, int(start) - 1, int(end)))
    return out


def presence_matrix(matrix, partitions):
    """bool [genomes x genes]: the genome has at least one residue in the gene's columns."""
    return np.stack([(matrix[:, s:e] != GAP).any(axis=1) for _, s, e in partitions], axis=1)


def occupancy(ids, matrix, partitions):
    """Per-genome table: genes present, gap fraction of the whole row."""
    present = presence_matrix(matrix, partitions)
    return {
        "genome": list(ids),
        "genes": present.sum(axis=1).tolist(),
        "gap_fraction": (matrix == GAP).mean(axis=1).tolist(),
    }, present


# ------------------------------------------------------------------------------ pairs
def label_vector(ids, taxonomy, rank):
    """Clade name per genome for ``rank`` ('' if unlabelled)."""
    return [taxonomy.get(i, {}).get(rank, "") for i in ids]


def sample_pairs(labels, n_pairs, rng):
    """(within, between): index pairs of labelled genomes of the same / different clade."""
    by_clade = defaultdict(list)
    for index, name in enumerate(labels):
        if name:
            by_clade[name].append(index)
    clades = [c for c, members in by_clade.items() if len(members) >= 2]
    labelled = [i for ls in by_clade.values() for i in ls]
    within, between = [], []
    if clades:
        weights = [len(by_clade[c]) for c in clades]
        for clade in rng.choices(clades, weights=weights, k=n_pairs):
            within.append(tuple(rng.sample(by_clade[clade], 2)))
    if len(by_clade) >= 2:
        guard = 0
        while len(between) < n_pairs and guard < n_pairs * 20:
            guard += 1
            i, j = rng.sample(labelled, 2)
            if labels[i] != labels[j]:
                between.append((i, j))
    return within, between


def pdist(a, b, min_sites=20):
    """Proportion of differing residues over columns where both rows have a residue (nan if few)."""
    shared = (a != GAP) & (b != GAP)
    n = int(shared.sum())
    if n < min_sites:
        return float("nan")
    return float(((a != b) & shared).sum() / n)


def pair_metrics(matrix, present, partitions, pairs, per_gene=True, min_sites=20):
    """Arrays for the pairs: overall p-distance, shared genes, gene Jaccard and (optionally)
    p-distance per gene [pairs x genes]."""
    n = len(pairs)
    overall = np.full(n, np.nan)
    shared_genes = np.zeros(n, dtype=int)
    jaccard = np.zeros(n)
    gene_dist = np.full((n, len(partitions)), np.nan) if per_gene else None
    for k, (i, j) in enumerate(pairs):
        overall[k] = pdist(matrix[i], matrix[j], min_sites)
        both = present[i] & present[j]
        either = present[i] | present[j]
        shared_genes[k] = both.sum()
        jaccard[k] = both.sum() / either.sum() if either.any() else 0.0
        if per_gene:
            for g, (_, s, e) in enumerate(partitions):
                if both[g]:
                    gene_dist[k, g] = pdist(matrix[i, s:e], matrix[j, s:e], min_sites)
    return {"pdist": overall, "shared_genes": shared_genes, "jaccard": jaccard, "gene_dist": gene_dist}


def separation(within, between):
    """P(within < between) over all (within, between) value pairs, ties counting half, nan ignored.
    1.0 = perfect separation, 0.5 = no information, <0.5 = inverted."""
    w = np.asarray(within, dtype=float)
    b = np.asarray(between, dtype=float)
    w, b = w[~np.isnan(w)], b[~np.isnan(b)]
    if len(w) == 0 or len(b) == 0:
        return float("nan")
    values = np.concatenate([w, b])
    order = values.argsort(kind="mergesort")
    ranks = np.empty(len(values))
    sorted_values = values[order]
    # average ranks for ties
    i = 0
    while i < len(sorted_values):
        j = i
        while j + 1 < len(sorted_values) and sorted_values[j + 1] == sorted_values[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2 + 1
        i = j + 1
    u = ranks[:len(w)].sum() - len(w) * (len(w) + 1) / 2   # P(w > b) * len(w)*len(b)
    return float(1 - u / (len(w) * len(b)))


def rank_report(ids, matrix, present, partitions, taxonomy, ranks=RANKS, n_pairs=2000, seed=0, per_gene=True):
    """{rank: {'within': metrics, 'between': metrics, 'summary': dict}} on sampled pairs."""
    rng = random.Random(seed)
    report = {}
    for rank in ranks:
        labels = label_vector(ids, taxonomy, rank)
        within, between = sample_pairs(labels, n_pairs, rng)
        if not within or not between:
            continue
        w = pair_metrics(matrix, present, partitions, within, per_gene)
        b = pair_metrics(matrix, present, partitions, between, per_gene)
        summary = {
            "rank": rank, "n_within": len(within), "n_between": len(between),
            "median_pdist_within": float(np.nanmedian(w["pdist"])),
            "median_pdist_between": float(np.nanmedian(b["pdist"])),
            "median_jaccard_within": float(np.median(w["jaccard"])),
            "median_jaccard_between": float(np.median(b["jaccard"])),
            "median_shared_genes_within": float(np.median(w["shared_genes"])),
            "median_shared_genes_between": float(np.median(b["shared_genes"])),
            "separation_pdist": separation(w["pdist"], b["pdist"]),
        }
        report[rank] = {"within": w, "between": b, "summary": summary}
    return report


def gene_separation(report, partitions, rank):
    """[(gene, separation, n_within_pairs_with_data)] for one rank, from the sampled pairs."""
    w, b = report[rank]["within"]["gene_dist"], report[rank]["between"]["gene_dist"]
    return [(partitions[g][0], separation(w[:, g], b[:, g]), int((~np.isnan(w[:, g])).sum()))
            for g in range(len(partitions))]


# ------------------------------------------------------------------------- multicopy
def kmers(sequence, k=3):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def jaccard_sets(a, b):
    return len(a & b) / len(a | b) if (a or b) else 0.0


def read_fasta(path):
    out, name, parts = {}, None, []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(parts)
                name, parts = line[1:].rstrip("\n"), []
            else:
                parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts)
    return out


def duplicated_cases(full_table):
    """{busco id: [row dicts, best score first]} for BUSCO genes with status Duplicated."""
    cases = defaultdict(list)
    columns = None
    with open(full_table, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("# Busco id"):
                columns = line.lstrip("# ").rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip() or columns is None:
                continue
            row = dict(zip(columns, line.rstrip("\n").split("\t")))
            if row.get("Status") == "Duplicated":
                try:
                    row["Score"] = float(row["Score"])
                except (KeyError, ValueError):
                    row["Score"] = float("-inf")
                cases[row["Busco id"]].append(row)
    for rows in cases.values():
        rows.sort(key=lambda r: -r["Score"])
    return dict(cases)


def copy_sequences(full_table, busco_id, rows):
    """Protein sequence of each Duplicated row (same matching rule as the aggregation)."""
    seq_file = Path(full_table).parent / "busco_sequences" / "multi_copy_busco_sequences" / f"{busco_id}.faa"
    if not seq_file.exists():
        return None
    records = read_fasta(seq_file)
    sequences, used = [], set()
    for row in rows:
        try:
            want = (str(row["Sequence"]), min(int(row["Gene Start"]), int(row["Gene End"])),
                    max(int(row["Gene Start"]), int(row["Gene End"])), str(row["Strand"]))
        except (KeyError, ValueError):
            sequences.append(None)
            continue
        found = None
        for header, seq in records.items():
            if header in used:
                continue
            m = re.search(r"\|([^|:]+):(\d+)-(\d+)\|([+-])(?:\s|$)", header)
            if m and (m.group(1), min(int(m.group(2)), int(m.group(3))),
                      max(int(m.group(2)), int(m.group(3))), m.group(4)) == want:
                found = header
                break
        if found is None:
            sequences.append(None)
        else:
            used.add(found)
            sequences.append(records[found].rstrip("*").replace("*", "X"))
    return sequences


def evaluate_copy_choice(copies, scores, peers, k=3):
    """For one Duplicated gene: how does the heuristic's choice (highest score) compare with the
    copy most similar to ``peers`` (sequences of the same gene in related genomes)?

    Returns None if there is nothing to compare, else a dict with the similarity of each copy
    (median k-mer Jaccard to the peers), the index chosen by score, and the index closest to peers.
    """
    if len(copies) < 2 or not peers or any(c is None for c in copies):
        return None
    peer_sets = [kmers(p, k) for p in peers]
    similarity = []
    for copy in copies:
        cset = kmers(copy, k)
        similarity.append(float(np.median([jaccard_sets(cset, ps) for ps in peer_sets])))
    chosen = int(np.argmax(scores))
    closest = int(np.argmax(similarity))
    ordered = sorted(scores, reverse=True)
    return {
        "n_copies": len(copies), "similarity": similarity, "chosen": chosen, "closest": closest,
        "agree": chosen == closest,
        "similarity_gap": similarity[closest] - similarity[chosen],
        "score_ratio": (ordered[1] / ordered[0]) if ordered[0] > 0 else float("nan"),
        "copy_identity": jaccard_sets(kmers(copies[0], k), kmers(copies[1], k)),
    }


def multicopy_analysis(busco_tables, taxonomy, merged_dir, max_cases=400, min_peers=5, n_peers=25, seed=0,
                       ranks=("family", "order", "class", "phylum")):
    """Evaluate the top-score heuristic on randomly chosen Duplicated genes.

    ``busco_tables``: [(genome accession, full_table path)]; ``merged_dir`` holds the per-gene
    FASTAs (``genome|busco|copy`` ids) the tree is built from.
    """
    rng = random.Random(seed)
    tables = list(busco_tables)
    rng.shuffle(tables)
    gene_cache, results = {}, []

    def genome_index(gene):
        if gene not in gene_cache:
            path = Path(merged_dir) / f"{gene}.faa"
            by_genome = {}
            if path.exists():
                for header, seq in read_fasta(path).items():
                    by_genome[header.split("|")[0]] = seq.rstrip("*").replace("*", "X")
            gene_cache[gene] = by_genome
        return gene_cache[gene]

    for genome, table in tables:
        if len(results) >= max_cases:
            break
        if not Path(table).exists():
            continue
        cases = duplicated_cases(table)
        if not cases:
            continue
        gene = rng.choice(sorted(cases))
        rows = cases[gene]
        copies = copy_sequences(table, gene, rows)
        if copies is None:
            continue
        index = genome_index(gene)
        lineage = taxonomy.get(genome, {})
        peers, used_rank = [], None
        for rank in ranks:
            name = lineage.get(rank)
            if not name:
                continue
            candidates = [g for g in index if g != genome and taxonomy.get(g, {}).get(rank) == name]
            if len(candidates) >= min_peers:
                peers = [index[g] for g in rng.sample(candidates, min(n_peers, len(candidates)))]
                used_rank = rank
                break
        outcome = evaluate_copy_choice(copies, [r["Score"] for r in rows], peers)
        if outcome:
            outcome.update({"genome": genome, "gene": gene, "peer_rank": used_rank, "n_peers": len(peers)})
            results.append(outcome)
    return results
