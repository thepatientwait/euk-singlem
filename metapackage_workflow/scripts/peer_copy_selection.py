"""Choose which copy of a Duplicated BUSCO gene represents a genome.

The top BUSCO score is not a good guide to which paralog is the orthologous one (it matches the copy
most similar to the genome's taxonomic peers only 73-77% of the time). For every genome with several
candidate copies of a gene, each copy is aligned (local, BLOSUM62) to Complete (single-copy) sequences
of the same gene from genomes of the same taxon, and the copy with the highest mean identity is kept.

Peers are the Complete sequences of other genomes sharing the candidate genome's lineage down to the
deepest rank (genus ... phylum) that has at least ``min_peers`` such genomes; at most ``max_peers`` are
used, picked evenly from the sorted genome ids so reruns choose the same ones. Genomes without a usable
taxonomy, or without any Complete peer even at phylum level, fall back to the top-scoring copy.

Ties: higher mean identity, then higher BUSCO score, then longer, then lower copy number.
"""
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

from joblib import Parallel, delayed

_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZX*")
RANK_LEVELS = (6, 5, 4, 3, 2)  # lineage prefix lengths tried: genus, family, order, class, phylum (with domain)

_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio import Align
        from Bio.Align import substitution_matrices
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        aligner.open_gap_score = -11
        aligner.extend_gap_score = -1
        _ALIGNER = aligner
    return _ALIGNER


def clean(seq: str) -> str:
    """Upper-case protein with anything BLOSUM62 does not know replaced by X, without stops."""
    seq = seq.upper().replace("*", "")
    return "".join(c if c in _ALPHABET else "X" for c in seq)


def identity(a: str, b: str) -> float:
    """Identical aligned residues / length of the shorter sequence (best local alignment)."""
    if not a or not b:
        return 0.0
    aligner = _aligner()
    best = aligner.align(a, b)[0]
    return best.counts().identities / min(len(a), len(b))


def mean_identity(candidate: str, peers: Sequence[str]) -> float:
    return sum(identity(candidate, peer) for peer in peers) / len(peers)


def protein_of(info: Dict) -> str:
    sequence = info["sequence"]
    return clean("".join(sequence.strip().split("\n")[1:]) if ">" in sequence else sequence)


def lineage_parts(taxonomy: Optional[str]) -> List[str]:
    return [part.strip() for part in (taxonomy or "").split(";")]


def lineage_key(parts: List[str], level: int) -> Optional[Tuple[str, ...]]:
    """Lineage prefix of ``level`` ranks, or None if that rank is missing/blank."""
    if len(parts) < level:
        return None
    prefix = parts[:level]
    if not prefix[-1] or prefix[-1].endswith("__"):
        return None
    return tuple(prefix)


def evenly_sampled(items: List, limit: int) -> List:
    if len(items) <= limit:
        return items
    step = len(items) / limit
    return [items[int(i * step)] for i in range(limit)]


def find_peers(parts: List[str], index: Dict[Tuple[str, ...], List[Tuple[str, str]]],
               min_peers: int, max_peers: int):
    """Peers (list of (genome, protein)) at the deepest level with >= min_peers, else the best
    available level; (level, peers). ``level`` is 0 when there are none."""
    fallback = (0, [])
    for level in RANK_LEVELS:
        key = lineage_key(parts, level)
        if key is None:
            continue
        peers = index.get(key, [])
        if len(peers) >= min_peers:
            return level, evenly_sampled(peers, max_peers)
        if peers and not fallback[1]:
            fallback = (level, peers)
    return fallback[0], evenly_sampled(fallback[1], max_peers)


def _score_task(candidates: List[str], peers: List[str]) -> List[float]:
    return [mean_identity(candidate, peers) for candidate in candidates]


def select_peer_copies(gene_sequences: Dict[str, Dict[str, Dict]], usable_taxonomy: Dict[str, str],
                       min_peers: int = 3, max_peers: int = 8, n_jobs: int = 1):
    """Reduce every genome with several candidate copies of a gene to one copy.

    gene_sequences: {busco_id: {sequence_id: info}} where info has genome_acc, status, score, length,
    sequence and sequence_id (``genome|busco|copy``). Returns (selected, report_rows):
    ``selected`` has the same shape with one entry per genome and gene (renamed to copy 1);
    report_rows describes every decision.
    """
    selected: Dict[str, Dict[str, Dict]] = {}
    tasks = []  # (busco_id, genome, [infos], level, n_peers, candidate proteins, peer proteins)
    fallback_rows = []

    for busco_id, seqs in gene_sequences.items():
        by_genome = defaultdict(list)
        for info in seqs.values():
            by_genome[info["genome_acc"]].append(info)

        selected[busco_id] = {}
        has_candidates = any(len(infos) > 1 for infos in by_genome.values())
        index: Dict[Tuple[str, ...], List[Tuple[str, str]]] = defaultdict(list)
        for genome, infos in by_genome.items():
            if len(infos) > 1:
                continue  # candidates, resolved below
            info = infos[0]
            selected[busco_id][info["sequence_id"]] = info
            if not has_candidates or info.get("status", "Complete") != "Complete":
                continue  # only single-copy sequences serve as peers
            parts = lineage_parts(usable_taxonomy.get(str(genome)))
            protein = None
            for level in RANK_LEVELS:
                key = lineage_key(parts, level)
                if key is not None:
                    if protein is None:
                        protein = protein_of(info)
                    index[key].append((genome, protein))
        for peers in index.values():
            peers.sort()

        for genome, infos in sorted(by_genome.items()):
            if len(infos) == 1:
                continue
            infos = sorted(infos, key=lambda i: i.get("copy_number", 1))
            parts = lineage_parts(usable_taxonomy.get(str(genome)))
            level, peers = find_peers(parts, index, min_peers, max_peers)
            if not peers:
                fallback_rows.append((busco_id, genome, infos))
                continue
            tasks.append((busco_id, genome, infos, level, len(peers),
                          [protein_of(i) for i in infos], [p for _, p in peers]))

    scores = Parallel(n_jobs=n_jobs, batch_size="auto")(
        delayed(_score_task)(task[5], task[6]) for task in tasks
    )

    report = []
    for task, identities in zip(tasks, scores):
        busco_id, genome, infos, level, n_peers, _, _ = task
        order = sorted(
            range(len(infos)),
            key=lambda i: (-identities[i], -float(infos[i].get("score") or 0.0),
                           -infos[i].get("length", 0), infos[i].get("copy_number", 1)),
        )
        best = order[0]
        top_score = max(range(len(infos)), key=lambda i: (float(infos[i].get("score") or 0.0),
                                                           infos[i].get("length", 0), -i))
        report.append({
            "busco_id": busco_id, "genome_acc": genome, "n_candidates": len(infos),
            "chosen_copy": infos[best].get("copy_number", best + 1),
            "chosen_identity": round(identities[best], 4),
            "top_score_copy": infos[top_score].get("copy_number", top_score + 1),
            "top_score_identity": round(identities[top_score], 4),
            "agrees_with_top_score": int(best == top_score),
            "peer_level": level, "n_peers": n_peers,
        })
        _store(selected, busco_id, infos[best])

    for busco_id, genome, infos in fallback_rows:
        top = max(infos, key=lambda i: (float(i.get("score") or 0.0), i.get("length", 0)))
        report.append({
            "busco_id": busco_id, "genome_acc": genome, "n_candidates": len(infos),
            "chosen_copy": top.get("copy_number", 1), "chosen_identity": "",
            "top_score_copy": top.get("copy_number", 1), "top_score_identity": "",
            "agrees_with_top_score": 1, "peer_level": 0, "n_peers": 0,
        })
        _store(selected, busco_id, top)

    return selected, report


def _store(selected, busco_id, info):
    """Keep the chosen copy under the genome's single-copy id (``genome|busco|1``)."""
    new_id = f"{info['genome_acc']}|{busco_id}|1"
    selected[busco_id][new_id] = {**info, "sequence_id": new_id, "copy_number": 1}


REPORT_COLUMNS = ["busco_id", "genome_acc", "n_candidates", "chosen_copy", "chosen_identity",
                  "top_score_copy", "top_score_identity", "agrees_with_top_score", "peer_level", "n_peers"]


def write_report(rows: List[Dict], path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[c]) for c in REPORT_COLUMNS) + "\n")
