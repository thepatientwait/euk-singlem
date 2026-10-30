"""Intron helpers: where BUSCO genes have introns and whether a SingleM window crosses one.

Why: SingleM windows are matched against metagenomic *DNA* reads. A window whose 60 nt contain
an intron in the genome cannot be recovered from reads of that genome, so windows (and the
transcript sequences used to build the nucleotide database) must not span introns.

Coordinates
-----------
* ``exons``  : genomic (start, end) pairs, 1-based inclusive, from the miniprot GFF that BUSCO
  writes next to every hit (``busco_sequences/*/<busco_id>.gff``, one ``CDS`` line per exon).
* An *intron* is the genomic gap between two consecutive exons in transcript order. Gaps shorter
  than ``min_intron`` are frameshift/indel artefacts and are not introns.
* ``offset``  : position of an intron on the spliced CDS, in nt from the first coding base, i.e.
  the intron sits between CDS nt ``offset - 1`` and ``offset`` (0-based). ``offset // 3`` is the
  protein residue containing/after the intron and ``offset % 3`` the phase.
* A window of ``W`` HMM match columns starting at column ``p`` covers nt ``[3a, 3b + 3)`` where
  ``a``/``b`` are the residues aligned to its first/last column. It *spans* an intron iff
  ``3a < offset < 3b + 3``.
"""
from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

DEFAULT_MIN_INTRON = 30
WINDOW_NT = 60  # SingleM's default window; `singlem seqs` is run without --window-size


# --------------------------------------------------------------------------------------
# GFF / exon structure
# --------------------------------------------------------------------------------------
@dataclass
class Model:
    """One miniprot gene model (an ``mRNA`` with its ``CDS`` exons)."""

    contig: str
    start: int
    end: int
    strand: str
    exons: List[Tuple[int, int]] = field(default_factory=list)
    # miniprot marks alignments that needed a frameshift or crossed an in-frame stop codon;
    # exon offsets after such a point do not correspond to protein residues
    irregular: bool = False

    def transcript_order(self) -> List[Tuple[int, int]]:
        return sorted(self.exons, reverse=(self.strand == "-"))


def parse_gff_text(lines) -> List[Model]:
    """Gene models from the lines of a miniprot GFF (``mRNA`` + ``CDS`` lines linked by ``Parent``)."""
    models: Dict[str, Model] = {}
    order: List[str] = []
    for line in lines:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        feature = parts[2]
        attrs = dict(kv.split("=", 1) for kv in parts[8].split(";") if "=" in kv)
        if feature == "mRNA":
            model_id = attrs.get("ID", f"m{len(order)}")
            models[model_id] = Model(parts[0], int(parts[3]), int(parts[4]), parts[6])
            order.append(model_id)
        elif feature == "CDS":
            parent = attrs.get("Parent")
            if parent in models:
                models[parent].exons.append((int(parts[3]), int(parts[4])))
                if "Frameshift" in attrs or "StopCodon" in attrs:
                    models[parent].irregular = True
    return [models[i] for i in order if models[i].exons]


def parse_gff(path: Path) -> List[Model]:
    """Read the gene models of a miniprot GFF file (see :func:`parse_gff_text`)."""
    return parse_gff_text(Path(path).read_text().splitlines())


def cds_matches_protein(model: Model, cds_len: int, protein_len: int) -> bool:
    """True if the exon offsets of ``model`` can be read as protein residue positions.

    Offsets are counted from the first coding base, so they are exact for the whole gene when
    there is no frameshift/stop-codon event, the CDS is a whole number of codons, and its length
    agrees with the protein to within one residue (a terminal stop codon, or one residue that
    miniprot reports at the 3' end beyond the last CDS line).
    """
    return (not model.irregular) and cds_len % 3 == 0 and abs(cds_len // 3 - protein_len) <= 1


def intron_gaps(model: Model, min_intron: int = DEFAULT_MIN_INTRON):
    """Introns of a model: list of ``(offset, genomic_gap_length)`` plus the CDS length.

    ``offset`` is the spliced-CDS position of the intron (see module docstring); the CDS length
    is the sum of exon lengths (it can include a terminal stop codon).
    """
    exons = model.transcript_order()
    introns, cum = [], 0
    for i, (s, e) in enumerate(exons):
        cum += e - s + 1
        if i + 1 < len(exons):
            ns, ne = exons[i + 1]
            gap = (ns - e - 1) if model.strand == "+" else (s - ne - 1)
            if gap >= min_intron:
                introns.append((cum, gap))
    return introns, cum


def exon_segments(model: Model, min_intron: int = DEFAULT_MIN_INTRON) -> List[Tuple[int, int]]:
    """Genomic segments (transcript order) that contain no intron.

    Consecutive exons separated by less than ``min_intron`` are merged (the gap is a
    frameshift/indel, not an intron), so each returned ``(start, end)`` can be read from the
    genome as one contiguous piece.
    """
    exons = model.transcript_order()
    segments: List[List[int]] = []
    for i, (s, e) in enumerate(exons):
        if segments:
            ps, pe = segments[-1]
            gap = (s - pe - 1) if model.strand == "+" else (ps - e - 1)
            if gap < min_intron:
                segments[-1] = [min(ps, s), max(pe, e)]
                continue
        segments.append([s, e])
    return [(s, e) for s, e in segments]


def select_model(models: Sequence[Model], contig: str, start: int, end: int,
                 strand: Optional[str] = None, tolerance: int = 6) -> Optional[Model]:
    """The model matching a BUSCO table row (contig and span; BUSCO's span is ~3 nt shorter
    than miniprot's because of the stop codon and a 0/1-based shift)."""
    lo, hi = min(start, end), max(start, end)
    for m in models:
        if m.contig != contig or (strand and m.strand != strand):
            continue
        if abs(m.start - lo) <= tolerance and abs(m.end - hi) <= tolerance:
            return m
    return None


# --------------------------------------------------------------------------------------
# Protein <-> HMM columns
# --------------------------------------------------------------------------------------
def hmm_length(hmm_path: Path) -> int:
    with open(hmm_path) as handle:
        for line in handle:
            if line.startswith("LENG"):
                return int(line.split()[1])
    raise ValueError(f"No LENG line in {hmm_path}")


_CLEAN = re.compile(r"[^A-Za-z]")


def clean_protein(seq: str) -> str:
    """Letters only (drops ``*`` and gaps); hmmalign copes with ambiguity codes."""
    return _CLEAN.sub("", seq)


def hmmalign_matrix(hmm_path: Path, proteins: Dict[str, str], hmmalign: str = "hmmalign"):
    """Align proteins to an HMM and return, per protein, the residue index at each match column.

    Returns ``(names, matrix)``; ``matrix[i, c]`` is the 0-based index (in the cleaned protein)
    of the residue aligned to HMM match column ``c`` of protein ``i``, or ``-1`` for a deletion.
    Insert residues (lowercase) are counted in the residue index but occupy no column.
    """
    names = list(proteins)
    if not names:
        return names, np.zeros((0, hmm_length(hmm_path)), dtype=np.int32)
    fasta = "".join(f">{n}\n{proteins[n]}\n" for n in names)
    result = subprocess.run([hmmalign, "--amino", "--outformat", "pfam", str(hmm_path), "/dev/stdin"],
                            input=fasta, capture_output=True, text=True, check=True)
    aligned: Dict[str, List[str]] = {}
    for line in result.stdout.splitlines():
        if not line or line.startswith(("#", "//")):
            continue
        name, chunk = line.split(None, 1)
        aligned.setdefault(name, []).append(chunk.strip())
    ncol = hmm_length(hmm_path)
    matrix = np.full((len(names), ncol), -1, dtype=np.int32)
    for i, name in enumerate(names):
        col = res = -1
        for ch in "".join(aligned.get(name, [])):
            if ch == "-":
                col += 1
            elif ch == ".":
                continue
            else:
                res += 1
                if ch.isupper():
                    col += 1
                    matrix[i, col] = res
    return names, matrix


def residue_columns(matrix_row: np.ndarray, n_residues: int) -> np.ndarray:
    """Inverse of one matrix row: per residue, its match column (``-1`` for inserts)."""
    cols = np.full(n_residues, -1, dtype=np.int32)
    present = np.nonzero(matrix_row >= 0)[0]
    cols[matrix_row[present]] = present
    return cols


def aligned_extent(matrix: np.ndarray):
    """First and last HMM match column each sequence is aligned to, and how many it occupies.

    Returns three int arrays (``-1, -1, 0`` for a sequence aligned to no column). Used to draw
    the alignment: a sequence spans ``first..last`` with ``n_cols`` of those columns occupied
    (the rest are deletions).
    """
    present = matrix >= 0
    n_cols = present.sum(axis=1).astype(np.int64)
    has = n_cols > 0
    first = np.where(has, present.argmax(axis=1), -1).astype(np.int64)
    last = np.where(has, matrix.shape[1] - 1 - present[:, ::-1].argmax(axis=1), -1).astype(np.int64)
    return first, last, n_cols


# --------------------------------------------------------------------------------------
# Windows
# --------------------------------------------------------------------------------------
def window_scan(matrix: np.ndarray, introns: Sequence[Sequence[int]], valid: np.ndarray,
                window_cols: int = WINDOW_NT // 3):
    """Score every window start for coverage and for intron spanning.

    Parameters
    ----------
    matrix : (n_genomes, n_cols) residue index per match column (``-1`` = not aligned).
    introns : per genome, spliced-CDS offsets of its introns.
    valid : per genome, False if its gene model could not be mapped onto the protein (those
        genomes still count for coverage but not for intron statistics).
    window_cols : window length in HMM columns (nt / 3).

    Returns
    -------
    dict of arrays indexed by window start ``p`` (0-based): ``n_cover`` genomes aligned at
    both ends, ``aligned_bases`` (SingleM's score: aligned positions inside the window summed
    over covering genomes), ``n_cover_valid`` and ``n_span`` (valid covering genomes whose
    window spans an intron), plus ``spans`` (per genome x start boolean) when requested.
    """
    n, ncol = matrix.shape
    n_starts = ncol - window_cols + 1
    if n_starts < 1:
        raise ValueError("HMM shorter than the window")
    aligned = matrix >= 0
    cum = np.concatenate([np.zeros((n, 1), dtype=np.int32), np.cumsum(aligned, axis=1, dtype=np.int32)], axis=1)
    covers = aligned[:, :n_starts] & aligned[:, window_cols - 1:window_cols - 1 + n_starts]
    aligned_in_window = cum[:, window_cols:window_cols + n_starts] - cum[:, :n_starts]
    aligned_bases = (aligned_in_window * covers).sum(axis=0)
    first = matrix[:, :n_starts]
    last = matrix[:, window_cols - 1:window_cols - 1 + n_starts]
    spans = np.zeros((n, n_starts), dtype=bool)
    for g in range(n):
        if not valid[g] or not introns[g]:
            continue
        for offset in introns[g]:
            spans[g] |= (first[g] <= (offset - 1) // 3) & (last[g] >= offset // 3)
    spans &= covers
    return {
        "n_cover": covers.sum(axis=0),
        "aligned_bases": aligned_bases,
        "n_cover_valid": (covers & valid[:, None]).sum(axis=0),
        "n_span": spans.sum(axis=0),
        "spans": spans,
        "covers": covers,
    }


def choose_window(scan: dict, max_span_fraction: float, score_tolerance: float = 0.005, min_valid: int = 20):
    """Choose the window start (0-based) in HMM columns, avoiding introns.

    ``singlem seqs`` picks the window covered by the most aligned bases (earliest on ties). This
    does the same in HMM coordinates, restricted to windows that cross an intron in at most
    ``max_span_fraction`` of the valid genomes covering them. Windows scoring within
    ``score_tolerance`` of the best allowed score count as equally well covered; among those the
    one crossing introns least wins (then the earliest). With no allowed window the least-spanning
    one wins (``violation=True``).

    Returns ``(start, reason, violation, best_coverage_start)`` where the last item is the
    unconstrained SingleM-style choice, for comparison.
    """
    score = scan["aligned_bases"].astype(float)
    valid = scan["n_cover_valid"]
    span_frac = np.where(valid >= min_valid, scan["n_span"] / np.maximum(valid, 1), 1.0)
    best_coverage = int(np.argmax(score))  # earliest maximum, like singlem seqs
    allowed = (span_frac <= max_span_fraction) & (valid >= min_valid)
    if allowed.any():
        top = score[allowed].max()
        near = np.nonzero(allowed & (score >= (1 - score_tolerance) * top))[0]
        start = int(near[np.lexsort((near, span_frac[near]))[0]])
        unconstrained_ok = bool(allowed[best_coverage])
        return start, ("best_coverage" if unconstrained_ok and start == best_coverage else
                       "moved_for_introns" if not unconstrained_ok else "best_coverage_fewer_introns"), False, best_coverage
    candidates = np.nonzero(valid >= min_valid)[0]
    if len(candidates) == 0:
        return best_coverage, "insufficient_intron_data", False, best_coverage
    order = np.lexsort((candidates, -score[candidates], span_frac[candidates]))
    return int(candidates[order[0]]), "no_window_under_threshold", True, best_coverage


def spans_intron(first_residue: int, last_residue: int, offsets: Iterable[int]) -> bool:
    """True if an intron lies strictly inside the nt range covered by the two residues."""
    return any(3 * first_residue < o < 3 * last_residue + 3 for o in offsets)
