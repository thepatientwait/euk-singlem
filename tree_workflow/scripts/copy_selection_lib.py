"""Computations behind the copy-selection notebook (notebooks/copy_selection.ipynb).

The aggregation writes ``multicopy_selection.tsv`` (scripts/peer_copy_selection.py) next to the per-gene
FASTAs: one row for every genome x BUSCO gene that had several candidate copies, with the copy chosen by
peer similarity, the copy the old "top BUSCO score" rule would have chosen, and the mean identity of both
to Complete copies of the same gene from the same taxon.

Questions answered:

* how many choices were made, and how many by peers (rather than the fallback to the top score)?
* how often does the peer rule disagree with the top score, and how much better is the chosen copy?
* where do the disagreements cluster (genes, genomes, ranks of the peer set)?
* which choices are weak (few peers, or no copy resembles its peers)?
"""
from pathlib import Path

import polars as pl

LEVEL_NAMES = {0: "no peers (top score)", 2: "phylum", 3: "class", 4: "order", 5: "family", 6: "genus"}
NUMERIC = {"n_candidates": pl.Int64, "chosen_copy": pl.Int64, "top_score_copy": pl.Int64,
           "agrees_with_top_score": pl.Int64, "peer_level": pl.Int64, "n_peers": pl.Int64,
           "chosen_identity": pl.Float64, "top_score_identity": pl.Float64}


def load_selection(path, source=None):
    """Read one multicopy_selection.tsv; adds ``source``, ``decided_by_peers`` and ``switched``."""
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)
    frame = frame.with_columns([pl.col(c).cast(t, strict=False) for c, t in NUMERIC.items() if c in frame.columns])
    return frame.with_columns(
        pl.lit(source or Path(path).parent.name).alias("source"),
        pl.col("chosen_identity").is_not_null().alias("decided_by_peers"),
        (pl.col("agrees_with_top_score") == 0).alias("switched"),
    )


def load_all(entries):
    """entries: [{"label", "path"}]; existing files only. Returns (frame or None, [missing labels])."""
    frames, missing = [], []
    for entry in entries:
        if Path(entry["path"]).is_file():
            frames.append(load_selection(entry["path"], entry["label"]))
        else:
            missing.append(entry["label"])
    return (pl.concat(frames, how="vertical_relaxed") if frames else None), missing


def overview(frame):
    """One row per source: decisions, share decided by peers, share switched away from the top score."""
    return (
        frame.group_by("source", maintain_order=True)
        .agg(
            pl.len().alias("decisions"),
            pl.col("genome_acc").n_unique().alias("genomes"),
            pl.col("busco_id").n_unique().alias("genes"),
            pl.col("decided_by_peers").sum().alias("by_peers"),
            (pl.col("decided_by_peers") & pl.col("switched")).sum().alias("switched"),
        )
        .with_columns(
            (pl.col("by_peers") / pl.col("decisions")).alias("share_by_peers"),
            (pl.col("switched") / pl.col("by_peers")).alias("share_switched"),
        )
    )


def by_candidates(frame):
    """Switch rate by number of candidate copies (peer-decided choices only)."""
    return (
        frame.filter(pl.col("decided_by_peers"))
        .group_by("n_candidates").agg(pl.len().alias("decisions"), pl.col("switched").mean().alias("share_switched"))
        .sort("n_candidates")
    )


def by_peer_level(frame):
    """Decisions and switch rate by the taxonomic level the peers came from."""
    return (
        frame.group_by("peer_level")
        .agg(pl.len().alias("decisions"), pl.col("switched").mean().alias("share_switched"),
             pl.col("n_peers").mean().alias("mean_peers"))
        .with_columns(pl.col("peer_level").replace_strict(LEVEL_NAMES, default="other").alias("level"))
        .sort("peer_level", descending=True)
    )


def gains(frame):
    """Identity gain (chosen - top-score copy) for switched, peer-decided choices."""
    return (
        frame.filter(pl.col("decided_by_peers") & pl.col("switched"))
        .with_columns((pl.col("chosen_identity") - pl.col("top_score_identity")).alias("gain"))
    )


def gene_table(frame, min_decisions=20):
    """Per BUSCO gene: how often the peer rule overrules the top score."""
    return (
        frame.filter(pl.col("decided_by_peers"))
        .group_by("busco_id")
        .agg(pl.len().alias("decisions"), pl.col("switched").mean().alias("share_switched"),
             pl.col("chosen_identity").mean().alias("mean_chosen_identity"))
        .filter(pl.col("decisions") >= min_decisions)
        .sort("share_switched", descending=True)
    )


def genome_table(frame, min_decisions=10):
    """Per genome: share of its multicopy genes switched. A very high share can point to an assembly
    problem (e.g. a duplicated or contaminated genome) rather than to the selection rule."""
    return (
        frame.filter(pl.col("decided_by_peers"))
        .group_by("genome_acc")
        .agg(pl.len().alias("decisions"), pl.col("switched").mean().alias("share_switched"),
             pl.col("chosen_identity").mean().alias("mean_chosen_identity"))
        .filter(pl.col("decisions") >= min_decisions)
        .sort("share_switched", descending=True)
    )


def weak_choices(frame, min_peers=3, low_identity=0.4):
    """Peer-decided choices that deserve less trust: thin peer set, or even the best copy resembles
    its peers poorly (no real ortholog among the candidates)."""
    peer = frame.filter(pl.col("decided_by_peers"))
    return {
        "decided by peers": peer.height,
        "fewer than %d peers" % min_peers: peer.filter(pl.col("n_peers") < min_peers).height,
        "best copy identity < %.1f" % low_identity: peer.filter(pl.col("chosen_identity") < low_identity).height,
        "switched and gain < 0.01": gains(frame).filter(pl.col("gain") < 0.01).height,
    }


def attach_phylum(frame, taxonomy):
    """Add a ``phylum`` column from {genome: 'd__..;p__..;...'}."""
    def phylum(acc):
        for part in (taxonomy.get(acc) or "").split(";"):
            part = part.strip()
            if part.startswith("p__") and len(part) > 3:
                return part[3:]
        return "unknown"
    return frame.with_columns(pl.col("genome_acc").map_elements(phylum, return_dtype=pl.Utf8).alias("phylum"))


def phylum_table(frame, min_decisions=50):
    return (
        frame.filter(pl.col("decided_by_peers"))
        .group_by("phylum").agg(pl.len().alias("decisions"), pl.col("switched").mean().alias("share_switched"))
        .filter(pl.col("decisions") >= min_decisions).sort("share_switched", descending=True)
    )
