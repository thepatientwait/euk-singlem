#!/usr/bin/env python3
"""Pick a small, reproducible subset of the genomes in the BUSCO QC manifest, for quick checks of the
metapackage rules and the intron survey (``singlem.introns.subset`` / ``--config subset=...``).

Modes
-----
* ``family5``          ``per_family`` random genomes from every family, all sources
* ``jgi_p10k``         every JGI (fungi, algae) and P10K genome
* ``jgi_p10k_family5`` ``per_family`` random JGI/P10K genomes per family

Only genomes with a usable taxonomy (``tax_status`` complete or partial) are considered; the family
modes also need a family. Selection is deterministic: genomes are ranked by an md5 of
``seed:accession``, so it does not depend on the order of the input.

The subset is taken from the QC manifest as it is on disk, without rerunning BUSCO, so run it after the
BUSCO results you want to test with exist. It is not tracked as a Snakemake input (that would pull in the
whole BUSCO phase); use ``--forcerun make_subset_manifest`` after the QC manifest changes.
"""
import csv
import hashlib
import sys
from pathlib import Path

JGI_P10K_GROUPS = {"fungi", "algae", "p10k"}
MODES = ("family5", "jgi_p10k", "jgi_p10k_family5")


def _rank(seed, accession):
    return hashlib.md5(f"{seed}:{accession}".encode()).hexdigest()


def select_subset(mode, rows, taxonomy, per_family=5, seed=1):
    """Sorted accessions to keep.

    Parameters
    ----------
    mode : str
        One of MODES.
    rows : iterable of (group, accession)
        The genomes of the QC manifest.
    taxonomy : dict
        ``{accession: (family, tax_status)}``.
    per_family : int, optional
        Genomes per family in the family modes.
    seed : int, optional
        Seed of the ranking.
    """
    if mode not in MODES:
        raise ValueError(f"unknown subset mode {mode!r}; choose from {MODES}")
    candidates = []
    for group, accession in rows:
        if mode.startswith("jgi_p10k") and group not in JGI_P10K_GROUPS:
            continue
        family, status = taxonomy.get(accession, ("", ""))
        if status not in ("complete", "partial"):
            continue
        candidates.append((accession, family or ""))
    if mode == "jgi_p10k":
        return sorted({a for a, _ in candidates})
    by_family = {}
    for accession, family in candidates:
        if family:
            by_family.setdefault(family, []).append(accession)
    chosen = []
    for family in sorted(by_family):
        chosen += sorted(by_family[family], key=lambda a: _rank(seed, a))[:per_family]
    return sorted(set(chosen))


def read_taxonomy(path):
    with open(path, "r", encoding="utf-8") as handle:
        return {row["acc"]: (row.get("family") or "", row.get("tax_status") or "")
                for row in csv.DictReader(handle, delimiter="\t")}


def subset_manifest(qc_manifest, taxonomy_path, out_path, mode, per_family=5, seed=1):
    """Write the QC manifest rows of the selected genomes; returns how many were kept."""
    with open(qc_manifest, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header, rows = reader.fieldnames, list(reader)
    chosen = set(select_subset(mode, [(r["group"], r["genome_acc"]) for r in rows], read_taxonomy(taxonomy_path), per_family, seed))
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            if row["genome_acc"] in chosen:
                out.write("\t".join(row[c] for c in header) + "\n")
                kept += 1
    return kept


if __name__ == "__main__":
    from log_config import configure_logger
    logger = configure_logger(__file__)
    p = snakemake.params
    if not Path(p.qc_manifest).exists():
        raise FileNotFoundError(f"{p.qc_manifest} not found: run the BUSCO QC step first")
    kept = subset_manifest(p.qc_manifest, snakemake.input.taxonomy, snakemake.output.manifest, p.mode, int(p.per_family), int(p.seed))
    logger.info("subset %s: %d genomes kept", p.mode, kept)
    if kept == 0:
        sys.exit("the subset is empty")
