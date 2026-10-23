#!/usr/bin/env python3
"""Pick a small, reproducible subset of genomes for quick checks of the tree workflow.

Modes (``tree.subset`` in config.yaml):

* ``family5``          ``per_family`` random genomes from every family, all sources
* ``jgi_p10k``         every JGI (fungi, algae) and P10K genome
* ``jgi_p10k_family5`` ``per_family`` random JGI/P10K genomes per family

Only genomes with a usable taxonomy (``tax_status`` complete/partial) are considered; the
family modes need a family. ``n_external`` external genomes are always added (chosen at random).
Selection is deterministic: genomes are ranked by an md5 of ``seed:accession``, so it does not
depend on input order.
"""
import csv
import hashlib
import sys
from pathlib import Path

JGI_P10K_GROUPS = {"fungi", "algae", "p10k"}
MODES = ("family5", "jgi_p10k", "jgi_p10k_family5")

try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:
    redirect_to_snakemake_log = None


def _rank(seed, accession):
    return hashlib.md5(f"{seed}:{accession}".encode()).hexdigest()


def select_subset(mode, base_rows, taxonomy, external, per_family=5, n_external=50, seed=1):
    """Returns (base_selected, external_selected) as sorted accession lists.

    base_rows: iterable of (group, accession); taxonomy: {accession: (family, tax_status)}.
    """
    if mode not in MODES:
        raise ValueError(f"unknown subset mode {mode!r}; choose from {MODES}")
    candidates = []
    for group, accession in base_rows:
        if mode.startswith("jgi_p10k") and group not in JGI_P10K_GROUPS:
            continue
        family, status = taxonomy.get(accession, ("", ""))
        if status not in ("complete", "partial"):
            continue
        candidates.append((group, accession, family or ""))

    if mode == "jgi_p10k":
        chosen = [a for _, a, _ in candidates]
    else:
        by_family = {}
        for _, accession, family in candidates:
            if family:
                by_family.setdefault(family, []).append(accession)
        chosen = []
        for family in sorted(by_family):
            chosen += sorted(by_family[family], key=lambda a: _rank(seed, a))[:per_family]

    ext = sorted(sorted(external, key=lambda a: _rank(seed, a))[:n_external])
    return sorted(set(chosen)), ext


def read_taxonomy(path):
    out = {}
    with open(path, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            out[row["acc"]] = (row.get("family") or "", row.get("tax_status") or "")
    return out


def read_manifest(path, base):
    """[(group, accession, absolute busco_result path, complete_count)] from a BUSCO QC manifest."""
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result = Path(row["busco_result"])
            if not result.is_absolute():
                result = Path(base) / result
            rows.append((row["group"], row["genome_acc"], str(result.resolve()), row.get("complete_count", "")))
    return rows


def main():
    if redirect_to_snakemake_log:
        redirect_to_snakemake_log()
    p = snakemake.params
    base = read_manifest(snakemake.input.base_qc, p.base_dir)
    external = [r[1] for r in read_manifest(snakemake.input.external_qc, p.external_dir)]
    taxonomy = read_taxonomy(snakemake.input.taxonomy)
    chosen, ext = select_subset(p.mode, [(g, a) for g, a, _, _ in base], taxonomy, external,
                                int(p.per_family), int(p.n_external), int(p.seed))
    keep = set(chosen)
    Path(snakemake.output.genomes).parent.mkdir(parents=True, exist_ok=True)
    with open(snakemake.output.genomes, "w", encoding="utf-8") as handle:
        handle.write("\n".join(chosen + ext) + "\n")
    with open(snakemake.output.manifest, "w", encoding="utf-8") as handle:
        handle.write("group\tgenome_acc\tbusco_result\tcomplete_count\n")
        for group, accession, result, count in base:
            if accession in keep:
                handle.write(f"{group}\t{accession}\t{result}\t{count}\n")
    print(f"subset {p.mode}: {len(chosen)} base genomes, {len(ext)} external genomes", file=sys.stderr)


if __name__ == "__main__":
    main()
