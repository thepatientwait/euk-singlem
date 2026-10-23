#!/usr/bin/env python3
"""Merge the metapackage and external per-gene protein FASTAs for one gene.

* Trailing ``*`` (stop) is removed and internal ``*`` become ``X``: alignment tools do not
  accept them, and ~half of the BUSCO proteins carry one.
* Duplicate sequence ids are an error: two inputs claiming the same ``genome|busco|copy`` would
  otherwise be silently concatenated and distort the alignment.
* ``copy``: ``multi`` keeps every sequence (Complete BUSCOs and the top copy of Duplicated ones);
  ``single`` keeps only Complete ones. The status comes from the ``<gene>_status.tsv`` sidecar
  written by the aggregation; inputs without a sidecar were aggregated single-copy, so all their
  sequences count as Complete (reported as ``no_sidecar``).
* ``keep`` (genome accessions; the dereplicated tree) restricts the base input only: external
  genomes are already clustered.
* ``subset`` (genome accessions; quick-check subsets) restricts base and external input.
"""
from pathlib import Path

try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:  # imported from tests
    redirect_to_snakemake_log = None


def read_records(path: Path):
    """Yield (id, sequence) from a FASTA file (missing/empty file: nothing)."""
    if not path.exists() or path.stat().st_size == 0:
        return
    header, parts = None, []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:].split()[0], []
            elif line:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def clean_sequence(sequence: str) -> str:
    return sequence.rstrip("*").replace("*", "X")


def read_status(fasta: Path):
    """{sequence id: status} from the sidecar next to ``fasta`` (None if there is no sidecar)."""
    sidecar = fasta.with_name(fasta.stem + "_status.tsv")
    if not sidecar.exists():
        return None
    status = {}
    with open(sidecar, "r", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                status[fields[0]] = fields[1]
    return status


def merge(base: Path, external: Path, output: Path, use_external: bool, keep=None, subset=None,
          copy: str = "multi") -> dict:
    if copy not in ("single", "multi"):
        raise ValueError(f"copy must be 'single' or 'multi', not {copy!r}")
    seen = set()
    stats = {"base": 0, "external": 0, "dropped_not_kept": 0, "dropped_not_in_subset": 0,
             "dropped_duplicated": 0, "no_sidecar": 0, "empty_after_clean": 0}
    output.parent.mkdir(parents=True, exist_ok=True)
    sources = [("base", base)] + ([("external", external)] if use_external else [])
    with open(output, "w", encoding="utf-8") as out:
        for label, path in sources:
            status = read_status(path) if copy == "single" else None
            if copy == "single" and status is None and path.exists():
                stats["no_sidecar"] += 1
            for record_id, sequence in read_records(path):
                genome = record_id.split("|")[0]
                if label == "base" and keep is not None and genome not in keep:
                    stats["dropped_not_kept"] += 1
                    continue
                if subset is not None and genome not in subset:
                    stats["dropped_not_in_subset"] += 1
                    continue
                if status is not None and status.get(record_id, "Complete") == "Duplicated":
                    stats["dropped_duplicated"] += 1
                    continue
                if record_id in seen:
                    raise ValueError(f"Duplicate sequence id {record_id!r} while merging {path}")
                seen.add(record_id)
                sequence = clean_sequence(sequence)
                if not sequence:
                    stats["empty_after_clean"] += 1
                    continue
                out.write(f">{record_id}\n{sequence}\n")
                stats[label] += 1
    return stats


def load_keep(path):
    """Genome accessions from a one-per-line file (or None for 'keep everything')."""
    if not path:
        return None
    with open(path, "r", encoding="utf-8") as handle:
        return {line.strip().split("\t")[0] for line in handle if line.strip()}


def main() -> None:
    redirect_to_snakemake_log()
    stats = merge(
        Path(str(snakemake.params.base_faa)),
        Path(str(snakemake.params.external_faa)),
        Path(str(snakemake.output.merged_faa)),
        bool(snakemake.params.use_external),
        load_keep(snakemake.input.get("keep")),
        load_keep(snakemake.input.get("subset_list")),
        str(snakemake.params.copy_mode),
    )
    print("merge_gene_fastas:", stats)


if __name__ == "__main__":
    main()
