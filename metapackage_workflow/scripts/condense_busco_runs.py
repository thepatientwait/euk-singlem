#!/usr/bin/env python3
"""Condense finished BUSCO result directories (see busco_run.py): ~540 files / 170 MB per genome become 3 small files.

    condense_busco_runs.py outputs/busco [outputs/busco ...]            # report what would be done (nothing is changed)
    condense_busco_runs.py outputs/busco --apply [--threads 8]          # do it

Every ``run_*`` directory found under the given folders that still has the original per-gene files is condensed: the
concatenated ``proteins.faa.gz`` / ``models.gff.gz`` / ``cds.fna.gz`` are written first and only then are the large
intermediates (``busco_sequences``, ``miniprot_output`` with the 167 MB ``ref.mpi``, ``hmmer_output``, ``metaeuk_output``) removed.
A directory is skipped, and reported, if it has no ``full_table.tsv`` or the condensed files could not be read back.
Repeatable and safe to interrupt: a directory is either untouched or complete. Downstream results do not change (the
readers give the same sequences and gene models for both layouts), so Snakemake sees no reason to rerun anything.
"""
import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from busco_run import BuscoRun, condense_run  # noqa: E402


def find_runs(roots):
    """``run_*`` directories containing a full_table.tsv, below the roots (depth <= 4, no deeper walk into runs)."""
    found = []
    for root in roots:
        for dirpath, dirnames, filenames in os.walk(root):
            name = os.path.basename(dirpath)
            if name.startswith("run_") and "full_table.tsv" in filenames:
                found.append(Path(dirpath))
                dirnames[:] = []
    return sorted(found)


def dir_size(path: Path) -> int:
    total = 0
    for sub, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(sub, f))
            except OSError:
                pass
    return total


def process(run_dir: Path, apply: bool):
    layout = BuscoRun(run_dir).layout
    if layout == "condensed":
        left = [d for d in ("busco_sequences", "miniprot_output", "hmmer_output", "metaeuk_output") if (run_dir / d).exists()]
        if left and apply:
            condense_run(run_dir)
        return run_dir, "condensed" if not left else "tidied", 0
    if layout == "none":
        return run_dir, "no_sequences", 0
    if not apply:
        return run_dir, "would_condense", dir_size(run_dir.parent)
    before = dir_size(run_dir.parent)
    condense_run(run_dir)
    check = BuscoRun(run_dir)
    if check.layout != "condensed" or not (run_dir / "proteins.faa.gz").exists():
        return run_dir, "FAILED", 0
    return run_dir, "condensed_now", before - dir_size(run_dir.parent)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("roots", nargs="+")
    parser.add_argument("--apply", action="store_true", help="really condense (default: only report)")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    runs = find_runs(args.roots)
    print(f"{len(runs)} BUSCO run directories found", flush=True)
    counts, freed = {}, 0
    with ThreadPoolExecutor(max(1, args.threads)) as pool:
        for i, (run_dir, status, size) in enumerate(pool.map(lambda r: process(r, args.apply), runs), 1):
            counts[status] = counts.get(status, 0) + 1
            freed += size
            if status == "FAILED":
                print(f"FAILED: {run_dir}", file=sys.stderr)
            if i % 500 == 0:
                print(f"  {i}/{len(runs)} {counts}", flush=True)
    verb = "freed" if args.apply else "could free"
    print(f"{counts}; {verb} {freed / 1e9:.1f} GB" + ("" if args.apply else " (run again with --apply)"))
    return 1 if counts.get("FAILED") else 0


if __name__ == "__main__":
    sys.exit(main())
