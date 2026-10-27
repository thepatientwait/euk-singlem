#!/usr/bin/env python3
"""Summarise which expected genomes are not (validly) on disk, and why.

Read-only apart from the report it writes. Run from metapackage_workflow/ at any time:

    python scripts/download_summary.py [--report PATH] [--max-list N]

It is also called when the download phase exits (success or failure), so it is clear what is
still outstanding. Categories:

  jgi_restore_in_progress  JGI file on tape, restore requested but not finished
  jgi_restore_refused      JGI refused/could not restore the file
  jgi_on_tape              JGI file on tape, no restore requested yet
  zero_byte_marker         failure marker (assembly permanently unavailable, e.g. suppressed)
  invalid_gzip             file on disk fails the gzip check (see failed_genome_cleanup_report.tsv)
  missing_other            no file and no known reason (failed or not yet run download job)
"""
import argparse
import collections
import csv
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DEFAULT_REPORT = "outputs/metadata/reports/download_summary.tsv"
REPORT_DIR = Path("outputs/metadata/reports")


def _read_tsv(path):
    path = Path(path)
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_config():
    import yaml
    config = {}
    for name in ("config.yaml", "paths.yaml"):
        with open(ROOT / "config" / name, encoding="utf-8") as handle:
            config.update(yaml.safe_load(handle) or {})
    return config


def collect(groups):
    """Return a list of (category, group, genome_acc, detail) for every outstanding genome."""
    from helpers import get_expected_genome_paths, get_jgi_not_ready_paths

    expected = get_expected_genome_paths(groups)
    try:
        jgi_not_ready = get_jgi_not_ready_paths(groups)
    except OSError:
        jgi_not_ready = {}

    # states from the restore step: genome_acc -> state (restore_refused / restore_in_progress)
    restore_state = {}
    for group in groups:
        for row in _read_tsv(REPORT_DIR / f"jgi_not_ready.{group}.tsv"):
            restore_state[(group, row["genome_acc"])] = row.get("state", "")

    invalid = {
        (row["group"], row["genome_acc"]): row["reason"]
        for row in _read_tsv(REPORT_DIR / "failed_genome_cleanup_report.tsv")
    }

    outstanding = []
    for path in expected:
        p = Path(path)
        group, acc = p.parent.name, p.name.replace(".fasta.gz", "")
        if p.exists():
            if p.stat().st_size == 0:
                outstanding.append(("zero_byte_marker", group, acc, "empty file"))
            elif (group, acc) in invalid:
                outstanding.append(("invalid_gzip", group, acc, invalid[(group, acc)]))
            continue
        if path in jgi_not_ready:
            state = restore_state.get((group, acc), "")
            if state == "restore_refused":
                category = "jgi_restore_refused"
            elif state == "restore_in_progress" or jgi_not_ready[path] == "RESTORE_REGISTERED":
                category = "jgi_restore_in_progress"
            else:
                category = "jgi_on_tape"
            outstanding.append((category, group, acc, jgi_not_ready[path]))
        else:
            outstanding.append(("missing_other", group, acc, ""))
    return len(expected), outstanding


def summarise(report_path=DEFAULT_REPORT, max_list=10, out=sys.stdout):
    config = _load_config()
    total, outstanding = collect(config["groups"])

    report = Path(report_path)
    report.parent.mkdir(parents=True, exist_ok=True)
    with open(report, "w", encoding="utf-8") as handle:
        handle.write("category\tgroup\tgenome_acc\tdetail\n")
        for row in sorted(outstanding):
            handle.write("\t".join(row) + "\n")

    by_cat = collections.defaultdict(list)
    for category, group, acc, _ in outstanding:
        by_cat[category].append((group, acc))

    bar = "=" * 70
    print(bar, file=out)
    print("DOWNLOAD SUMMARY", file=out)
    print(bar, file=out)
    print(f"Expected genomes: {total}; outstanding: {len(outstanding)}", file=out)
    if not outstanding:
        print("All expected genomes are on disk.", file=out)
    for category in sorted(by_cat):
        items = by_cat[category]
        per_group = collections.Counter(group for group, _ in items)
        groups_text = ", ".join(f"{g}={n}" for g, n in sorted(per_group.items()))
        print(f"\n{category}: {len(items)} ({groups_text})", file=out)
        shown = ", ".join(acc for _, acc in sorted(items)[:max_list])
        more = f", ... (+{len(items) - max_list} more)" if len(items) > max_list else ""
        print(f"  {shown}{more}", file=out)
    print(f"\nFull list: {report}", file=out)
    if any(c.startswith("jgi_restore") or c == "jgi_on_tape" for c in by_cat):
        print("JGI genomes still on tape: refresh the token in config/secrets.yaml and re-run the "
              "download phase (jobs poll until restored).", file=out)
    print(bar, file=out)
    return total, outstanding


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--report", default=DEFAULT_REPORT)
    parser.add_argument("--max-list", type=int, default=10)
    args = parser.parse_args()
    os.chdir(ROOT)
    summarise(args.report, args.max_list)


if __name__ == "__main__":
    main()
