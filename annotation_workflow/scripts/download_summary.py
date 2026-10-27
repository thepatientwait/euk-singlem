#!/usr/bin/env python3
"""Print and write what of the configured OrthoDB files is (not) on disk. Read-only apart from the report.

    python scripts/download_summary.py          (run from annotation_workflow/)
Also called when the workflow exits (success or failure), like metapackage_workflow/scripts/download_summary.py.
"""
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent


def collect(config):
    raw = Path(config["paths"]["orthodb_raw"])
    rows = []
    for key, entry in config["orthodb"]["files"].items():
        path = raw / entry["name"]
        partial = path.with_name(path.name + ".partial")
        if path.exists() and path.stat().st_size >= int(config.get("min_download_bytes", 1)):
            rows.append((key, entry["name"], "present", path.stat().st_size))
        elif path.exists():
            rows.append((key, entry["name"], "too_small", path.stat().st_size))
        elif partial.exists():
            rows.append((key, entry["name"], "partial", partial.stat().st_size))
        else:
            rows.append((key, entry["name"], "missing", 0))
    return rows


def summarise(config=None, out=sys.stdout):
    if config is None:
        config = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())
    rows = collect(config)
    report = Path(config["paths"]["summary"])
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("key\tfile\tstatus\tbytes\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    bar = "=" * 70
    print(bar + "\nORTHODB DOWNLOAD SUMMARY\n" + bar, file=out)
    local = config["orthodb"].get("local_dir")
    if local:
        print(f"Using local copy: {local}", file=out)
    for key, name, status, size in rows:
        print(f"  {status:9} {name} ({size / 1e9:.2f} GB)", file=out)
    bad = [r for r in rows if r[2] != "present"]
    print(f"\n{len(rows) - len(bad)} of {len(rows)} files present. Full list: {report}", file=out)
    if bad:
        print("Re-run the workflow to resume; if OrthoDB is unreachable from the cluster, download the files "
              "elsewhere and set orthodb.local_dir.", file=out)
    print(bar, file=out)
    return rows


if __name__ == "__main__":
    import os
    os.chdir(ROOT)
    summarise()
