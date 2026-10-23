#!/usr/bin/env python3
"""Copy the final trees and their metadata into a new date-time stamped directory.

    <publish_dir>/<YYYY-MM-DD_HHMMSS>[_subset_<name>]/
        <derep>/<copies>/tree.nwk, tree.rerooted.nwk, red_decorated.tree, red_cutoffs.dict,
                         tree_metadata.tsv, taxonomy.revised_reference.tsv,
                         target_novelty_taxonomy.tsv, notebooks/*.ipynb
        conservative_taxonomy.tsv, comparison/tree_comparison.ipynb
        MANIFEST.tsv   (every file with size and md5)
        RUN_INFO.txt   (date, git commit, settings)

Each call makes a new directory, so earlier results are never overwritten. ``<publish_dir>/latest``
is pointed at the newest one when the file system allows symlinks.
"""
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:  # imported from tests
    redirect_to_snakemake_log = None


def md5(path, chunk=1 << 20):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def git_state(repo):
    """(commit, 'clean' | 'dirty') of the repository containing ``repo``; ('unknown', '') if not a repo."""
    try:
        commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return commit, "dirty" if dirty else "clean"
    except (OSError, subprocess.CalledProcessError):
        return "unknown", ""


def snapshot_name(subset="none", now=None):
    stamp = time.strftime("%Y-%m-%d_%H%M%S", now or time.localtime())
    return stamp if subset in ("none", "", None) else f"{stamp}_subset_{subset}"


def publish(files, publish_dir, subset="none", info=None, repo=None, now=None):
    """Copy ``files`` ([(source, relative destination)]) into a new snapshot directory.

    Returns the snapshot directory. Missing sources are an error: a half-published snapshot would
    look complete.
    """
    missing = [str(src) for src, _ in files if not Path(src).is_file()]
    if missing:
        raise FileNotFoundError("cannot publish, missing files: " + ", ".join(missing))

    root = Path(publish_dir).expanduser()
    root.mkdir(parents=True, exist_ok=True)
    snapshot = root / snapshot_name(subset, now)
    suffix = 1
    while snapshot.exists():  # two snapshots in the same second
        suffix += 1
        snapshot = root / f"{snapshot_name(subset, now)}_{suffix}"
    snapshot.mkdir()

    manifest = []
    for source, relative in files:
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        manifest.append((relative, target.stat().st_size, md5(target)))
    with open(snapshot / "MANIFEST.tsv", "w", encoding="utf-8") as handle:
        handle.write("file\tbytes\tmd5\n")
        for relative, size, digest in sorted(manifest):
            handle.write(f"{relative}\t{size}\t{digest}\n")

    commit, state = git_state(repo or Path.cwd())
    lines = [f"published   : {time.strftime('%Y-%m-%d %H:%M:%S', now or time.localtime())}",
             f"subset      : {subset}",
             f"git commit  : {commit} ({state})",
             f"files       : {len(manifest)}"]
    for key, value in (info or {}).items():
        lines.append(f"{key:12}: {json.dumps(value) if not isinstance(value, str) else value}")
    (snapshot / "RUN_INFO.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

    latest = root / "latest"
    try:
        if latest.is_symlink() or latest.exists():
            latest.unlink()
        latest.symlink_to(snapshot.name)
    except OSError:
        pass  # file system without symlinks: the dated directory is still there
    return snapshot


def main():
    if redirect_to_snakemake_log:
        redirect_to_snakemake_log()
    p = snakemake.params
    files = [(src, rel) for src, rel in p.files]
    snapshot = publish(files, p.publish_dir, p.subset, dict(p.info), repo=p.repo)
    Path(snakemake.output.marker).write_text(f"{snapshot}\n", encoding="utf-8")
    print(f"published {len(files)} files to {snapshot}")


if __name__ == "__main__":
    main()
