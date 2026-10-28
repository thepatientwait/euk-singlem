#!/usr/bin/env python3
"""Snapshot of the intron / SingleM-window survey, to repeat and compare after the metapackage changes.

    snapshot_intron_survey.py create  [--name NAME] [--note TEXT] [--subset NAME]   freeze the current survey
    snapshot_intron_survey.py compare OLD NEW                      what changed between two snapshots
    snapshot_intron_survey.py render  SNAPSHOT                     re-run the notebook on a snapshot

A snapshot lives in ``snapshots/intron_survey/<name>/`` (default name: the date and time) and holds

* ``data/``    - the per-gene tables of rule ``intron_window`` (gzipped) and the window files of both variants,
                 the taxonomy mapping, the protein names and the BUSCO manifest the survey used, so it can be
                 re-rendered on its own;
* ``notebook/``- the executed notebook, its HTML export and its tables/figures, if present;
* ``config/``  - ``config.yaml`` and ``paths.yaml``; ``code/`` - the scripts and notebook generator as they were (the working tree is
                 often uncommitted);
* ``inputs.tsv`` - size and md5 of every input (gene FASTAs, HMMs, manifests) so a later run can tell
                 exactly which inputs changed;
* ``RUN_INFO.json`` - when, git commit, tool versions and settings; ``MANIFEST.tsv`` - md5 of every file.

To replicate after the metapackage is updated: re-run the workflow (rule ``intron_window``), take a new
snapshot with ``create``, then ``compare`` the two. Snapshots are large (hundreds of MB) and are not tracked by git.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parent.parent
SNAPSHOT_ROOT = ROOT / "snapshots" / "intron_survey"
SOURCE_DIR = ROOT / "outputs" / "singlem" / "introns"
GENE_INPUT_DIR = ROOT / "outputs" / "singlem" / "input"
HMM_DIR = ROOT / "busco_downloads" / "lineages" / "eukaryota_odb12" / "hmms"
BUSCO_MANIFEST = ROOT / "outputs" / "metadata" / "qc" / "busco_quality_control_results.tsv"
TAXONOMY = ROOT / "outputs" / "metadata" / "taxonomy" / "taxonomy_mapping.tsv"
GENE_NAMES = HMM_DIR.parent / "links_to_ODB12.txt"
SUMMARY = ROOT / "outputs" / "metadata" / "reports" / "window_intron_summary.tsv"
NOTEBOOK_DIR = ROOT / "notebooks" / "rendered"
TABLE_SUFFIXES = ("introns", "genomes", "window_scan", "window_choice")
CODE_FILES = ("scripts/subset_busco_manifest.py", "scripts/intron_utils.py", "scripts/intron_window.py", "scripts/window_intron_summary.py", "scripts/extract_busco_transcripts.py",
              "scripts/snapshot_intron_survey.py", "notebooks/build_intron_survey.py", "notebooks/intron_survey.ipynb")


def use_subset(subset: str) -> None:
    """Point the module at the outputs of a subset run (``outputs/subset_<name>/``)."""
    global SOURCE_DIR, GENE_INPUT_DIR, BUSCO_MANIFEST, SUMMARY, NOTEBOOK_DIR
    root = ROOT / "outputs" / f"subset_{subset}"
    SOURCE_DIR = root / "singlem" / "introns"
    GENE_INPUT_DIR = root / "singlem" / "input"
    BUSCO_MANIFEST = root / "metadata" / "busco_quality_control_results.tsv"
    SUMMARY = root / "metadata" / "reports" / "window_intron_summary.tsv"
    NOTEBOOK_DIR = ROOT / "notebooks" / "rendered" / f"subset_{subset}"


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gzip_copy(source: Path, target: Path) -> None:
    """Gzip ``source`` to ``target`` (deterministic: no timestamp in the header)."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with open(source, "rb") as src, open(target, "wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=6) as out:
        shutil.copyfileobj(src, out, 1 << 20)


def git(*args) -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


def hmmalign_version() -> str:
    """HMMER version line of the hmmalign that the survey uses (the singlem pixi env, else PATH)."""
    candidates = [ROOT / ".pixi" / "envs" / "singlem" / "bin" / "hmmalign", shutil.which("hmmalign")]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            result = subprocess.run([str(candidate), "-h"], capture_output=True, text=True)
            for line in (result.stdout + result.stderr).splitlines():
                if "HMMER" in line:
                    return line.strip("# ").strip()
    return "unavailable"


def read_config() -> dict:
    import yaml
    return yaml.safe_load((ROOT / "config" / "config.yaml").read_text())


def busco_versions(manifest: Path) -> dict:
    """BUSCO versions found in the first line of every full table the survey used (version -> count)."""
    versions = {}
    for path in pl.read_csv(manifest, separator="\t")["busco_result"].to_list():
        try:
            with open(ROOT / path) as handle:
                line = handle.readline().strip().lstrip("# ")
        except OSError:
            line = "missing table"
        versions[line] = versions.get(line, 0) + 1
    return versions


def create(args) -> Path:
    if getattr(args, "subset", None) and args.subset != "none":
        use_subset(args.subset)
    genes = sorted(p.name.split(".")[0] for p in SOURCE_DIR.glob("*.window_choice.tsv"))
    if not genes:
        sys.exit(f"no survey found in {SOURCE_DIR}; run rule intron_window first")
    name = args.name or datetime.datetime.now().strftime("%Y-%m-%d_%H%M%S")
    out = SNAPSHOT_ROOT / name
    if out.exists():
        sys.exit(f"{out} exists; choose another --name")
    data = out / "data"
    print(f"snapshot {out}: {len(genes)} genes")

    # --- survey tables ---------------------------------------------------------------
    for gene in genes:
        for suffix in TABLE_SUFFIXES:
            source = SOURCE_DIR / f"{gene}.{suffix}.tsv"
            if not source.exists():
                sys.exit(f"missing {source}")
            gzip_copy(source, data / f"{gene}.{suffix}.tsv.gz")
        for variant in ("naive", "corrected"):
            shutil.copy2(SOURCE_DIR / f"{gene}.window_position.{variant}.txt", data / f"{gene}.window_position.{variant}.txt")
    if SUMMARY.exists():
        shutil.copy2(SUMMARY, data / SUMMARY.name)

    # --- what the survey was computed from ----------------------------------------------
    gzip_copy(TAXONOMY, data / "taxonomy_mapping.tsv.gz")
    if GENE_NAMES.exists():
        gzip_copy(GENE_NAMES, data / "gene_names.tsv.gz")   # BUSCO id -> protein name, for plot labels
    gzip_copy(BUSCO_MANIFEST, data / "busco_quality_control_results.tsv.gz")
    (out / "config").mkdir(parents=True)
    for config in ("config.yaml", "paths.yaml"):
        shutil.copy2(ROOT / "config" / config, out / "config" / config)
    inputs = []
    for label, path in [("busco_manifest", BUSCO_MANIFEST), ("taxonomy_mapping", TAXONOMY)] + \
                       [("gene_fasta", GENE_INPUT_DIR / f"{g}.faa") for g in genes] + \
                       [("hmm", HMM_DIR / f"{g}.hmm") for g in genes]:
        if path.exists():
            inputs.append({"kind": label, "path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "md5": md5(path)})
    pl.DataFrame(inputs).write_csv(out / "inputs.tsv", separator="\t")

    # --- notebook ----------------------------------------------------------------------
    if NOTEBOOK_DIR.exists():
        shutil.copytree(NOTEBOOK_DIR, out / "notebook")

    # --- the code that produced the survey (the working tree may not be committed) ------------
    for name_ in CODE_FILES:
        source = ROOT / name_
        if source.exists():
            (out / "code" / Path(name_).parent).mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, out / "code" / name_)

    # --- run info ----------------------------------------------------------------------
    choices = pl.concat([pl.read_csv(data / f"{g}.window_choice.tsv.gz", separator="\t") for g in genes])
    n_genomes = pl.read_csv(data / "busco_quality_control_results.tsv.gz", separator="\t")["genome_acc"].n_unique()
    config = read_config()
    info = {
        "name": name, "note": args.note or "",
        "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "git_commit": git("rev-parse", "HEAD"), "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain", "--", ".")),
        "n_genes": len(genes), "n_genomes_in_busco_manifest": n_genomes,
        "subset": getattr(args, "subset", None) or "none",
        "mean_naive_span_fraction": float(choices["naive_span_fraction"].mean()),
        "mean_corrected_span_fraction": float(choices["corrected_span_fraction"].mean()),
        "genes_corrected_over_limit": int(choices["violation"].sum()),
        "intron_settings": (config.get("singlem") or {}).get("introns"),
        "busco_versions": busco_versions(BUSCO_MANIFEST), "hmmer": hmmalign_version(),
        "python": sys.version.split()[0], "polars": pl.__version__,
        "command": " ".join(sys.argv),
        "replicate": [
            "snakemake -s subworkflows/generate_metapackage.smk --use-conda --forcerun intron_window  # after updating the metapackage (add --config subset=family5 for a subset)",
            "python scripts/snapshot_intron_survey.py create --name <new name>",
            f"python scripts/snapshot_intron_survey.py compare {name} <new name>",
        ],
    }
    (out / "RUN_INFO.json").write_text(json.dumps(info, indent=1))

    # --- manifest, readme, latest ---------------------------------------------------------
    files = [p for p in sorted(out.rglob("*")) if p.is_file() and p.name != "MANIFEST.tsv"]
    pl.DataFrame({"path": [str(p.relative_to(out)) for p in files], "bytes": [p.stat().st_size for p in files],
                  "md5": [md5(p) for p in files]}).write_csv(out / "MANIFEST.tsv", separator="\t")
    (out / "README.md").write_text(f"""# Intron survey snapshot {name}

{args.note or ''}

Created {info['created']} at git commit {info['git_commit'][:10]}{' (working tree had uncommitted changes)' if info['git_dirty'] else ''}.
{len(genes)} genes, {n_genomes} genomes in the BUSCO manifest.

* Re-render the notebook on this data:  `python scripts/snapshot_intron_survey.py render {name}`
* Compare with another snapshot:        `python scripts/snapshot_intron_survey.py compare {name} <other>`
* What the survey used: `inputs.tsv` (md5 of gene FASTAs, HMMs, manifests), `config/`, `RUN_INFO.json`.
* Tables are gzipped TSV in `data/`; `data/<gene>.window_position.{{naive,corrected}}.txt` are the windows of the two variants.
""")
    link = SNAPSHOT_ROOT / "latest"
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(name)
    size = sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) / 1e6
    print(f"done: {len(files)} files, {size:.0f} MB; 'latest' -> {name}")
    return out


def resolve(name: str) -> Path:
    path = Path(name)
    if not path.exists():
        path = SNAPSHOT_ROOT / name
    if not (path / "RUN_INFO.json").exists():
        sys.exit(f"{name}: not a snapshot (no RUN_INFO.json)")
    return path.resolve()


# names used by snapshots taken before legacy -> naive and final -> corrected (best_coverage is the naive window)
OLD_COLUMNS = {"best_coverage_start": "naive_start", "best_coverage_span_fraction": "naive_span_fraction",
               "final_start": "corrected_start", "final_span_fraction": "corrected_span_fraction",
               "final_score_ratio": "corrected_score_ratio"}


def read_choices(snapshot: Path) -> pl.DataFrame:
    frame = pl.concat([pl.read_csv(p, separator="\t") for p in sorted((snapshot / "data").glob("*.window_choice.tsv.gz"))])
    return frame.rename({old: new for old, new in OLD_COLUMNS.items() if old in frame.columns and new not in frame.columns})


def compare(args) -> None:
    old, new = resolve(args.old), resolve(args.new)
    a, b = read_choices(old), read_choices(new)
    info_a, info_b = json.loads((old / "RUN_INFO.json").read_text()), json.loads((new / "RUN_INFO.json").read_text())
    print(f"old {old.name}  ({info_a['created']}, {info_a['n_genes']} genes, {info_a['n_genomes_in_busco_manifest']} genomes)")
    print(f"new {new.name}  ({info_b['created']}, {info_b['n_genes']} genes, {info_b['n_genomes_in_busco_manifest']} genomes)")

    # keyed by file name, so a subset run (other directory) compares with a full run
    inputs_a, inputs_b = (pl.read_csv(s / "inputs.tsv", separator="\t").with_columns(pl.col("path").str.split("/").list.last().alias("file"))
                          for s in (old, new))
    joined = inputs_a.join(inputs_b, on=["kind", "file"], how="full", suffix="_new", coalesce=True)
    changed = joined.filter((pl.col("md5") != pl.col("md5_new")) | pl.col("md5").is_null() | pl.col("md5_new").is_null())
    print("\ninputs that changed:", changed.group_by("kind").len().sort("kind").to_dicts() or "none")

    both = a.join(b, on="busco_id", suffix="_new", how="inner")
    only_old, only_new = sorted(set(a["busco_id"]) - set(b["busco_id"])), sorted(set(b["busco_id"]) - set(a["busco_id"]))
    print(f"\ngenes: {both.height} in both, {len(only_old)} only in old, {len(only_new)} only in new")
    summary = [
        ("corrected window moved", int((both["corrected_start"] != both["corrected_start_new"]).sum())),
        ("naive window moved", int((both["naive_start"] != both["naive_start_new"]).sum())),
        ("genes over the limit (corrected): old", int(a["violation"].sum())), ("genes over the limit (corrected): new", int(b["violation"].sum())),
        ("mean corrected span fraction: old", round(float(a["corrected_span_fraction"].mean()), 5)),
        ("mean corrected span fraction: new", round(float(b["corrected_span_fraction"].mean()), 5)),
        ("mean naive span fraction: old", round(float(a["naive_span_fraction"].mean()), 5)),
        ("mean naive span fraction: new", round(float(b["naive_span_fraction"].mean()), 5)),
    ]
    for key, value in summary:
        print(f"  {key:44s} {value}")

    diff = (both.with_columns((pl.col("corrected_span_fraction_new") - pl.col("corrected_span_fraction")).alias("corrected_span_change"),
                              (pl.col("n_mapped_new") - pl.col("n_mapped")).alias("genomes_mapped_change"))
            .select("busco_id", "corrected_start", "corrected_start_new", "corrected_span_fraction", "corrected_span_fraction_new", "corrected_span_change",
                    "naive_start", "naive_start_new", "n_mapped", "n_mapped_new", "genomes_mapped_change"))
    target = Path(args.out) if args.out else new / f"comparison_with_{old.name}.tsv"
    diff.write_csv(target, separator="\t")
    with pl.Config(tbl_rows=15, tbl_cols=12, tbl_width_chars=200):
        print("\nlargest changes in corrected window span fraction:")
        print(diff.sort(pl.col("corrected_span_change").abs(), descending=True).head(10))
    print(f"\nper-gene comparison written to {target}")


def render(args) -> None:
    snapshot = resolve(args.snapshot)
    out = snapshot / "notebook_rendered"
    out.mkdir(exist_ok=True)
    info = json.loads((snapshot / "RUN_INFO.json").read_text())
    params = {"intron_dir": str(snapshot / "data"), "taxonomy": str(snapshot / "data" / "taxonomy_mapping.tsv.gz"),
              "gene_names": str(snapshot / "data" / "gene_names.tsv.gz"),
              "max_span_fraction": (info.get("intron_settings") or {}).get("max_span_fraction", 0.05),
              "build_naive": bool((info.get("intron_settings") or {}).get("build_naive", False)),
              "subset": info.get("subset", "none"), "out_dir": str(out / "tables")}
    # the notebook as it was when the snapshot was taken (older snapshots use other column names)
    template = snapshot / "code" / "notebooks" / "intron_survey.ipynb"
    if not template.exists():
        template = ROOT / "notebooks" / "intron_survey.ipynb"
    params_file = out / "intron_survey.params.json"
    params_file.write_text(json.dumps(params, indent=1))
    env = dict(os.environ, NB_PARAMS=str(params_file), MPLBACKEND="Agg")
    cmd = [sys.executable, "-m", "nbconvert", "--to", "notebook", "--execute", "--ExecutePreprocessor.timeout=-1",
           str(template), "--output", "intron_survey.ipynb", "--output-dir", str(out)]
    subprocess.run(cmd, check=True, env=env)
    subprocess.run([sys.executable, "-m", "nbconvert", "--to", "html", str(out / "intron_survey.ipynb")], check=True)
    print(f"rendered: {out / 'intron_survey.html'}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("create", help="freeze the current survey")
    p.add_argument("--name"); p.add_argument("--note")
    p.add_argument("--subset", default="none", help="snapshot a subset run (outputs/subset_<name>/) instead of the full survey")
    p.set_defaults(func=create)
    p = sub.add_parser("compare", help="compare two snapshots")
    p.add_argument("old"); p.add_argument("new"); p.add_argument("--out")
    p.set_defaults(func=compare)
    p = sub.add_parser("render", help="re-run the notebook on a snapshot")
    p.add_argument("snapshot")
    p.set_defaults(func=render)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
