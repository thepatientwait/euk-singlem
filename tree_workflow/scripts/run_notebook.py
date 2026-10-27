#!/usr/bin/env python3
"""Execute a notebook template and write the executed copy (``.ipynb``, outputs embedded).

The template in ``notebooks/`` is never modified, so results do not end up in git. Parameters
(``params.nb_params``, a dict with absolute paths) are written to ``<notebook>.params.json`` next
to the executed notebook and handed to it through the NB_PARAMS environment variable.

Cells tagged ``findings-top`` (the key findings, which can only be computed after the analyses) are
moved to just below the title cell of the executed copy, so the verdict is the first thing read.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:
    redirect_to_snakemake_log = None


def move_tagged_to_top(notebook_path, tag="findings-top"):
    """Move the cells tagged ``tag`` (key findings, computed late) to just after the title cell."""
    path = Path(notebook_path)
    nb = json.loads(path.read_text())
    cells = nb["cells"]
    tagged = [c for c in cells if tag in c.get("metadata", {}).get("tags", [])]
    if not tagged:
        return 0
    rest = [c for c in cells if not any(c is t for t in tagged)]
    nb["cells"] = rest[:1] + tagged + rest[1:]
    path.write_text(json.dumps(nb, indent=1))
    return len(tagged)


def main():
    if redirect_to_snakemake_log:
        redirect_to_snakemake_log()
    executed = Path(snakemake.output.notebook).resolve()
    executed.parent.mkdir(parents=True, exist_ok=True)
    params = dict(snakemake.params.nb_params)
    params.setdefault("out_dir", str(executed.parent))
    params_file = executed.with_suffix(".params.json")
    params_file.write_text(json.dumps(params, indent=1))
    env = dict(os.environ, NB_PARAMS=str(params_file), MPLBACKEND="Agg")
    cmd = [sys.executable, "-m", "nbconvert", "--to", "notebook", "--execute",
           "--ExecutePreprocessor.timeout=-1", "--ExecutePreprocessor.startup_timeout=300", str(Path(snakemake.params.notebook).resolve()),
           "--output", executed.name, "--output-dir", str(executed.parent)]
    print("running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)
    moved = move_tagged_to_top(executed)
    if moved:
        print(f"moved {moved} key-findings cells to the top of {executed.name}", flush=True)


if __name__ == "__main__":
    main()
