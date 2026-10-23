#!/usr/bin/env python3
"""Execute a notebook with nbconvert. Its parameters are ``params.nb_params`` (a dict with absolute
paths), written to ``<output dir>/<name>.params.json`` and handed over through NB_PARAMS."""
import json
import os
import subprocess
import sys
from pathlib import Path

try:
    from log_redirect import redirect_to_snakemake_log
except ImportError:
    redirect_to_snakemake_log = None


def main():
    if redirect_to_snakemake_log:
        redirect_to_snakemake_log()
    html = Path(snakemake.output.html).resolve()
    html.parent.mkdir(parents=True, exist_ok=True)
    params = dict(snakemake.params.nb_params)
    params.setdefault("out_dir", str(html.parent))
    params_file = html.with_suffix(".params.json")
    params_file.write_text(json.dumps(params, indent=1))
    env = dict(os.environ, NB_PARAMS=str(params_file), MPLBACKEND="Agg")
    cmd = [sys.executable, "-m", "nbconvert", "--to", "html", "--execute",
           "--ExecutePreprocessor.timeout=-1", str(Path(snakemake.params.notebook).resolve()),
           "--output", html.name, "--output-dir", str(html.parent)]
    print("running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)


if __name__ == "__main__":
    main()
