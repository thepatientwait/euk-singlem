#!/usr/bin/env python3
"""Execute notebooks/tree_diagnostics.ipynb with nbconvert for one {derep} value.

The notebook reads its parameters from the JSON file named in TREE_DIAG_PARAMS (written here).
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


def main():
    if redirect_to_snakemake_log:
        redirect_to_snakemake_log()
    html = Path(snakemake.output.html).resolve()
    out_dir = html.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    params = {
        "scripts_dir": str(Path(snakemake.params.scripts_dir).resolve()),
        "alignment": str(Path(snakemake.input.alignment).resolve()),
        "partitions": str(Path(snakemake.input.partitions).resolve()),
        "taxonomy": str(Path(snakemake.input.taxonomy).resolve()),
        "trees": [str(Path(t).resolve()) for t in snakemake.params.trees],
        "merged_dir": str(Path(snakemake.params.merged_dir).resolve()),
        "busco_manifests": [{"manifest": str(Path(m).resolve()), "base": str(Path(b).resolve())}
                            for m, b in snakemake.params.busco_manifests],
        "out_dir": str(out_dir),
        "n_pairs": int(snakemake.params.n_pairs),
        "max_multicopy_cases": int(snakemake.params.max_multicopy_cases),
    }
    params_file = out_dir / "params.json"
    params_file.write_text(json.dumps(params, indent=1))
    env = dict(os.environ, TREE_DIAG_PARAMS=str(params_file), MPLBACKEND="Agg")
    cmd = [sys.executable, "-m", "nbconvert", "--to", "html", "--execute",
           "--ExecutePreprocessor.timeout=-1", str(Path(snakemake.params.notebook).resolve()),
           "--output", html.name, "--output-dir", str(out_dir)]
    print("running:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, env=env)


if __name__ == "__main__":
    main()
