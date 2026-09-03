#!/usr/bin/env python3

import os
import shlex
import shutil
import subprocess
from pathlib import Path


def main() -> None:
    genome = Path(str(snakemake.input.genome))
    done_file = Path(str(snakemake.output.done))

    outdir = Path(str(snakemake.params.outdir))
    command_template = str(snakemake.params.command_template)
    db_path = str(getattr(snakemake.params, "db_path", "") or "").strip()
    threads = int(snakemake.threads)

    # Always start from a clean per-genome output directory.
    # Partial epa-ng/MetaEuk artifacts can cause deterministic failures.
    if outdir.exists():
        shutil.rmtree(outdir, ignore_errors=True)
    outdir.mkdir(parents=True, exist_ok=True)
    done_file.parent.mkdir(parents=True, exist_ok=True)

    command = command_template.format(
        genome=str(genome),
        outdir=str(outdir),
        db=db_path,
        prefix=snakemake.wildcards.genome_acc,
        threads=threads,
    )

    env = dict(os.environ)
    if db_path:
        env["EUKCC2_DB"] = db_path

    if "--db" not in command and not env.get("EUKCC2_DB"):
        raise ValueError(
            "EukCC database is not configured. Set eukcc.db in tree_workflow/config.yaml or include --db {db} in eukcc.command_template."
        )

    try:
        subprocess.run(shlex.split(command), check=True, env=env)
    except subprocess.CalledProcessError:
        # EukCC can leave partial epa-ng/MetaEuk outputs that break retries.
        # Remove the per-genome output directory so the next attempt starts clean.
        if outdir.exists():
            shutil.rmtree(outdir, ignore_errors=True)
        raise

    done_file.touch()


if __name__ == "__main__":
    main()
