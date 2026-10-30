#!/usr/bin/env python3

import os
import subprocess
import shutil
import tempfile
from pathlib import Path

from log_config import configure_logger
from busco_run import condense_run


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)


def finalize(result_dir: Path, final_dir: Path, condense: bool = True) -> dict:
    """Move a finished BUSCO result into its place, condensed (see busco_run.py) unless ``condense`` is False.

    BUSCO is run in a scratch folder (the node's TMPDIR) so its large intermediates (a 167 MB miniprot index per genome,
    hundreds of per-gene files) never reach the shared disk. ``final_dir`` is replaced as a whole.
    """
    removed = {"removed_files": 0, "removed_bytes": 0}
    if condense:
        for run_dir in sorted(result_dir.glob("run_*")):
            stats = condense_run(run_dir)
            removed["removed_files"] += stats["removed_files"]
            removed["removed_bytes"] += stats["removed_bytes"]
    if final_dir.exists():
        shutil.rmtree(final_dir)
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(result_dir), str(final_dir))
    return removed


def main() -> None:
    genome_file = Path(snakemake.input.genome)
    output_file = Path(snakemake.output.full_table)

    lineage = snakemake.params.lineage
    mode = snakemake.params.mode
    out_path_dir = snakemake.params.out_path_dir
    download_path = getattr(snakemake.params, "download_path", None)
    offline = bool(getattr(snakemake.params, "offline", False))
    condense = bool(getattr(snakemake.params, "condense", True))

    group = snakemake.wildcards.group
    genome_acc = snakemake.wildcards.genome_acc
    threads = snakemake.threads

    output_file.parent.mkdir(parents=True, exist_ok=True)

    if shutil.which("busco") is None:
        raise FileNotFoundError(
            "BUSCO executable not found on PATH. Ensure this rule runs with a conda environment that includes 'busco'."
        )

    if not genome_file.exists() or genome_file.stat().st_size == 0:
        logger.warning(
            "Skipping BUSCO for %s/%s because genome input is missing or empty: %s",
            group,
            genome_acc,
            genome_file,
        )
        output_file.touch()
        return

    scratch = Path(tempfile.mkdtemp(prefix=f"busco_{genome_acc}_"))   # TMPDIR of the job (node-local on the cluster)
    cmd = [
        "busco",
        "-i",
        str(genome_file),
        "-l",
        lineage,
        "-m",
        mode,
        "-o",
        genome_acc,
        "--cpu",
        str(threads),
        "--out_path",
        str(scratch),
        "-f",
        "-q",
    ]

    if download_path:
        cmd.extend(["--download_path", str(download_path)])

    if offline:
        cmd.append("--offline")

    logger.info("Running BUSCO for %s/%s in %s", group, genome_acc, scratch)
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error("BUSCO failed for %s/%s", group, genome_acc)
            if result.stdout:
                logger.error("BUSCO stdout:\n%s", result.stdout)
            if result.stderr:
                logger.error("BUSCO stderr:\n%s", result.stderr)
            raise RuntimeError(f"BUSCO failed for {group}/{genome_acc}")

        produced = scratch / genome_acc
        if not produced.is_dir():
            raise RuntimeError(f"BUSCO reported success but wrote no {produced}")
        removed = finalize(produced, Path(out_path_dir) / genome_acc, condense)
        logger.info("%s/%s: result moved to %s%s", group, genome_acc, Path(out_path_dir) / genome_acc,
                    f" (condensed: {removed['removed_files']} files, {removed['removed_bytes'] / 1e6:.0f} MB not kept)"
                    if condense else "")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    if not output_file.exists() or output_file.stat().st_size == 0:
        logger.warning(
            "BUSCO finished but expected output is missing or empty for %s/%s; touching %s",
            group,
            genome_acc,
            output_file,
        )
        output_file.touch()


if __name__ == "__main__":
    main()
