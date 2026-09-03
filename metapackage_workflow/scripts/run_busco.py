#!/usr/bin/env python3

import os
import subprocess
import shutil
from pathlib import Path

from log_config import configure_logger


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)


def main() -> None:
    genome_file = Path(snakemake.input.genome)
    output_file = Path(snakemake.output.full_table)

    lineage = snakemake.params.lineage
    mode = snakemake.params.mode
    out_path_dir = snakemake.params.out_path_dir
    download_path = getattr(snakemake.params, "download_path", None)
    offline = bool(getattr(snakemake.params, "offline", False))

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
        out_path_dir,
        "-f",
        "-q",
    ]

    if download_path:
        cmd.extend(["--download_path", str(download_path)])

    if offline:
        cmd.append("--offline")

    logger.info("Running BUSCO for %s/%s", group, genome_acc)
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error("BUSCO failed for %s/%s", group, genome_acc)
        if result.stdout:
            logger.error("BUSCO stdout:\n%s", result.stdout)
        if result.stderr:
            logger.error("BUSCO stderr:\n%s", result.stderr)
        raise RuntimeError(f"BUSCO failed for {group}/{genome_acc}")

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
