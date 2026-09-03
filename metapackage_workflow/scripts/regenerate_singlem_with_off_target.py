#!/usr/bin/env python3

import shutil
import subprocess
from pathlib import Path
from typing import Any

from log_config import configure_logger


logger = configure_logger(__file__)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake: Any = None


def _count_rows(tsv_path: Path) -> int:
    if not tsv_path.exists():
        return 0
    with open(tsv_path, "r", encoding="utf-8") as handle:
        lines = [line for line in handle if line.strip()]
    if len(lines) == 0:
        return 0
    if lines[0].startswith("sequence_id\t"):
        return max(0, len(lines) - 1)
    return len(lines)


def main() -> None:
    spkg = Path(str(snakemake.input.singlem_spkg))
    seqs = Path(str(snakemake.input.seqs))
    taxonomy = Path(str(snakemake.input.taxonomy))
    off_target_seqs = Path(str(snakemake.input.off_target_seqs))
    off_target_taxonomy = Path(str(snakemake.input.off_target_taxonomy))
    output_spkg = Path(str(snakemake.output.spkg))
    sequence_prefix = str(snakemake.params.sequence_prefix)

    if output_spkg.exists():
        shutil.rmtree(output_spkg)

    off_target_count = _count_rows(off_target_taxonomy)

    cmd = [
        "singlem",
        "regenerate",
        "--input-singlem-package",
        str(spkg),
        "--output-singlem-package",
        str(output_spkg),
        "--min-aligned-percent",
        "0",
        "--input-sequences",
        str(seqs),
        "--input-taxonomy",
        str(taxonomy),
        "--sequence-prefix",
        sequence_prefix,
    ]

    if off_target_count > 0 and off_target_seqs.exists():
        cmd.extend(
            [
                "--candidate-decoy-sequences",
                str(off_target_seqs),
                "--candidate-decoy-taxonomy",
                str(off_target_taxonomy),
            ]
        )
        logger.info("Regenerating with %d off-target decoy sequences", off_target_count)
    else:
        cmd.append("--no-further-euks")
        logger.info("Regenerating with no off-target decoys for %s", spkg)

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
