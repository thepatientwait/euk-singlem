#!/usr/bin/env python3

import os
from pathlib import Path
from typing import Optional, Tuple

from log_config import configure_logger


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)


def _extract_group_and_genome_acc(path: Path) -> Tuple[str, str]:
    path_parts = path.parts
    if len(path_parts) >= 5:
        return path_parts[-4], path_parts[-3]
    return "unknown", "unknown"


def _count_complete_markers(busco_result: Path) -> Optional[int]:
    try:
        complete_count = 0
        with open(busco_result, "r", encoding="utf-8") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 1 and parts[1] == "Complete":
                    complete_count += 1
        return complete_count
    except Exception as error:
        logger.warning("Could not parse BUSCO result %s: %s", busco_result, error)
        return None


def main() -> None:
    busco_results = [Path(p) for p in snakemake.input.busco_results]
    output_file = Path(snakemake.output.manifest)
    min_complete_markers = int(snakemake.params.min_complete_markers)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Running BUSCO quality control for %d candidates with threshold > %d Complete markers",
        len(busco_results),
        min_complete_markers,
    )

    quality_controlled = []
    parse_failures = 0
    missing_files = 0
    for busco_result in busco_results:
        if not busco_result.exists() or busco_result.stat().st_size == 0:
            missing_files += 1
            continue

        complete_count = _count_complete_markers(busco_result)
        if complete_count is None:
            parse_failures += 1
            continue

        if complete_count > min_complete_markers:
            group, genome_acc = _extract_group_and_genome_acc(busco_result)
            quality_controlled.append((group, genome_acc, str(busco_result), complete_count))

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write("group\tgenome_acc\tbusco_result\tcomplete_count\n")
        for group, genome_acc, busco_result, complete_count in sorted(quality_controlled):
            handle.write(f"{group}\t{genome_acc}\t{busco_result}\t{complete_count}\n")

    logger.info(
        "Wrote %d quality-controlled BUSCO results to %s (missing_or_empty=%d, parse_failures=%d)",
        len(quality_controlled),
        output_file,
        missing_files,
        parse_failures,
    )


if __name__ == "__main__":
    main()
