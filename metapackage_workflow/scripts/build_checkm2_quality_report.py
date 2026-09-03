#!/usr/bin/env python3
"""
Generate a checkm2-compatible quality report from BUSCO results.

Inputs BUSCO quality control manifest + short_summary JSON files, outputs:
    - TSV with headers [Name | Completeness_General | Contamination |
        Completeness_Specific | Completeness_Model_Used]

Where:
    - Name = genome_acc.fasta (must match CoverM's file_stem for .fasta.gz files)
    - Completeness_General = (complete_count / n_markers) * 100
    - Contamination = 0
    - Completeness_Specific = 0
    - Completeness_Model_Used = "Gradient Boost (General Model)"

Missing or non-numeric values are coerced to 0 so records are retained.
"""

import json
import os
from pathlib import Path
from typing import Dict, Optional, Tuple

from log_config import configure_logger


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake = None  # type: ignore[assignment]


def _extract_group_and_genome_acc(busco_path: Path) -> Tuple[str, str]:
    """Extract group and genome accession from BUSCO result path."""
    # outputs/busco/{group}/{genome_acc}/run_eukaryota_odb12/full_table.tsv
    path_parts = busco_path.parts
    if len(path_parts) >= 4:
        return path_parts[-4], path_parts[-3]
    raise ValueError(f"Unexpected BUSCO path format: {busco_path}")


def _get_short_summary_path(group: str, genome_acc: str) -> Path:
    """Get path to short_summary JSON for given genome."""
    return Path(
        f"outputs/busco/{group}/{genome_acc}/"
        f"short_summary.specific.eukaryota_odb12.{genome_acc}.json"
    )


def _read_short_summary(json_path: Path) -> Optional[Dict]:
    """Read and parse short_summary JSON, return results dict or None on failure."""
    try:
        if not json_path.exists():
            logger.warning("Short summary JSON not found: %s", json_path)
            return None

        with open(json_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)

        return data.get("results", {})
    except Exception as error:
        logger.warning("Could not parse short_summary JSON %s: %s", json_path, error)
        return None


def main() -> None:
    if snakemake is None:
        raise RuntimeError("This script is intended to be run by Snakemake")

    busco_manifest = Path(snakemake.input.busco_manifest)
    output_file = Path(snakemake.output.report)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Building checkm2-compatible quality report from BUSCO results")

    # Read BUSCO quality control manifest
    # Expected format: group\tgenome_acc\tbusco_result\tcomplete_count
    busco_records = {}
    with open(busco_manifest, "r", encoding="utf-8") as handle:
        header = handle.readline().strip()
        logger.debug("BUSCO manifest header: %s", header)

        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            group, genome_acc, busco_result_path, complete_count_str = parts[0:4]
            try:
                complete_count = int(complete_count_str)
            except ValueError:
                logger.warning(
                    "Could not parse complete_count for %s/%s (value=%r); using 0",
                    group,
                    genome_acc,
                    complete_count_str,
                )
                complete_count = 0

            busco_records[(group, genome_acc)] = (busco_result_path, complete_count)

    logger.info(
        "Read %d BUSCO records from manifest",
        len(busco_records),
    )

    # Build quality report: read short_summary for each genome, extract n_markers.
    # Keep all genomes and coerce missing/invalid values to 0.
    quality_records = []
    missing_summaries = 0
    parse_failures = 0

    for (group, genome_acc), (busco_result_path, complete_count) in busco_records.items():
        short_summary_path = _get_short_summary_path(group, genome_acc)
        results = _read_short_summary(short_summary_path)

        completeness = 0.0
        contamination = 0.0

        if results is None:
            missing_summaries += 1
        else:
            n_markers_raw = results.get("n_markers", 0)
            try:
                n_markers = int(n_markers_raw) if n_markers_raw is not None else 0
            except (ValueError, TypeError):
                logger.warning(
                    "Could not parse n_markers for %s/%s (value=%r); using 0",
                    group,
                    genome_acc,
                    n_markers_raw,
                )
                parse_failures += 1
                n_markers = 0

            if n_markers > 0:
                completeness = (complete_count / n_markers) * 100.0
            else:
                logger.warning(
                    "n_markers missing/zero for %s/%s; using completeness=0",
                    group,
                    genome_acc,
                )

        quality_records.append((f"{genome_acc}.fasta", completeness, contamination))

    logger.info(
        "Processed %d genomes: %d quality records, "
        "%d missing summaries, %d parse failures",
        len(busco_records),
        len(quality_records),
        missing_summaries,
        parse_failures,
    )

    # Write checkm2-compatible TSV.
    # Use the 5-column format to avoid a known parser panic in coverm/checkm
    # when a 3-column file is provided.
    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write(
            "Name\tCompleteness_General\tContamination\t"
            "Completeness_Specific\tCompleteness_Model_Used\n"
        )
        for name, completeness, contamination in sorted(quality_records):
            handle.write(
                f"{name}\t{completeness:.6f}\t{contamination:.6f}\t"
                "0.000000\tGradient Boost (General Model)\n"
            )

    logger.info(
        "Wrote %d quality records to %s",
        len(quality_records),
        output_file,
    )


if __name__ == "__main__":
    main()
