#!/usr/bin/env python3

"""Filter SingleM packages with invalid stored window positions.

A package is valid when:
    singlem_hmm_position + (singlem_window_size / 3) <= HMM LENG

The script writes:
1. A newline-delimited manifest of valid package paths.
2. A TSV report containing validation status for all packages.
"""

import json
from pathlib import Path

import sys
sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger

logger = configure_logger(__file__)


def read_hmm_length(spkg_path: Path, graftm_package_path: str) -> int:
    hmm_files = list(spkg_path.glob("*/*.hmm"))
    if len(hmm_files) == 0:
        raise ValueError(f"No HMM files found in {spkg_path}")

    # Prefer the marker alignment HMM, e.g. <marker>/<marker>.hmm.
    expected_hmm = spkg_path / graftm_package_path / f"{graftm_package_path}.hmm"
    if expected_hmm.exists():
        hmm_path = expected_hmm
    else:
        # Fallback: choose a non-search helper HMM if present.
        non_search_hmms = [h for h in hmm_files if "search" not in h.name]
        if len(non_search_hmms) == 1:
            hmm_path = non_search_hmms[0]
        elif len(hmm_files) == 1:
            hmm_path = hmm_files[0]
        else:
            raise ValueError(
                "Could not uniquely determine marker HMM in {}. Found: {}".format(
                    spkg_path,
                    ", ".join(sorted(str(h.relative_to(spkg_path)) for h in hmm_files)),
                )
            )

    with open(hmm_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("LENG"):
                return int(line.split()[1])

    raise ValueError(f"Could not find LENG in HMM file: {hmm_path}")


def validate_package(spkg_path: Path) -> dict:
    contents_path = spkg_path / "CONTENTS.json"
    with open(contents_path, "r", encoding="utf-8") as handle:
        contents = json.load(handle)

    start = int(contents["singlem_hmm_position"])
    window_nt = int(contents["singlem_window_size"])

    if window_nt % 3 != 0:
        raise ValueError(
            f"Window size is not divisible by 3 for {spkg_path}: {window_nt}"
        )

    stretch_aa = window_nt // 3
    graftm_package_path = str(contents["graftm_package_path"])
    hmm_len = read_hmm_length(spkg_path, graftm_package_path)
    max_start = hmm_len - stretch_aa
    margin = hmm_len - (start + stretch_aa)
    is_valid = start <= max_start

    return {
        "spkg_path": str(spkg_path),
        "singlem_hmm_position": start,
        "singlem_window_size_nt": window_nt,
        "stretch_aa": stretch_aa,
        "hmm_length": hmm_len,
        "max_valid_start": max_start,
        "margin": margin,
        "status": "valid" if is_valid else "invalid",
    }


def main():
    package_paths = [Path(p) for p in snakemake.input.final_packages]
    valid_manifest_path = Path(str(snakemake.output.valid_packages_manifest))
    invalid_report_path = Path(str(snakemake.output.invalid_packages_report))

    logger.info("Validating %d SingleM packages", len(package_paths))

    results = []
    for spkg_path in package_paths:
        if not spkg_path.exists():
            raise FileNotFoundError(f"Expected package does not exist: {spkg_path}")
        results.append(validate_package(spkg_path))

    valid_packages = [row["spkg_path"] for row in results if row["status"] == "valid"]
    invalid_packages = [row for row in results if row["status"] == "invalid"]

    valid_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    invalid_report_path.parent.mkdir(parents=True, exist_ok=True)

    with open(valid_manifest_path, "w", encoding="utf-8") as handle:
        for pkg in valid_packages:
            handle.write(f"{pkg}\n")

    with open(invalid_report_path, "w", encoding="utf-8") as handle:
        handle.write(
            "spkg_path\tstatus\tsinglem_hmm_position\tsinglem_window_size_nt\t"
            "stretch_aa\thmm_length\tmax_valid_start\tmargin\n"
        )
        for row in results:
            handle.write(
                f"{row['spkg_path']}\t{row['status']}\t"
                f"{row['singlem_hmm_position']}\t{row['singlem_window_size_nt']}\t"
                f"{row['stretch_aa']}\t{row['hmm_length']}\t"
                f"{row['max_valid_start']}\t{row['margin']}\n"
            )

    logger.info("Valid packages: %d", len(valid_packages))
    logger.info("Invalid packages filtered out: %d", len(invalid_packages))

    if invalid_packages:
        logger.warning("Invalid packages written to report: %s", invalid_report_path)

    if len(valid_packages) == 0:
        raise ValueError("No valid packages remain after filtering")


if __name__ == "__main__":
    main()
