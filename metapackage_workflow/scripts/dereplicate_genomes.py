#!/usr/bin/env python3
"""Dereplicate genomes and transcriptomes with galah (skani ANI) on their marker-gene CDS.

Each genome is represented by the spliced CDS of its BUSCO marker genes (``extract_marker_cds.py``), so a
transcriptome and a genome of the same species compare on the same sequence. Representatives are chosen by the
quality report (BUSCO completeness; transcriptomes are penalised so a genome wins a shared cluster).
Outputs keep the format coverm produced: the cluster definition (representative path, member path),
a list of representative accessions and the BUSCO manifest filtered to them.
"""

import subprocess
from pathlib import Path
from typing import Any

from log_config import configure_logger


logger = configure_logger(__file__)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake: Any = None


MIN_CDS_FILE_BYTES = 200  # an empty gzip FASTA (no marker CDS) is about 20 bytes


def _read_manifest_marker_cds(manifest_file: Path, template: str) -> list[str]:
    """Marker-CDS FASTA paths of the genomes in the BUSCO QC manifest (columns ``group``, ``genome_acc``).

    Genomes without marker CDS (missing or empty file) cannot be clustered and are left out, with a warning.
    """
    paths: list[str] = []
    empty = 0
    with open(manifest_file, "r", encoding="utf-8") as handle:
        _header = handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            path = Path(template.format(group=parts[0], genome_acc=parts[1]))
            if path.exists() and path.stat().st_size >= MIN_CDS_FILE_BYTES:
                paths.append(str(path))
            else:
                empty += 1
    if empty:
        logger.warning("%d genomes have no marker CDS and are not dereplicated", empty)
    return sorted(set(paths))


def _write_genome_list(genome_paths: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for genome_path in genome_paths:
            handle.write(f"{genome_path}\n")


def _run_galah_cluster(
    genome_list_file: Path,
    quality_report: Path,
    clusters_out: Path,
    precluster_ani: Any,
    ani: Any,
    threads: int,
    min_aligned_fraction: Any,
    precluster_method: str,
    cluster_method: str,
) -> None:
    cmd = [
        "galah",
        "cluster",
        "--genome-fasta-list",
        str(genome_list_file),
        "--checkm2-quality-report",
        str(quality_report),
        "--precluster-ani",
        str(precluster_ani),
        "--ani",
        str(ani),
        "--min-aligned-fraction",
        str(min_aligned_fraction),
        "--precluster-method",
        str(precluster_method),
        "--cluster-method",
        str(cluster_method),
        "--output-cluster-definition",
        str(clusters_out),
        "--threads",
        str(threads),
    ]
    logger.info("Running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def _extract_representatives(clusters_file: Path) -> list[str]:
    representatives = set()
    with open(clusters_file, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rep_path = line.split("\t")[0].split()[0]
            rep_name = Path(rep_path).name.replace(".fasta.gz", "").replace(".gz", "")  # marker CDS file = <acc>.fasta.gz
            representatives.add(rep_name)
    return sorted(representatives)


def _write_representatives(representatives: list[str], output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w", encoding="utf-8") as handle:
        for rep in representatives:
            handle.write(f"{rep}\n")


def _filter_busco_manifest(
    input_manifest: Path,
    output_manifest: Path,
    representatives: set[str],
) -> None:
    output_manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(input_manifest, "r", encoding="utf-8") as infile, open(
        output_manifest,
        "w",
        encoding="utf-8",
    ) as outfile:
        header = infile.readline()
        outfile.write(header)
        for line in infile:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if parts[1] in representatives:
                outfile.write(line)


def main() -> None:
    if snakemake is None:
        raise RuntimeError("This script must be executed by Snakemake")

    quality_report = Path(str(snakemake.input.quality_report))
    busco_manifest = Path(str(snakemake.input.busco_manifest))

    clusters_out = Path(str(snakemake.output.clusters))
    representatives_out = Path(str(snakemake.output.representatives))
    filtered_manifest_out = Path(str(snakemake.output.filtered_manifest))

    precluster_ani = snakemake.params.precluster_ani
    ani = snakemake.params.ani
    min_aligned_fraction = snakemake.params.min_aligned_fraction
    precluster_method = snakemake.params.precluster_method
    cluster_method = snakemake.params.cluster_method
    cds_template = snakemake.params.cds_template
    threads = int(snakemake.threads)

    clusters_out.parent.mkdir(parents=True, exist_ok=True)

    genome_list_file = clusters_out.parent / "genome_list.tmp"

    genome_paths = _read_manifest_marker_cds(busco_manifest, cds_template)
    if len(genome_paths) == 0:
        raise ValueError("No genomes found from BUSCO manifest for dereplication")

    logger.info("Preparing to dereplicate %d genomes", len(genome_paths))
    _write_genome_list(genome_paths, genome_list_file)

    try:
        _run_galah_cluster(
            genome_list_file=genome_list_file,
            quality_report=quality_report,
            clusters_out=clusters_out,
            precluster_ani=precluster_ani,
            ani=ani,
            threads=threads,
            min_aligned_fraction=min_aligned_fraction,
            precluster_method=precluster_method,
            cluster_method=cluster_method,
        )
    finally:
        if genome_list_file.exists():
            genome_list_file.unlink()

    representatives = _extract_representatives(clusters_out)
    _write_representatives(representatives, representatives_out)

    _filter_busco_manifest(
        input_manifest=busco_manifest,
        output_manifest=filtered_manifest_out,
        representatives=set(representatives),
    )

    logger.info(
        "Dereplication complete: %d representatives selected (clusters=%s)",
        len(representatives),
        clusters_out,
    )


if __name__ == "__main__":
    main()
