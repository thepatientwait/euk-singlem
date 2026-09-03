#!/usr/bin/env python3

import csv
import gzip
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any

from log_config import configure_logger


logger = configure_logger(__file__)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake: Any = None


def read_hmm_gene_map(ids_path: Path) -> dict[str, str]:
    hmm_to_gene = {}
    with open(ids_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["gene"].strip()
            hmm_name = row["hmm_name"].strip()
            if gene and hmm_name:
                hmm_to_gene[hmm_name] = gene
    return hmm_to_gene


def parse_fasta(path: Path) -> dict[str, str]:
    seqs = {}
    current_id = None
    chunks = []
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        seqs[current_id] = "".join(chunks)
    return seqs


def materialize_genome_file(source_path: Path, destination_dir: Path) -> Path:
    destination_dir.mkdir(parents=True, exist_ok=True)
    if source_path.suffix == ".gz":
        output_path = destination_dir / source_path.name[:-3]
        with gzip.open(source_path, "rt", encoding="utf-8", errors="ignore") as src, open(
            output_path,
            "w",
            encoding="utf-8",
        ) as dst:
            shutil.copyfileobj(src, dst)
        return output_path

    output_path = destination_dir / source_path.name
    shutil.copy2(source_path, output_path)
    return output_path


def parse_tblout(tblout_path: Path, hmm_to_gene: dict[str, str]) -> dict[str, set[str]]:
    gene_hits: dict[str, set[str]] = defaultdict(set)
    with open(tblout_path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            parts = raw.split()
            if len(parts) < 3:
                continue
            sequence_id = parts[0]
            hmm_name = parts[2]
            gene = hmm_to_gene.get(hmm_name)
            if gene is not None:
                gene_hits[gene].add(sequence_id)
    return gene_hits


def write_result(path: Path, rows: list[tuple[str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", "sequence_id", "sequence"])
        for gene, sequence_id, sequence in rows:
            writer.writerow([gene, sequence_id, sequence])


def main() -> None:
    hmm_file = Path(str(snakemake.input.hmm_file))
    hmm_ids = Path(str(snakemake.input.hmm_ids))
    genome_path = Path(str(snakemake.input.genome))
    result_path = Path(str(snakemake.output.result))

    hmm_evalue = str(snakemake.params.hmm_evalue)
    off_target_id = str(snakemake.wildcards.off_target_id)

    hmm_to_gene = read_hmm_gene_map(hmm_ids)

    if not genome_path.exists():
        logger.warning("Skipping missing off-target genome file: %s", genome_path)
        write_result(result_path, [])
        return

    logger.info("Running off-target hmmsearch for %s", off_target_id)

    with tempfile.TemporaryDirectory(prefix="off_target_hmmsearch_") as tmpdir:
        tmp = Path(tmpdir)
        search_path = materialize_genome_file(genome_path, tmp / "decompressed")

        tblout_path = tmp / f"{off_target_id}.tblout"
        cmd = [
            "hmmsearch",
            "-E",
            hmm_evalue,
            "--tblout",
            str(tblout_path),
            str(hmm_file),
            str(search_path),
        ]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        gene_hits = parse_tblout(tblout_path, hmm_to_gene)
        if len(gene_hits) == 0:
            write_result(result_path, [])
            return

        seqs = parse_fasta(search_path)
        rows = []
        for gene in sorted(gene_hits.keys()):
            for sequence_id in sorted(gene_hits[gene]):
                seq = seqs.get(sequence_id)
                if seq:
                    prefixed_id = f"{off_target_id}|{sequence_id}"
                    rows.append((gene, prefixed_id, seq))

        write_result(result_path, rows)

    logger.info("Completed off-target search for %s", off_target_id)


if __name__ == "__main__":
    main()
