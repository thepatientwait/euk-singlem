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


def read_hmm_gene_map(ids_path: Path) -> tuple[dict[str, str], list[str]]:
    hmm_to_gene = {}
    genes = []
    with open(ids_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["gene"].strip()
            hmm_name = row["hmm_name"].strip()
            if gene and hmm_name:
                hmm_to_gene[hmm_name] = gene
                genes.append(gene)
    return hmm_to_gene, sorted(set(genes))


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
    """Copy or decompress a genome file into the temp directory for searching."""
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


def write_outputs(
    genes: list[str],
    sequences_by_gene: dict[str, list[tuple[str, str]]],
    aggregate_faa_pattern: str,
    aggregate_tax_pattern: str,
    default_taxonomy: str,
) -> None:
    for gene in genes:
        faa_path = Path(aggregate_faa_pattern.format(busco_gene=gene))
        tax_path = Path(aggregate_tax_pattern.format(busco_gene=gene))
        faa_path.parent.mkdir(parents=True, exist_ok=True)
        tax_path.parent.mkdir(parents=True, exist_ok=True)

        with open(faa_path, "w", encoding="utf-8") as faa_out:
            for seq_id, seq in sequences_by_gene.get(gene, []):
                faa_out.write(f">{seq_id}\n{seq}\n")

        with open(tax_path, "w", encoding="utf-8") as tax_out:
            tax_out.write("sequence_id\ttaxonomy\n")
            for seq_id, _ in sequences_by_gene.get(gene, []):
                tax_out.write(f"{seq_id}\t{default_taxonomy}\n")


def main() -> None:
    hmm_file = Path(str(snakemake.input.hmm_file))
    hmm_ids = Path(str(snakemake.input.hmm_ids))
    off_target_genomes = [Path(str(path)) for path in snakemake.input.off_target_genomes]

    aggregate_faa_pattern = str(snakemake.params.aggregate_faa_pattern)
    aggregate_tax_pattern = str(snakemake.params.aggregate_tax_pattern)
    hmm_evalue = str(snakemake.params.hmm_evalue)
    default_taxonomy = str(snakemake.params.default_taxonomy)

    done_path = Path(str(snakemake.output.done))
    done_path.parent.mkdir(parents=True, exist_ok=True)

    hmm_to_gene, genes = read_hmm_gene_map(hmm_ids)
    logger.info("Preparing off-target inputs for %d genes from %d genome files", len(genes), len(off_target_genomes))

    sequences_by_gene: dict[str, list[tuple[str, str]]] = defaultdict(list)

    with tempfile.TemporaryDirectory(prefix="off_target_hmmsearch_") as tmpdir:
        tmp = Path(tmpdir)
        decompressed_dir = tmp / "decompressed"

        for genome_path in off_target_genomes:
            if not genome_path.exists():
                logger.warning("Skipping missing off-target genome file: %s", genome_path)
                continue

            search_path = materialize_genome_file(genome_path, decompressed_dir)

            tblout_path = tmp / f"{search_path.stem}.tblout"
            cmd = [
                "hmmsearch",
                "-E",
                hmm_evalue,
                "--tblout",
                str(tblout_path),
                str(hmm_file),
                str(search_path),
            ]

            logger.info("Running hmmsearch for %s (search path: %s)", genome_path, search_path)
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            gene_hits = parse_tblout(tblout_path, hmm_to_gene)
            if len(gene_hits) == 0:
                continue

            seqs = parse_fasta(search_path)
            for gene, sequence_ids in gene_hits.items():
                for sequence_id in sorted(sequence_ids):
                    seq = seqs.get(sequence_id)
                    if seq:
                        prefixed_id = f"{search_path.stem}|{sequence_id}"
                        sequences_by_gene[gene].append((prefixed_id, seq))

    write_outputs(
        genes=genes,
        sequences_by_gene=sequences_by_gene,
        aggregate_faa_pattern=aggregate_faa_pattern,
        aggregate_tax_pattern=aggregate_tax_pattern,
        default_taxonomy=default_taxonomy,
    )

    done_path.touch()
    logger.info("Completed off-target input preparation")


if __name__ == "__main__":
    main()
