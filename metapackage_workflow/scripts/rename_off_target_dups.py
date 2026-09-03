#!/usr/bin/env python3

from collections import defaultdict
from pathlib import Path
from typing import Any

from log_config import configure_logger


logger = configure_logger(__file__)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake: Any = None


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    current_id = None
    chunks = []
    if not path.exists():
        return records

    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks)))
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        records.append((current_id, "".join(chunks)))

    return records


def read_taxonomy(path: Path) -> dict[str, str]:
    mapping = {}
    if not path.exists():
        return mapping

    with open(path, "r", encoding="utf-8") as handle:
        header = next(handle, "")
        if "sequence_id" not in header:
            parts = header.rstrip("\n").split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
        for raw in handle:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]

    return mapping


def main() -> None:
    aggregate_dir = Path(str(snakemake.input.aggregate_dir))
    output_seqs = Path(str(snakemake.output.seqs))
    output_taxonomy = Path(str(snakemake.output.taxonomy))
    busco_gene = str(snakemake.wildcards.busco_gene)

    input_seqs = aggregate_dir / f"{busco_gene}.faa"
    input_taxonomy = aggregate_dir / f"{busco_gene}_taxonomy.tsv"

    output_seqs.parent.mkdir(parents=True, exist_ok=True)
    output_taxonomy.parent.mkdir(parents=True, exist_ok=True)

    records = read_fasta(input_seqs)
    taxonomy_map = read_taxonomy(input_taxonomy)

    counts = defaultdict(int)
    renamed_records = []

    for seq_id, seq in records:
        counts[seq_id] += 1
        if counts[seq_id] == 1:
            new_id = seq_id
        else:
            new_id = f"{seq_id}--{counts[seq_id]}"
        renamed_records.append((new_id, seq, taxonomy_map.get(seq_id, "d__Bacteria")))

    with open(output_seqs, "w", encoding="utf-8") as seq_out:
        for seq_id, seq, _ in renamed_records:
            seq_out.write(f">{seq_id}\n{seq}\n")

    with open(output_taxonomy, "w", encoding="utf-8") as tax_out:
        tax_out.write("sequence_id\ttaxonomy\n")
        for seq_id, _, taxonomy in renamed_records:
            tax_out.write(f"{seq_id}\t{taxonomy}\n")

    logger.info("Renamed %d records (%d unique original IDs)", len(renamed_records), len(counts))


if __name__ == "__main__":
    main()
