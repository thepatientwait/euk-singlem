#!/usr/bin/env python3
"""Spliced coding sequences of one genome's BUSCO marker genes, for dereplication.

Whole-genome ANI is a poor way to compare a transcriptome with a genome (or two assemblies of very
different completeness), so genomes are dereplicated on the genes they share: the marker CDS.

Output (``<acc>.fasta.gz``): one record per marker gene (``>busco_id``), the spliced CDS in gene
orientation, for every Complete gene and the top-scoring copy of every Duplicated gene.

* genomes: the exons of the miniprot gene model BUSCO wrote (``busco_sequences/*/<gene>.gff``) are cut out
  of the genome and joined (introns removed);
* transcriptomes (BUSCO transcriptome mode, group ``tsa``): the MetaEuk CDS (``<gene>.fna``).

Genes without a usable model or sequence are skipped (counted in the log).
"""
import gzip
import sys
from pathlib import Path
from typing import Dict, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger
import intron_utils as iu
from busco_table import match_record, read_full_table
from extract_busco_transcripts import load_genome_sequences

logger = configure_logger(__file__)

MIN_CDS_NT = 90


def top_rows(full_table: Path):
    """One row per gene: Complete, or the best-scoring copy of a Duplicated gene (ties: longer, then position)."""
    import polars as pl
    df = read_full_table(full_table).filter(pl.col("Status").is_in(["Complete", "Duplicated"]))
    if df.height == 0:
        return df
    df = df.with_columns(pl.col("Score").cast(pl.Float64, strict=False))
    sort_cols = ["Busco id", "Score"] + (["Length"] if "Length" in df.columns else []) + ["Sequence", "Gene Start"]
    desc = [False, True] + ([True] if "Length" in df.columns else []) + [False, False]
    return df.sort(sort_cols, descending=desc, nulls_last=True).unique(subset=["Busco id"], keep="first", maintain_order=True)


def spliced_cds(model: iu.Model, contig: SeqRecord) -> Optional[Seq]:
    """Exons joined in transcript order, each read in the gene's orientation."""
    parts = []
    for start, end in model.transcript_order():
        if start < 1 or end > len(contig.seq):
            return None
        piece = contig.seq[start - 1:end]
        parts.append(piece.reverse_complement() if model.strand == "-" else piece)
    return sum(parts[1:], parts[0]) if parts else None


def genome_cds(row: dict, run_dir: Path, genome: Dict[str, SeqRecord]) -> Optional[Seq]:
    subdir = "single_copy_busco_sequences" if row["Status"] == "Complete" else "multi_copy_busco_sequences"
    gff = run_dir / "busco_sequences" / subdir / f"{row['Busco id']}.gff"
    models = iu.parse_gff(gff) if gff.exists() else []
    model = iu.select_model(models, str(row["Sequence"]), int(row["Gene Start"]), int(row["Gene End"]), str(row["Strand"]))
    if model is None and len(models) == 1:
        model = models[0]
    contig = genome.get(model.contig) if model is not None else None
    return spliced_cds(model, contig) if contig is not None else None


def transcript_cds(row: dict, run_dir: Path) -> Optional[Seq]:
    subdir = "single_copy_busco_sequences" if row["Status"] == "Complete" else "multi_copy_busco_sequences"
    fna = run_dir / "busco_sequences" / subdir / f"{row['Busco id']}.fna"
    if not fna.exists():
        return None
    with open(fna) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    record = match_record(records, str(row["Sequence"]), row["Gene Start"], row["Gene End"], row.get("Strand"), set())
    if record is None and len(records) == 1:
        record = records[0]
    return record.seq if record is not None else None


def extract_marker_cds(genome_file, full_table_file, output_file) -> int:
    full_table = Path(full_table_file)
    run_dir = full_table.parent
    rows = top_rows(full_table)
    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    records, skipped = [], 0
    if rows.height:
        transcriptome = bool(rows["Transcriptome"][0])
        genome = None if transcriptome else load_genome_sequences(genome_file)
        for row in rows.sort("Busco id").iter_rows(named=True):
            seq = transcript_cds(row, run_dir) if transcriptome else genome_cds(row, run_dir, genome)
            if seq is None or len(seq) < MIN_CDS_NT:
                skipped += 1
                continue
            records.append(SeqRecord(seq.upper(), id=str(row["Busco id"]), description=""))
    with gzip.open(out, "wt", compresslevel=6) as handle:
        SeqIO.write(records, handle, "fasta")
    logger.info("%s: %d marker CDS written, %d skipped", full_table.parts[-3], len(records), skipped)
    return len(records)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]
    except NameError:
        sys.exit("run through snakemake")
    extract_marker_cds(snakemake.input.genome, snakemake.input.full_table, snakemake.output.cds)  # type: ignore[name-defined]
