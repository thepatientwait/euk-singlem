#!/usr/bin/env python3

"""
Extract nucleotide sequences of a genome's complete BUSCO genes (mode ``exon_blocks`` or ``gene_span``).

These sequences are searched with ``singlem pipe`` to build the nucleotide database of the
metapackage, which is later queried with metagenomic *DNA* reads. A read cannot cross an intron,
so a window that spans one is never observed in reads. Two modes (wildcard ``introns``):

* ``exon_blocks`` (corrected): the genes are cut at every intron (using the exon structure of the
  miniprot gene model BUSCO writes next to each hit) and each intron-free block is written as its
  own record; no ORF, and so no window, can bridge two exons;
* ``gene_span`` (naive): the whole genomic span of each BUSCO hit, introns included, as one record.

Exons separated by less than ``min_intron`` nt are frameshift/indel artefacts rather than
introns and stay together. Blocks shorter than ``min_exon_nt`` cannot hold a window and are
dropped.
"""

import gzip
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append(str(Path(__file__).parent))
from log_config import configure_logger
import intron_utils as iu

logger = configure_logger(__file__)

DEFAULT_MIN_EXON_NT = 60


def parse_busco_file(busco_file: Path) -> Tuple[str, str, pl.DataFrame]:
    """Group, genome accession and the Complete rows of a BUSCO ``full_table.tsv``.

    Path format: outputs/busco/{group}/{genome_acc}/run_eukaryota_odb12/full_table.tsv
    """
    group = busco_file.parts[-4]
    genome_acc = busco_file.parts[-3]
    logger.info(f"Parsing BUSCO file for {group}/{genome_acc}")
    busco_df = pl.read_csv(busco_file, separator='\t', skip_rows=2).rename({'# Busco id': 'Busco id'})
    complete_genes = busco_df.filter(pl.col('Status') == 'Complete')
    logger.info(f"Found {len(complete_genes)} complete genes for {group}/{genome_acc}")
    return group, genome_acc, complete_genes


def load_genome_sequences(genome_file) -> Dict[str, SeqRecord]:
    """Genome FASTA (optionally gzipped) as ``{contig: record}``."""
    genome_file = str(genome_file)
    opener = gzip.open if genome_file.endswith('.gz') else open
    with opener(genome_file, 'rt') as handle:
        sequences = {record.id: record for record in SeqIO.parse(handle, 'fasta')}
    logger.info(f"Loaded {len(sequences)} sequences from {genome_file}")
    return sequences


def extract_segment(contig: SeqRecord, start: int, end: int, strand: str) -> Optional[Seq]:
    """Bases ``start..end`` (1-based, inclusive) of a contig in the gene's orientation."""
    if start < 1 or end > len(contig.seq):
        return None
    seq = contig.seq[start - 1:end]
    return seq.reverse_complement() if strand == '-' else seq


def gene_records(busco_id: str, model: iu.Model, genome: Dict[str, SeqRecord], min_intron: int, min_exon_nt: int):
    """Intron-free records of one gene: ``(records, n_introns, n_short_dropped)``."""
    contig = genome.get(model.contig)
    if contig is None:
        logger.warning(f"Contig {model.contig} not found in genome for {busco_id}")
        return [], 0, 0
    segments = iu.exon_segments(model, min_intron)
    records, short = [], 0
    for index, (start, end) in enumerate(segments, start=1):
        seq = extract_segment(contig, start, end, model.strand)
        if seq is None:
            logger.warning(f"Coordinates out of bounds for {busco_id} on {model.contig}: {start}-{end}")
            continue
        if len(seq) < min_exon_nt:
            short += 1
            continue
        records.append(SeqRecord(
            seq,
            id=f"{busco_id}|exon{index}of{len(segments)}",
            description=f"{model.contig}:{start}-{end} strand {model.strand}",
        ))
    return records, len(segments) - 1, short


def gene_span_record(busco_id: str, row: dict, genome: Dict[str, SeqRecord]) -> Optional[SeqRecord]:
    """The whole genomic span (introns included) of one BUSCO hit, in the gene's orientation."""
    contig = genome.get(str(row['Sequence']))
    if contig is None:
        return None
    start, end = int(row['Gene Start']), int(row['Gene End'])
    lo, hi = min(start, end), max(start, end)
    seq = extract_segment(contig, lo, hi, str(row['Strand']))
    if seq is None:
        return None
    return SeqRecord(seq, id=busco_id, description=f"{row['Sequence']}:{lo}-{hi} strand {row['Strand']} (gene span)")


def extract_busco_transcripts(genome_file, full_table_file, output_file,
                              min_intron: int = iu.DEFAULT_MIN_INTRON,
                              min_exon_nt: int = DEFAULT_MIN_EXON_NT,
                              mode: str = "exon_blocks") -> int:
    """Write the sequences of every Complete BUSCO gene (see module doc for ``mode``); returns the record count."""
    if mode not in ("exon_blocks", "gene_span"):
        raise ValueError(f"unknown mode {mode!r}")
    full_table = Path(full_table_file)
    group, genome_acc, complete_genes = parse_busco_file(full_table)
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if len(complete_genes) == 0:
        logger.warning("No complete BUSCO genes found in full table")
        output_path.touch()
        return 0

    genome = load_genome_sequences(genome_file)
    gff_dir = full_table.parent / "busco_sequences" / "single_copy_busco_sequences"
    records, no_model, with_introns, short_total = [], 0, 0, 0

    for row in complete_genes.iter_rows(named=True):
        busco_id = str(row['Busco id'])
        if mode == "gene_span":
            record = gene_span_record(busco_id, row, genome)
            if record is None:
                no_model += 1
            else:
                records.append(record)
            continue
        gff = gff_dir / f"{busco_id}.gff"
        models = iu.parse_gff(gff) if gff.exists() else []
        model = iu.select_model(models, str(row['Sequence']), int(row['Gene Start']), int(row['Gene End']),
                                str(row['Strand'])) or (models[0] if len(models) == 1 else None)
        if model is None:
            no_model += 1
            continue
        gene_recs, n_introns, short = gene_records(busco_id, model, genome, min_intron, min_exon_nt)
        records.extend(gene_recs)
        with_introns += n_introns > 0
        short_total += short

    if no_model:
        logger.warning(f"{no_model} of {len(complete_genes)} complete genes had no usable gene model "
                       f"(GFF missing or not matching the table) and were skipped")
    logger.info(f"{group}/{genome_acc}: {len(complete_genes)} complete genes, {with_introns} with introns, "
                f"{len(records)} intron-free records written, {short_total} blocks shorter than {min_exon_nt} nt dropped")

    with open(output_path, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
    return len(records)


def main():
    """Entry point when run by Snakemake."""
    try:
        count = extract_busco_transcripts(
            snakemake.input.genome, snakemake.input.full_table, snakemake.output.transcripts,
            min_intron=int(snakemake.params.get('min_intron', iu.DEFAULT_MIN_INTRON)),
            min_exon_nt=int(snakemake.params.get('min_exon_nt', DEFAULT_MIN_EXON_NT)),
            mode=str(snakemake.params.get('mode', 'exon_blocks')),
        )
        if count == 0:
            logger.warning("No sequences extracted")
            Path(snakemake.output.transcripts).touch()
    except Exception as e:
        logger.error(f"Error extracting transcripts: {e}")
        raise


if __name__ == "__main__":
    main()
