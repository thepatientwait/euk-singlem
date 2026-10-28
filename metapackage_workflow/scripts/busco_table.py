"""Read BUSCO ``full_table.tsv`` files of both modes with one set of columns.

Genome mode (miniprot) columns: ``Busco id, Status, Sequence, Gene Start, Gene End, Strand, Score, Length, ...``
with ``Sequence`` = contig.  Transcriptome mode (MetaEuk) has no coordinate columns: ``Sequence`` is
``contig:start-end`` and the strand is only in the headers of ``busco_sequences/*/<gene>.{faa,fna}``
(``>gene_..|contig:start-end|+``).

``read_full_table`` returns the genome-mode columns for either mode, plus ``Transcriptome`` (bool).  For
transcriptome tables ``Sequence`` becomes the bare contig and ``Strand`` is the empty string (unknown);
record matching (``match_record``) then ignores the strand.
"""
import re
from pathlib import Path
from typing import Iterable, Optional

import polars as pl

_SEQ_COORDS = re.compile(r"^(?P<contig>.+):(?P<start>\d+)-(?P<end>\d+)$")
_HEADER_COORDS = re.compile(r"\|([^|:]+):(\d+)-(\d+)\|([+-])(?:\s|$)")


def normalise_full_table(df: pl.DataFrame) -> pl.DataFrame:
    """Add ``Gene Start``/``Gene End``/``Strand`` to a transcriptome-mode table (no-op for genome mode)."""
    if "Gene Start" in df.columns:
        return df.with_columns(pl.lit(False).alias("Transcriptome"))
    parts = df["Sequence"].cast(pl.String).str.extract_groups(_SEQ_COORDS.pattern)
    return df.with_columns(
        pl.when(parts.struct.field("contig").is_not_null()).then(parts.struct.field("contig"))
        .otherwise(pl.col("Sequence")).alias("Sequence"),
        parts.struct.field("start").cast(pl.Int64, strict=False).alias("Gene Start"),
        parts.struct.field("end").cast(pl.Int64, strict=False).alias("Gene End"),
        pl.lit("").alias("Strand"),
        pl.lit(True).alias("Transcriptome"),
    )


def read_full_table(path) -> pl.DataFrame:
    """The BUSCO full table with normalised columns (see module doc); ``Busco id`` is renamed from ``# Busco id``."""
    df = pl.read_csv(path, separator="\t", skip_rows=2,
                     schema_overrides={"Sequence": pl.String}).rename({"# Busco id": "Busco id"})
    return normalise_full_table(df)


def header_coords(description: str):
    """(contig, lo, hi, strand) from a BUSCO sequence header, or None."""
    match = _HEADER_COORDS.search(description)
    if match is None:
        return None
    a, b = int(match.group(2)), int(match.group(3))
    return match.group(1), min(a, b), max(a, b), match.group(4)


def match_record(records: Iterable, contig: str, start, end, strand: Optional[str], used_ids: set):
    """The FASTA record whose header coordinates equal the table row (strand ignored when unknown/empty)."""
    try:
        lo, hi = min(int(start), int(end)), max(int(start), int(end))
    except (TypeError, ValueError):
        return None
    strand = (strand or "").strip()
    if strand in ("None", "nan"):
        strand = ""
    for record in records:
        if record.id in used_ids:
            continue
        key = header_coords(record.description)
        if key is None:
            continue
        if key[:3] == (str(contig), lo, hi) and (not strand or key[3] == strand):
            return record
    return None
