#!/usr/bin/env python3
"""Reduce OrthoDB's big tab-separated tables to one clade, and turn the OG sequence table into a FASTA.

    orthodb_tables.py subset IN.gz OUT.gz --column N --level-taxid 2759
    orthodb_tables.py fasta  IN.gz OUT.faa --og-column 0 --gene-column 1 --seq-column 2 --level-taxid 2759

OrthoDB names an orthologous group ``<number>at<taxid of its level>`` (BUSCO uses the same ids), so the OGs of
a clade are recognised by the ``at<taxid>`` suffix in the OG column, without needing the taxonomy tables.
Everything streams: the input tables are tens of gigabytes.

The FASTA headers are ``<gene>|<og>`` so a DIAMOND hit names both the gene and the group it belongs to.
Column positions are configuration (config/config.yaml): they are taken from OrthoDB's documentation and
must be confirmed on the real files (the pilot does that).
"""
import argparse
import gzip
import io
import re
import sys
from contextlib import contextmanager


@contextmanager
def open_text(path, mode):
    if str(path).endswith(".gz"):
        raw = gzip.open(path, mode + "b")
        handle = io.TextIOWrapper(raw, encoding="utf-8", newline="\n")
    else:
        handle = open(path, mode, encoding="utf-8", newline="\n")
    try:
        yield handle
    finally:
        handle.close()


def og_matcher(level_taxid):
    return re.compile(rf"^\d+at{int(level_taxid)}$")


def subset(src, dst, column, level_taxid):
    """Copy the lines whose ``column`` is an OG of ``level_taxid``. Returns (kept, total)."""
    match = og_matcher(level_taxid)
    kept = total = 0
    with open_text(src, "r") as fin, open_text(dst, "w") as fout:
        for line in fin:
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) > column and match.search(fields[column]):
                fout.write(line)
                kept += 1
    if total and not kept:
        raise ValueError(f"no line of {src} has an OG of level {level_taxid} in column {column}: "
                         "is the column number right for this OrthoDB release?")
    return kept, total


def to_fasta(src, dst, og_column, gene_column, seq_column, level_taxid):
    """Write ``>gene|og`` records for the OGs of ``level_taxid``. Returns (records, total)."""
    match = og_matcher(level_taxid)
    need = max(og_column, gene_column, seq_column)
    records = total = 0
    with open_text(src, "r") as fin, open(dst, "w", encoding="utf-8") as fout:
        for line in fin:
            total += 1
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= need or not match.search(fields[og_column]):
                continue
            sequence = fields[seq_column].strip().replace("*", "").replace(" ", "")
            if not sequence:
                continue
            fout.write(f">{fields[gene_column]}|{fields[og_column]}\n{sequence}\n")
            records += 1
    if total and not records:
        raise ValueError(f"no usable record in {src}: check og/gene/sequence columns ({og_column}/{gene_column}/{seq_column})")
    return records, total


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    a = sub.add_parser("subset")
    a.add_argument("src"); a.add_argument("dst")
    a.add_argument("--column", type=int, required=True)
    a.add_argument("--level-taxid", type=int, default=2759)
    b = sub.add_parser("fasta")
    b.add_argument("src"); b.add_argument("dst")
    b.add_argument("--og-column", type=int, default=0)
    b.add_argument("--gene-column", type=int, default=1)
    b.add_argument("--seq-column", type=int, default=2)
    b.add_argument("--level-taxid", type=int, default=2759)
    args = parser.parse_args(argv)
    if args.command == "subset":
        kept, total = subset(args.src, args.dst, args.column, args.level_taxid)
    else:
        kept, total = to_fasta(args.src, args.dst, args.og_column, args.gene_column, args.seq_column, args.level_taxid)
    print(f"{args.command}: kept {kept:,} of {total:,} lines")


if __name__ == "__main__":
    try:
        main()
    except ValueError as err:
        sys.exit(str(err))
