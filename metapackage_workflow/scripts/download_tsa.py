#!/usr/bin/env python3
"""Download one NCBI TSA project (all its bulk part files) into outputs/genomes/tsa/<acc>.fasta.gz.

The project's file names come from the metadata written by discover_tsa.py (column ``parts``).
NCBI publishes no checksums for these files, so each part is verified by its Content-Length and
by reading the whole gzip stream.  Only contigs of the requested assembly version
(``<PREFIX><version:02d>...``) are kept, in case a file holds several versions.
A permanently missing file (HTTP 404) leaves a zero-byte marker, like the other download rules.
"""
import csv
import gzip
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import requests

from log_config import configure_logger

logger = configure_logger(os.path.splitext(os.path.basename(__file__))[0])

FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genbank/tsa"
ATTEMPTS = 4


def project_row(metadata: Path, genome_acc: str) -> dict:
    with open(metadata, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["genome_acc"] == genome_acc:
                return row
    raise KeyError(f"{genome_acc} is not in {metadata}")


def fetch_part(url: str, dest: Path) -> None:
    """Stream ``url`` to ``dest``; raises on a short or corrupt file (retried by the caller)."""
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        expected = int(response.headers.get("Content-Length") or 0)
        with open(dest, "wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                handle.write(chunk)
    if expected and dest.stat().st_size != expected:
        raise IOError(f"short download {dest.stat().st_size} of {expected} bytes: {url}")
    with gzip.open(dest, "rb") as handle:
        while handle.read(1024 * 1024):
            pass


def filter_version(sources, out_path: Path, prefix: str, version: int) -> int:
    """Concatenate the gzip parts, keeping records whose id starts with ``<prefix><version:02d>``."""
    wanted = f">{prefix}{version:02d}"
    kept = 0
    keep = False
    with gzip.open(out_path, "wt", compresslevel=6) as out:
        for src in sources:
            with gzip.open(src, "rt") as handle:
                for line in handle:
                    if line.startswith(">"):
                        keep = line.startswith(wanted)
                        kept += keep
                    if keep:
                        out.write(line)
    return kept


def download(genome_acc: str, metadata: Path, output: Path) -> None:
    if output.exists() and output.stat().st_size > 0:
        logger.info("%s already downloaded", genome_acc)
        return
    row = project_row(metadata, genome_acc)
    prefix, version = row["prefix"], int(row["version"])
    letter = prefix[0]
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=f".{genome_acc}.") as tmp:
        parts = []
        for name in row["parts"].split(";"):
            dest = Path(tmp) / name
            url = f"{FTP_ROOT}/{letter}/{name}"
            for attempt in range(ATTEMPTS):
                try:
                    fetch_part(url, dest)
                    break
                except requests.HTTPError as err:
                    if err.response is not None and err.response.status_code == 404:
                        logger.error("%s: %s is gone (404); leaving a zero-byte marker", genome_acc, url)
                        output.touch()
                        return
                    logger.warning("%s failed (%s), attempt %d", url, err, attempt + 1)
                except (requests.RequestException, IOError, EOFError, OSError) as err:
                    logger.warning("%s failed (%s), attempt %d", url, err, attempt + 1)
                time.sleep(2 ** attempt + 1)
            else:
                raise RuntimeError(f"could not download {url}")
            parts.append(dest)
        partial = output.with_suffix(output.suffix + ".partial")
        kept = filter_version(parts, partial, prefix, version)
        if kept == 0:
            partial.unlink(missing_ok=True)
            raise RuntimeError(f"{genome_acc}: no contigs with id {prefix}{version:02d}* in {row['parts']}")
        shutil.move(str(partial), str(output))
    logger.info("%s: %d contigs from %d part file(s)", genome_acc, kept, len(parts))


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]
    except NameError:
        sys.exit("run through snakemake")
    download(snakemake.wildcards.genome_acc, Path(snakemake.input.metadata),  # type: ignore[name-defined]
             Path(snakemake.output.genome))
