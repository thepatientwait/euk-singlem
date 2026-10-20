#!/usr/bin/env python3
"""Download one genome from NCBI, P10K, or JGI into the legacy output path."""
import gzip
import json
import shutil
import subprocess
import tempfile
import zipfile
import os
from pathlib import Path

import requests


from log_config import configure_logger


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)



def validate_gzip_file(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except (EOFError, OSError, gzip.BadGzipFile):
        return False


def download_p10k(url: str, output: Path) -> bool:
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()
    with output.open("wb") as handle:
        for chunk in response.iter_content(1024 * 1024):
            if chunk:
                handle.write(chunk)
    return True


def download_ncbi(accession: str, output: Path) -> bool:
    with tempfile.TemporaryDirectory(prefix="ncbi_download_") as temporary_dir:
        command = [
            "datasets",
            "download",
            "genome",
            "accession",
            accession,
            "--include",
            "genome",
            "--filename",
            str(Path(temporary_dir) / "genome.zip"),
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(result.stderr.strip() or "NCBI Datasets download failed")

        archive = Path(temporary_dir) / "genome.zip"
        extract_dir = Path(temporary_dir) / "extracted"
        with zipfile.ZipFile(archive) as zip_handle:
            zip_handle.extractall(extract_dir)

        fasta_files = sorted(extract_dir.rglob("*.fna"))
        if not fasta_files:
            report_file = extract_dir / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
            if report_file.exists():
                with report_file.open("r") as f:
                    for line in f:
                        if not line.strip():
                            continue
                        try:
                            data = json.loads(line)
                            status = data.get("assemblyInfo", {}).get("assemblyStatus")
                            if status == "suppressed":
                                return False
                        except json.JSONDecodeError:
                            continue

            raise FileNotFoundError(f"No FASTA file found in NCBI archive for {accession}")

        with fasta_files[0].open("rb") as source, gzip.open(output, "wb") as destination:
            shutil.copyfileobj(source, destination)

        return True


def download_jgi(accession: str, name: str, output: Path) -> bool:
    output_dir = output.parent
    downloaded_file = output_dir / name
    command = [
        "python",
        str(Path(__file__).with_name("jgi-query.py")),
        accession,
        "-r",
        f"^{name}$",
    ]
    result = subprocess.run(command, cwd=output_dir, capture_output=True, text=True)
    known_missing = "AttributeError: 'NoneType' object has no attribute 'values'" in result.stderr
    if result.returncode != 0 and not known_missing:
        raise RuntimeError(result.stderr.strip() or "JGI download failed")
    if not downloaded_file.exists() or downloaded_file.stat().st_size == 0:
        return False
    downloaded_file.replace(output)
    return True


def main() -> None:
    accession = str(snakemake.wildcards.genome_acc)
    source = str(snakemake.params.source).upper()
    output = Path(snakemake.output.genome)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output.with_suffix(output.suffix + ".partial")
    temporary_output.unlink(missing_ok=True)

    try:
        if source == "JGI":
            downloaded = download_jgi(accession, str(snakemake.params.name), temporary_output)
        elif source == "P10K":
            downloaded = download_p10k(str(snakemake.params.url), temporary_output)
        elif source == "NCBI":
            downloaded = download_ncbi(accession, temporary_output)
        else:
            raise ValueError(f"Unsupported download source: {source}")

        if not downloaded:
            output.touch()
            return

        if not validate_gzip_file(temporary_output):
            raise RuntimeError(f"Downloaded file is not a valid gzip stream: {accession}")
        temporary_output.replace(output)
    except Exception:
        temporary_output.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    main()