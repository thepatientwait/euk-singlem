#!/usr/bin/env python3
"""Download one genome from NCBI, P10K, or JGI into the legacy output path."""
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import time
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


NCBI_FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
FTP_ATTEMPTS = 3
GZIP_LEVEL = 6  # python's default of 9 is several times slower for little gain


def ftp_parent_url(accession: str) -> str | None:
    """Directory listing URL holding this accession's assembly folder, e.g. .../GCA/000/149/185/."""
    match = re.fullmatch(r"(GC[AF])_(\d{9})\.\d+", accession)
    if not match:
        return None
    prefix, digits = match.group(1), match.group(2)
    return f"{NCBI_FTP_ROOT}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"


def find_ftp_assembly_dir(accession: str, session: requests.Session) -> str | None:
    """URL of the assembly folder (``<acc>_<asm name>/``) on the NCBI FTP, or None if absent."""
    parent = ftp_parent_url(accession)
    if parent is None:
        return None
    response = session.get(parent, timeout=60)
    if response.status_code == 404:
        return None
    response.raise_for_status()
    match = re.search(rf'href="({re.escape(accession)}_[^"/]+)/"', response.text)
    return f"{parent}{match.group(1)}/" if match else None


def expected_md5(session: requests.Session, directory_url: str, file_name: str) -> str | None:
    response = session.get(f"{directory_url}md5checksums.txt", timeout=60)
    if response.status_code != 200:
        return None
    for line in response.text.splitlines():
        parts = line.split()
        if len(parts) == 2 and Path(parts[1]).name == file_name:
            return parts[0]
    return None


def download_ncbi_ftp(accession: str, output: Path) -> bool:
    """Fetch the already-gzipped ``*_genomic.fna.gz`` from the NCBI FTP and check its md5.

    Returns False if the assembly has no FTP folder / file (suppressed, replaced, too new), so
    the caller can fall back to NCBI Datasets. Transient errors are retried, then raised.
    """
    session = requests.Session()
    for attempt in range(1, FTP_ATTEMPTS + 1):
        try:
            directory = find_ftp_assembly_dir(accession, session)
            if directory is None:
                logger.info("%s: no NCBI FTP folder", accession)
                return False
            file_name = f"{directory.rstrip('/').rsplit('/', 1)[1]}_genomic.fna.gz"
            md5 = expected_md5(session, directory, file_name)
            with session.get(directory + file_name, stream=True, timeout=(30, 120)) as response:
                if response.status_code == 404:
                    logger.info("%s: %s not on the FTP", accession, file_name)
                    return False
                response.raise_for_status()
                digest = hashlib.md5()
                with output.open("wb") as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        handle.write(chunk)
                        digest.update(chunk)
            if md5 is not None and digest.hexdigest() != md5:
                raise IOError(f"md5 mismatch for {file_name}: {digest.hexdigest()} != {md5}")
            logger.info("%s: downloaded %s from the NCBI FTP (%d bytes)", accession, file_name, output.stat().st_size)
            return True
        except (requests.RequestException, IOError) as exc:
            output.unlink(missing_ok=True)
            if attempt == FTP_ATTEMPTS:
                raise
            delay = 15 * attempt
            logger.warning("%s: FTP download failed (attempt %d/%d): %s. Retrying in %ds",
                           accession, attempt, FTP_ATTEMPTS, exc, delay)
            time.sleep(delay)
    return False


def gzip_stream(source, output: Path, threads: int) -> None:
    """Compress a binary stream into ``output``; multi-threaded with pigz when it is installed."""
    pigz = shutil.which("pigz")
    if pigz:
        with output.open("wb") as destination:
            process = subprocess.Popen([pigz, "-c", "-p", str(max(1, threads)), f"-{GZIP_LEVEL}"],
                                       stdin=subprocess.PIPE, stdout=destination)
            try:
                shutil.copyfileobj(source, process.stdin, 1024 * 1024)
            finally:
                process.stdin.close()
            if process.wait() != 0:
                raise RuntimeError("pigz failed")
    else:
        with gzip.open(output, "wb", compresslevel=GZIP_LEVEL) as destination:
            shutil.copyfileobj(source, destination, 1024 * 1024)


def download_ncbi_datasets(accession: str, output: Path, threads: int = 1) -> bool:
    """Fallback: NCBI Datasets zip (uncompressed FASTA), streamed straight into the compressor."""
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
        with zipfile.ZipFile(archive) as zip_handle:
            fasta_members = sorted(n for n in zip_handle.namelist() if n.endswith(".fna"))
            if not fasta_members:
                report = "ncbi_dataset/data/assembly_data_report.jsonl"
                if report in zip_handle.namelist():
                    with zip_handle.open(report) as handle:
                        for line in handle:
                            if not line.strip():
                                continue
                            try:
                                data = json.loads(line)
                            except json.JSONDecodeError:
                                continue
                            if data.get("assemblyInfo", {}).get("assemblyStatus") == "suppressed":
                                return False

                raise FileNotFoundError(f"No FASTA file found in NCBI archive for {accession}")

            with zip_handle.open(fasta_members[0]) as source:
                gzip_stream(source, output, threads)

        return True


def download_ncbi(accession: str, output: Path, threads: int = 1) -> bool:
    """FTP first (no recompression, md5-checked); NCBI Datasets when the FTP has no file."""
    try:
        if download_ncbi_ftp(accession, output):
            return True
    except Exception as exc:  # FTP unreachable / corrupt after retries: Datasets is the safety net
        logger.warning("%s: FTP download failed (%s: %s)", accession, type(exc).__name__, exc)
    logger.info("%s: falling back to NCBI Datasets", accession)
    output.unlink(missing_ok=True)
    return download_ncbi_datasets(accession, output, threads)


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

    # Reuse existing downloads: discovery regeneration must not trigger
    # re-downloads of genomes already on disk. Zero-byte files are failure
    # markers from cleanup_failed_genomes and remain eligible for retry.
    if output.exists() and output.stat().st_size > 0:
        logger.info("Genome already present, skipping download: %s", output)
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output.with_suffix(output.suffix + ".partial")
    temporary_output.unlink(missing_ok=True)

    try:
        if source == "JGI":
            downloaded = download_jgi(accession, str(snakemake.params.name), temporary_output)
        elif source == "P10K":
            downloaded = download_p10k(str(snakemake.params.url), temporary_output)
        elif source == "NCBI":
            downloaded = download_ncbi(accession, temporary_output, int(snakemake.threads))
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