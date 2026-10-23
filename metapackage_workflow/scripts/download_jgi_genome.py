#!/usr/bin/env python3
"""Download one JGI genome that is already on disk at JGI (RESTORED / BACKUP_COMPLETE).

Files that are still on tape (PURGED) or mid-restore (RESTORE_REGISTERED) are waited for when
``params.wait`` is set: the job re-queries the file status every ``params.poll_seconds`` and
downloads as soon as it is RESTORED. The wait has no limit of its own; the scheduler's walltime
ends the attempt and Snakemake retries it with a longer one. Without ``wait`` such a file fails at
once. The token is re-read from the secrets file before every authenticated call, and a
rejected/expired token fails immediately.
"""
import gzip
import os
import random
import shutil
import time
import zipfile
from pathlib import Path

from jgi_common import (
    API_TIMEOUT,
    DOWNLOADABLE_STATUSES,
    JDP_DOWNLOAD_URL,
    JDP_RESTORE_URL,
    JDP_SEARCH_URL,
    JgiAuthError,
    JgiNotRestoredError,
    check_token,
    make_session,
    request_with_retries,
    read_token,
    restore_entry,
    WAITING_STATUSES,
)
from log_config import configure_logger

logger = configure_logger(os.path.splitext(os.path.basename(__file__))[0])


def validate_gzip_file(file_path):
    """Ensure the gzip stream can be read to EOF."""
    try:
        with gzip.open(file_path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except (EOFError, OSError, gzip.BadGzipFile) as exc:
        logger.error("Gzip validation failed: %s: %s", type(exc).__name__, exc)
        return False


def find_file(organisms, genome_acc, target_filename):
    """Exact file name first, then any assembly file mentioning the accession."""
    for org in organisms:
        for f in org.get("files", []):
            if f.get("file_name") == target_filename and f.get("_id"):
                return org, f
    for org in organisms:
        for f in org.get("files", []):
            fname = f.get("file_name", "")
            if genome_acc in fname and fname.endswith((".fa.gz", ".fasta.gz")) and f.get("_id"):
                return org, f
    return None, None


def refused_restore(report, genome_acc):
    """True if jgi_request_restore recorded that JGI cannot restore this genome (nothing to wait for)."""
    if not report or not os.path.exists(report):
        return False
    with open(report, "r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == genome_acc:
                return len(fields) > 2 and fields[2] == "restore_refused"
    return False


def query_file(genome_acc, target_filename, anonymous, authed):
    """Current (org, file) from the JDP file list. Anonymous first (like the group list); if JGI
    insists on a token, fall back to the authenticated session."""
    query = dict(params={"organism": genome_acc, "api_version": 2, "a": "false", "h": "false",
                         "d": "asc", "p": 1, "x": 50, "t": "simple"}, timeout=API_TIMEOUT)
    try:
        response = request_with_retries(anonymous, "GET", JDP_SEARCH_URL, logger, **query)
    except JgiAuthError:
        response = request_with_retries(authed(), "GET", JDP_SEARCH_URL, logger, **query)
    if response.status_code != 200:
        raise RuntimeError(f"Mycocosm file list query failed: HTTP {response.status_code}: {response.text[:300]}")
    return find_file(response.json().get("organisms", []), genome_acc, target_filename)


def main():
    import requests

    genome_acc = snakemake.wildcards.genome_acc
    output_file = Path(snakemake.output.genome)
    target_filename = snakemake.params.name
    initial_token = snakemake.params.token
    secrets_path = snakemake.params.secrets
    wait = bool(snakemake.params.wait)
    poll_seconds = float(snakemake.params.poll_seconds)
    request_restore = bool(snakemake.params.request_restore)
    not_ready_report = snakemake.input.get("not_ready")

    if output_file.exists() and output_file.stat().st_size > 0:
        logger.info("Genome already present, skipping download: %s", output_file)
        return

    def fresh_token():
        token = read_token(secrets_path, initial_token)
        check_token(token, min_seconds=300)  # the call itself only needs a few minutes
        return token

    def authed():
        return make_session(fresh_token())

    if not wait:
        check_token(read_token(secrets_path, initial_token))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    downloaded_file = output_file.parent / target_filename
    temp_zip_file = output_file.parent / f"{genome_acc}_download.zip"
    anonymous = requests.Session()

    restore_sent = False
    waited = 0
    while True:
        org, file_obj = query_file(genome_acc, target_filename, anonymous, authed)
        if file_obj is None:
            logger.warning("No file matching %s / %s on JDP; writing an empty failure marker.", genome_acc, target_filename)
            output_file.touch()
            return

        status = file_obj.get("file_status", "UNKNOWN")
        logger.info("Target file %s status: %s", file_obj["_id"], status)
        if status in DOWNLOADABLE_STATUSES:
            break
        if not wait or status not in WAITING_STATUSES:
            # A download attempt would only answer HTTP 409; fail at once and do not retry
            raise JgiNotRestoredError(
                f"{genome_acc} is not downloadable (status {status}) and waiting is "
                f"{'not enabled' if not wait else 'pointless for this status'}."
            )
        if refused_restore(not_ready_report, genome_acc):
            raise JgiNotRestoredError(f"{genome_acc}: JGI refused the restore request; nothing to wait for.")
        if status == "PURGED" and request_restore and not restore_sent:
            org_key, sub = restore_entry(org, file_obj)
            payload = {"ids": {org_key: sub}, "send_mail": False, "api_version": "2"}
            resp = request_with_retries(authed(), "POST", JDP_RESTORE_URL, logger, json=payload, timeout=60)
            if resp.status_code >= 400:
                raise JgiNotRestoredError(f"{genome_acc}: restore request refused: HTTP {resp.status_code} {resp.text[:200]}")
            logger.info("Restore requested: HTTP %d %s", resp.status_code, resp.text[:200].replace("\n", " "))
            restore_sent = True
        delay = poll_seconds * random.uniform(0.9, 1.1)
        logger.info("Waiting for JGI to restore the file (%.0f min so far); next check in %.0f min",
                    waited / 60, delay / 60)
        time.sleep(delay)
        waited += delay

    # authenticated from here on, with a token read right now (it may have been replaced meanwhile)
    session = authed()

    org_key, sub = restore_entry(org, file_obj)
    payload = {"ids": {org_key: sub}, "send_mail": False, "api_version": "2"}
    logger.info("Requesting download of %s", target_filename)
    with request_with_retries(session, "POST", JDP_DOWNLOAD_URL, logger, json=payload,
                              stream=True, timeout=180) as r:
        if r.status_code != 200:
            raise RuntimeError(f"Download failed: HTTP {r.status_code}: {r.text[:300]}")
        written = 0
        with open(temp_zip_file, "wb") as handle:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    handle.write(chunk)
                    written += len(chunk)
    logger.info("Archive downloaded: %d bytes", written)

    try:
        with zipfile.ZipFile(temp_zip_file) as zf:
            names = zf.namelist()
            member = next((m for m in names if Path(m).name == target_filename), None)
            if member is None and len(names) == 1:
                member = names[0]
            if member is None:
                raise RuntimeError(f"{target_filename} not found in zip contents: {names}")
            with zf.open(member) as src, open(downloaded_file, "wb") as dst:
                shutil.copyfileobj(src, dst, 1024 * 1024)
    finally:
        temp_zip_file.unlink(missing_ok=True)

    if not downloaded_file.exists() or downloaded_file.stat().st_size == 0:
        logger.error("Extracted file is missing or 0 bytes; writing an empty failure marker.")
        downloaded_file.unlink(missing_ok=True)
        output_file.touch()
        return

    if not validate_gzip_file(downloaded_file):
        downloaded_file.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed for {genome_acc}: gzip verification failed.")

    downloaded_file.replace(output_file)
    logger.info("Downloaded, verified and placed genome for %s", genome_acc)


if __name__ == "__main__":
    try:
        main()
    except (JgiAuthError, JgiNotRestoredError) as exc:
        logger.error("%s", exc)
        raise SystemExit(1)
