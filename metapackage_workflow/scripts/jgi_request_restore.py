#!/usr/bin/env python3
"""Ask JGI to restore every missing, PURGED (tape) genome of one group in a few batched requests.

Does not wait: restores finish on JGI's side. Refresh the metadata and run the download rule
again later; genomes whose status has become RESTORED are then downloaded.
"""
import csv
import os
import time
from pathlib import Path

from jgi_common import (
    JDP_RESTORE_URL,
    JgiAuthError,
    check_token,
    largest_file,
    make_session,
    request_with_retries,
    restore_entry,
)
from log_config import configure_logger

logger = configure_logger(os.path.splitext(os.path.basename(__file__))[0])


def main():
    import json

    group = snakemake.wildcards.group
    token = snakemake.params.token
    batch_size = int(snakemake.params.batch_size)
    limit = int(snakemake.params.limit or 0)
    check_token(token)

    with open(snakemake.input.json, "r", encoding="utf-8") as handle:
        organisms = json.load(handle)["organisms"]

    pending = []
    for org in organisms:
        if not org.get("files"):
            continue
        f = largest_file(org)
        acc = org.get("portal_detail_id")
        path = Path(f"outputs/genomes/{group}/{acc}.fasta.gz")
        have = path.exists() and path.stat().st_size > 0
        if not have and f.get("file_status") == "PURGED" and f.get("_id"):
            pending.append((acc, org, f))
    if limit:
        pending = pending[:limit]
    logger.info("%s: %d missing PURGED genomes to restore (batch size %d)", group, len(pending), batch_size)

    session = make_session(token)
    rows = []
    failed = []

    def send(batch):
        """POST one batch; on a client error split it in halves to isolate the bad entries."""
        ids = {}
        for _, org, f in batch:
            key, sub = restore_entry(org, f)
            if key in ids:
                ids[key]["file_ids"].extend(sub["file_ids"])
            else:
                ids[key] = sub
        payload = {"ids": ids, "send_mail": False, "api_version": "2"}
        resp = request_with_retries(session, "POST", JDP_RESTORE_URL, logger, json=payload, timeout=60)
        if resp.status_code < 400:
            try:
                request_id = resp.json().get("request_id", "")
            except ValueError:
                request_id = ""
            logger.info("%d files: HTTP %d %s", len(batch), resp.status_code, resp.text[:200].replace("\n", " "))
            stamp = time.strftime("%Y-%m-%d %H:%M:%S")
            rows.extend((acc, f["_id"], stamp, resp.status_code, request_id) for acc, _, f in batch)
        elif resp.status_code < 500 and len(batch) > 1:
            mid = len(batch) // 2
            send(batch[:mid])
            send(batch[mid:])
        elif resp.status_code < 500:
            acc, _, f = batch[0]
            logger.warning("JGI cannot restore %s (file %s, %s): HTTP %d", acc, f["_id"], f.get("file_name"), resp.status_code)
            failed.append(acc)
        else:
            raise RuntimeError(f"Restore request failed: HTTP {resp.status_code}: {resp.text[:500]}")

    for start in range(0, len(pending), batch_size):
        send(pending[start:start + batch_size])
        time.sleep(1)

    if failed:
        logger.warning("%d genomes could not be restored: %s", len(failed), ", ".join(failed))

    out = Path(snakemake.output.report)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome_acc", "file_id", "requested_at", "http_status", "request_id"])
        writer.writerows(rows)
    logger.info("Requested restore of %d genomes; report: %s", len(rows), out)


if __name__ == "__main__":
    try:
        main()
    except JgiAuthError as exc:
        logger.error("%s", exc)
        raise SystemExit(1)
