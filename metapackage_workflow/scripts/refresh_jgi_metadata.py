#!/usr/bin/env python3
"""Refresh a JGI group file list (``{group}.json``) so tape-restore statuses stay current.

File statuses (PURGED -> RESTORE_REGISTERED -> RESTORED) change on JGI's side while a restore
runs. The workflow calls :func:`refresh_group` once per run, before the DAG is built.

The JSON is only replaced when the organisms actually changed, so an unchanged catalogue keeps
its mtime and triggers no reruns. Any failure keeps the existing file: a flaky API must
never break (or truncate) the metadata a run depends on.
"""
import collections
import json
import os
from pathlib import Path

try:
    from jgi_common import fetch_all_pages
except ImportError:  # imported as scripts.refresh_jgi_metadata
    from scripts.jgi_common import fetch_all_pages

# Refuse a refreshed list that lost more than this share of the organisms (truncated answer).
MIN_SIZE_RATIO = 0.9


def _status_counts(organisms):
    counts = collections.Counter()
    for org in organisms:
        files = org.get("files") or []
        if not files:
            counts["NOFILES"] += 1
            continue
        biggest = max(files, key=lambda f: int(f.get("file_size") or 0))
        counts[biggest.get("file_status", "UNKNOWN")] += 1
    return counts


def _key(organisms):
    return json.dumps(organisms, sort_keys=True, separators=(",", ":"))


def refresh_group(group, json_url, json_path, logger):
    """Re-download ``json_url`` into ``json_path`` if it changed. Returns True if replaced."""
    json_path = Path(json_path)
    if not json_path.exists():
        return False  # the download_jgi_metadata checkpoint creates it
    try:
        new = fetch_all_pages(json_url, logger, label=f"for {group}")
        with open(json_path, "r", encoding="utf-8") as handle:
            old = json.load(handle)
    except Exception as exc:  # network, JSON or HTTP errors: keep what we have
        logger.warning("JGI %s: metadata refresh failed (%s: %s); keeping existing %s",
                       group, type(exc).__name__, exc, json_path)
        return False

    old_orgs, new_orgs = old.get("organisms", []), new.get("organisms", [])
    if old_orgs and len(new_orgs) < MIN_SIZE_RATIO * len(old_orgs):
        logger.warning("JGI %s: refreshed list has %d organisms vs %d; refusing it, keeping %s",
                       group, len(new_orgs), len(old_orgs), json_path)
        return False
    if _key(old_orgs) == _key(new_orgs):
        logger.info("JGI %s: metadata unchanged (%d organisms)", group, len(new_orgs))
        return False

    tmp = json_path.with_name(json_path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(new, handle)
    os.replace(tmp, json_path)
    logger.info("JGI %s: metadata updated, %d -> %d organisms; best-file status %s -> %s",
                group, len(old_orgs), len(new_orgs),
                dict(_status_counts(old_orgs)), dict(_status_counts(new_orgs)))
    return True
