#!/usr/bin/env python3
"""Write a JGI group file list (``{group}.json``) only when its content changed.

Used by the ``download_jgi_metadata`` rule. The rule is not re-run automatically (its outputs
exist and it has no inputs); redo it on purpose with ``--forcerun download_jgi_metadata`` or by
deleting the output. Rewriting only on change keeps the mtime of an unchanged catalogue, so a
forced refresh of an identical list triggers no reruns downstream.
"""
import json
import os
from pathlib import Path

# Refuse a refreshed list that lost more than this share of the organisms (truncated answer).
MIN_SIZE_RATIO = 0.9


def _key(organisms):
    return json.dumps(organisms, sort_keys=True, separators=(",", ":"))


def write_json_if_changed(json_path, new, logger, label=""):
    """Replace ``json_path`` with ``new`` if the organisms differ (or it is missing).

    Returns True if the file was written. Raises if the new list looks truncated.
    """
    json_path = Path(json_path)
    new_orgs = new.get("organisms", [])
    old_orgs = []
    if json_path.exists():
        with open(json_path, "r", encoding="utf-8") as handle:
            old_orgs = json.load(handle).get("organisms", [])
        if old_orgs and len(new_orgs) < MIN_SIZE_RATIO * len(old_orgs):
            raise RuntimeError(f"JGI {label}: refreshed list has {len(new_orgs)} organisms vs {len(old_orgs)}; "
                               f"refusing it, keeping {json_path}")
        if _key(old_orgs) == _key(new_orgs):
            logger.info("JGI %s: metadata unchanged (%d organisms), %s left as is", label, len(new_orgs), json_path)
            return False

    json_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = json_path.with_name(json_path.name + ".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(new, handle)
    os.replace(tmp, json_path)
    logger.info("JGI %s: metadata written, %d -> %d organisms", label, len(old_orgs), len(new_orgs))
    return True


def write_bytes_if_changed(path, content):
    path = Path(path)
    if path.exists() and path.read_bytes() == content:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return True
