"""Shared helpers for the JGI (JDP / MycoCosm) download scripts."""
import base64
import json
import random
import re
import time
from concurrent.futures import ThreadPoolExecutor

import requests

JDP_SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
JDP_RESTORE_URL = "https://files.jgi.doe.gov/request_archived_files/"
JDP_DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/"
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
AUTH_STATUS_CODES = {401, 403}
MIN_TOKEN_SECONDS = 3600

# file_status values from which a download works right now; everything else (PURGED = on
# tape, RESTORE_REGISTERED = restore requested but not finished) answers HTTP 409.
DOWNLOADABLE_STATUSES = frozenset({"RESTORED", "BACKUP_COMPLETE"})

# Not downloadable yet but will become so: on tape, or restore requested and still running.
WAITING_STATUSES = frozenset({"PURGED", "RESTORE_REGISTERED"})

# Generous because the JGI API is slow when busy; a too-short timeout only adds more load
# through retries.
API_TIMEOUT = 120
BACKOFF_BASE_SECONDS = 10
BACKOFF_MAX_SECONDS = 300


class JgiAuthError(RuntimeError):
    """The JGI token is missing, expired or rejected. Retrying cannot help."""


class JgiNotRestoredError(RuntimeError):
    """The file is still on tape; request a restore and download it on a later pass."""


def normalize_token(token):
    """The bare JWT: any number of leading ``Bearer`` prefixes (a double paste is easy to make)
    and surrounding whitespace removed. Every use of the token goes through this."""
    token = (token or "").strip()
    while token.lower().startswith("bearer "):
        token = token[7:].strip()
    return token


def token_seconds_left(token):
    """Seconds until the JWT expires (negative if expired); None if it has no readable ``exp``."""
    try:
        payload = normalize_token(token).split(".")[1]
        data = json.loads(base64.urlsafe_b64decode(payload + "=" * (-len(payload) % 4)))
        return float(data["exp"]) - time.time()
    except Exception:
        return None


def check_token(token, min_seconds=MIN_TOKEN_SECONDS):
    """Raise JgiAuthError if the token is missing or (nearly) expired. Never prints the token."""
    token = normalize_token(token)
    if not token:
        raise JgiAuthError("JGI token missing: set jgi_token in config/secrets.yaml.")
    if token.count(".") != 2:
        raise JgiAuthError("jgi_token in config/secrets.yaml is not a JWT (expected three dot-separated parts); "
                           "copy it again from the JGI token page.")
    left = token_seconds_left(token)
    if left is None:
        return None  # not a readable JWT; let the server decide
    if left <= 0:
        raise JgiAuthError(
            f"JGI token expired {abs(left) / 3600:.1f} h ago: update jgi_token in config/secrets.yaml."
        )
    if left < min_seconds:
        raise JgiAuthError(
            f"JGI token expires in {left / 60:.0f} min (< {min_seconds // 60} min): refresh jgi_token in "
            "config/secrets.yaml before starting."
        )
    return left


def read_token(secrets_path, fallback=None):
    """Current ``jgi_token`` from the secrets file, so a token pasted in while a job waits is
    picked up; falls back to the value the workflow was started with."""
    try:
        import yaml
        with open(secrets_path, "r", encoding="utf-8") as handle:
            token = (yaml.safe_load(handle) or {}).get("jgi_token")
        return normalize_token(token) or fallback
    except Exception:
        return fallback


def make_session(token):
    session = requests.Session()
    session.headers.update({
        "Authorization": f"Bearer {normalize_token(token)}",
        "Accept": "application/json, application/zip, application/octet-stream",
        "Content-Type": "application/json",
    })
    return session


def backoff_delay(attempt, response=None):
    """Seconds to wait after failed attempt number ``attempt`` (1-based).

    Exponential with full-range jitter, so hundreds of jobs that failed together do not retry
    together. A numeric ``Retry-After`` header from the server wins (capped).
    """
    if response is not None:
        retry_after = response.headers.get("Retry-After", "")
        if retry_after.strip().isdigit():
            return min(int(retry_after), 600) + random.uniform(0, 5)
    base = min(BACKOFF_BASE_SECONDS * 2 ** (attempt - 1), BACKOFF_MAX_SECONDS)
    return base * random.uniform(0.5, 1.5)


def request_with_retries(session, method, url, logger, max_attempts=6, **kwargs):
    """API call with bounded retries for transient failures; 401/403 raise JgiAuthError at once."""
    for attempt in range(1, max_attempts + 1):
        try:
            response = session.request(method, url, **kwargs)
        except requests.RequestException as exc:
            if attempt == max_attempts:
                raise
            delay = backoff_delay(attempt)
            logger.warning("Request failed (attempt %d/%d): %s. Retrying in %.0fs", attempt, max_attempts, exc, delay)
            time.sleep(delay)
            continue

        if response.status_code in AUTH_STATUS_CODES:
            response.close()
            raise JgiAuthError(
                f"JGI rejected the token (HTTP {response.status_code}) for {url}: "
                "update jgi_token in config/secrets.yaml."
            )
        if response.status_code not in RETRYABLE_STATUS_CODES or attempt == max_attempts:
            return response
        delay = backoff_delay(attempt, response)
        logger.warning("HTTP %d (attempt %d/%d). Retrying in %.0fs", response.status_code, attempt, max_attempts, delay)
        response.close()
        time.sleep(delay)
    raise RuntimeError(f"API request failed after {max_attempts} attempts")


def fetch_all_pages(base_url, logger, label="", pause=0.5, workers=4):
    """Download every page of a paged ``mycocosm_file_list`` query; returns the merged JSON.

    Unauthenticated, retried with backoff. The server caps pages at 50 organisms and a large
    group takes ~10 s per page, so after page 1 the remaining pages (count taken from the
    response's ``total``) are fetched by a few parallel workers; without a ``total`` the pages are
    walked one by one. Raises if any page fails or the organism count does not match ``total``,
    so callers never see a silently truncated list.
    """
    def get_page(page):
        response = request_with_retries(requests.Session(), "GET", f"{base_url}&p={page}", logger,
                                        timeout=API_TIMEOUT)
        response.raise_for_status()
        logger.info("Fetched page %s %s", page, label)
        return response.json()

    first = get_page(1)
    data, organisms = first, list(first.get("organisms", []))
    total = first.get("total")
    page_size = len(organisms)

    if first.get("next_page") and isinstance(total, int) and page_size:
        last_page = -(-total // page_size)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for page_data in pool.map(get_page, range(2, last_page + 1)):
                organisms.extend(page_data.get("organisms", []))
                data = page_data
    else:
        page = first.get("next_page")
        while page:
            time.sleep(pause)
            data = get_page(page)
            organisms.extend(data.get("organisms", []))
            page = data.get("next_page")

    if isinstance(total, int) and len(organisms) != total:
        raise RuntimeError(f"Expected {total} organisms {label} but received {len(organisms)}")
    data["organisms"] = organisms
    return data


def restore_entry(org, file_obj):
    """(org_key, sub_dict) for the ``ids`` payload of request_archived_files / download_files."""
    org_id = org.get("id") or org.get("_id")
    portal_id = org.get("mycocosm_portal_id") or org.get("portal_id") or org.get("name")
    top_hit = (org.get("top_hit") or {}).get("_id")
    sub = {"file_ids": [file_obj["_id"]]}
    if top_hit:
        sub["top_hit"] = top_hit
    if portal_id:
        sub["mycocosm_portal_id"] = portal_id
    return org_id, sub


def largest_file(org):
    return max(org["files"], key=lambda f: int(f.get("file_size") or 0))
