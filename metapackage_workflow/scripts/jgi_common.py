"""Shared helpers for the JGI (JDP / MycoCosm) download scripts."""
import base64
import json
import re
import time

import requests

JDP_SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
JDP_RESTORE_URL = "https://files.jgi.doe.gov/request_archived_files/"
JDP_DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/"
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
AUTH_STATUS_CODES = {401, 403}
MIN_TOKEN_SECONDS = 3600


class JgiAuthError(RuntimeError):
    """The JGI token is missing, expired or rejected. Retrying cannot help."""


class JgiNotRestoredError(RuntimeError):
    """The file is still on tape; request a restore and download it on a later pass."""


def token_seconds_left(token):
    """Seconds until the JWT expires (negative if expired); None if it has no readable ``exp``."""
    try:
        payload = token.replace("Bearer ", "").strip().split(".")[1]
        data = json.loads(base64.urlsafe_b64decode(payload + "=" * (-len(payload) % 4)))
        return float(data["exp"]) - time.time()
    except Exception:
        return None


def check_token(token, min_seconds=MIN_TOKEN_SECONDS):
    """Raise JgiAuthError if the token is missing or (nearly) expired. Never prints the token."""
    if not token:
        raise JgiAuthError("JGI token missing: set jgi_token in config/secrets.yaml.")
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


def make_session(token):
    session = requests.Session()
    session.headers.update({
        "Authorization": token if token.startswith("Bearer ") else f"Bearer {token}",
        "Accept": "application/json, application/zip, application/octet-stream",
        "Content-Type": "application/json",
    })
    return session


def request_with_retries(session, method, url, logger, max_attempts=5, **kwargs):
    """API call with bounded retries for transient failures; 401/403 raise JgiAuthError at once."""
    for attempt in range(1, max_attempts + 1):
        try:
            response = session.request(method, url, **kwargs)
        except requests.RequestException as exc:
            if attempt == max_attempts:
                raise
            delay = min(5 * 2 ** (attempt - 1), 60)
            logger.warning("Request failed (attempt %d/%d): %s. Retrying in %ds", attempt, max_attempts, exc, delay)
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
        delay = min(5 * 2 ** (attempt - 1), 60)
        logger.warning("HTTP %d (attempt %d/%d). Retrying in %ds", response.status_code, attempt, max_attempts, delay)
        response.close()
        time.sleep(delay)
    raise RuntimeError(f"API request failed after {max_attempts} attempts")


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
