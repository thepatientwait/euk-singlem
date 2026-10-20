#!/usr/bin/env python3

import os
import gzip
import json
import time
import zipfile
import requests
from pathlib import Path
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

JDP_SEARCH_URL = "https://files.jgi.doe.gov/mycocosm_file_list/"
JDP_RESTORE_URL = "https://files.jgi.doe.gov/request_archived_files/"
JDP_DOWNLOAD_URL = "https://files-download.jgi.doe.gov/download_files/"
MAX_API_ATTEMPTS = 10
RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
RETRY_SLEEP_SECONDS = 5


def request_with_retries(session, method, url, **kwargs):
    """Make an API request with bounded retries for transient failures."""
    max_attempts = kwargs.pop("max_attempts", MAX_API_ATTEMPTS)
    for attempt in range(1, max_attempts + 1):
        try:
            response = session.request(method, url, **kwargs)
        except requests.RequestException as e:
            if attempt == max_attempts:
                raise
            delay = min(RETRY_SLEEP_SECONDS * 2 ** (attempt - 1), 60)
            logger.warning(
                f"API request failed on attempt {attempt}/{max_attempts}: {e}. "
                f"Retrying in {delay}s..."
            )
            time.sleep(delay)
            continue

        if response.status_code not in RETRYABLE_STATUS_CODES:
            return response

        if attempt == max_attempts:
            return response

        delay = min(RETRY_SLEEP_SECONDS * 2 ** (attempt - 1), 60)
        logger.warning(
            f"API request returned HTTP {response.status_code} on attempt "
            f"{attempt}/{max_attempts}. Retrying in {delay}s..."
        )
        response.close()
        time.sleep(delay)

    raise RuntimeError(f"API request failed after {max_attempts} attempts")


def validate_gzip_file(file_path):
    """Ensure the gzip stream can be read to EOF."""
    logger.info(f"Starting gzip integrity validation for: {file_path}")
    try:
        with gzip.open(file_path, 'rb') as handle:
            bytes_read = 0
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                bytes_read += len(chunk)
        logger.info(f"Gzip validation passed successfully. Total uncompressed bytes read: {bytes_read}")
        return True
    except (EOFError, OSError, gzip.BadGzipFile) as e:
        logger.error(f"Gzip validation failed with exception: {type(e).__name__}: {e}")
        return False


def main():
    logger.info("Entering script main execution block...")
    
    # 1. Extract Snakemake variables
    try:
        genome_acc = snakemake.wildcards.genome_acc
        output_file = snakemake.output.genome
        target_filename = snakemake.params.name
        token = snakemake.params.token
    except Exception as e:
        logger.error(f"Failed to pull variables from snakemake object: {e}")
        raise e

    if not token:
        raise ValueError("JGI token missing from Snakemake config/secrets.yaml.")

    # 2. Setup output directories
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    downloaded_file = output_dir / target_filename
    temp_zip_file = output_dir / f"{genome_acc}_download.zip"

    # 3. Configure HTTP session and headers
    session = requests.Session()
    auth_header = token if token.startswith("Bearer ") else f"Bearer {token}"
    session.headers.update({
        "Authorization": auth_header,
        "Accept": "application/json, application/zip, application/octet-stream",
        "Content-Type": "application/json"
    })

    # 4. Prepare Mycocosm file list query parameters
    params = {
        "organism": genome_acc,
        "api_version": 2,
        "a": "false",
        "h": "false",
        "d": "asc",
        "p": 1,
        "x": 50,
        "t": "simple"
    }

    req_prepped = requests.Request('GET', JDP_SEARCH_URL, params=params).prepare()

    # 5. Execute API Request
    response = request_with_retries(
        session,
        req_prepped.method,
        req_prepped.url,
        timeout=30,
    )

    if response.status_code in (401, 403):
        raise RuntimeError("JGI Authentication failed (401/403). Your session token has likely expired.")
    
    if response.status_code != 200:
        logger.error(f"Non-200 Response Text Output:\n{response.text}")
        raise RuntimeError(f"Mycocosm file list query failed with status code {response.status_code}")

    try:
        data = response.json()
    except Exception as e:
        logger.error(f"Failed to decode JSON. Full raw response text output follows:\n{response.text}")
        raise e

    organisms = data.get("organisms", [])
    target_file_obj = None
    target_org_key = None
    target_portal_id = None
    target_top_hit = None
    
    # 6. Iterate and match exact filename, recording proper dictionary grouping keys
    for org in organisms:
        org_id = org.get("id") or org.get("_id")
        portal_id = org.get("mycocosm_portal_id") or org.get("portal_id") or org.get("name")
        top_hit = org.get("top_hit", {}).get("_id")
        
        files = org.get("files", [])
        for f in files:
            fname = f.get("file_name")
            jamo_id = f.get("_id")
            if fname == target_filename and jamo_id:
                target_file_obj = f
                target_org_key = org_id or f"Mycocosm_{genome_acc}"
                target_portal_id = portal_id
                target_top_hit = top_hit
                break
        if target_file_obj:
            break

    # 7. Fallback pattern match if exact filename isn't hit
    if not target_file_obj:
        for org in organisms:
            org_id = org.get("id") or org.get("_id")
            portal_id = org.get("mycocosm_portal_id") or org.get("portal_id") or org.get("name")
            top_hit = org.get("top_hit", {}).get("_id")
            
            files = org.get("files", [])
            for f in files:
                fname = f.get("file_name", "")
                jamo_id = f.get("_id")
                if genome_acc in fname and (fname.endswith(".fa.gz") or fname.endswith(".fasta.gz")) and jamo_id:
                    target_file_obj = f
                    target_org_key = org_id or f"Mycocosm_{genome_acc}"
                    target_portal_id = portal_id
                    target_top_hit = top_hit
                    break
            if target_file_obj:
                break

    if not target_file_obj:
        logger.warning(f"No matching files found on JDP matching organism '{genome_acc}' or target '{target_filename}'.")
        Path(output_file).touch()
        return

    jamo_id = target_file_obj.get("_id")
    file_status = target_file_obj.get("file_status", "UNKNOWN")
    logger.info(f"Target file resolved. ID: {jamo_id}, Status: {file_status}")

    # 8. Handle PURGED (tape-archived) files via restoration request & polling
    if file_status == "PURGED":
        logger.info("File is PURGED on tape archive. Initiating restoration request...")
        
        sub_dict = {"file_ids": [jamo_id]}
        if target_top_hit:
            sub_dict["top_hit"] = target_top_hit
        if target_portal_id:
            sub_dict["mycocosm_portal_id"] = target_portal_id

        restore_payload = {
            "ids": {
                target_org_key: sub_dict
            },
            "send_mail": False,
            "api_version": "2"
        }
        
        try:
            restore_resp = request_with_retries(
                session,
                "POST",
                JDP_RESTORE_URL,
                json=restore_payload,
                timeout=30,
            )
            logger.info(f"Restoration response code: {restore_resp.status_code}")
            logger.debug(f"Restoration response text: {restore_resp.text}")
            restore_resp.raise_for_status()
        except Exception as e:
            logger.error(f"Failed to submit restoration request: {e}")
            raise e

        # Poll until RESTORED
        poll_attempts = 0
        max_poll_attempts = 120
        while poll_attempts < max_poll_attempts:
            poll_attempts += 1
            logger.info(f"Requerying file status (Attempt {poll_attempts}/{max_poll_attempts})...")
            
            try:
                poll_resp = request_with_retries(
                    session,
                    req_prepped.method,
                    req_prepped.url,
                    timeout=30,
                )
                poll_resp.raise_for_status()
                poll_data = poll_resp.json()
            except Exception as e:
                logger.warning(f"Error requerying status: {e}. Retrying in 60s...")
                time.sleep(180)
                continue

            current_status = "UNKNOWN"
            for org in poll_data.get("organisms", []):
                for f in org.get("files", []):
                    if f.get("_id") == jamo_id:
                        current_status = f.get("file_status")
                        break
                if current_status != "UNKNOWN":
                    break

            logger.info(f"Current file status from requery: {current_status}")

            if current_status == "RESTORED":
                logger.info("Files are successfully staged and RESTORED on disk!")
                time.sleep(5)
                break
            elif current_status in ["EXPIRED", "ERROR"]:
                raise RuntimeError(f"JGI file restoration failed with status: {current_status}")

            time.sleep(60)
        else:
            raise TimeoutError("Timed out waiting for JGI tape archive file restoration to complete.")

    # 9. Build the precise documentation-compliant download payload structure
    file_sub_dict = {
        "file_ids": [jamo_id]
    }
    if target_top_hit:
        file_sub_dict["top_hit"] = target_top_hit
    if target_portal_id:
        file_sub_dict["mycocosm_portal_id"] = target_portal_id

    download_payload = {
        "ids": {
            target_org_key: file_sub_dict
        },
        "send_mail": False,
        "api_version": "2"
    }

    logger.info(f"Sending POST request to download archive with proper payload schema: {JDP_DOWNLOAD_URL}")
    for download_attempt in range(1, MAX_API_ATTEMPTS + 1):
        try:
            with request_with_retries(
                session,
                "POST",
                JDP_DOWNLOAD_URL,
                json=download_payload,
                stream=True,
                timeout=180,
                max_attempts=1,
            ) as r:
                logger.info(f"Download Response Status Code: {r.status_code}")
                if r.status_code != 200:
                    logger.error(f"Error response body: {r.text}")
                r.raise_for_status()

                downloaded_bytes = 0
                with open(temp_zip_file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                            downloaded_bytes += len(chunk)
                logger.info(f"Zip archive download complete. Total bytes written: {downloaded_bytes}")
                break
        except (requests.RequestException, OSError) as e:
            if download_attempt == MAX_API_ATTEMPTS:
                logger.error(f"Exception occurred during download: {e}")
                raise
            delay = min(RETRY_SLEEP_SECONDS * 2 ** (download_attempt - 1), 60)
            logger.warning(
                f"Download failed on attempt {download_attempt}/{MAX_API_ATTEMPTS}: {e}. "
                f"Retrying in {delay}s..."
            )
            time.sleep(delay)

    # 10. Extract target file from the downloaded zip archive
    logger.info(f"Extracting target file '{target_filename}' from zip archive...")
    extracted_successfully = False
    try:
        with zipfile.ZipFile(temp_zip_file, 'r') as zf:
            namelist = zf.namelist()
            logger.debug(f"Files inside zip archive: {namelist}")
            
            matching_member = None
            for member in namelist:
                if Path(member).name == target_filename:
                    matching_member = member
                    break
            
            if not matching_member and len(namelist) == 1:
                matching_member = namelist[0]

            if matching_member:
                source_in_zip = zf.open(matching_member)
                with open(downloaded_file, "wb") as out_f:
                    out_f.write(source_in_zip.read())
                extracted_successfully = True
            else:
                raise RuntimeError(f"Target filename '{target_filename}' not found in zip archive contents: {namelist}")
    except Exception as e:
        logger.error(f"Failed to extract from zip archive: {e}")
        raise e
    finally:
        if temp_zip_file.exists():
            temp_zip_file.unlink()

    if not extracted_successfully or not downloaded_file.exists() or downloaded_file.stat().st_size == 0:
        logger.error("Extracted file is missing or 0 bytes.")
        Path(output_file).touch()
        return

    # 11. Validate Gzip integrity and finalize
    if validate_gzip_file(downloaded_file):
        logger.info(f"Renaming staging file to final output: '{output_file}'")
        downloaded_file.rename(output_file)
        logger.info(f"Successfully downloaded, verified, and placed genome for {genome_acc}")
        return

    try:
        downloaded_file.unlink()
    except FileNotFoundError:
        pass
    raise Exception(f"Download failed for {genome_acc}: gzip verification failed.")


if __name__ == "__main__":
    main()