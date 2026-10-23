"""Fetch the JGI file list (JSON) and catalogue (CSV) of one group.

Not re-run automatically: redo it with ``--forcerun download_jgi_metadata`` or by deleting the
outputs. Files are only rewritten when their content changed.
"""
import os

import requests

from jgi_common import API_TIMEOUT, fetch_all_pages, request_with_retries
from log_config import configure_logger
from refresh_jgi_metadata import write_bytes_if_changed, write_json_if_changed

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

group = snakemake.wildcards.group

json_data = fetch_all_pages(snakemake.params.json_url, logger, label=f"for {group}")
write_json_if_changed(snakemake.output.json, json_data, logger, label=group)

response = request_with_retries(requests.Session(), "GET", snakemake.params.meta_url, logger, timeout=API_TIMEOUT)
response.raise_for_status()
# replace '##' and '#' with 'number' in the header
content = response.content.decode('latin-1').replace('##', 'number').replace('#', 'number')
write_bytes_if_changed(snakemake.output.meta, content.encode('latin-1'))

logger.info("Fetched %s total organisms", len(json_data['organisms']))
