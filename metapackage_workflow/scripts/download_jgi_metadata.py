import json
import os
import time

import requests

from jgi_common import API_TIMEOUT, fetch_all_pages, request_with_retries
from log_config import configure_logger


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# Download JSON with pagination
json_data = fetch_all_pages(snakemake.params.json_url, logger, label=f"for {snakemake.wildcards.group}")
with open(snakemake.output.json, 'w') as f:
    json.dump(json_data, f)

# Download metadata CSV
response = request_with_retries(requests.Session(), "GET", snakemake.params.meta_url, logger, timeout=API_TIMEOUT)
response.raise_for_status()
with open(snakemake.output.meta, 'wb') as f:
    # replace '##' and '#' with 'number' in the header
    content = response.content.decode('latin-1').replace('##', 'number').replace('#', 'number')
    f.write(content.encode('latin-1')) 

logger.info("Downloaded %s total organisms", len(json_data['organisms']))