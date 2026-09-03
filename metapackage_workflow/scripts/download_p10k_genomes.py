# init
# ----
import os
from pathlib import Path

import requests

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
GENOME_URL          = snakemake.params.genome_url

OUTPUT_FILE         = snakemake.output.genome

# fucntions
# ---------

def download_genome(
    url:    str, 
    output: str
    ) ->    None:
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Check if the request was successful
        with open(output, 'wb') as f:
            f.write(response.content)
        logger.info(f"Successfully downloaded genome from {url}")
    except requests.exceptions.HTTPError as http_err:
        logger.warning(f"HTTP error occurred while downloading genome from {url}: {http_err}")
        Path(output).touch()  # Create an empty file to indicate no genome available
    except Exception as err:
        logger.error(f"An error occurred while downloading genome from {url}: {err}")
        raise

def main():
    download_genome(GENOME_URL, OUTPUT_FILE)


if __name__ == '__main__':
    main()