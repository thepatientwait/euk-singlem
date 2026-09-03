# init
# ----
import os

import requests
import json
import zipfile
import io
import csv

# logging init
# ------------
from log_config import configure_logger

logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

# snakefile parameters
# --------------------
URL     = snakemake.params.base_url

OUTPUT  = snakemake.output.metadata


# fucntions
# ---------

def download_metadata(
    url:    str, 
    output: str
    ) ->    None:

    payload = {
        "recordNumber": 3734,
        "lastInput": json.dumps({
            "draw": 2,
            "columns": [
                {
                    "data": "scientificName",
                    "name": "scientificName",
                    "searchable": True,
                    "orderable": True,
                    "search": {"value": "", "regex": False},
                },
                {
                    "data": "group",
                    "name": "group",
                    "searchable": True,
                    "orderable": True,
                    "search": {"value": "Protists", "regex": True},
                },
                {
                    "data": "source",
                    "name": "source",
                    "searchable": True,
                    "orderable": True,
                    "search": {"value": "NCBI", "regex": True},
                },
            ],
            "order": [{"column": 1, "dir": "asc"}],
            "start": 0,
            "length": 10,
            "search": {"value": "", "regex": False},
        })
    }

    logger.info("Fetching data from GWH...")
    response = requests.post(url, data=payload)

    if response.status_code == 200:
        # Open the ZIP in memory
        with zipfile.ZipFile(io.BytesIO(response.content)) as z:
            # Get the name of the internal file (whatever the server named it)
            internal_filename = z.namelist()[0]
            
            # Read the content of that internal file and write it to our custom name
            with z.open(internal_filename) as source:
                content = source.read().decode("utf-8")

            with open(output, "w", newline="", encoding="utf-8") as target:
                writer = csv.writer(target)
                reader = csv.reader(content.splitlines(), delimiter="\t")
                for row in reader:
                    writer.writerow(row)
                
        logger.info(f"Success! File saved as {output}")
    else:
        logger.info(f"Failed. Server returned status code: {response.status_code}")

def main():
    download_metadata(URL, OUTPUT)

if __name__ == "__main__":
    main()