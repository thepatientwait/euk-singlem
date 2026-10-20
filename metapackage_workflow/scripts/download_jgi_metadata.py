import json
import os
import time

import requests

from log_config import configure_logger


logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger = configure_logger(logger_name)

def download_all_pages(base_url):
    all_organisms = []
    page = 1
    
    while True:
        url = base_url + f'&p={page}'
        
        logger.info("Fetching page %s for %s", page, snakemake.wildcards.group)
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        all_organisms.extend(data.get('organisms', []))
        
        if not data.get('next_page', False):
            break
            
        page = data.get('next_page', False)
        time.sleep(0.5)
    
    data['organisms'] = all_organisms
    return data

# Download JSON with pagination
json_data = download_all_pages(snakemake.params.json_url)
with open(snakemake.output.json, 'w') as f:
    json.dump(json_data, f)

# Download metadata CSV
response = requests.get(snakemake.params.meta_url)
response.raise_for_status()
with open(snakemake.output.meta, 'wb') as f:
    # replace '##' and '#' with 'number' in the header
    content = response.content.decode('latin-1').replace('##', 'number').replace('#', 'number')
    f.write(content.encode('latin-1')) 

logger.info("Downloaded %s total organisms", len(json_data['organisms']))