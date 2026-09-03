#!/usr/bin/env python3

import subprocess
import gzip
from pathlib import Path


def validate_gzip_file(file_path):
    """
    Ensure the gzip stream can be read to EOF.
    Returns True when valid, False when truncated/corrupt.
    """
    try:
        with gzip.open(file_path, 'rb') as handle:
            while handle.read(1024 * 1024):
                pass
        return True
    except (EOFError, OSError, gzip.BadGzipFile):
        return False

def main():
    genome_acc = snakemake.wildcards.genome_acc
    group = snakemake.wildcards.group
    genome_name = snakemake.params.name
    output_file = snakemake.output.genome
    
    # Change to output directory
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Try to download the genome
    cmd = [
        'python', 
        '/home/n9985158/phd/euk-singlem/metapackage_workflow/scripts/jgi-query.py', 
        genome_acc, 
        '-r', 
        f'^{genome_name}$'
    ]

    downloaded_file = output_dir / genome_name

    def is_known_no_file_error(stderr_text):
        return "AttributeError: 'NoneType' object has no attribute 'values'" in stderr_text

    def run_download_once():
        return subprocess.run(cmd, cwd=output_dir, capture_output=True, text=True)

    result = run_download_once()
    if result.returncode != 0:
        if is_known_no_file_error(result.stderr):
            Path(output_file).touch()
            print(f"No files available for download for {genome_acc} (NoneType error)")
            return
        raise Exception(f"Download failed for {genome_acc}: {result.stderr}")

    if not (downloaded_file.exists() and downloaded_file.stat().st_size > 0):
        Path(output_file).touch()
        print(f"No files available for download for {genome_acc}")
        return

    if validate_gzip_file(downloaded_file):
        downloaded_file.rename(output_file)
        print(f"Successfully downloaded genome for {genome_acc}")
        return

    print(f"Detected corrupt or truncated gzip for {genome_acc}; deleting file")
    try:
        downloaded_file.unlink()
    except FileNotFoundError:
        pass
    raise Exception(f"Download failed for {genome_acc}: gzip integrity check failed")

if __name__ == "__main__":
    main()