#!/usr/bin/env python3
"""Download one file robustly: resume, ``.partial`` then atomic rename, retries, optional checksum.

    fetch_file.py URL OUTPUT [--md5 HEX] [--retries N] [--min-bytes N]

* writes to ``OUTPUT.partial`` (``curl -C -`` resumes it) and renames only after the download and checks passed,
  so a half-written file is never mistaken for a finished one;
* retries with exponential backoff and jitter (servers answer 429/5xx when hit by many jobs);
* refuses an empty or implausibly small file (error pages) and a checksum mismatch (the partial file is
  removed in that case, otherwise a corrupt resume would repeat forever);
* works for ``file://`` URLs as well, which is how the tests run it.
"""
import argparse
import hashlib
import os
import random
import subprocess
import sys
import time
from pathlib import Path


def md5sum(path, chunk=1 << 22):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, output, md5=None, retries=4, min_bytes=1, curl="curl", sleep=time.sleep):
    """Fetch ``url`` to ``output``. Returns the number of attempts used; raises RuntimeError on failure."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".partial")
    last = "no attempt"
    for attempt in range(1, retries + 1):
        cmd = [curl, "-fsSL", "-C", "-", "--retry", "0", "-o", str(partial), url]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 33 and partial.exists():
            # server refused the range request (partial file already complete, or no range support): start over
            partial.unlink()
            result = subprocess.run([c for c in cmd if c not in ("-C", "-")], capture_output=True, text=True)
        if result.returncode == 0:
            size = partial.stat().st_size if partial.exists() else 0
            if size < min_bytes:
                last = f"downloaded file is only {size} bytes"
                partial.unlink(missing_ok=True)
            elif md5 and md5sum(partial) != md5.lower():
                last = f"md5 mismatch (expected {md5})"
                partial.unlink(missing_ok=True)
            else:
                os.replace(partial, output)
                return attempt
        else:
            last = f"curl exit {result.returncode}: {result.stderr.strip()[:200]}"
        if attempt < retries:
            sleep(min(60, 2 ** attempt) + random.uniform(0, 3))
    raise RuntimeError(f"could not download {url} after {retries} attempts: {last}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("url")
    parser.add_argument("output")
    parser.add_argument("--md5")
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--min-bytes", type=int, default=1)
    args = parser.parse_args()
    attempts = download(args.url, args.output, args.md5, args.retries, args.min_bytes)
    print(f"{args.output}: {Path(args.output).stat().st_size} bytes after {attempts} attempt(s)")


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as err:
        sys.exit(str(err))
