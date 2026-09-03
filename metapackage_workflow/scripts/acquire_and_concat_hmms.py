#!/usr/bin/env python3

import json
from pathlib import Path
from typing import Any

from log_config import configure_logger


logger = configure_logger(__file__)

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake: Any = None


def _extract_hmm_name(hmm_path: Path) -> str:
    with open(hmm_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("NAME"):
                parts = line.strip().split()
                if len(parts) >= 2:
                    return parts[1]
    return hmm_path.stem


def main() -> None:
    gpkgs = [Path(p) for p in snakemake.input.gpkgs]
    output_hmm = Path(str(snakemake.output.hmm_file))
    output_ids = Path(str(snakemake.output.hmm_ids))
    output_done = Path(str(snakemake.output.done))

    output_hmm.parent.mkdir(parents=True, exist_ok=True)
    output_ids.parent.mkdir(parents=True, exist_ok=True)
    output_done.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(output_hmm, "w", encoding="utf-8") as hmm_out:
        for gpkg in sorted(gpkgs):
            contents_path = gpkg / "CONTENTS.json"
            if not contents_path.exists():
                logger.warning("Skipping %s, CONTENTS.json not found", gpkg)
                continue

            with open(contents_path, "r", encoding="utf-8") as handle:
                contents = json.load(handle)

            search_hmms = contents.get("search_hmms", [])
            if len(search_hmms) == 0:
                logger.warning("Skipping %s, no search_hmms entry", gpkg)
                continue

            hmm_rel_path = Path(search_hmms[0])
            hmm_path = gpkg / hmm_rel_path
            if not hmm_path.exists():
                logger.warning("Skipping %s, HMM path missing: %s", gpkg, hmm_path)
                continue

            with open(hmm_path, "r", encoding="utf-8") as hmm_in:
                hmm_out.write(hmm_in.read())
                if not str(hmm_path).endswith("\n"):
                    hmm_out.write("\n")

            busco_gene = gpkg.name.replace(".gpkg", "")
            hmm_name = _extract_hmm_name(hmm_path)
            rows.append((busco_gene, hmm_name, str(hmm_path)))

    with open(output_ids, "w", encoding="utf-8") as ids_out:
        ids_out.write("gene\thmm_name\thmm_filepath\n")
        for gene, hmm_name, hmm_path in rows:
            ids_out.write(f"{gene}\t{hmm_name}\t{hmm_path}\n")

    output_done.touch()
    logger.info("Wrote %d HMM mappings", len(rows))


if __name__ == "__main__":
    main()
