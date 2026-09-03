#!/usr/bin/env python3

from pathlib import Path
import polars as pl


def _pick_taxonomy_from_dir(eukcc_dir: Path, hint_file: str) -> str | None:
    candidate_paths = []

    hinted = eukcc_dir / hint_file
    if hinted.exists():
        candidate_paths.append(hinted)

    candidate_paths.extend(sorted(eukcc_dir.glob("*taxonomy*.tsv")))
    candidate_paths.extend(sorted(eukcc_dir.glob("*lineage*.tsv")))
    candidate_paths.extend(sorted(eukcc_dir.glob("*taxonomy*.txt")))
    candidate_paths.extend(sorted(eukcc_dir.glob("*lineage*.txt")))

    for candidate in candidate_paths:
        try:
            text = candidate.read_text(encoding="utf-8", errors="ignore").strip()
            if not text:
                continue

            lines = [line.strip() for line in text.splitlines() if line.strip()]
            for line in reversed(lines):
                if "__" in line or ";" in line:
                    return line
        except OSError:
            continue

    return None


def main() -> None:
    manifest = Path(str(snakemake.input.external_manifest_table))
    output = Path(str(snakemake.output.summary))
    hint_file = str(snakemake.params.taxonomy_hint_file)

    output.parent.mkdir(parents=True, exist_ok=True)

    manifest_df = pl.read_csv(manifest, separator="\t")
    rows = []

    for row in manifest_df.iter_rows(named=True):
        genome_acc = row["genome_acc"]
        eukcc_dir = Path(str(snakemake.params.eukcc_dir_pattern).format(genome_acc=genome_acc))

        taxonomy = _pick_taxonomy_from_dir(eukcc_dir, hint_file)
        rows.append(
            {
                "acc": genome_acc,
                "eukcc_taxonomy": taxonomy,
                "eukcc_dir": str(eukcc_dir),
            }
        )

    pl.DataFrame(rows).write_csv(output, separator="\t")


if __name__ == "__main__":
    main()
