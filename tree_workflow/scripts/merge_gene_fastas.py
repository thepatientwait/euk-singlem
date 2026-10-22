#!/usr/bin/env python3

from pathlib import Path

from log_redirect import redirect_to_snakemake_log

# send this rule's stdout/stderr to its log: file
redirect_to_snakemake_log()


def append_if_exists(path: Path, handle) -> None:
    if not path.exists() or path.stat().st_size == 0:
        return
    with open(path, "r", encoding="utf-8") as input_handle:
        for line in input_handle:
            handle.write(line)


def main() -> None:
    base_faa = Path(str(snakemake.input.base_faa))
    external_faa = Path(str(snakemake.params.external_faa))
    output_faa = Path(str(snakemake.output.merged_faa))

    use_external = snakemake.params.use_external
    if not use_external:   
        with open(output_faa, "w", encoding="utf-8") as output_handle:
            append_if_exists(base_faa, output_handle)
        return

    output_faa.parent.mkdir(parents=True, exist_ok=True)

    with open(output_faa, "w", encoding="utf-8") as output_handle:
        append_if_exists(base_faa, output_handle)
        append_if_exists(external_faa, output_handle)


if __name__ == "__main__":
    main()
