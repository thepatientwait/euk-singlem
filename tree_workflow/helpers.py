from pathlib import Path
from typing import Callable, Optional

import polars as pl


class _SafeFormatDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


def create_dirs(path_nodes, local_root: Path, substitutions: dict | None = None) -> None:
    def iter_path_templates(node):
        if isinstance(node, dict):
            for value in node.values():
                yield from iter_path_templates(value)
        elif isinstance(node, str):
            yield node

    substitutions = substitutions or {}

    for path_template in iter_path_templates(path_nodes):
        rendered = path_template.format_map(_SafeFormatDict(substitutions))
        path_str = rendered.split("{", 1)[0].rstrip("/") if "{" in rendered else rendered
        path = Path(path_str)
        if not path.is_absolute():
            path = (local_root / path).resolve()
        directory = path.parent if path.suffix else path
        directory.mkdir(parents=True, exist_ok=True)


def read_nonempty_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip() and not line.startswith("#")]


def genome_acc_from_path(path_str: str) -> str:
    path = Path(path_str)
    name = path.name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def paths_from_cluster_table(cluster_table: Path, base_dir: str) -> list[str]:
    if not cluster_table.exists() or cluster_table.stat().st_size == 0:
        return []

    cluster_df = pl.read_csv(cluster_table, separator="\t", has_header=False)
    first_col = cluster_df.columns[0]

    normalized = (
        cluster_df.select(pl.col(first_col).cast(pl.Utf8).str.strip_chars().alias("raw"))
        .filter(pl.col("raw").is_not_null() & (pl.col("raw") != ""))
        .with_columns(
            pl.when(pl.col("raw").str.starts_with("."))
            .then(pl.lit(base_dir) + pl.col("raw").str.slice(1))
            .otherwise(pl.col("raw"))
            .alias("genome_path")
        )
        .select("genome_path")
        .unique()
        .sort("genome_path")
    )

    return normalized["genome_path"].to_list()


def load_valid_genes(valid_packages_manifest: Path) -> list[str]:
    lines = read_nonempty_lines(valid_packages_manifest)
    genes = []
    for line in lines:
        basename = Path(line).name
        if basename.endswith(".spkg"):
            genes.append(basename[:-5])
    return sorted(set(genes))


def load_external_genome_map(
    *,
    external_enabled: bool,
    external_manifest: Path,
    external_config: dict,
) -> dict[str, str]:
    if not external_enabled:
        return {}

    genome_paths = []
    cluster_table = Path(str(external_config.get("cluster_table", ""))).expanduser()
    cluster_base_dir = str(external_config.get("cluster_base_dir", "")).rstrip("/")

    if cluster_table.exists() and cluster_base_dir:
        genome_paths.extend(paths_from_cluster_table(cluster_table, cluster_base_dir))

    genome_paths.extend(read_nonempty_lines(external_manifest))

    mapping = {}
    for genome_path in genome_paths:
        mapping[genome_acc_from_path(genome_path)] = genome_path

    return dict(sorted(mapping.items()))


def build_external_busco_targets_factory(external_paths: dict, external_group: str, external_accessions: list[str], use_external: bool):
    def external_busco_targets(_wildcards):
        if not use_external:
            return []
        return [
            external_paths["busco"]["full_table"].format(group=external_group, genome_acc=genome_acc)
            for genome_acc in external_accessions
        ]

    return external_busco_targets


def build_external_eukcc_done_targets_factory(external_paths: dict, external_accessions: list[str], use_external: bool):
    def external_eukcc_done_targets(_wildcards):
        if not use_external:
            return []
        return [
            external_paths["eukcc"]["per_genome_done"].format(genome_acc=genome_acc)
            for genome_acc in external_accessions
        ]

    return external_eukcc_done_targets


def build_external_gene_path_factory(external_paths: dict, checkpoints, use_external: bool):
    def external_gene_path(gene: str) -> Path:
        if not use_external:
            return Path(external_paths["aggregate"]["sequences_dir"]) / f"{gene}.faa"
        checkpoint_output = checkpoints.aggregate_external_busco_sequences.get().output.sequences_dir
        return Path(str(checkpoint_output)) / f"{gene}.faa"

    return external_gene_path


def build_merged_genes_factory(valid_genes: list[str], checkpoints, use_external: bool):
    def merged_genes(_wildcards):
        genes = set(valid_genes)
        if use_external:
            try:
                seq_dir = Path(str(checkpoints.aggregate_external_busco_sequences.get().output.sequences_dir))
                if seq_dir.exists() and seq_dir.is_dir():
                    external_genes = [path.stem for path in seq_dir.glob("*.faa")]
                    if external_genes:
                        genes.update(external_genes)
            except Exception:
                # Checkpoint not available during DAG construction; rely on valid_genes only
                pass
        return sorted(genes)

    return merged_genes


def build_merged_bmge_targets_factory(tree_paths: dict, merged_genes_func):
    def merged_bmge_targets(wildcards):
        return [tree_paths["bmge"].format(gene=gene) for gene in merged_genes_func(wildcards)]

    return merged_bmge_targets


def build_merged_taxonomy_input_factory(external_paths: dict, metapackage_taxonomy: Path, use_external: bool):
    def merged_taxonomy_input(_wildcards):
        if use_external:
            return external_paths["metadata"]["ncbi_rebuilt_taxonomy"]
        return str(metapackage_taxonomy)

    return merged_taxonomy_input


def create_resource_profile_factory(
    base_mem_mb: int,
    base_runtime: str,
    total_attempts: int,
) -> Callable:
    def _build_resources(
        mem_scale: float = 1.0,
        runtime: Optional[str] = None,
        busco: bool = False,
    ) -> dict:
        scaled_mem_mb = max(1, int(base_mem_mb * mem_scale))
        runtime_base = runtime or base_runtime
        mem_fn = (
            get_busco_mem_mb(scaled_mem_mb, total_attempts)
            if busco
            else get_mem_mb(scaled_mem_mb, total_attempts)
        )
        return {
            "mem_mb": mem_fn,
            "runtime": get_runtime(runtime_base, total_attempts),
        }

    return _build_resources


def get_mem_mb(base_mem_mb: int, total_attempts: int) -> Callable:
    def _get_mem_mb(wildcards, attempt):
        return base_mem_mb * (2 ** (attempt - 1))

    return _get_mem_mb


def get_runtime(base_runtime: str, total_attempts: int) -> Callable:
    def _get_runtime(wildcards, attempt):
        time = int("".join([char for char in base_runtime if char.isdigit()]))
        unit = "".join([char for char in base_runtime if not char.isdigit()])
        return f"{int(time * attempt)}{unit}"

    return _get_runtime


def get_busco_mem_mb(base_mem_mb: int, total_attempts: int) -> Callable:
    def _get_busco_mem_mb(wildcards, attempt, input):
        try:
            input_size_mb = input.size_mb
        except Exception:
            input_size_mb = 100

        scale_factor = max(1, input_size_mb / 100)
        scaled_mem = int(base_mem_mb * scale_factor)

        if attempt == total_attempts:
            return scaled_mem
        return scaled_mem * (2 ** (attempt - 1))

    return _get_busco_mem_mb
