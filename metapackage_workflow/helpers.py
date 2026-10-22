import json
from functools import lru_cache
from pathlib import Path
from typing import Callable, List, Optional, Tuple

import jq
import polars as pl

from scripts.utils import CacheManager


cache_manager = CacheManager()


def pixi_env(workflow_dir, name, envs_dir=None, check=True):
    """Absolute path of pixi environment `name` for use in a rule's `conda:` directive.

    Snakemake treats a directory path as an existing conda prefix and activates it
    without creating or modifying it. If the directory does not exist it silently falls
    back to treating the path as an environment *name*, which fails later with a cryptic
    "Environment names cannot contain path separators"; so with `check=True` a missing
    environment raises immediately with instructions.

    `envs_dir` overrides the default `<workflow_dir>/.pixi/envs` (e.g. when pixi's
    `detached-environments` is set).
    """
    base = Path(envs_dir) if envs_dir else Path(workflow_dir) / ".pixi" / "envs"
    path = (base / name).resolve()
    if check and not (path / "conda-meta").is_dir():
        raise FileNotFoundError(
            f"pixi environment '{name}' not found at {path}.\n"
            f"Run `pixi install --all` in {workflow_dir} first (needs network). If pixi keeps its "
            f"environments elsewhere, pass --config pixi_envs_dir=/path/to/envs. "
            f"Dry runs (-n) skip this check."
        )
    return path.as_posix()


def _file_signature(path: Path) -> Tuple[bool, int, int]:
    if not path.exists():
        return (False, 0, 0)
    stat = path.stat()
    return (True, stat.st_mtime_ns, stat.st_size)


def _dir_signature(path: Path) -> Tuple[bool, int]:
    if not path.exists():
        return (False, 0)
    return (True, path.stat().st_mtime_ns)


@cache_manager.cache_results
def _read_genome_manifest_cached(
    manifest_path: str,
    signature: Tuple[bool, int, int],
) -> List[str]:
    manifest = Path(manifest_path)
    with open(manifest, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    return sorted(set(lines))


@cache_manager.cache_results
def _build_busco_targets_cached(
    manifest_path: str,
    signature: Tuple[bool, int, int],
) -> List[str]:
    targets = []
    for genome_path in _read_genome_manifest_cached(manifest_path, signature):
        group, genome_acc = genome_path_to_group_genome_acc(genome_path)
        targets.append(
            f"outputs/busco/{group}/{genome_acc}/run_eukaryota_odb12/full_table.tsv"
        )
    return sorted(targets)


@cache_manager.cache_results
def _get_genomes_with_transcripts_cached(
    manifest_file: str,
    signature: Tuple[bool, int, int],
) -> List[str]:
    path = Path(manifest_file)

    if not path.exists() or path.stat().st_size == 0:
        return []

    manifest_df = pl.read_csv(path, separator="\t")
    required_cols = {"group", "genome_acc", "busco_result"}
    if not required_cols.issubset(set(manifest_df.columns)):
        return []

    unique_pairs = manifest_df.select(["group", "genome_acc"]).unique()
    return sorted(
        [
            f"outputs/transcripts/taxonomy/{row['group']}/{row['genome_acc']}.otu_table.tsv"
            for row in unique_pairs.iter_rows(named=True)
        ]
    )


@cache_manager.cache_results
def _get_discovered_busco_genes_cached(
    singlem_input_dir: str,
    signature: Tuple[bool, int],
) -> List[str]:
    base = Path(singlem_input_dir)
    if not base.exists():
        return []

    busco_genes = []
    for faa_file in base.glob("*.faa"):
        busco_gene = faa_file.stem
        tax_file = base / f"{busco_gene}_taxonomy.tsv"
        if tax_file.exists():
            busco_genes.append(busco_gene)

    return sorted(busco_genes)


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
    """Double runtime on each retry, then cap at 1h on final attempt."""

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


@lru_cache(maxsize=128)
@cache_manager.cache_results
def get_jgi_genome_accs(group: str) -> Tuple[List[str], pl.DataFrame]:
    """Parse JGI data and return list of genome accessions for a group."""
    with open(f"outputs/metadata/source/jgi/{group}.json", "r", encoding="utf-8") as handle:
        data = json.load(handle)

    query = jq.compile(
        """
        .organisms[]
        | .files
        | max_by(.file_size)
        | {
            name: .file_name,
            file_id: .file_id,
            genome_acc: .portal_detail_id,
          }
    """
    )

    file_id = pl.DataFrame(query.input(data).all()).unique()
    ncbi_id = (
        pl.read_csv(
            f"outputs/metadata/source/jgi/{group}_meta.csv",
            encoding="latin-1",
            ignore_errors=True,
        )
        .select(["name", "NCBI Taxon", "portal"])
    )
    df = file_id.join(ncbi_id, left_on="genome_acc", right_on="portal", how="left")

    return df["genome_acc"].to_list(), df


def get_genome_name(group: str, genome_acc: str) -> str:
    _, df = get_jgi_genome_accs(group)
    return df.filter(pl.col("genome_acc") == genome_acc)["name"][0]


@lru_cache(maxsize=1)
def get_p10k_genome_accs() -> List[str]:  # type: ignore[misc]
    df = pl.read_csv("outputs/metadata/source/p10k/p10k_meta.csv")

    # filter bad urls
    df = df.filter(pl.col("DNA").str.contains(r"https?://"))

    return tuple(df["Accession"].to_list())

_get_p10k_genome_accs_cached = get_p10k_genome_accs
def get_p10k_genome_accs_list() -> List[str]:  # type: ignore[no-redef]
    return list(_get_p10k_genome_accs_cached())


@lru_cache(maxsize=1)
def _get_p10k_genome_urls():  # type: ignore[misc]
    df = pl.read_csv("outputs/metadata/source/p10k/p10k_meta.csv")
    return df

def get_p10k_genome_url(genome_acc: str) -> str:
    df = _get_p10k_genome_urls()
    return df.filter(pl.col("Accession") == genome_acc)["DNA"][0]


@lru_cache(maxsize=1)
def get_ncbi_genome_accs() -> Tuple[str, ...]:
    metadata = Path("outputs/metadata/source/ncbi/ncbi_genomes.tsv")
    if not metadata.exists():
        return tuple()
    return tuple(
        pl.read_csv(metadata, separator="\t")["genome_acc"].unique().sort().to_list()
    )


@lru_cache(maxsize=128)
def get_successful_genome_accs(group: str) -> Tuple[str, ...]:
    """Get sorted list of successfully downloaded genomes for a group. Cached for DAG speed."""
    successful_ids = []
    genome_accs, _ = get_jgi_genome_accs(group)

    for genome_acc in genome_accs:
        genome_file = Path(f"outputs/genomes/{group}/{genome_acc}.fasta.gz")
        if genome_file.exists() and genome_file.stat().st_size > 0:
            successful_ids.append(genome_acc)

    return tuple(sorted(successful_ids))

def get_successful_genome_accs_list(group: str) -> List[str]:  # type: ignore[no-redef]
    """Return successful genomes as list (public API)."""
    return list(get_successful_genome_accs(group))


@lru_cache(maxsize=1)
def get_successful_p10k_genome_accs() -> Tuple[str, ...]:
    """Get sorted list of successfully downloaded p10k genomes. Cached for DAG speed."""
    successful_ids = []

    for genome_acc in sorted(get_p10k_genome_accs()):
        genome_file = Path(f"outputs/genomes/p10k/{genome_acc}.fasta.gz")
        if genome_file.exists() and genome_file.stat().st_size > 0:
            successful_ids.append(genome_acc)

    return tuple(sorted(successful_ids))

def get_successful_p10k_genome_accs_list() -> List[str]:  # type: ignore[no-redef]
    """Return successful p10k genomes as list (public API)."""
    return list(get_successful_p10k_genome_accs())


def get_expected_genome_paths(
    groups: List[str],
    jgi_checkpoint_getter: Optional[Callable] = None,
    p10k_checkpoint_getter: Optional[Callable] = None,
    ncbi_checkpoint_getter: Optional[Callable] = None,
) -> List[str]:
    """All genome paths listed in the metadata (whether or not they exist on disk).

    Checkpoint getters are only needed inside the DAG, to wait for the metadata files.
    """
    paths = []

    for group in groups:
        if jgi_checkpoint_getter:
            jgi_checkpoint_getter(group)
        genome_accs, _ = get_jgi_genome_accs(group)
        paths.extend(f"outputs/genomes/{group}/{acc}.fasta.gz" for acc in sorted(genome_accs))

    if p10k_checkpoint_getter:
        p10k_checkpoint_getter()
    paths.extend(f"outputs/genomes/p10k/{acc}.fasta.gz" for acc in sorted(set(get_p10k_genome_accs())))

    if ncbi_checkpoint_getter:
        ncbi_checkpoint_getter()
    paths.extend(f"outputs/genomes/ncbi/{acc}.fasta.gz" for acc in get_ncbi_genome_accs())

    return paths


def needs_download(path: str, retry_failed: bool = False) -> bool:
    """True if the genome file is absent, or (with retry_failed) a zero-byte failure marker.

    Zero-byte files are written by the download scripts for assemblies that are
    permanently unavailable (e.g. suppressed). They are left in place and not re-requested
    unless retry_failed is set.
    """
    p = Path(path)
    if not p.exists():
        return True
    return retry_failed and p.stat().st_size == 0


def build_all_genomes_to_download(
    groups: List[str],
    jgi_checkpoint_getter: Callable,
    p10k_checkpoint_getter: Callable,
    ncbi_checkpoint_getter: Callable,
    retry_failed: bool = False,
) -> Callable:
    """Input function returning only genomes that still need a download job.

    Genomes already on disk get no job at all (so no queue time and no rerun-trigger
    re-downloads). The checkpoint getters are still called so the DAG waits for metadata.
    """
    def _get_all_genomes_to_download(wildcards):
        expected = get_expected_genome_paths(
            groups, jgi_checkpoint_getter, p10k_checkpoint_getter, ncbi_checkpoint_getter
        )
        if retry_failed:
            # An existing zero-byte file counts as up to date for Snakemake, so remove the
            # failure marker to make the download rule run again (idempotent).
            for path in expected:
                p = Path(path)
                if p.exists() and p.stat().st_size == 0:
                    p.unlink()
        return [path for path in expected if needs_download(path, retry_failed)]

    return _get_all_genomes_to_download


def get_on_disk_genome_paths(groups: List[str]) -> List[str]:
    """Metadata-listed genomes present on disk, including zero-byte failure markers."""
    return sorted(p for p in get_expected_genome_paths(groups) if Path(p).exists())


def get_existing_genome_paths_from_metadata(groups: List[str]) -> List[str]:
    genome_paths = []

    for group in groups:
        genome_accs, _ = get_jgi_genome_accs(group)
        for genome_acc in sorted(genome_accs):
            path = Path(f"outputs/genomes/{group}/{genome_acc}.fasta.gz")
            if path.exists() and path.stat().st_size > 0:
                genome_paths.append(str(path))

    print(f"Found {len(genome_paths)} existing genomes from JGI metadata.")
    jgi_count = len(genome_paths)

    for genome_acc in sorted(get_p10k_genome_accs()):
        path = Path(f"outputs/genomes/p10k/{genome_acc}.fasta.gz")
        if path.exists() and path.stat().st_size > 0:
            genome_paths.append(str(path))

    print(f"Found {len(genome_paths) - jgi_count} existing genomes from P10K metadata.")
    p10k_count = len(genome_paths) - jgi_count

    for genome_acc in get_ncbi_genome_accs():
        path = Path(f"outputs/genomes/ncbi/{genome_acc}.fasta.gz")
        if path.exists() and path.stat().st_size > 0:
            genome_paths.append(str(path))

    print(f"Found {len(genome_paths) - jgi_count - p10k_count} existing genomes from NCBI metadata.")

    return sorted(genome_paths)


def write_genome_manifest_nls(manifest_path: str, genome_paths: List[str]) -> None:
    out = Path(manifest_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    with open(out, "w", encoding="utf-8") as handle:
        for genome_path in sorted(set(genome_paths)):
            handle.write(f"{genome_path}\n")


def read_genome_manifest_nls(manifest_path: str) -> List[str]:
    manifest = Path(manifest_path)
    if not manifest.exists():
        raise FileNotFoundError(
            f"Genome manifest not found at {manifest_path}. Run phase 1 first."
        )

    return _read_genome_manifest_cached(str(manifest), _file_signature(manifest))


def genome_path_to_group_genome_acc(genome_path: str) -> Tuple[str, str]:
    path = Path(genome_path)
    parts = path.parts

    # outputs/genomes/{group}/{genome_acc}.fasta.gz
    if len(parts) < 4:
        raise ValueError(f"Unexpected genome path format: {genome_path}")

    group = parts[-2]
    genome_acc = path.name.replace(".fasta.gz", "")
    return group, genome_acc


def build_busco_targets_from_manifest(manifest_path: str) -> Callable:
    def _build_busco_targets(wildcards):
        manifest = Path(manifest_path)
        return _build_busco_targets_cached(str(manifest), _file_signature(manifest))

    return _build_busco_targets


def get_genomes_with_transcripts_factory(manifest_file: str) -> Callable:
    def _get_genomes_with_transcripts(wildcards):
        path = Path(manifest_file)
        return _get_genomes_with_transcripts_cached(str(path), _file_signature(path))

    return _get_genomes_with_transcripts


def get_discovered_busco_genes(
    checkpoint_getter: Optional[Callable] = None,
    singlem_input_dir: str = "outputs/singlem/input",
) -> List[str]:
    if checkpoint_getter is not None:
        checkpoint_getter()

    base = Path(singlem_input_dir)
    return _get_discovered_busco_genes_cached(str(base), _dir_signature(base))


def collect_final_packages_factory(
    checkpoint_getter: Optional[Callable] = None,
    final_packages_dir: str = "outputs/singlem/packages/final",
) -> Callable:
    def _collect_final_packages(wildcards):
        busco_genes = get_discovered_busco_genes(checkpoint_getter=checkpoint_getter)
        return [f"{final_packages_dir}/{busco_gene}.spkg" for busco_gene in busco_genes]

    return _collect_final_packages


def collect_package_paths_factory(
    package_dir: str,
    checkpoint_getter: Optional[Callable] = None,
) -> Callable:
    def _collect_package_paths(wildcards):
        if checkpoint_getter is not None:
            checkpoint_getter()

        base = Path(package_dir)
        if not base.exists():
            return []

        return [str(path) for path in sorted(base.glob("*.spkg"))]

    return _collect_package_paths


@cache_manager.cache_results
def get_off_target_representative_accessions(
    archaea_metadata_path: str,
    bacteria_metadata_path: str,
) -> List[str]:
    off_target_reps = (
        pl.scan_csv(
            [archaea_metadata_path, bacteria_metadata_path],
            separator="\t",
            has_header=True,
            ignore_errors=True,
        )
        .select(
            [
                "accession",
                "gtdb_taxonomy",
                "checkm2_completeness",
                "checkm2_contamination",
                "gtdb_representative",
            ]
        )
        .filter(pl.col("gtdb_representative") == "t")
        .with_columns(
            pl.col("gtdb_taxonomy").str.split(";").list.get(-3).alias("family")
        )
        .sort(
            ["family", "checkm2_completeness", "checkm2_contamination"],
            descending=[False, True, False],
        )
        .group_by("family")
        .head(1)
        .collect()
        ["accession"]
        .to_list()
    )

    return sorted(off_target_reps)


def build_off_target_genome_targets_factory(
    archaea_metadata_path: str,
    bacteria_metadata_path: str,
    genome_template: Optional[str] = None,
    gtdb_genomes_database_dir: Optional[str] = None,
) -> Callable:
    def _build_off_target_genome_targets(wildcards):
        accessions = get_off_target_representative_accessions(
            archaea_metadata_path=archaea_metadata_path,
            bacteria_metadata_path=bacteria_metadata_path,
        )

        if gtdb_genomes_database_dir:
            return [
                gtdb_genome_path_from_accession(acc, gtdb_genomes_database_dir)
                for acc in accessions
            ]

        if genome_template:
            return [genome_template.format(accession=acc) for acc in accessions]

        raise ValueError(
            "No off-target input path mode configured. "
            "Set either off_target.input.gtdb_genomes_database_dir "
            "or off_target.input.genome_by_accession."
        )

    return _build_off_target_genome_targets


def normalize_gtdb_accession(accession: str) -> str:
    """
    Normalize GTDB accessions for path construction.

    GTDB metadata commonly prefixes accessions with RS_ or GB_.
    For filesystem paths under the GTDB database tree, we use the bare
    RefSeq/GenBank accession, e.g.:
      RS_GCF_000006685.1 -> GCF_000006685.1
      GB_GCA_000001405.28 -> GCA_000001405.28
    """
    normalized = accession.strip()
    normalized = normalized.replace("RS_", "")
    normalized = normalized.replace("GB_", "")
    return normalized


def gtdb_genome_path_from_accession(accession: str, database_dir: str) -> str:
    """
    Build a GTDB genome path from accession using nested numeric directories.

    Example:
      accession: GCF_000006685.1
      output:    {database_dir}/GCF/000/006/685/GCF_000006685.1_genomic.fna.gz
    """
    normalized = normalize_gtdb_accession(accession)

    # Expected format: GCF_000006685.1 or GCA_000006685.1
    if "_" not in normalized:
        raise ValueError(f"Unexpected accession format: {accession}")

    prefix, suffix = normalized.split("_", 1)
    digits = suffix.split(".", 1)[0]

    if len(digits) < 9 or not digits[:9].isdigit():
        raise ValueError(f"Unexpected accession numeric component: {accession}")

    dir1 = digits[0:3]
    dir2 = digits[3:6]
    dir3 = digits[6:9]

    return str(
        Path(database_dir)
        / prefix
        / dir1
        / dir2
        / dir3
        / f"{normalized}_genomic.fna.gz"
    )


def build_dereplicated_genome_list_factory(
    dereplicated_representatives_file: str,
    checkpoint_getter: Optional[Callable] = None,
) -> Callable:
    """
    Factory to collect dereplicated genome paths from checkpoint output.
    
    Reads dereplicated_genome_representatives.txt (accessions, one per line)
    and expands to full genome paths for downstream rules.
    """
    def _build_dereplicated_genome_list(wildcards):
        if checkpoint_getter is not None:
            checkpoint_getter()

        reps_file = Path(dereplicated_representatives_file)
        if not reps_file.exists():
            raise FileNotFoundError(
                f"Dereplicated representatives file not found: {reps_file}"
            )

        dereplicated_genomes = []
        with open(reps_file, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                
                # Parse genome_acc from path or accession
                # Path format: outputs/genomes/{group}/{genome_acc}.fasta.gz
                # Or just accession: {genome_acc}
                if "/" in line:
                    # Full path
                    dereplicated_genomes.append(line)
                else:
                    # Just accession; need to infer group from filesystem
                    # Try to find the genome in outputs/genomes
                    genome_accs = []
                    for group in ["fungi", "algae", "p10k"]:
                        genome_path = Path(f"outputs/genomes/{group}/{line}.fasta.gz")
                        if genome_path.exists():
                            dereplicated_genomes.append(str(genome_path))
                            break

        if not dereplicated_genomes:
            raise ValueError(
                f"No dereplicated genomes found at {reps_file}"
            )

        return sorted(dereplicated_genomes)

    return _build_dereplicated_genome_list
