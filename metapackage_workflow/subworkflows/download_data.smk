import sys
from pathlib import Path

from snakemake_interface_executor_plugins.settings import ExecMode

ROOT = Path(workflow.basedir).resolve().parent
sys.path.insert(0, str(ROOT))

# config
# ------
CONFIG_DIR          = ROOT / "config"
CONFIG_PATH         = CONFIG_DIR / "config.yaml"
PATHS_CONFIG_PATH   = CONFIG_DIR / "paths.yaml"
SECRETS_CONFIG_PATH = CONFIG_DIR / "secrets.yaml"

configfile: str(CONFIG_PATH)
configfile: str(PATHS_CONFIG_PATH)
configfile: str(SECRETS_CONFIG_PATH)

SCRIPTS_DIR         = ROOT / "scripts"

from helpers import (
    build_all_genomes_to_download,
    create_resource_profile_factory,
    get_existing_genome_paths_from_metadata,
    get_on_disk_genome_paths,
    get_genome_name,
    get_p10k_genome_url,
    get_threads,
    log_path,
    redirect_to_log,
    write_genome_manifest_nls,
)

# workflow
GROUPS       = config["groups"]  # TODO: rename to JGI_GROUPS
GROUP_CONFIG = config["jgi_config"]
DOWNLOAD     = config["download"]
METADATA     = DOWNLOAD["metadata"]
GENOMES      = DOWNLOAD["genomes"]
CLEANUP      = DOWNLOAD["cleanup"]
MANIFEST     = DOWNLOAD["manifest"]

# resources
BASE_MEM_MB    = config["resources"]["base_mem_mb"]
BASE_RUNTIME   = config["resources"]["base_runtime"]
TOTAL_ATTEMPTS = config["resources"]["total_attempts"]

RESOURCE_PROFILE   = create_resource_profile_factory(BASE_MEM_MB, BASE_RUNTIME, TOTAL_ATTEMPTS)
DEFAULT_RESOURCES  = RESOURCE_PROFILE()

# Refresh the JGI file lists (restore statuses change on JGI's side) before the DAG is built.
# Only in the main process of a real run: not in dry runs (no side effects) and not when a
# cluster job re-parses this file. The JSON is replaced only if its content changed, so an
# unchanged catalogue keeps its mtime and triggers no downloads or reruns.
def _flag(name, default=True):
    """Boolean config value; `--config name=false` arrives as the string "false"."""
    return str(config.get(name, default)).strip().lower() not in ("false", "0", "no", "n", "")


class _FormattingLogger:
    """snakemake's logger does not %-format its arguments."""
    def __init__(self, base):
        self._base = base
    def info(self, msg, *args):
        self._base.info(msg % args if args else msg)
    def warning(self, msg, *args):
        self._base.warning(msg % args if args else msg)


if (_flag("refresh_jgi_metadata")
        and workflow.exec_mode == ExecMode.DEFAULT
        and not workflow.dryrun
        and not any(flag in sys.argv for flag in ("-n", "--dry-run", "--dryrun", "--touch", "-t"))):
    from scripts.refresh_jgi_metadata import refresh_group
    from snakemake.logging import logger as _smk_logger
    for _group in GROUPS:
        refresh_group(_group, GROUP_CONFIG[_group]["json_url"],
                      METADATA["jgi"]["json"].format(group=_group), _FormattingLogger(_smk_logger))

# Token preflight: a missing/expired JGI token otherwise fails every JGI job separately, after
# queueing. Real runs stop here; dry runs only warn. Skipped with request_jgi_restores: false
# (NCBI-only work).
if workflow.exec_mode == ExecMode.DEFAULT and _flag("request_jgi_restores"):
    from scripts.jgi_common import JgiAuthError, check_token
    try:
        _token_left = check_token(config.get("jgi_token"))
    except JgiAuthError as _err:
        if workflow.dryrun or "--dry-run" in sys.argv or "-n" in sys.argv:
            from snakemake.logging import logger as _smk_logger
            _smk_logger.warning(f"JGI token check failed (ignored in a dry run): {_err}")
        else:
            raise WorkflowError(f"{_err} (set request_jgi_restores: false to skip JGI.)")

ALL_GENOME_DOWNLOAD_TARGETS = build_all_genomes_to_download(
    groups=GROUPS,
    jgi_checkpoint_getter=lambda group: checkpoints.download_jgi_metadata.get(group=group),
    p10k_checkpoint_getter=lambda: checkpoints.download_p10k_metadata.get(),
    ncbi_checkpoint_getter=lambda: checkpoints.discover_ncbi_genomes.get(),
    # Genomes already on disk (including zero-byte failure markers) get no download job.
    # Set retry_failed_downloads: true in config.yaml to re-request the zero-byte ones.
    retry_failed=bool(config.get("retry_failed_downloads", False)),
    skip_jgi_unrestored=_flag("jgi_skip_unrestored"),
)

localrules: download_jgi_metadata, jgi_request_restore, download_p10k_metadata, discover_ncbi_genomes, build_unified_manifest, fetch_ncbi_taxonomy, cleanup_failed_genomes, write_genome_manifest,


rule all:
    input:
        MANIFEST["genome"],
        CLEANUP["report"],
        CLEANUP["done_flag"],
        MANIFEST["unified"],
        METADATA["taxonomy"]["mapping"]


checkpoint download_jgi_metadata:
    output:
        json = METADATA["jgi"]["json"],
        meta = METADATA["jgi"]["csv"]
    params:
        json_url = lambda wildcards: GROUP_CONFIG[wildcards.group]["json_url"],
        meta_url = lambda wildcards: GROUP_CONFIG[wildcards.group]["meta_url"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("download_jgi_metadata", "group")
    script:
        str(SCRIPTS_DIR / "download_jgi_metadata.py")


rule jgi_request_restore:
    """Ask JGI to restore this group's missing tape (PURGED) genomes; does not wait for them.

    Re-runs whenever the group metadata changes (it is refreshed at the start of each run), which
    re-requests anything still on tape. Also reports every missing genome and why it is missing.
    """
    input:
        json = METADATA["jgi"]["json"]
    output:
        report = METADATA["jgi"]["restore_report"],
        not_ready = METADATA["jgi"]["not_ready_report"]
    params:
        genome_pattern = lambda wildcards: GENOMES["by_group"],
        token = config.get("jgi_token"),
        batch_size = config.get("jgi_restore_batch_size", 50),
        limit = config.get("jgi_restore_limit", 0)
    threads: 1
    resources: **DEFAULT_RESOURCES
    log:
        log_path("jgi_request_restore", "group")
    script:
        str(SCRIPTS_DIR / "jgi_request_restore.py")


rule download_jgi_genome:
    input:
        json = METADATA["jgi"]["json"],
        meta = METADATA["jgi"]["csv"],
        # manifest = MANIFEST["unified"]
    output:
        genome = GENOMES["by_group"]
    params:
        name=lambda wildcards: get_genome_name(wildcards.group, wildcards.genome_acc),
        source="JGI",
        token=config.get("jgi_token")
    threads: 1
    # downloads only files already on disk at JGI (no tape wait): small and short; an auth
    # error fails immediately and a transient error gets one retry
    retries: 1
    # jgi_api: JGI's API answered HTTP 500/429/timeouts when ~600 jobs hit it at once. The
    # orchestrator passes --resources jgi_api=<jgi_max_concurrent>; set that when running
    # this file directly.
    resources:
        jgi_api = 1,
        **RESOURCE_PROFILE(mem_scale=0.25, runtime="1h", fixed_mem=True)
    wildcard_constraints:
        group = "|".join(GROUPS)
    log:
        log_path("download_jgi_genome", "group", "genome_acc")
    script:
        str(SCRIPTS_DIR / "download_jgi_genome.py")


checkpoint download_p10k_metadata:
    output:
        metadata = METADATA["p10k"]["csv"]
    params:
        base_url = "https://ngdc.cncb.ac.cn/gwh/api/esdt/gwh_assemblies/download_table_metadata"
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("download_p10k_metadata")
    script:
        str(SCRIPTS_DIR / "download_p10k_metadata.py")


checkpoint discover_ncbi_genomes:
    input:
        # ancient(): refreshing JGI statuses (restores) must not re-run NCBI discovery and
        # everything downstream of it; the JGI/P10K catalogues only rarely change
        p10k_metadata = ancient(METADATA["p10k"]["csv"]),
        fungi_metadata = ancient(METADATA["jgi"]["csv"].format(group="fungi")),
        algae_metadata = ancient(METADATA["jgi"]["csv"].format(group="algae")),
    output:
        metadata = METADATA["ncbi"]["metadata"]
    params:
        tax_id = config["ncbi"]["tax_id"],
        reference_only = config["ncbi"]["reference_only"],
        include_metazoa_viridiplantae_refs = config["ncbi"].get("include_metazoa_viridiplantae_refs", False),
        email = config["ncbi"]["email"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("discover_ncbi_genomes")
    script:
        str(SCRIPTS_DIR / "discover_ncbi_genomes.py")

# TODO: REMOVE
rule build_unified_manifest:
    """Consolidate NCBI, P10K, and JGI metadata into unified download manifest."""
    input:
        jgi_csvs   = [METADATA["jgi"]["csv"].format(group=group) for group in GROUPS],
        # ancient(): a restore-status refresh must not rebuild the manifest and re-query NCBI taxonomy
        jgi_jsons  = [ancient(METADATA["jgi"]["json"].format(group=group)) for group in GROUPS],
        p10k_csv   = METADATA["p10k"]["csv"],
        ncbi_tsv   = METADATA["ncbi"]["metadata"]
    output:
        manifest = MANIFEST["unified"]
    params:
        ncbi_email = config["ncbi"]["email"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("build_unified_manifest")
    script:
        str(SCRIPTS_DIR / "build_unified_manifest.py")


rule download_ncbi_genome:
    input:
        metadata = METADATA["p10k"]["csv"],
        # manifest = MANIFEST["unified"]
    output:
        genome = GENOMES["p10k"]
    params:
        source="NCBI"
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("download_ncbi_genome", "genome_acc")
    script:
        str(SCRIPTS_DIR / "download_genome.py")


rule download_discovered_ncbi_genome:
    input:
        # No metadata input: a regenerated ncbi_genomes.tsv must not invalidate
        # already-downloaded genomes (download_genome.py skips existing outputs).
        # New accessions still get jobs via the checkpoint-derived target list.
        # metadata = METADATA["ncbi"]["metadata"],
        # manifest = MANIFEST["unified"]
    output:
        genome = GENOMES["ncbi"]
    params:
        source = "NCBI"
    # CPUs double on every retry (only used for the NCBI Datasets fallback with pigz); memory
    # stays fixed because the download is streamed
    threads: get_threads(1)
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=0.5, fixed_mem=True)
    log:
        log_path("download_discovered_ncbi_genome", "genome_acc")
    script:
        str(SCRIPTS_DIR / "download_genome.py")


rule fetch_ncbi_taxonomy:
    input:
        unified_manifest = MANIFEST["unified"]
    output:
        taxonomy = METADATA["taxonomy"]["mapping"]
    params:
        email = config["ncbi"]["email"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("fetch_ncbi_taxonomy")
    script:
        str(SCRIPTS_DIR / "fetch_ncbi_taxonomy.py")


# rule cleanup_failed_genomes:
#     input:
#         genomes = ALL_GENOME_DOWNLOAD_TARGETS
#     output:
#         report = CLEANUP["report"],
#         done   = CLEANUP["done_flag"]
#     retries: TOTAL_ATTEMPTS
#     resources: **DEFAULT_RESOURCES
#     run:
#         import gzip
        
#         def is_valid_gzip(file_path, chunk_size=1024 * 1024):
#             """Validate gzip by fully decompressing it so EOF/trailer errors are detected."""
#             try:
#                 with gzip.open(file_path, "rb") as f:
#                     while f.read(chunk_size):
#                         pass
#                 return True
#             except (EOFError, OSError, gzip.BadGzipFile):
#                 # EOFError: truncated stream / missing trailer
#                 # BadGzipFile/OSError: invalid gzip structure or corrupt data
#                 return False
        
#         report_path = Path(output.report)
#         report_path.parent.mkdir(parents=True, exist_ok=True)

#         removed_rows = []
#         for genome_path in sorted([Path(p) for p in input.genomes]):
#             if not genome_path.exists():
#                 continue
            
#             file_size = genome_path.stat().st_size
#             reason = None
            
#             # Check for empty file
#             if file_size == 0:
#                 reason = "empty_file"
#             # Check for corrupt gzip
#             elif file_size > 0 and not is_valid_gzip(genome_path):
#                 reason = "corrupt_gzip"
            
#             if reason:
#                 genome_acc = genome_path.name.replace(".fasta.gz", "")
#                 group = genome_path.parent.name
#                 genome_path.unlink()
#                 removed_rows.append((group, genome_acc, str(genome_path), reason))

#         with open(report_path, "w", encoding="utf-8") as handle:
#             handle.write("group\tgenome_acc\tgenome_path\treason\n")
#             for group, genome_acc, genome_path, reason in removed_rows:
#                 handle.write(f"{group}\t{genome_acc}\t{genome_path}\t{reason}\n")

#         Path(output.done).touch()


rule cleanup_failed_genomes:
    input:
        # Only genomes still needing a download; being listed here is what triggers downloads.
        # Validation below covers every genome on disk, not just these.
        genomes = ALL_GENOME_DOWNLOAD_TARGETS,
        # ask JGI to restore tape genomes (set request_jgi_restores: false to skip, e.g. when
        # only NCBI genomes are wanted and the JGI token is not current)
        restore_requests = (
            [METADATA["jgi"][key].format(group=group) for group in GROUPS
             for key in ("restore_report", "not_ready_report")]
            if _flag("request_jgi_restores") else []
        )
    output:
        report = CLEANUP["report"],
        done   = CLEANUP["done_flag"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("cleanup_failed_genomes")
    run:
        with redirect_to_log(log[0]):
            import struct
            from pathlib import Path
            from tqdm import tqdm

            from scripts.log_config import configure_logger

            logger_name = os.path.splitext(os.path.basename(__file__))[0]
            logger = configure_logger(logger_name)



            def has_valid_gzip_eof(genome_path):
                if not genome_path.exists():
                    return None
            
                file_size = genome_path.stat().st_size
                reason = None
            
                # Check for empty file or a file too small to contain a gzip header (10 bytes) + trailer (8 bytes)
                if file_size < 18:
                    reason = "empty_or_truncated_file"
                else:
                    try:
                        with open(genome_path, "rb") as f:
                            print(f"Checking {genome_path} (size: {file_size} bytes)")
                            # 1. Check magic header (first 2 bytes)
                            if f.read(2) != b"\x1f\x8b":
                                reason = "invalid_header"
                            else:
                                # 2. Check EOF trailer (last 8 bytes)
                                f.seek(-8, 2)
                                trailer = f.read(8)
                                if len(trailer) < 8:
                                    reason = "truncated_trailer"
                                else:
                                    # Unpack CRC32 and original size (little-endian unsigned 32-bit ints)
                                    crc32, uncompressed_size = struct.unpack("<II", trailer)
                    except (OSError, IOError):
                        reason = "read_error"
            
                if reason:
                    genome_acc = genome_path.name.replace(".fasta.gz", "")
                    group = genome_path.parent.name
                    return (group, genome_acc, str(genome_path), reason)
            
                return None

            report_path = Path(output.report)
            report_path.parent.mkdir(parents=True, exist_ok=True)
            # Validate everything on disk (already-present genomes are not rule inputs).
            # Failures are only reported, never deleted: zero-byte files stay in place as
            # markers so they are not re-requested on every run.
            genome_paths_to_check = get_on_disk_genome_paths(GROUPS)
            logger.info(f"[cleanup_failed_genomes] Checking {len(genome_paths_to_check)} genome files for validity...")

            removed_rows = []
            # for genome_path in tqdm(sorted([Path(p) for p in input.genomes]), desc="Checking genome files"):
            for genome_path in sorted(Path(p) for p in genome_paths_to_check):
                res = has_valid_gzip_eof(genome_path)
                if res:
                    group, genome_acc, genome_path_str, reason = res
                    # genome_path.unlink()
                    removed_rows.append((group, genome_acc, genome_path_str, reason))

            with open(report_path, "w", encoding="utf-8") as handle:
                handle.write("group\tgenome_acc\tgenome_path\treason\n")
                for group, genome_acc, genome_path, reason in removed_rows:
                    handle.write(f"{group}\t{genome_acc}\t{genome_path}\t{reason}\n")

            Path(output.done).touch()


rule write_genome_manifest:
    input:
        cleanup_done = CLEANUP["done_flag"],
        taxonomy     = METADATA["taxonomy"]["mapping"]
    output:
        manifest = MANIFEST["genome"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("write_genome_manifest")
    run:
        with redirect_to_log(log[0]):
            import polars as pl

            genome_paths = get_existing_genome_paths_from_metadata(GROUPS)

            taxonomy_df = pl.read_csv(input.taxonomy, separator='\t')
            if "acc" not in taxonomy_df.columns or "tax_status" not in taxonomy_df.columns:
                raise ValueError(
                    f"Taxonomy mapping {input.taxonomy} lacks 'acc'/'tax_status' columns; "
                    "regenerate it with the fetch_ncbi_taxonomy rule"
                )

            taxonomy_exact = set(
                taxonomy_df.filter(pl.col("tax_status").is_in(["complete", "partial"]))["acc"].to_list()
            )

            valid_genome_paths = []
            missing_taxonomy = []

            for genome_path in genome_paths:
                genome_acc = Path(genome_path).name.replace(".fasta.gz", "")

                if genome_acc in taxonomy_exact:
                    valid_genome_paths.append(genome_path)
                else:
                    missing_taxonomy.append(genome_acc)

            if missing_taxonomy:
                print(
                    "[write_genome_manifest] Excluding {} genomes without taxonomy mapping (examples: {}).".format(
                        len(missing_taxonomy),
                        ", ".join(sorted(missing_taxonomy)[:10]),
                    )
                )

            if len(valid_genome_paths) == 0:
                raise ValueError("No genomes have taxonomy mappings; cannot write manifest")

            write_genome_manifest_nls(output.manifest, valid_genome_paths)
