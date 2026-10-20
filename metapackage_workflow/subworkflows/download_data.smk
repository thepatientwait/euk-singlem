import sys
from pathlib import Path

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
    get_genome_name,
    get_p10k_genome_url,
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

ALL_GENOME_DOWNLOAD_TARGETS = build_all_genomes_to_download(
    groups=GROUPS,
    jgi_checkpoint_getter=lambda group: checkpoints.download_jgi_metadata.get(group=group),
    p10k_checkpoint_getter=lambda: checkpoints.download_p10k_metadata.get(),
    ncbi_checkpoint_getter=lambda: checkpoints.discover_ncbi_genomes.get(),
)

localrules: download_jgi_metadata, download_p10k_metadata, discover_ncbi_genomes, build_unified_manifest, fetch_ncbi_taxonomy, cleanup_failed_genomes, write_genome_manifest,


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
    script:
        str(SCRIPTS_DIR / "download_jgi_metadata.py")


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
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    wildcard_constraints:
        group = "|".join(GROUPS)
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
    script:
        str(SCRIPTS_DIR / "download_p10k_metadata.py")


checkpoint discover_ncbi_genomes:
    input:
        p10k_metadata = METADATA["p10k"]["csv"],
        fungi_metadata = METADATA["jgi"]["csv"].format(group="fungi"),
        algae_metadata = METADATA["jgi"]["csv"].format(group="algae"),
    output:
        metadata = METADATA["ncbi"]["metadata"]
    params:
        tax_id = config["ncbi"]["tax_id"],
        reference_only = config["ncbi"]["reference_only"],
        email = config["ncbi"]["email"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    script:
        str(SCRIPTS_DIR / "discover_ncbi_genomes.py")

# TODO: REMOVE
rule build_unified_manifest:
    """Consolidate NCBI, P10K, and JGI metadata into unified download manifest."""
    input:
        jgi_csvs   = [METADATA["jgi"]["csv"].format(group=group) for group in GROUPS],
        jgi_jsons  = [METADATA["jgi"]["json"].format(group=group) for group in GROUPS],
        p10k_csv   = METADATA["p10k"]["csv"],
        ncbi_tsv   = METADATA["ncbi"]["metadata"]
    output:
        manifest = MANIFEST["unified"]
    params:
        ncbi_email = config["ncbi"]["email"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
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
    script:
        str(SCRIPTS_DIR / "download_genome.py")


rule download_discovered_ncbi_genome:
    input:
        metadata = METADATA["ncbi"]["metadata"],
        # manifest = MANIFEST["unified"]
    output:
        genome = GENOMES["ncbi"]
    params:
        source = "NCBI"
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
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
    script:
        str(SCRIPTS_DIR / "fetch_ncbi_taxonomy_v2.py")


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
        genomes = ALL_GENOME_DOWNLOAD_TARGETS
    output:
        report = CLEANUP["report"],
        done   = CLEANUP["done_flag"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    run:
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
        logger.info(f"[cleanup_failed_genomes] Checking {len(input.genomes)} genome files for validity...")

        removed_rows = []
        # for genome_path in tqdm(sorted([Path(p) for p in input.genomes]), desc="Checking genome files"):
        for genome_path in sorted([Path(p) for p in input.genomes]):
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
    run:
        import polars as pl

        genome_paths = get_existing_genome_paths_from_metadata(GROUPS)

        taxonomy_df = pl.read_csv(input.taxonomy, separator='\t')
        if "acc" not in taxonomy_df.columns or "parsed_taxon_string" not in taxonomy_df.columns:
            raise ValueError(
                f"Taxonomy mapping file missing required columns: {input.taxonomy}"
            )

        taxonomy_exact = set(
            taxonomy_df
            .filter(pl.col("parsed_taxon_string").is_not_null() & (pl.col("parsed_taxon_string") != ""))
            ["acc"]
            .to_list()
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
