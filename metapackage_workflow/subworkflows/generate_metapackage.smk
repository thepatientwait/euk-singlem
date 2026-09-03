import sys
from pathlib import Path

ROOT = Path(workflow.basedir).resolve().parent
sys.path.insert(0, str(ROOT))

# config
# ------
CONFIG_PATH       = ROOT / "config.yaml"
PATHS_CONFIG_PATH = ROOT / "paths.yaml"
SCRIPTS_DIR       = ROOT / "scripts"

from helpers import (
    build_off_target_genome_targets_factory,
    build_busco_targets_from_manifest,
    build_dereplicated_genome_list_factory,
    collect_final_packages_factory,
    create_resource_profile_factory,
    get_genomes_with_transcripts_factory,
)


configfile: str(CONFIG_PATH)
configfile: str(PATHS_CONFIG_PATH)

# workflow
GROUPS       = config["groups"]
BUSCO_GROUPS = GROUPS + ["p10k"]
DOWNLOAD     = config["download"]
BUSCO        = config["generate"]["busco"]
SINGLEM      = config["generate"]["singlem"]
GRAFTM       = config["generate"]["graftm"]
TRANSCRIPTS  = config["generate"]["transcripts"]
COMMON       = config["common"]
OFF_TARGET   = SINGLEM.get("off_target", {})
OFF_TARGET_CONFIG = config["singlem"].get("off_target", {})
OFF_TARGET_ENABLED = OFF_TARGET_CONFIG.get("enabled", False)
OFF_TARGET_INPUT = OFF_TARGET["input"]
SINGLEM_ENV  = (ROOT / COMMON["environments"]["singlem"]).resolve()

# resources
BASE_MEM_MB    = config["resources"]["base_mem_mb"]
BASE_RUNTIME   = config["resources"]["base_runtime"]
TOTAL_ATTEMPTS = config["resources"]["total_attempts"]


RESOURCE_PROFILE   = create_resource_profile_factory(BASE_MEM_MB, BASE_RUNTIME, TOTAL_ATTEMPTS)
DEFAULT_RESOURCES  = RESOURCE_PROFILE()
BUSCO_RESOURCES    = RESOURCE_PROFILE(busco=True)
BUSCO_DOWNLOAD_PATH = str((ROOT / "busco_downloads").resolve())
BUSCO_LINEAGE_READY = str((ROOT / "outputs" / "metadata" / f".busco_lineage_{config['busco']['lineage']}.ready").resolve())

BUSCO_RESULTS_FROM_MANIFEST = build_busco_targets_from_manifest(DOWNLOAD["manifest"]["genome"])
GENOMES_WITH_TRANSCRIPTS = get_genomes_with_transcripts_factory(
    BUSCO["quality_control_results"]
)
PACKAGE_DIR = (
    str(Path(SINGLEM["packages"]["regenerated_spkg"]).parent)
    if OFF_TARGET_ENABLED
    else str(Path(SINGLEM["packages"]["final_spkg"]).parent)
)
COLLECT_FINAL_PACKAGES = collect_final_packages_factory(
    checkpoint_getter = lambda: checkpoints.aggregate_busco_sequences.get(),
    final_packages_dir = PACKAGE_DIR,
)
COLLECT_OFF_TARGET_GENOMES = build_off_target_genome_targets_factory(
    archaea_metadata_path=OFF_TARGET_INPUT["gtdb_metadata"]["archaea"],
    bacteria_metadata_path=OFF_TARGET_INPUT["gtdb_metadata"]["bacteria"],
    genome_template=OFF_TARGET_INPUT.get("genome_by_accession"),
    gtdb_genomes_database_dir=OFF_TARGET_INPUT.get("gtdb_genomes_database_dir"),
)

DEREP_ENABLED = config["singlem"]["derep"].get("enabled", True)


def BUSCO_MANIFEST_FOR_AGGREGATION(wildcards):
    if DEREP_ENABLED:
        return checkpoints.dereplicate_genomes.get().output.filtered_manifest
    return BUSCO["quality_control_results"]


def COLLECT_GRAFTM_PACKAGES(wildcards):
    busco_genes = checkpoints.aggregate_busco_sequences.get().output.sequences_dir
    busco_gene_dir = Path(str(busco_genes))
    genes = sorted(p.stem for p in busco_gene_dir.glob("*.faa"))
    return [GRAFTM["package"].format(busco_gene=gene) for gene in genes]


def COLLECT_VALID_FINAL_PACKAGES(wildcards):
    checkpoint_output = checkpoints.filter_invalid_singlem_packages.get().output
    valid_manifest = Path(str(checkpoint_output.valid_packages_manifest))

    if not valid_manifest.exists():
        raise FileNotFoundError(f"Valid package manifest not found: {valid_manifest}")

    with open(valid_manifest, "r", encoding="utf-8") as handle:
        packages = [line.strip() for line in handle if line.strip()]

    if len(packages) == 0:
        raise ValueError("No valid SingleM packages were found after filtering")

    return packages


def off_target_id_from_path(path_str):
    name = Path(path_str).name
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name.endswith("_genomic"):
        name = name[: -len("_genomic")]
    return name


OFF_TARGET_GENOME_PATHS = sorted(set(COLLECT_OFF_TARGET_GENOMES(None))) if OFF_TARGET_ENABLED else []
OFF_TARGET_GENOME_MAP = {
    off_target_id_from_path(path): path
    for path in OFF_TARGET_GENOME_PATHS
}
OFF_TARGET_IDS = sorted(OFF_TARGET_GENOME_MAP.keys())


def COLLECT_PER_GENOME_OFF_TARGET_RESULTS(wildcards):
    return [
        SINGLEM["off_target"]["per_genome"]["result"].format(off_target_id=off_target_id)
        for off_target_id in OFF_TARGET_IDS
    ]


rule all:
    input:
        SINGLEM["metapackages"]["final"]


rule init_busco_lineage_cache:
    output:
        ready_marker = BUSCO_LINEAGE_READY
    params:
        lineage = config["busco"]["lineage"],
        download_path = BUSCO_DOWNLOAD_PATH
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    shell:
        """
        mkdir -p {params.download_path}
        mkdir -p $(dirname {output.ready_marker})
        if [ ! -s {params.download_path}/file_versions.tsv ] || [ ! -d {params.download_path}/lineages/{params.lineage} ]; then
            busco --download {params.lineage} --download_path {params.download_path}
        fi
        test -s {params.download_path}/file_versions.tsv
        test -d {params.download_path}/lineages/{params.lineage}
        touch {output.ready_marker}
        """


rule run_busco:
    input:
        genome = DOWNLOAD["genomes"]["by_group"],
        lineage_ready = BUSCO_LINEAGE_READY
    output:
        full_table = BUSCO["full_table"]
    params:
        lineage      = config["busco"]["lineage"],
        mode         = config["busco"]["mode"],
        out_path_dir = lambda wildcards: BUSCO["output_dir"].format(group=wildcards.group),
        download_path = BUSCO_DOWNLOAD_PATH,
        offline       = True,
    threads: config["resources"]["busco_threads"]
    retries: TOTAL_ATTEMPTS
    resources: **BUSCO_RESOURCES
    wildcard_constraints:
        group = "|".join(BUSCO_GROUPS)
    script:
        str(SCRIPTS_DIR / "run_busco.py")


rule busco_results_quality_control:
    """Create a quality-controlled BUSCO manifest for downstream processing."""
    input:
        busco_results = BUSCO_RESULTS_FROM_MANIFEST
    output:
        manifest = BUSCO["quality_control_results"]
    params:
        min_complete_markers = config["busco"]["quality_control"]["min_complete_markers"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    script:
        str(SCRIPTS_DIR / "busco_results_quality_control.py")


def get_mem_mb(base_mem_mb: int) -> Callable:
    def _get_mem_mb(wildcards, attempt):
        return base_mem_mb * (2 ** (attempt - 1))

    return _get_mem_mb


def get_runtime(base_runtime: str) -> Callable:
    """Double runtime on each retry."""

    def _get_runtime(wildcards, attempt):
        time = int("".join([char for char in base_runtime if char.isdigit()]))
        unit = "".join([char for char in base_runtime if not char.isdigit()])
        return f"{int(time * attempt)}{unit}"

    return _get_runtime

rule build_checkm2_quality_report:
    """Build a checkm2-compatible quality report from BUSCO results."""
    input:
        busco_manifest = BUSCO["quality_control_results"]
    output:
        report = config["generate"]["quality_report"]["checkm2_style"]
    retries: TOTAL_ATTEMPTS
    # resources: **DEFAULT_RESOURCES
    resources:
        mem_mb = get_mem_mb(BASE_MEM_MB),
        runtime = get_runtime(BASE_RUNTIME),
    script:
        str(SCRIPTS_DIR / "build_checkm2_quality_report.py")


if DEREP_ENABLED:
    checkpoint dereplicate_genomes:
        """Dereplicate genomes using coverm cluster and emit a filtered BUSCO manifest."""
        input:
            quality_report = config["generate"]["quality_report"]["checkm2_style"],
            busco_manifest = BUSCO["quality_control_results"]
        output:
            clusters = config["generate"]["derep"]["cluster_definitions"],
            representatives = config["generate"]["derep"]["representatives"],
            filtered_manifest = config["generate"]["derep_manifest"]["busco_filtered"]
        params:
            precluster_ani = config["singlem"]["derep"]["precluster_ani"],
            ani = config["singlem"]["derep"]["ani"]
        threads: config["resources"]["coverm_cluster_threads"]
        retries: TOTAL_ATTEMPTS
        resources: **RESOURCE_PROFILE(mem_scale=50, runtime="4h")
        conda:
            str(SINGLEM_ENV)
        script:
            str(SCRIPTS_DIR / "dereplicate_genomes.py")


def BUSCO_MANIFEST_FOR_AGGREGATION(wildcards):
    if DEREP_ENABLED:
        return checkpoints.dereplicate_genomes.get().output.filtered_manifest
    return BUSCO["quality_control_results"]


checkpoint aggregate_busco_sequences:
    input:
        quality_report = config["generate"]["quality_report"]["checkm2_style"],
        metadata_files = [DOWNLOAD["metadata"]["jgi"]["csv"].format(group=group) for group in GROUPS] + [DOWNLOAD["metadata"]["p10k"]["csv"]],
        busco_manifest = BUSCO_MANIFEST_FOR_AGGREGATION,
        taxonomy      = DOWNLOAD["metadata"]["taxonomy"]["mapping"]
    output:
        sequences_dir        = directory(SINGLEM["gene_input"]["sequences_dir"]),
        genome_input_genomes = SINGLEM["genome_input"]["genome_protein_files"],
        genome_input_hmms    = SINGLEM["genome_input"]["hmms_and_names"]
    params:
        min_sequences_per_gene = config["singlem"]["min_sequences_per_gene"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    threads: 16
    script:
        str(SCRIPTS_DIR / "aggregate_busco_sequences.py")


rule graftm_create:
    """Create GraftM packages from BUSCO sequences"""
    input:
        aa_sequences_file = SINGLEM["gene_input"]["gene"]["fasta"],
        taxonomy_file     = SINGLEM["gene_input"]["gene"]["taxonomy"]
    output:
        gpkg      = directory(GRAFTM["package"]),
        alignment = GRAFTM["alignment"]
    params:
        hmm_path = lambda wildcards: COMMON["busco_hmm"]["file"].format(busco_gene=wildcards.busco_gene)
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=2.0)
    conda:
        str(SINGLEM_ENV)
    shell:
        "graftM create --min_aligned_percent 0 --force "
        "--sequences {input.aa_sequences_file} "
        "--hmm {params.hmm_path} "
        "--output {output.gpkg} "
        "--threads {threads} "
        "--no_tree "
        "--taxonomy {input.taxonomy_file}"


rule singlem_seqs:
    """Generate window positions for SingleM packages"""
    input:
        alignment = GRAFTM["alignment"]
    output:
        window_position_file = SINGLEM["windows"]["positions"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem seqs "
        "--alignment {input.alignment} "
        "--alignment-type aa > {output.window_position_file}"


rule singlem_create:
    """Create SingleM packages from GraftM packages"""
    input:
        gpkg                 = GRAFTM["package"],
        window_position_file = SINGLEM["windows"]["positions"]
    output:
        spkg = directory(SINGLEM["packages"]["initial_spkg"])
    params:
        taxonomy_file  = SINGLEM["gene_input"]["gene"]["taxonomy"],
        target_domains = config["singlem"]["target_domains"]
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=0.5, runtime="1h")
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem create "
        "--input-graftm-package {input.gpkg} "
        "--output {output.spkg} "
        "--target-domains {params.target_domains} "
        "--hmm-position $(cat {input.window_position_file}) "
        "--gene-description {wildcards.busco_gene} "
        "--force "
        "--input-taxonomy {params.taxonomy_file}"


rule singlem_chainsaw:
    """Optimize SingleM packages with chainsaw"""
    input:
        spkg = SINGLEM["packages"]["initial_spkg"]
    output:
        chainsaw_spkg = directory(SINGLEM["packages"]["final_spkg"])
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem chainsaw "
        "--input-singlem-package {input.spkg} "
        "--output-singlem-package {output.chainsaw_spkg}"


rule acquire_and_concat_hmms:
    """Collect search HMMs from per-gene GraftM packages for off-target scanning."""
    input:
        gpkgs = COLLECT_GRAFTM_PACKAGES
    output:
        hmm_file = SINGLEM["off_target"]["hmms"]["concatenated"],
        hmm_ids = SINGLEM["off_target"]["hmms"]["ids"],
        done = SINGLEM["off_target"]["hmms"]["done"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    conda:
        str(SINGLEM_ENV)
    script:
        str(SCRIPTS_DIR / "acquire_and_concat_hmms.py")


rule prepare_off_target_per_genome:
    """Run off-target search for one genome and emit normalized per-genome matches."""
    input:
        hmms_done = SINGLEM["off_target"]["hmms"]["done"],
        hmm_file = SINGLEM["off_target"]["hmms"]["concatenated"],
        hmm_ids = SINGLEM["off_target"]["hmms"]["ids"],
        genome = lambda wildcards: OFF_TARGET_GENOME_MAP[wildcards.off_target_id]
    output:
        result = SINGLEM["off_target"]["per_genome"]["result"]
    params:
        hmm_evalue = OFF_TARGET_CONFIG.get("hmm_evalue", 1e-5),
        default_taxonomy = OFF_TARGET_CONFIG.get("default_taxonomy", "d__Bacteria")
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=1.0, runtime="2h")
    wildcard_constraints:
        off_target_id = "[A-Za-z0-9._-]+"
    conda:
        str(SINGLEM_ENV)
    script:
        str(SCRIPTS_DIR / "prepare_off_target_per_genome.py")


rule gather_off_target_sequences:
    """Gather per-genome off-target matches into per-gene FASTA/taxonomy files."""
    input:
        hmms_done = SINGLEM["off_target"]["hmms"]["done"],
        hmm_ids = SINGLEM["off_target"]["hmms"]["ids"],
        per_genome_results = COLLECT_PER_GENOME_OFF_TARGET_RESULTS
    output:
        aggregate_dir = directory(SINGLEM["off_target"]["aggregate"]["dir"]),
        done = SINGLEM["off_target"]["aggregate"]["done"]
    params:
        aggregate_faa_pattern = lambda wildcards: SINGLEM["off_target"]["aggregate"]["sequences"],
        aggregate_tax_pattern = lambda wildcards: SINGLEM["off_target"]["aggregate"]["taxonomy"],
        default_taxonomy = OFF_TARGET_CONFIG.get("default_taxonomy", "d__Bacteria")
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=1.5, runtime="2h")
    conda:
        str(SINGLEM_ENV)
    script:
        str(SCRIPTS_DIR / "gather_off_target_sequences.py")


rule off_target_dup_rename:
    """Ensure duplicate off-target sequence names are unique per BUSCO gene."""
    input:
        prepare_done = SINGLEM["off_target"]["aggregate"]["done"],
        aggregate_dir = SINGLEM["off_target"]["aggregate"]["dir"]
    output:
        seqs = SINGLEM["off_target"]["prepare"]["sequences"],
        taxonomy = SINGLEM["off_target"]["prepare"]["taxonomy"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    conda:
        str(SINGLEM_ENV)
    script:
        str(SCRIPTS_DIR / "rename_off_target_dups.py")


rule singlem_regenerate:
    """Regenerate SingleM packages with optional off-target decoy sequences."""
    input:
        singlem_spkg = SINGLEM["packages"]["initial_spkg"],
        seqs = SINGLEM["gene_input"]["gene"]["fasta"],
        taxonomy = SINGLEM["gene_input"]["gene"]["taxonomy"],
        off_target_seqs = SINGLEM["off_target"]["prepare"]["sequences"],
        off_target_taxonomy = SINGLEM["off_target"]["prepare"]["taxonomy"]
    output:
        spkg = directory(SINGLEM["packages"]["regenerated_spkg"])
    params:
        sequence_prefix = "{busco_gene}~"
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=4.0, runtime="2h")
    conda:
        str(SINGLEM_ENV)
    script:
        str(SCRIPTS_DIR / "regenerate_singlem_with_off_target.py")


checkpoint filter_invalid_singlem_packages:
    """Filter out SingleM packages with invalid window positions."""
    input:
        final_packages = COLLECT_FINAL_PACKAGES
    output:
        valid_packages_manifest = SINGLEM["metapackages"]["valid_packages_manifest"],
        invalid_packages_report = SINGLEM["metapackages"]["invalid_packages_report"]
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    script:
        str(SCRIPTS_DIR / "filter_invalid_singlem_packages.py")


rule create_draft_metapackage:
    """Create draft metapackage from protein packages for transcript processing"""
    input:
        final_packages = COLLECT_VALID_FINAL_PACKAGES,
        invalid_packages_report = lambda wildcards: checkpoints.filter_invalid_singlem_packages.get().output.invalid_packages_report
    output:
        metapackage = directory(SINGLEM["metapackages"]["draft"])
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=4.0, runtime="8h")
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem metapackage "
        "--singlem-packages {input.final_packages} "
        "--no-nucleotide-sdb "
        "--no-taxon-genome-lengths "
        "--metapackage {output.metapackage} "
        "--threads {threads} "


rule extract_transcripts:
    """Extract transcript sequences from genome using BUSCO annotations"""
    input:
        genome     = DOWNLOAD["genomes"]["by_group"],
        full_table = BUSCO["full_table"]
    output:
        transcripts = TRANSCRIPTS["extracted"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    script:
        str(SCRIPTS_DIR / "extract_busco_transcripts.py")


rule eukaryotic_transcripts:
    """Run SingleM on transcript sequences using draft metapackage"""
    input:
        transcripts       = TRANSCRIPTS["extracted"],
        draft_metapackage = SINGLEM["metapackages"]["draft"]
    output:
        otu_table = TRANSCRIPTS["singlem_otu"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem pipe "
        "--sequences {input.transcripts} "
        "--metapackage {input.draft_metapackage} "
        "--no-assign-taxonomy "
        "--otu-table {output.otu_table} "
        "--threads {threads}"


rule assign_eukaryotic_taxonomy:
    """Assign eukaryotic taxonomy to transcript hits"""
    input:
        otu_table        = TRANSCRIPTS["singlem_otu"],
        taxonomy_mapping = DOWNLOAD["metadata"]["taxonomy"]["mapping"]
    output:
        taxonomy_otu = TRANSCRIPTS["taxonomy_otu"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    script:
        str(SCRIPTS_DIR / "assign_eukaryotic_taxonomy.py")


rule aggregate_transcript_taxonomy:
    """Aggregate all per-taxon transcript taxonomy results from QC-passing BUSCO genomes"""
    input:
        busco_manifest  = BUSCO["quality_control_results"],
        taxonomy_tables = GENOMES_WITH_TRANSCRIPTS
    output:
        combined_otu = TRANSCRIPTS["aggregate"]["combined_otu"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    shell:
        """
        tables=({input.taxonomy_tables})
        if [ ${{#tables[@]}} -eq 0 ]; then
            echo "No transcript taxonomy tables discovered from QC-passing BUSCO genomes." >&2
            exit 1
        fi

        head -1 "${{tables[0]}}" > {output.combined_otu}
        for file in "${{tables[@]}}"; do
            tail -n +2 "$file" >> {output.combined_otu}
        done
        """


rule make_nucleotide_sdb:
    """Create nucleotide sequence database from transcript OTU table"""
    input:
        otu_table = TRANSCRIPTS["aggregate"]["combined_otu"]
    output:
        sdb = directory(TRANSCRIPTS["aggregate"]["sdb"])
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=2.0, runtime="2h")
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem makedb "
        "--otu-table {input.otu_table} "
        "--db {output.sdb} "
        "--threads {threads}"


rule create_final_metapackage:
    """Create final metapackage from all SingleM packages with nucleotide SDB"""
    input:
        final_packages = COLLECT_VALID_FINAL_PACKAGES,
        sdb           = TRANSCRIPTS["aggregate"]["sdb"]
    output:
        metapackage = directory(SINGLEM["metapackages"]["final"])
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=8.0, runtime="8h")
    conda:
        str(SINGLEM_ENV)
    shell:
        "singlem metapackage "
        "--singlem-packages {input.final_packages} "
        "--nucleotide-sdb {input.sdb} "
        "--no-taxon-genome-lengths "
        "--calculate-average-num-genes-per-species "
        "--metapackage {output.metapackage} "
        "--threads {threads}"
