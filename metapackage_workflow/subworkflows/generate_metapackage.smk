import sys
from pathlib import Path

ROOT = Path(workflow.basedir).resolve().parent
sys.path.insert(0, str(ROOT))

# config
# ------
SCRIPTS_DIR         = ROOT / "scripts"
CONFIG_DIR          = ROOT / "config"
CONFIG_PATH         = CONFIG_DIR / "config.yaml"
PATHS_CONFIG_PATH   = CONFIG_DIR / "paths.yaml"

from helpers import (
    build_off_target_genome_targets_factory,
    build_busco_targets_from_manifest,
    build_dereplicated_genome_list_factory,
    collect_final_packages_factory,
    create_resource_profile_factory,
    get_genomes_with_transcripts_factory,
    log_path,
    pixi_env,
)


configfile: str(CONFIG_PATH)
configfile: str(PATHS_CONFIG_PATH)

# workflow
GROUPS       = config["groups"]
BUSCO_GROUPS = GROUPS + ["p10k", "ncbi", "tsa"]
# Groups that are transcriptomes: BUSCO runs in transcriptome mode (MetaEuk), whose full table has no gene
# coordinates (see scripts/busco_table.py), and they have no introns to survey.
TRANSCRIPTOME_GROUPS = {"tsa"}
DOWNLOAD     = config["download"]
BUSCO        = config["generate"]["busco"]
SINGLEM      = config["generate"]["singlem"]
GRAFTM       = config["generate"]["graftm"]
TRANSCRIPTS  = config["generate"]["transcripts"]
COMMON       = config["common"]
INTRONS      = config["singlem"]["introns"]

# Intron variants (wildcard `introns`): corrected is always built, naive (the negative control) only on request
# (singlem.introns.build_naive, or `--config build_naive=true`).
def _as_bool(value):
    """`--config build_naive=false` can arrive as the string "false"."""
    return str(value).strip().lower() not in ("false", "0", "no", "n", "", "none")


BUILD_NAIVE = _as_bool(config.get("build_naive", INTRONS.get("build_naive", False)))
INTRON_VARIANTS = ["corrected"] + (["naive"] if BUILD_NAIVE else [])

# Quick-check subset (singlem.introns.subset, or `--config subset=family5`): the genomes of the QC manifest are
# reduced to a small set and every output moves to outputs/subset_<name>/. Off-target decoys are skipped.
SUBSET = str(config.get("subset", INTRONS.get("subset", "none")))
SUBSET_ROOT = f"outputs/subset_{SUBSET}"


def _relocate(node):
    """Move a subset run's outputs next to, not over, the real ones."""
    if isinstance(node, dict):
        return {key: _relocate(value) for key, value in node.items()}
    if isinstance(node, str):
        for old, new in (("outputs/singlem/", f"{SUBSET_ROOT}/singlem/"), ("outputs/transcripts/", f"{SUBSET_ROOT}/transcripts/"),
                         ("outputs/metadata/qc/valid_", f"{SUBSET_ROOT}/metadata/qc/valid_"),
                         ("outputs/metadata/reports/", f"{SUBSET_ROOT}/metadata/reports/")):
            if node.startswith(old):
                return new + node[len(old):]
    return node


if SUBSET != "none":
    SINGLEM, GRAFTM, TRANSCRIPTS = _relocate(SINGLEM), _relocate(GRAFTM), _relocate(TRANSCRIPTS)
    SINGLEM["windows"]["notebook"] = f"notebooks/rendered/subset_{SUBSET}/intron_survey.ipynb"
SUBSET_MANIFEST = f"{SUBSET_ROOT}/metadata/busco_quality_control_results.tsv"
# the genomes the survey, the aggregation and the transcript steps work on
SURVEY_MANIFEST = BUSCO["quality_control_results"] if SUBSET == "none" else SUBSET_MANIFEST

OFF_TARGET   = SINGLEM.get("off_target", {})
OFF_TARGET_CONFIG = config["singlem"].get("off_target", {})
OFF_TARGET_ENABLED = OFF_TARGET_CONFIG.get("enabled", False) and SUBSET == "none"
OFF_TARGET_INPUT = OFF_TARGET["input"]

wildcard_constraints:
    introns = "corrected|naive"

# one pixi environment per group of similar rules (see pixi.toml); override the
# location with `--config pixi_envs_dir=/path` if pixi uses detached environments
PIXI_ENVS_USED = {}


def env(name):
    PIXI_ENVS_USED[name] = pixi_env(ROOT, name, config.get("pixi_envs_dir"), check=False)
    return PIXI_ENVS_USED[name]


onstart:
    # a missing env directory would otherwise be mistaken for an env *name* by snakemake
    if not workflow.dryrun:
        for _name in PIXI_ENVS_USED:
            pixi_env(ROOT, _name, config.get("pixi_envs_dir"))


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
MARKER_CDS = config["generate"]["derep"]["marker_cds"]
# marker CDS of every genome in the QC manifest, for dereplicate_genomes
MARKER_CDS_FROM_MANIFEST = get_genomes_with_transcripts_factory(
    BUSCO["quality_control_results"], MARKER_CDS, genome_manifest=DOWNLOAD["manifest"]["genome"])
GENOMES_WITH_TRANSCRIPTS = get_genomes_with_transcripts_factory(
    SURVEY_MANIFEST, TRANSCRIPTS["taxonomy_otu"],
    genome_manifest=DOWNLOAD["manifest"]["genome"] if SUBSET == "none" else None)
COLLECT_FINAL_PACKAGES = collect_final_packages_factory(
    checkpoint_getter = lambda: checkpoints.aggregate_busco_sequences.get(),
    package_template = SINGLEM["packages"]["regenerated_spkg" if OFF_TARGET_ENABLED else "final_spkg"],
    input_dir = SINGLEM["gene_input"]["sequences_dir"],
)
COLLECT_OFF_TARGET_GENOMES = build_off_target_genome_targets_factory(
    archaea_metadata_path=OFF_TARGET_INPUT["gtdb_metadata"]["archaea"],
    bacteria_metadata_path=OFF_TARGET_INPUT["gtdb_metadata"]["bacteria"],
    genome_template=OFF_TARGET_INPUT.get("genome_by_accession"),
    gtdb_genomes_database_dir=OFF_TARGET_INPUT.get("gtdb_genomes_database_dir"),
)

# Whether the metapackage itself is built from dereplicated genomes. The dereplication rule
# below always exists: the tree workflow builds its "derep" tree from the representatives list
# regardless of this switch (run `snakemake outputs/metadata/qc/dereplicated_genome_representatives.txt`).
DEREP_ENABLED = bool(config["singlem"]["derep"].get("enabled", False))


def COLLECT_GRAFTM_PACKAGES(wildcards):
    busco_genes = checkpoints.aggregate_busco_sequences.get().output.sequences_dir
    busco_gene_dir = Path(str(busco_genes))
    genes = sorted(p.stem for p in busco_gene_dir.glob("*.faa"))
    return [GRAFTM["package"].format(busco_gene=gene) for gene in genes]


def COLLECT_INTRON_FILES(key):
    """Per-gene output ``key`` of rule intron_window for every aggregated BUSCO gene."""
    def _collect(wildcards):
        genes = sorted(p.stem for p in Path(str(checkpoints.aggregate_busco_sequences.get().output.sequences_dir)).glob("*.faa"))
        return [SINGLEM["windows"][key].format(busco_gene=gene) for gene in genes]
    return _collect


def intron_window_path(variant):
    """Window position file of one variant with the gene still a wildcard (an output of rule intron_window)."""
    return SINGLEM["windows"]["positions"].replace("{introns}", variant)


def COLLECT_VALID_FINAL_PACKAGES(wildcards):
    checkpoint_output = checkpoints.filter_invalid_singlem_packages.get(introns=wildcards.introns).output
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
        expand(SINGLEM["metapackages"]["final"], introns=INTRON_VARIANTS),
        SINGLEM["windows"]["notebook"],
        [] if SUBSET != "none" else BUSCO_RESULTS_FROM_MANIFEST


rule init_busco_lineage_cache:
    output:
        ready_marker = BUSCO_LINEAGE_READY
    params:
        lineage = config["busco"]["lineage"],
        download_path = BUSCO_DOWNLOAD_PATH
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("init_busco_lineage_cache")
    conda:
        env("busco")
    shell:
        """
        exec > {log} 2>&1
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
        lineage_ready = ancient(BUSCO_LINEAGE_READY)
    output:
        full_table = BUSCO["full_table"]
    params:
        lineage      = config["busco"]["lineage"],
        mode         = lambda wildcards: "transcriptome" if wildcards.group in TRANSCRIPTOME_GROUPS else config["busco"]["mode"],
        out_path_dir = lambda wildcards: BUSCO["output_dir"].format(group=wildcards.group),
        download_path = BUSCO_DOWNLOAD_PATH,
        offline       = True,
        condense      = bool(config["busco"].get("condense", True)),
    threads: config["resources"]["busco_threads"]
    retries: TOTAL_ATTEMPTS
    resources: **BUSCO_RESOURCES
    wildcard_constraints:
        group = "|".join(BUSCO_GROUPS)
    log:
        log_path("run_busco", "group", "genome_acc")
    conda:
        env("busco")
    script:
        str(SCRIPTS_DIR / "run_busco.py")


rule busco_results_quality_control:
    """Create a quality-controlled BUSCO manifest for downstream processing."""
    input:
        busco_results = BUSCO_RESULTS_FROM_MANIFEST
    output:
        manifest = BUSCO["quality_control_results"]
    params:
        min_complete_markers = config["busco"]["quality_control"]["min_complete_markers"],
        use_multi_copy = config["busco"]["quality_control"].get("use_multi_copy", False)
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("busco_results_quality_control")
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
    params:
        tsa_penalty = config["singlem"]["derep"].get("tsa_quality_penalty", 0)
    retries: TOTAL_ATTEMPTS
    # resources: **DEFAULT_RESOURCES
    resources:
        mem_mb = get_mem_mb(BASE_MEM_MB),
        runtime = get_runtime(BASE_RUNTIME),
    log:
        log_path("build_checkm2_quality_report")
    script:
        str(SCRIPTS_DIR / "build_checkm2_quality_report.py")


rule extract_marker_cds:
    """Spliced CDS of a genome's (or transcriptome's) BUSCO marker genes: the sequence it is dereplicated on."""
    input:
        genome     = DOWNLOAD["genomes"]["by_group"],
        full_table = BUSCO["full_table"]
    output:
        cds = MARKER_CDS
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    wildcard_constraints:
        group = "|".join(BUSCO_GROUPS)
    log:
        log_path("extract_marker_cds", "group", "genome_acc")
    script:
        str(SCRIPTS_DIR / "extract_marker_cds.py")


checkpoint dereplicate_genomes:
    """Dereplicate genomes and transcriptomes with galah on marker CDS and emit a filtered BUSCO manifest."""
    input:
        quality_report = config["generate"]["quality_report"]["checkm2_style"],
        busco_manifest = BUSCO["quality_control_results"],
        marker_cds = MARKER_CDS_FROM_MANIFEST
    output:
        clusters = config["generate"]["derep"]["cluster_definitions"],
        representatives = config["generate"]["derep"]["representatives"],
        filtered_manifest = config["generate"]["derep_manifest"]["busco_filtered"]
    params:
        precluster_ani = config["singlem"]["derep"]["precluster_ani"],
        ani = config["singlem"]["derep"]["ani"],
        min_aligned_fraction = config["singlem"]["derep"]["min_aligned_fraction"],
        precluster_method = config["singlem"]["derep"]["precluster_method"],
        cluster_method = config["singlem"]["derep"]["cluster_method"],
        cds_template = lambda wildcards: MARKER_CDS
    threads: config["resources"]["coverm_cluster_threads"]
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=50, runtime="4h")
    log:
        log_path("dereplicate_genomes")
    conda:
        env("galah")
    script:
        str(SCRIPTS_DIR / "dereplicate_genomes.py")


def BUSCO_MANIFEST_FOR_AGGREGATION(wildcards):
    if SUBSET != "none":
        return SUBSET_MANIFEST
    if DEREP_ENABLED:
        return checkpoints.dereplicate_genomes.get().output.filtered_manifest
    return BUSCO["quality_control_results"]


if SUBSET != "none":
    rule make_subset_manifest:
        """Quick-check subset of the genomes in the BUSCO QC manifest (singlem.introns.subset)."""
        input:
            taxonomy = DOWNLOAD["metadata"]["taxonomy"]["mapping"]
        output:
            manifest = SUBSET_MANIFEST
        params:
            mode = SUBSET,
            per_family = int(INTRONS.get("subset_per_family", 5)),
            seed = int(INTRONS.get("subset_seed", 1)),
            qc_manifest = BUSCO["quality_control_results"]  # read, not tracked: see scripts/subset_busco_manifest.py
        log:
            log_path("make_subset_manifest")
        script:
            str(SCRIPTS_DIR / "subset_busco_manifest.py")

# TODO: ADD NCBI GENOMES...

checkpoint aggregate_busco_sequences:
    input:
        quality_report = config["generate"]["quality_report"]["checkm2_style"] if SUBSET == "none" else SUBSET_MANIFEST,
        busco_manifest = BUSCO_MANIFEST_FOR_AGGREGATION,
        taxonomy      = DOWNLOAD["metadata"]["taxonomy"]["mapping"]
    output:
        sequences_dir        = directory(SINGLEM["gene_input"]["sequences_dir"]),
        genome_input_genomes = SINGLEM["genome_input"]["genome_protein_files"],
        genome_input_hmms    = SINGLEM["genome_input"]["hmms_and_names"]
    params:
        min_sequences_per_gene = config["singlem"]["min_sequences_per_gene"],
        use_multi_copy = config["busco"]["quality_control"].get("use_multi_copy", False),
        multi_copy_selection = config["busco"]["quality_control"].get("multi_copy_selection", "top_score"),
        multi_copy_min_gain = config["busco"]["quality_control"].get("multi_copy_min_gain", 0.02),
        # Keep sequences/genomes lacking a usable taxonomy (default: exclude them)
        allow_missing_taxonomy = config["singlem"].get("allow_missing_taxonomy", False)
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    threads: 16
    log:
        log_path("aggregate_busco_sequences")
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
    log:
        log_path("graftm_create", "busco_gene")
    conda:
        env("graftm")
    shell:
        "exec > {log} 2>&1; "
        "graftM create --min_aligned_percent 0 --force "
        "--sequences {input.aa_sequences_file} "
        "--hmm {params.hmm_path} "
        "--output {output.gpkg} "
        "--threads {threads} "
        "--no_tree "
        "--taxonomy {input.taxonomy_file}; "
        # the aligned FASTA must have exactly one column per HMM match state: SingleM window positions are
        # HMM match columns (stop symbols in the input used to add extra columns)
        "hmm_len=$(awk '/^LENG/{{print $2; exit}}' {params.hmm_path}); "
        "aln_len=$(awk '!/^>/{{print length($0); exit}}' {output.alignment}); "
        "if [ \"$hmm_len\" != \"$aln_len\" ]; then "
        "echo \"ERROR: alignment width $aln_len != HMM length $hmm_len\" >&2; exit 1; fi"


rule intron_window:
    """Survey introns of one BUSCO gene and choose the SingleM window of each variant.

    SingleM windows are matched by metagenomic reads, which cannot cross an intron. This maps the introns of
    every genome (BUSCO's miniprot gene models) onto the HMM and writes two windows, in HMM match columns (what
    `singlem create --hmm-position` and `singlem pipe` use): naive, the best-covered window, and corrected, the
    best-covered window that crosses an intron in at most max_span_fraction of the genomes covering it.
    The per-genome and per-window tables feed the intron notebook.
    """
    input:
        sequences      = SINGLEM["gene_input"]["gene"]["fasta"],
        hmm            = lambda wildcards: COMMON["busco_hmm"]["file"].format(busco_gene=wildcards.busco_gene),
        busco_manifest = SURVEY_MANIFEST
    output:
        position_naive     = intron_window_path("naive"),
        position_corrected = intron_window_path("corrected"),
        introns  = SINGLEM["windows"]["introns"],
        genomes  = SINGLEM["windows"]["genomes"],
        scan     = SINGLEM["windows"]["scan"],
        choice   = SINGLEM["windows"]["choice"]
    params:
        min_intron        = INTRONS["min_intron"],
        max_span_fraction = INTRONS["max_span_fraction"],
        score_tolerance   = INTRONS["score_tolerance"],
        min_valid         = INTRONS["min_valid_genomes"],
        window_nt         = 60,  # SingleM's default window
        use_multi_copy    = config["busco"]["quality_control"].get("use_multi_copy", False)
    threads: 2
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=2.0, runtime="2h")
    log:
        log_path("intron_window", "busco_gene")
    conda:
        env("singlem")
    script:
        str(SCRIPTS_DIR / "intron_window.py")


rule window_intron_summary:
    """Collect the per-gene window decisions; with introns.strict, fail if any window still spans introns."""
    input:
        choices = COLLECT_INTRON_FILES("choice")
    output:
        summary = SINGLEM["windows"]["summary"]
    params:
        strict = bool(INTRONS.get("strict", False))
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("window_intron_summary")
    script:
        str(SCRIPTS_DIR / "window_intron_summary.py")


rule intron_survey_notebook:
    """Executed notebook: intron presence, conservation, and where the SingleM windows fall."""
    input:
        template = str(ROOT / "notebooks" / "intron_survey.ipynb"),
        choices  = COLLECT_INTRON_FILES("choice"),
        genomes  = COLLECT_INTRON_FILES("genomes"),
        introns  = COLLECT_INTRON_FILES("introns"),
        scans    = COLLECT_INTRON_FILES("scan"),
        taxonomy = DOWNLOAD["metadata"]["taxonomy"]["mapping"],
        gene_names = str(Path(COMMON["busco_hmm"]["dir"]).parent / "links_to_ODB12.txt")  # BUSCO id -> protein name
    output:
        notebook = SINGLEM["windows"]["notebook"]
    params:
        notebook = str(ROOT / "notebooks" / "intron_survey.ipynb"),
        nb_params = lambda wildcards, input: {
            "intron_dir": str(Path(SINGLEM["windows"]["choice"].format(busco_gene="x")).parent.resolve()),
            "taxonomy": str(Path(input.taxonomy).resolve()),
            "gene_names": str(Path(input.gene_names).resolve()),
            "max_span_fraction": INTRONS["max_span_fraction"],
            "build_naive": BUILD_NAIVE,
            "subset": SUBSET,
            "out_dir": str((Path(SINGLEM["windows"]["notebook"]).parent / "tables").resolve()),
        }
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=4)
    log:
        log_path("intron_survey_notebook")
    conda:
        env("notebook")
    script:
        str(SCRIPTS_DIR / "run_notebook.py")


rule singlem_create:
    """Create SingleM packages from GraftM packages"""
    input:
        gpkg                 = GRAFTM["package"],
        window_position_file = lambda wildcards: SINGLEM["windows"]["positions"].format(busco_gene=wildcards.busco_gene, introns=wildcards.introns)
    output:
        spkg = directory(SINGLEM["packages"]["initial_spkg"])
    params:
        taxonomy_file  = SINGLEM["gene_input"]["gene"]["taxonomy"],
        target_domains = config["singlem"]["target_domains"]
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=0.5, runtime="1h")
    log:
        log_path("singlem_create", "introns", "busco_gene")
    conda:
        env("singlem")
    shell:
        "exec > {log} 2>&1; "
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
    log:
        log_path("singlem_chainsaw", "introns", "busco_gene")
    conda:
        env("singlem")
    shell:
        "exec > {log} 2>&1; "
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
    log:
        log_path("acquire_and_concat_hmms")
    conda:
        env("singlem")
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
    log:
        log_path("prepare_off_target_per_genome", "off_target_id")
    conda:
        env("singlem")
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
    log:
        log_path("gather_off_target_sequences")
    conda:
        env("singlem")
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
    log:
        log_path("off_target_dup_rename", "busco_gene")
    conda:
        env("singlem")
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
    log:
        log_path("singlem_regenerate", "introns", "busco_gene")
    conda:
        env("singlem")
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
    log:
        log_path("filter_invalid_singlem_packages", "introns")
    script:
        str(SCRIPTS_DIR / "filter_invalid_singlem_packages.py")


rule create_draft_metapackage:
    """Create draft metapackage from protein packages for transcript processing"""
    input:
        final_packages = COLLECT_VALID_FINAL_PACKAGES,
        invalid_packages_report = lambda wildcards: checkpoints.filter_invalid_singlem_packages.get(introns=wildcards.introns).output.invalid_packages_report,
        # window check: stops here with introns.strict (the naive variant is a control and is not gated)
        window_summary = lambda wildcards: SINGLEM["windows"]["summary"] if wildcards.introns == "corrected" else []
    output:
        metapackage = directory(SINGLEM["metapackages"]["draft"])
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **RESOURCE_PROFILE(mem_scale=4.0, runtime="8h")
    log:
        log_path("create_draft_metapackage", "introns")
    conda:
        env("singlem")
    shell:
        "exec > {log} 2>&1; "
        "singlem metapackage "
        "--singlem-packages {input.final_packages} "
        "--no-nucleotide-sdb "
        "--no-taxon-genome-lengths "
        "--metapackage {output.metapackage} "
        "--threads {threads} "


rule extract_transcripts:
    """Extract BUSCO gene sequences from the genome: corrected = cut at introns (one record per exon block), naive = whole gene span"""
    input:
        genome     = DOWNLOAD["genomes"]["by_group"],
        full_table = BUSCO["full_table"]
    output:
        transcripts = TRANSCRIPTS["extracted"]
    params:
        min_intron  = INTRONS["min_intron"],
        min_exon_nt = INTRONS["min_exon_nt"],
        mode        = lambda wildcards: "exon_blocks" if wildcards.introns == "corrected" else "gene_span"
    threads: 1
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        log_path("extract_transcripts", "introns", "group", "genome_acc")
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
    log:
        log_path("eukaryotic_transcripts", "introns", "group", "genome_acc")
    conda:
        env("singlem")
    shell:
        "exec > {log} 2>&1; "
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
    log:
        log_path("assign_eukaryotic_taxonomy", "introns", "group", "genome_acc")
    script:
        str(SCRIPTS_DIR / "assign_eukaryotic_taxonomy.py")


rule aggregate_transcript_taxonomy:
    """Aggregate all per-taxon transcript taxonomy results from QC-passing BUSCO genomes"""
    input:
        busco_manifest  = SURVEY_MANIFEST,
        taxonomy_tables = GENOMES_WITH_TRANSCRIPTS
    output:
        combined_otu = TRANSCRIPTS["aggregate"]["combined_otu"]
    threads: 1
    retries: TOTAL_ATTEMPTS
    log:
        log_path("aggregate_transcript_taxonomy", "introns")
    conda:
        env("singlem")
    resources: **DEFAULT_RESOURCES
    shell:
        """
        exec > {log} 2>&1
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
    log:
        log_path("make_nucleotide_sdb", "introns")
    conda:
        env("singlem")
    shell:
        "exec > {log} 2>&1; "
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
    log:
        log_path("create_final_metapackage", "introns")
    conda:
        env("singlem")
    shell:
        "exec > {log} 2>&1; "
        "singlem metapackage "
        "--singlem-packages {input.final_packages} "
        "--nucleotide-sdb {input.sdb} "
        "--no-taxon-genome-lengths "
        "--calculate-average-num-genes-per-species "
        "--metapackage {output.metapackage} "
        "--threads {threads}"
