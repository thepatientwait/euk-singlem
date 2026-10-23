# euk-singlem

Snakemake pipelines for building and quality-controlling **eukaryotic SingleM
metapackages**, and for building a phylogenetic tree with taxonomy assigned to
the genomes that go into them.

Remote: <https://github.com/thepatientwait/euk-singlem>

## Workflows

| Directory | Purpose |
| --- | --- |
| `metapackage_workflow/` | Download genomes (JGI MycoCosm/PhycoCosm, NCBI, P10K), QC them with BUSCO/CheckM2, run SingleM + GraftM, screen off-target (bacterial/archaeal) sequences, and assemble draft/final metapackages. |
| `tree_workflow/` | Build a genome tree from metapackage gene families (FAMSA + BMGE + tree building), decorate with PhyloRank, and assign GTDB-rank taxonomy (`assign_taxonomy_v5.py`), including optional outgroup and external (non-metapackage) genomes. |

Output locations are documented in `metapackage_workflow/OUTPUT_LAYOUT.md`;
`tree_workflow` writes under `tree_workflow/output/` (see `tree_workflow/paths.yaml`).

### metapackage_workflow

Two subworkflows, run separately (config lives in `metapackage_workflow/config/`):

```bash
cd metapackage_workflow
pixi install --all   # once: the snakemake driver (default) + per-rule envs (singlem, graftm, busco, coverm)

# phase 1: metadata + genome download (JGI / NCBI / P10K -> unified manifest)
pixi run snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda \
    -s subworkflows/download_data.smk

# phase 2: BUSCO/CheckM2 QC, SingleM packages, off-target screen, metapackage
pixi run snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda \
    -s subworkflows/generate_metapackage.smk
```

Software environments are managed with pixi (`pixi.toml`). Snakemake has no native
pixi support, but its `conda:` directive accepts the path of an existing environment,
so rules use `conda: env("<name>")`, which resolves to `.pixi/envs/<name>`. If pixi
keeps its environments elsewhere, pass `--config pixi_envs_dir=/path/to/envs`.

Notes:

- `config/config.yaml`, `config/paths.yaml`, `config/dryrun.yaml` are versioned;
  `config/secrets.yaml` (JGI API token) is **gitignored** — create it locally:
  ```yaml
  jgi_token: "Bearer <token>"
  ```
- `metapackage_workflow/snakefile` is an orchestrator that runs the two phases above
  in order (`pixi run snakemake --profile ... --use-conda`), each as a nested
  `snakemake` call. `--config key=value` overrides are forwarded to both phases. Each
  phase is skipped once its done-flag (`outputs/orchestration/*.done`) exists; delete a
  flag to force that phase to re-run.
- `--config dryrun=true` runs both phases as `snakemake --dry-run` instead (the same
  as `config/dryrun.yaml`). Each phase's plan is written to
  `logs/orchestration/{download,generate}.dryrun.log`, and the real done-flags are not
  touched. Phase 2 is planned from the current state of the outputs, so it does not
  include jobs that phase 1 would create.
- Every rule writes its stdout/stderr (including subprocess output) to a log in the
  workflow directory via `log:`: `logs/<rule>/<wildcards>.log`, or `logs/<rule>.log`
  for rules without wildcards. Logs are overwritten each time a job runs. The raw PBS
  `.OU`/`.ER` files still go to `~/qsub_logs` (check those for jobs killed before they
  wrote anything, e.g. out of memory).

- **JGI downloads.** Files on JGI tape (`PURGED`) must be restored by JGI first (up to ~24 h).
  `jgi_request_restore` sends one batched request per group; every missing genome then gets
  a `download_jgi_genome` job that polls JGI every `jgi_poll_minutes` (15) and downloads as
  soon as the file is restored. Attempts run 2 h, 4 h and 6 h (~12 h, about the life of the
  token); it re-reads `config/secrets.yaml` before each authenticated call, so a fresh token
  can be pasted in while jobs wait. Downloads run at most `jgi_max_concurrent` (10) at a
  time (`--resources jgi_api=N`, set by the orchestrator); waiting jobs do not count.
  `outputs/metadata/reports/jgi_not_ready.{group}.tsv` lists what is still missing and why.
  The JGI file lists are **not** re-fetched automatically: refresh them with
  `snakemake --forcerun download_jgi_metadata` (or delete `outputs/metadata/source/jgi/*.json`);
  they are only rewritten if their content changed.
- NCBI genomes are fetched as the already-gzipped `*_genomic.fna.gz` from the NCBI FTP
  (md5-checked) and fall back to NCBI Datasets (compressed with pigz, CPUs double per retry).
- `busco.quality_control.use_multi_copy: true` keeps Duplicated BUSCO genes, each represented
  by its top-scoring copy (one sequence per gene per genome, never several copies).
- The NCBI/P10K/JGI sources are consolidated by `scripts/build_unified_manifest.py`;
  `scripts/ncbi_taxid_resolver.py` and `scripts/fetch_ncbi_taxonomy.py` resolve
  taxonomy (Datasets API with Entrez fallback for replaced/suppressed assemblies).

### tree_workflow

```bash
cd tree_workflow
pixi install --all   # once: default (driver) + busco, eukcc, align, taxonomy envs
pixi run snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda -n
```

Trees are built once per `{derep}` x `{copies}` value (`tree.derep_modes`, `tree.copy_modes`
in `config.yaml`), i.e. up to four trees:

- `{derep}`: `all` (every genome of the metapackage input plus the external genomes) or `derep`
  (only the representatives in `metapackage_workflow/outputs/metadata/qc/dereplicated_genome_representatives.txt`,
  made by the metapackage rule `dereplicate_genomes`; external genomes are always kept).
- `{copies}`: `single` (Complete BUSCO genes only) or `multi` (plus the top-scoring copy of
  Duplicated genes). Both come from one aggregation: `<gene>_status.tsv` next to each per-gene
  FASTA says which sequences are Duplicated, and `merge_gene_fastas` filters on it.

Outputs are under `output/{derep}/{copies}/`. FastTree runs with the protein model
(`tree.fasttree_args`, default `-lg -gamma`). Two **global notebooks** (one executed copy each in
`output/notebooks/`, covering every `{derep}/{copies}` variant side by side, mostly as figures; the
numbers and PNG copies of each figure are in `output/notebooks/tables/<notebook>/`):

- `tree_diagnostics` (`notebooks/tree_diagnostics.ipynb`, alignment side): gene occupancy, divergence and
  gene overlap within vs between clades, per-gene separation, single vs multi copy, multicopy top-hit check.
- `tree_accuracy` (`notebooks/tree_accuracy.ipynb`, tree side): rooting-independent clade fragmentation,
  support of the fragmenting splits, tips that look misplaced (sister-group swaps excluded) versus BUSCO
  completeness/duplication and taxonomy status, RED statistics, and Robinson-Foulds distances between variants.

The notebooks in `notebooks/` are templates; the rules execute them and write the copy with outputs
and figures. `scripts/tree_fragmentation.py` can also compare finished trees from the command line.

`publish_trees` (part of `rule all`) copies the final trees and their metadata into a new
`<tree.publish_dir>/<YYYY-MM-DD_HHMMSS>[_subset_<name>]/` directory (default `~/trees`; override with
`--config publish_dir=...`): per `{derep}/{copies}` the `tree.nwk`, `tree.rerooted.nwk`, the RED-decorated
tree and cutoffs, `tree_metadata.tsv` and the reference/target taxonomy, plus the two global notebooks,
`MANIFEST.tsv` (size, md5), `RUN_INFO.txt` (git commit, settings) and a `latest` symlink. Every publish
gets its own directory; it re-runs only when a published file changed.

**Quick check on a subset.** `tree.subset` (or `--config subset=family5`) runs the whole workflow
on a small genome set and writes to `output_subset_<name>/` so the full run is untouched:
`family5` = 5 random genomes per family (about 2.2k) plus 50 external genomes, `jgi_p10k`,
`jgi_p10k_family5`. The subset's BUSCO genes are aggregated on their own
(`aggregate_subset_busco_sequences`), so no full re-aggregation is needed.

PhyloRank's RED is a median over rootings, which makes some branch lengths negative and lets a
child sit below its parent; `assign_taxonomy` ignores those lengths, raises each node's RED to at
least its parent's, and logs how many nodes it changed.

Config: `config.yaml` (metapackage inputs, `external` genomes, `outgroup`,
BUSCO, EukCC, alignment, taxonomy) and `paths.yaml`. Example inputs are
committed as `tree_metadata_example.tsv` and `trusted_taxa_example.txt`.

## Related repositories (managed separately, gitignored here)

- `singlem/`, `singlem-test/` — SingleM checkouts used by the workflows.
- `benchmarking_workflow/` — metapackage benchmarking.

## Version-control notes

- Generated outputs (`output/`, `outputs/`), caches, databases, logs and large
  downloaded test references are excluded via `.gitignore`.
- `metapackage_workflow/test/test.fasta` (119 MB) and
  `metapackage_workflow/test/Coolia_malayensis_transcripts.fasta` (55 MB) are
  downloaded test references, not versioned.
- `main` tracks the stable state; `feature/v2-ncbi-hybrid` carries the NCBI/P10K
  hybrid download + taxonomy rework (metapackage) and the tree-workflow taxonomy
  v5 / outgroup changes.
