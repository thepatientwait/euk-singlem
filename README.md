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
pixi install --all   # once: the snakemake driver (default) + per-rule envs (singlem, graftm, busco, galah, notebook)

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

- **TSA transcriptomes.** NCBI Transcriptome Shotgun Assembly projects are added as group `tsa` (config `tsa:`).
  `discover_tsa` (checkpoint, `scripts/discover_tsa.py`) queries Entrez for non-Metazoa/Viridiplantae `tsa-master`
  records (lineage re-checked), drops taxa a genome already covers (`only_novel_taxa`, judged by species name or tax ID,
  `novel_level`), looks up NCBI's bulk files (`genbank/tsa/<letter>/tsa.<PREFIX>.<n>.fsa_nt.gz`), keeps the largest
  project per tax ID above `min_size_mb`, and writes `outputs/metadata/source/tsa/tsa_projects.tsv`.
  `download_tsa` fetches the parts (Content-Length + gzip check; NCBI publishes no md5) into
  `outputs/genomes/tsa/<master accession>.fasta.gz`. TSA enters the unified manifest, taxonomy fetch and
  `tax_status` gate like any genome. The taxonomy mapping is read by path (it depends on this rule), so run
  `--forcerun discover_tsa` after it changes. In phase 2 BUSCO runs in `transcriptome` mode for group `tsa`
  (MetaEuk; its `full_table.tsv` has no `Gene Start/End/Strand`, handled by `scripts/busco_table.py`);
  `extract_transcripts` writes the MetaEuk CDS (no introns to cut), and TSA is excluded from the intron survey.
  TSA proteins reach the per-gene FASTAs, and so the tree, when their taxonomy is `complete` or `partial`.
- **Dereplication** (`dereplicate_genomes`, used by the tree's `derep` mode and optionally the metapackage) runs
  `galah cluster` (skani) on each genome's **marker CDS** (rule `extract_marker_cds`: spliced CDS of the Complete and
  best Duplicated BUSCO genes), so genomes and transcriptomes compare on the same sequence. Marker-CDS ANI runs higher
  than whole-genome ANI: `ani: 97` reproduces the old coverm 95% clusters (ARI 0.986 on 1,303 genomes; 590 vs 600
  clusters). Transcriptomes lose `tsa_quality_penalty` completeness points in the quality report, so a genome is the
  representative whenever one shares the cluster.

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
- **Introns and SingleM windows.** SingleM windows (20 aa / 60 nt) are matched against
  metagenomic DNA reads, which cannot cross an intron, so a window that contains an intron in
  some genome cannot be recovered from that genome. Introns are read from the exon structure of
  the miniprot gene models BUSCO writes (`busco_sequences/*/<gene>.gff`; a gap between exons of
  at least `singlem.introns.min_intron` nt is an intron). The metapackage is built in up to two
  **intron variants** (wildcard `introns`, in the paths of packages, metapackages and transcripts):
  - `corrected` (always built): the window is the best-covered one that crosses an intron in at most
    `max_span_fraction` (5%) of the genomes covering it (windows within `score_tolerance` of the best
    score count as equally covered; the one with fewest introns wins), and the sequences of the
    nucleotide database (rule `extract_transcripts`) are cut at the introns: one record per
    intron-free block (`<gene>|exon<i>of<n>`), blocks shorter than `min_exon_nt` dropped.
  - `naive` (the negative control, **off by default**: `singlem.introns.build_naive: true` or
    `--config build_naive=true`): the best-covered window, introns ignored, and whole-gene spans
    (introns included) in the nucleotide database. It doubles the packaging and per-genome transcript
    work, which is why it is optional. The survey always computes both windows, so the notebook compares
    them whether or not the naive metapackage was built.
  - `intron_window` (per gene, shared by both variants) aligns every genome's protein to the BUSCO HMM, maps its
    introns onto the HMM columns, scans every 20-column window start and writes
    `outputs/singlem/introns/<gene>.window_position.{naive,corrected}.txt` plus per-intron, per-genome and
    per-window tables. Windows are numbered in HMM match columns, which is what `singlem create --hmm-position`
    and `singlem pipe` use. Genes with a frameshift or a CDS that does not match the protein are not used for the
    intron statistics (still counted for coverage).
  - `window_intron_summary` collects the decisions (`outputs/metadata/reports/window_intron_summary.tsv`)
    and gates `create_draft_metapackage` of the corrected variant: with `introns.strict: true` the run stops
    if any gene has no corrected window under the limit.
  - **Stop symbols.** BUSCO proteins contain `*`. `graftM create --hmm` (unlike its mafft path) does not strip them;
    `hmmalign` writes them in upper case, GraftM only deletes lower-case insert columns, and the columns holding a `*`
    stay in the alignment as nearly empty columns (later rewritten to `X`). The alignment is then wider than the HMM and
    `singlem seqs` window numbers (alignment columns) no longer match the HMM columns that `create`/`pipe` use.
    `aggregate_busco_sequences.py` therefore strips `*`, and `graftm_create` fails if the alignment width differs from
    the HMM length. (`singlem seqs` is no longer used: windows come from the survey.)
  - **Quick-check subset.** `singlem.introns.subset` (or `--config subset=family5`; also `jgi_p10k`, `jgi_p10k_family5`)
    reduces the genomes of the existing BUSCO QC manifest to `subset_per_family` random genomes per family
    (`subset_seed`) and writes everything to `outputs/subset_<name>/`. Off-target decoys are skipped. The subset
    is read from the QC manifest on disk without rerunning BUSCO (`--forcerun make_subset_manifest` after it changes).
  - `intron_survey_notebook` writes `notebooks/rendered/intron_survey.ipynb`
    (`notebooks/rendered/subset_<name>/` for a subset run) with tables and PNGs in `tables/`. The notebook
    (`notebooks/intron_survey.ipynb`) is generated: edit `notebooks/build_intron_survey.py` and run
    `.pixi/envs/notebook/bin/python notebooks/build_intron_survey.py`. Besides the summary figures it draws the
    alignment of the most covered, cleanest and worst genes as one row per genome (coloured by phylum) with the
    introns and the naive (red) and corrected (magenta) windows. Executed notebooks live in
    `<workflow>/notebooks/rendered/` (gitignored), next to their templates.
  - **Snapshots.** `scripts/snapshot_intron_survey.py` freezes the survey so it can be repeated after the
    metapackage changes (use the `notebook` pixi env):
    `python scripts/snapshot_intron_survey.py create --name <name> [--subset family5]` copies the per-gene tables (gzipped),
    the windows of both variants, the taxonomy, protein names and BUSCO manifest, the config, the code and the notebook into
    `snapshots/intron_survey/<name>/`, with `inputs.tsv` (md5 of every gene FASTA/HMM/manifest), `RUN_INFO.json` (git commit,
    tool versions, settings) and `MANIFEST.tsv`. `compare <old> <new>` reports which inputs changed and how windows and span
    fractions moved per gene (snapshots from before the naive/corrected rename still compare); `render <name>` re-runs the
    notebook on a snapshot alone. Snapshots are not tracked by git.
  - To inspect one gene by hand: `.pixi/envs/singlem/bin/python scripts/intron_window.py --help`.
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
`notebooks/rendered/`, covering every `{derep}/{copies}` variant side by side, mostly as figures; the
numbers and PNG copies of each figure are in `notebooks/rendered/tables/<notebook>/`):

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

PhyloRank's RED is a median over rootings, which can make some branch lengths negative and let a
child sit below its parent. `phylorank_red` therefore post-processes the RED trees with
`scripts/monotone_red_tree.py` (RED made non-decreasing from root to tip, branch lengths recomputed, raw trees kept as
`*.phylorank_raw.tree`); `assign_taxonomy` keeps its own guard and warns if it ever has to raise a node.

**Source groups.** `tree_metadata.tsv` has a `source_group` column (fungi, algae, p10k, ncbi, tsa, external; `unknown` if a genome
is in no BUSCO QC manifest), and `tree_accuracy` / `tree_diagnostics` break misplaced tips and genes per genome down by it, so
transcriptome (`tsa`) tips can be judged separately from genome tips (`scripts/source_groups.py`).

**Unresolved genomes as targets.** With `external.include_unresolved: true` the metapackage genomes whose `tax_status` is
`unresolved` (no usable lineage, e.g. 263 p10k "Eukaryota sp.") are added to the external genomes (`helpers.unresolved_genome_paths`,
read from `genomes_to_process.paths` and `taxonomy_mapping.tsv`; transcriptomes are skipped because the external route runs BUSCO in
genome mode). They are BUSCO'd again here, aggregated with `allow_missing_taxonomy` and enter the tree as targets named from their
neighbours; the metapackage still drops their sequences. When checking the DAG with `-n`, include the `input` rerun trigger
(`--rerun-triggers input mtime`): with `mtime` alone Snakemake does not cascade the new BUSCO jobs to their consumers.

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
