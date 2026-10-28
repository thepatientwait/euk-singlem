# ARCHITECTURE: where things are

Use this to find code without searching. Commands and config details are in `README.md`; state in `MEMORY.md`;
open work in `TODO.md`. Paths are relative to the repo root. Tests are `unittest`-style: run
`<workflow>/.pixi/envs/default/bin/python -m unittest discover -s tests -q` (metapackage default env has Biopython/polars). Tree tests need two envs: `default` (alignment, diagnostics, accuracy; no dendropy) and `taxonomy` (assign_taxonomy, label masking, monotone RED; no Biopython); each errors only on the other's modules.

## Data flow

```
JGI / NCBI / P10K  --download_data.smk-->  outputs/genomes/{fungi,algae,ncbi,p10k}/*.fasta.gz
        + taxonomy_mapping.tsv (tax_status)      + genomes_to_process.paths (manifest)
   --generate_metapackage.smk-->  BUSCO per genome -> QC -> (derep) -> per-gene FASTAs (aggregate_busco_sequences)
        -> graftm_create -> intron_window (per gene) -> singlem_create -> off-target screen -> draft metapackage
        -> extract_transcripts (exon blocks) -> nucleotide DB -> final metapackage
   per-gene FASTAs  --tree_workflow-->  alignment -> FastTree -> reroot -> PhyloRank RED -> monotone RED
        -> assign_taxonomy (RED-based names) -> tree_metadata -> notebooks -> publish_trees (~/trees/<date>/)
   genomes  --annotation_workflow (planned)-->  OrthoDB / own OGs to test an own marker set
```

## metapackage_workflow/

| Path | What |
|---|---|
| `snakefile` | Orchestrator: runs phase 1 then phase 2 as nested snakemake calls; done-flags in `outputs/orchestration/*.done` (delete a flag to force a phase) |
| `subworkflows/download_data.smk` | Phase 1. Rules: `download_jgi_metadata`, `jgi_request_restore`, `download_jgi_genome`, `download_p10k_metadata`, `discover_ncbi_genomes`, `build_unified_manifest`, `download_ncbi_genome`, `fetch_ncbi_taxonomy`, `cleanup_failed_genomes`, `write_genome_manifest`, plus `discover_tsa` (checkpoint) and `download_tsa` (NCBI TSA transcriptomes). `onsuccess/onerror` print the download summary |
| `subworkflows/generate_metapackage.smk` | Phase 2 (~25 rules): `run_busco`, `busco_results_quality_control`, `extract_marker_cds`, `dereplicate_genomes` (galah on marker CDS), `aggregate_busco_sequences`, `graftm_create`, `intron_window`, `window_intron_summary`, `intron_survey_notebook`, `singlem_create`, off-target rules, `create_draft_metapackage`, `extract_transcripts`, `eukaryotic_transcripts`, `aggregate_transcript_taxonomy`, `make_nucleotide_sdb`, `create_final_metapackage` |
| `config/config.yaml` | Behaviour (multicopy: `busco.quality_control.multi_copy_selection`/`multi_copy_min_gain`, `singlem.introns.*`, BUSCO, derep, download limits). `paths.yaml` = every output path. `secrets.yaml` = JGI token (gitignored). `dryrun.yaml` |
| `helpers.py` | Config/paths loading, DAG filters (skip PURGED, present genomes), `get_genomes_with_transcripts_factory`, bool parsing of `--config` overrides |
| `scripts/` | One script per rule. Key ones below |
| `notebooks/` | `build_intron_survey.py` generates `intron_survey.ipynb`; executed copies in `notebooks/rendered/` (gitignored) |
| `snapshots/intron_survey/` | Frozen survey results (gitignored); `scripts/snapshot_intron_survey.py create/compare/render` |
| `tests/` | downloads, taxonomy, multicopy, peer copy selection, introns, transcripts, snapshot, tsa (discovery, download filter, transcriptome BUSCO tables, marker CDS, quality penalty) |
| `OUTPUT_LAYOUT.md` | Output directory map (partly pre-`{introns}`) |
| `sandbox.ipynb` | Scratch; its TSA discovery cells were ported to `scripts/discover_tsa.py`. `nuccore_tsa_metadata.tsv` = old output |

Key scripts: `fetch_ncbi_taxonomy.py` (single taxonomy path, `tax_status`), `build_unified_manifest.py`,
`download_genome.py` / `download_jgi_genome.py` / `jgi_request_restore.py` / `jgi_common.py`,
`download_summary.py` (classifies undownloaded genomes), `run_busco.py`, `aggregate_busco_sequences.py`
(strips `*`, top-scoring or peer copy, `<gene>_status.tsv` sidecar), `peer_copy_selection.py`,
`intron_utils.py` + `intron_window.py` + `window_intron_summary.py` (intron-aware windows),
`extract_busco_transcripts.py` (exon blocks), `subset_busco_manifest.py` (subset mode), `run_notebook.py`,
`discover_tsa.py` + `download_tsa.py` (TSA transcriptomes), `busco_table.py` (genome- and transcriptome-mode full tables with one column set), `extract_marker_cds.py` (spliced marker CDS per genome), `dereplicate_genomes.py` (galah/skani on marker CDS), `build_checkm2_quality_report.py` (BUSCO completeness, TSA penalty), `log_config.py` (per-rule logging).

Output layout (under `outputs/`, gitignored): `genomes/<group>/` (groups fungi, algae, p10k, ncbi, tsa), `derep/marker_cds/<group>/`, `metadata/source/tsa/tsa_projects.tsv`, `metadata/{source,manifests,taxonomy,qc,reports}/`,
`busco/<group>/<acc>/`, `singlem/{input,genome_input,packages,metapackages,introns,off_target}/`,
`transcripts/<introns>/{extracted,...}`, `subset_<name>/` (subset runs), `orchestration/`.
The `{introns}` wildcard (`corrected` | `naive`) is in package, metapackage and transcript paths.

## tree_workflow/

| Path | What |
|---|---|
| `snakefile` | All rules: external genomes BUSCO, `aggregate_*`, `merge_gene_fastas`, `famsa`, `bmge`, `generate_alignment`, `build_tree`, `reroot_tree`, `phylorank_decorate`, `phylorank_red`, `prepare_taxonomy`, `assign_taxonomy`, `generate_metadata`, `label_masking`, notebook rules, `publish_trees` |
| `config.yaml`, `paths.yaml` | Tree settings (`tree.derep_modes`, `copy_modes`, `fasttree_args`, `outgroup`, `subset`); paths incl. notebook dir |
| `scripts/assign_taxonomy_v5.py` | RED-based naming; pure tree logic separate from I/O. MIXED-node majority naming, name-change assertion |
| `scripts/monotone_red_tree.py` | Makes RED non-decreasing root to tip, checks cutoff clearance; raw trees kept as `*.phylorank_raw.tree` |
| `scripts/label_masking.py` | Hide 15% of reference labels, re-assign, report recovery |
| `scripts/copy_selection_lib.py`, `tree_accuracy_lib.py`, `tree_diagnostics_lib.py`, `tree_fragmentation.py` | Notebook computations |
| `scripts/notebook_narrative.py`, `run_notebook.py` | Generated key findings/verdicts; nbconvert runner (moves tagged `findings-top` cells up) |
| `scripts/publish_trees.py` | Copies finals to `~/trees/<date>/` with MANIFEST/RUN_INFO |
| `notebooks/` | `tree_diagnostics`, `tree_accuracy`, `copy_selection` templates; `rendered/` executed copies |
| `output/{derep}/{copies}/`, `output_subset_<name>/` | Results; four trees = `{all,derep}` x `{single,multi}` |
| `tests/` | assign_taxonomy, monotone RED, label masking, copy selection, alignment, diagnostics, accuracy, publish, subset, narrative |
| `quick_dist.py`, `test_trimmed_alignment.faa`, `*_example.*` | Misc/example inputs |

## annotation_workflow/ (new, skeleton)

`snakefile` (rules `fetch_orthodb_file`/`link_orthodb_file`, `orthodb_subset`, `orthodb_fasta`, `orthodb_diamond_index`), `config/config.yaml` (OrthoDB version/URLs/files; unverified),
`scripts/fetch_file.py` (resumable atomic download), `orthodb_tables.py`, `download_summary.py`, `tests/`.
Not yet: gene calling, pilot, marker-set analysis.

## benchmarking_workflow/

Separate nested git repo (gitignored here). Early stage; `plan.md` holds the expansion plan. Needs `final.smpkg`.

## Conventions

- Environments: pixi per workflow (`pixi.toml`), rules use `conda: env("<name>")` -> `.pixi/envs/<name>`.
- Logs: every rule writes `logs/<rule>/<wildcards>.log`.
- Re-running: use `--rerun-triggers mtime` (or the orchestrator's `DOWNLOAD_FLAGS`) so code edits do not rerun BUSCO.
- galah's tool check calls `which -- <tool>`; RHEL's exported `which` bash function rejects `--`. `unset -f which` (or use the pixi `galah` env on the compute nodes) when running galah by hand.
- Heavy runs go to the batch queue (`mqsub` / `snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua`); the sandbox cannot submit via `snakemake_mqsub`.
- Never rewrite reference database taxonomy names (asserted in `assign_taxonomy`).
- Docs to keep fresh: `MEMORY.md` (state/decisions), `HISTORY.md` (append per session), `TODO.md`, this file.
