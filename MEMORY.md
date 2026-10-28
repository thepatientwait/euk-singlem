# MEMORY: current state, decisions and gotchas

Last updated: 2026-10-28 (branches: `transcripts` = TSA + galah derep, `annotation` = OrthoDB/own markers; base `feature/v2-ncbi-hybrid`). Update the status table and decisions when they change; put narrative in `HISTORY.md`,
open work in `TODO.md`, code locations in `ARCHITECTURE.md`.

## Project

Eukaryotic SingleM metapackage (`metapackage_workflow`), a genome tree with RED-based taxonomy (`tree_workflow`),
an annotation/marker-set exploration (`annotation_workflow`, skeleton) and a benchmark (`benchmarking_workflow`, early).
Branch `feature/v2-ncbi-hybrid`. Remote `github.com/thepatientwait/euk-singlem`. Run on PBS (aqua); the sandbox cannot
submit via `snakemake_mqsub`, the user starts big runs.

## Status

| Stage | Status |
|---|---|
| Metadata and taxonomy (39,251 genomes) | Done (`taxonomy_mapping.tsv` in new schema; 9,056 complete, 1,116 partial, 263 unresolved) |
| Genome downloads | **Done** (Oct 28). 39,244 expected; 7 p10k zero-byte markers accepted as lost. On disk: ncbi 33,178, p10k 2,751, algae 260, fungi 3,861 |
| BUSCO/QC on the new manifest | ~25% (8,931 of ~39k). Dry run of phase 2 plans 29,426 `run_busco` jobs |
| Intron-aware windows, exon-block transcripts, `*` fix | Implemented and committed; verified on `family5` subset only. Full survey needs phase 2 |
| Per-gene packages (129 genes) | Must be regenerated (new `{introns}` paths, `*` fix) |
| Final metapackage | Not built |
| TSA transcriptomes | Code written on branch `transcripts`: discovery (377 novel projects, ~7 GB), download, manifest, BUSCO transcriptome mode, aggregation, CDS extraction, intron-survey exclusion. Coolia BUSCO trial OK. Nothing downloaded into `outputs/genomes/tsa/` yet; BUSCO on TSA not run |
| Derep | galah on marker CDS implemented (`extract_marker_cds`, `dereplicate_genomes`); validated on 1,303 genomes. Not run on the full set (needs phase-2 BUSCO) |
| tree_workflow | RED monotone fix, MIXED majority naming, peer copy selection, label masking all done; tested on `family5` only. Full tree not run |
| annotation_workflow | OrthoDB download rules written; file names/columns unverified (403 from sandbox) |
| benchmarking_workflow | Early stage; needs `final.smpkg`; `plan.md` uncommitted in nested repo |

## Decisions in force

- Tree variants: `{derep: all|derep} x {copies: single|multi}`. `single` = Complete BUSCOs only; `multi` adds one Duplicated copy.
- Duplicated copy choice: peer similarity (`multi_copy_selection: peer`) with `min_gain` 0.02; falls back to top BUSCO score.
- MIXED nodes above their rank's RED cutoff take the majority reference name (listed in `mixed_nodes.tsv`). Polyphyletic clade splitting stays **shelved**.
- Reference database taxonomy names are never changed (asserted). Placeholders are `Eukarya_pN`.
- 263 unresolved p10k genomes: tree targets only; excluded from the metapackage by `tax_status`.
- Midpoint rooting (`outgroup.enabled: false`); Asgard outgroup caused the PhyloRank crash.
- Windows: best-covered window crossing an intron in <=5% of covering genomes (`max_span_fraction`); `naive` variant is a negative control, off by default.
- Genomic data only for intron analysis; TSA transcripts are spliced, so intron status for them comes from genome relatives or is unknown.
- Transcriptomes: only taxa no genome covers (genomes are trusted more); they go into the tree too (complete/partial taxonomy). BUSCO transcriptome mode (MetaEuk). Derep: galah on marker CDS, `ani: 97`, transcriptomes lose 50 completeness points so genomes win clusters.
- Token rotation is handled by the user.

## Gotchas

- `--config key=false` arrives as the string `"false"`; use the helper that parses booleans.
- Use `--rerun-triggers mtime` (or the orchestrator's flags) or code edits rerun hundreds of BUSCO jobs. `ancient()` guards `run_external_busco` inputs.
- `build_tree` ran out of memory at 7 GB on `family5`; raise to ~24 GB for the full run.
- The QC manifest `busco_quality_control_results.tsv` goes stale until BUSCO is rerun; `extract_transcripts` now intersects it with `genomes_to_process.paths`.
- Orchestrator skips a phase if its `outputs/orchestration/phase*.done` flag exists. Delete the flag to force a rerun.
- Dry run of the full phase 2 takes ~2 minutes (exceeds the sandbox tool limit; run it in the background).
- Two BUSCO versions are mixed in results (6.0.0: 5,286 genomes; 6.1.0: 3,645).
- `.pixi/envs/notebook` lacks Biopython, so metapackage tests need the `default` env.
- BUSCO transcriptome mode has no `Gene Start/End/Strand` columns (Sequence = `contig:start-end`, strand only in FASTA headers); `scripts/busco_table.py` normalises both modes.
- Running a BUSCO job from the tool shell with `nohup &` dies when the tool call times out; use `setsid nohup ... < /dev/null &`.
- The tree reads `metapackage_workflow/outputs/singlem/input`, which only holds genomes with `complete`/`partial` taxonomy. The 263 `unresolved` p10k genomes therefore do NOT reach the tree (the Oct 27 decision 'make them tree targets' is only half done: they pass BUSCO but not aggregation).
- NCBI needs an explicit CA bundle (`/etc/pki/tls/certs/ca-bundle.crt`) in the sandbox. OrthoDB returns 403 there.
- The sandbox cannot see `/work/microbiome/users/n9985158/` (the user's conda envs).
- Sandbox `/tmp/pbs.*/opencode` is the writable scratch space; its job id changes per session.
- Old Sep-29 `merged_input` had copy numbers up to 23; stale outputs in `output*` are a common source of confusing results.

## Facts worth remembering

- 129 marker genes (`*at2759`, BUSCO `eukaryota_odb12`). Existing legacy windows are wrong-coordinate for 105 genes (see `HISTORY.md` section 8).
- `family5` subset = 5 random genomes per family (2,209) + 50 external; trees: `all` 2,259 tips, `derep` 936.
- Label masking on `family5` all/single: phylum 99.4%, class 95.1%, order 88.6%, family 74.0%, genus 64.3% recovered.
- Fragmented clades: 77 of 132 extra order-level pieces are single tips; suspects are label errors and contamination (e.g. `Umbisa1`, `Umbsp_AD052_1`, `Umbelo1` flagged in every variant).
- Intron-survey snapshots: `metapackage_workflow/snapshots/intron_survey/{baseline_2026-10-26,subset_family5_2026-10-27}`.
- Session history is queryable: `~/.local/share/opencode/opencode.db` (read-only SQLite; tables `session`, `message`, `part`, `todo`). The saved `todo` table is unreliable; use `TODO.md`.
