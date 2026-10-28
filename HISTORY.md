# HISTORY: what was done, in order

Append-only log. Current state lives in `MEMORY.md`, open work in `TODO.md`, code layout in `ARCHITECTURE.md`.
Sessions are stored in `~/.local/share/opencode/opencode.db` (SQLite, read-only queries work), but this file is
the cheap way to recover context. Add a section at the end of each working session.

| Session | Dates | Topic |
|---|---|---|
| Greeting and conversation start | Oct 19-20 | `graftm_create` crash, taxonomy leaks, non-reference NCBI genomes |
| Singlem not found in PATH | Oct 20 | Broken conda env (stopped early) |
| PhyloRank sample population error | Oct 20-21 | Tree rooting, `assign_taxonomy` rewrite, taxonomy query consolidation |
| Genomes missing taxonomy mapping | Oct 22 | Download DAG, pixi conversion, per-rule logs, JGI tape restore |
| Debugging download failures | Oct 23 | Download hardening, tree diagnostics, 4-tree design, subset mode, notebooks |
| Consolidate sessions into one | Oct 26-27 | Introns and SingleM windows, `*` fix, `{introns}` variants, snapshots |
| Previous session todo list | Oct 27 | Decisions, RED/MIXED/copy-selection work, annotation workflow |
| Context check | Oct 28 | Recovered context, refreshed these docs, commits |

---

## 1. Taxonomy and manifest fixes (Oct 19)

**Trigger:** `graftm_create` failed with "Unable to find sequence ... in the taxonomy definition" for 48 `GCA_*` sequences.

**Root cause:** 137 p10k genomes had an empty `tax_id`. `write_genome_manifest` had its taxonomy filter commented out, so they entered BUSCO and gene aggregation. `aggregate_busco_sequences.py` wrote all sequences to the FASTA but only taxonomy-resolved ones to the TSV.

**Why the tax_ids were missing (confirmed against the NCBI API):**
- 126 had strain-level tax IDs that the old Entrez lookup never matched to the manifest's accession version.
- 11 were replaced or suppressed assemblies. Datasets drops them, and Entrez reports the paired `GCF_` accession.

**Fixes:**
- Restored the taxonomy gate in `write_genome_manifest`.
- `aggregate_busco_sequences.py` now writes the FASTA and TSV from the same resolved set.
- New `scripts/ncbi_taxid_resolver.py`: Datasets REST in batches of 20, with an Entrez fallback and GCA/GCF matching. Wired into `build_unified_manifest.py`.
- Patched the live data: all 137 recovered. `taxonomy_mapping.tsv` went from 10,298 to 10,435 rows.

**Environment notes:**
- The user's `euk-singlem` and `tree` conda envs live under `/work/microbiome/users/n9985158/`, which is not mounted in the sandbox.
- NCBI network access needs an explicit CA bundle (`/etc/pki/tls/certs/ca-bundle.crt`).

## 2. NCBI scope change (Oct 19-20)

- **Counts for Eukaryota (2759):** 74,299 total (25,515 reference, 48,784 non-reference).
- **With the Metazoa/Viridiplantae exclusion:** 30,328 total (6,670 reference, 23,658 non-reference).
- **After the p10k/JGI anti-joins:** 9,944 net-new non-reference genomes.
- **Stock-take of `outputs/genomes/ncbi/` (24,658 files):**
  - 17,762 were off-target Metazoa/Viridiplantae references.
  - 6,588 were on-target references.
  - 13 were non-reference.
  - 295 were replaced or unknown.
- **Implemented (user's two requests):**
  - `discover_ncbi_genomes.py` now runs two queries. References come without subclade exclusions (`include_metazoa_viridiplantae_refs: true`), and non-references are always subclade-filtered.
  - `reference_only: false`.
  - `download_genome.py` skips genomes already on disk.
  - Removed `ncbi_genomes.tsv` as a download-rule input, so regenerating discovery doesn't force re-downloads.
- **Result:** 32,576 assemblies planned, 66% already cached, 11,077 to download (~2.5 TB).

## 3. Tree rooting and PhyloRank crash (Oct 20)

- **Crash:** `phylorank outliers` raised `ValueError: Sample larger than population`.
- **Cause (corrected mid-session):** the Asgard-outgroup reroot (`outgroup.enabled: true`) produced a 3299/440 root split. `p__Ascomycota` landed on the root label, so every leaf looked like Ascomycota. The first diagnosis blamed midpoint rooting, which was wrong. The midpoint fallback never ran.
- **Actions:**
  - Set `outgroup.enabled: false` and kept midpoint rooting.
  - `generate_alignment` now drops the Asgard accessions when the outgroup is off.
  - `helpers.py` no longer swallows `IncompleteCheckpointException`. That swallowing had produced an empty gene list and the polars "expected at least 1 input" error.
  - Added `use_multi_copy` params to the tree BUSCO QC rule.
  - Added an `allow_missing_taxonomy` flag to `aggregate_busco_sequences.py`. It defaults to false in the metapackage and is true in the tree workflow, because 636 external genomes have no taxonomy.
- **Re-running:** use `--rerun-triggers mtime` so rule-code and env changes don't rerun 640 BUSCO jobs. `--rerun-incomplete` clears the stale marker on `concatinated_alignment.fasta`.

## 4. `assign_taxonomy_v5` rewrite and taxonomy consolidation (Oct 20-21)

**Six bugs found and fixed:**
- Singletons not counted as clades.
- Conflicts not propagating up and down.
- Lower ranks leaking from the parent.
- Blank-taxonomy genomes classed as reference.
- No truncation below the phylum RED cutoff.
- Exact purity causing `Name_N` suffixes.

**Result of the rewrite:**
- Clade splitting is shelved.
- Database labels are kept below the phylum cutoff.
- The tree logic is separate from file I/O.
- 341 s became 0.9 s.
- 19 unit tests.
- The rule now writes headers, and `generate_tree_metadata.py` reads both formats.
- On the real tree, 0 reference labels changed, 27 genomes moved from reference to target, and distinct phyla went from 46 to 68 because of new `Eukarya_pN` placeholders.

**Taxonomy query consolidated:**
- `metapackage_workflow/scripts/fetch_ncbi_taxonomy.py` is now the single query and cleaning path.
- Outputs: explicit rank columns, `tax_status` (complete, partial or unresolved) and `parsed_taxon_string`.
- Domain is always `Eukarya`, and `k__` is dropped.
- Same-named different species get `_1`, `_2`.
- `prepare_taxonomy` now just splits the mapping by `tax_status`.
- Deleted dead scripts: old fetch v1, `assign_taxonomy_v4`, `summarize_eukcc_results`, `aggregate_gene_sequences`, `create_gene_taxonomy_files`.
- Real-data result: 10,435 rows (9,056 complete, 1,116 partial, 263 unresolved).
- The 263 unresolved are p10k "Eukaryota sp." / "uncultured eukaryote" genomes. They are dropped by the manifest gate (undecided whether to give them a placeholder).

## 5. Download DAG, pixi, logging, orchestrator (Oct 22)

- **Missing genomes:**
  - Many genomes were missing: 11,077 NCBI, ~600 JGI, and 7 zero-byte p10k files.
  - A `--touch` on Oct 20 had marked everything done, so the dry run said "nothing to do".
  - The DAG now schedules only genomes that are not on disk (`build_all_genomes_to_download`).
  - Zero-byte markers are left in place. `retry_failed_downloads: false` is the opt-in to retry them.
  - `download_jgi_genome.py` got a skip-if-exists guard.
- **Pixi conversion (commits `7274535`, `f8458be`, `8c7384e`, `bdf95ec`):**
  - Both workflows now follow the `env_predictors` layout: a default env plus per-rule envs (`singlem`, `graftm`, `busco`, `coverm`, and tree `busco`/`eukcc`/`taxonomy`/`align`).
  - Rules use `conda: env("<name>")`.
  - Old `envs/*.yml` were deleted.
  - Pins to keep: tree `busco>=5` (the solver otherwise picked `4.0.beta1`).
- **Orchestrator (`fc62f8c`, `f65bfb2`):**
  - Fixed config paths, secrets loading and `sys.executable -m snakemake`.
  - `--config dryrun=true` dry-runs both phases (download 11,691 jobs, generate 42,188).
- **Per-rule logs (`4e6ca7a`):**
  - Every rule writes `logs/<rule>/<wildcards>.log`, overwritten on each run.
  - `configure_logger()` dup2s fds 1 and 2.
  - `shell:` rules use `exec > {log} 2>&1`.
  - `run:` rules use `redirect_to_log`.
- **JGI tape-restore redesign (`55a70fa`):**
  - `jgi_request_restore` sends batched restore requests.
  - Batches that fail with a 4xx are split to isolate bad entries. `Aspparac1` can't be restored (404, its file is named `Asppha1_...`).
  - Fail-fast on 401/403.
  - 586 fungi and 15 algae genomes were requested.
- **Cleanup:** 859 queued `metapackage_wor*` jobs from an earlier run were `qdel`'d.

## 6. Download hardening (Oct 23)

**Failure analysis of the first big download:**
- **JGI 409 (261):** files in `RESTORE_REGISTERED` were not excluded from the DAG.
- **JGI 500/timeouts:** ~600 jobs hit the API at once, with no concurrency cap and no jitter.
- **Stale metadata:** the JGI JSONs never refreshed (the API said RESTORED for 295 genomes the JSON listed as registered).
- **NCBI:** 3 huge genomes (12/30/11 GB) hit the 4 h walltime while gzipping.
- **Token:** a later failure was a double `Bearer Bearer` prefix in `secrets.yaml`.

**Fixes (commits `9d829e4`, `80cc344`, `d4b3393`, `54dc9e7`, `468a6f0`, `233a976`):**
- NCBI downloads fetch the pre-gzipped `*_genomic.fna.gz` from FTP with an md5 check. They fall back to `datasets` plus pigz, and CPUs double on retry.
- JGI: `jgi_api` resource cap (10), jitter and `Retry-After`, a 120 s file-list timeout, and a token preflight that rejects non-JWTs. A double `Bearer` prefix is now tolerated.
- Waiting jobs: genomes still on tape get jobs that poll every 15 min with no limit of their own. Runtimes are 2/4/6 h across three attempts. They re-read `secrets.yaml` before each authenticated call.
- `download_jgi_metadata` is a normal rule again. Force it with `--forcerun`, or delete the JSONs.
- Deleted 2,569 orphan genomes (134 GB): 2,557 NCBI, 7 fungi, 5 algae. The list is in `outputs/metadata/reports/deleted_orphan_genomes.tsv`.
- Removed 3 stale `.partial` files and the 7 empty p10k files (all available from the NCBI FTP).
- Memory fixes: capped polars/rayon threads for 190-CPU nodes (a JGI job was OOM-killed at 1 GB while parsing the Snakefile), and bumped JGI jobs to 2 GB. My 2 GB setting for the big NCBI genomes was a bug and is now 8 GB.
- Final state: NCBI 33,178, p10k 2,751, algae 250, fungi 3,631 on disk. 68 JGI genomes were still restoring and 1 (`Aspparac1`) was refused.
- Run with `--keep-going --resources jgi_api=10`. The orchestrator adds these, but a direct `snakemake -s subworkflows/download_data.smk` call does not.

## 7. Tree quality investigation (Oct 23)

**User complaint:** negative branches, Bigyra and others fragmented, `_N` suffixes on clade names, and a remembered better tree.

**Findings:**
- **Negative branches:** these come from PhyloRank's RED rescaling (median over rootings, then child minus parent). They are not from FastTree. Only 3 nodes in the old tree, worst -0.027. `assign_taxonomy` now clamps child RED to at least the parent's. The branch lengths are not fixed.
- **`_N` suffixes:** stale outputs from the old script (Oct 20 19:40). The current script does no splitting.
- **`fasttree -nt` on a protein alignment** was in the first commit. Switching to `-lg -gamma` (now the default) explains only ~10% of the fragmentation.
- **Fragmented clades:** 77 of 132 extra order-level pieces are single tips. Support is high (82-94% have SH-like support of at least 0.95). Labels, contamination and NCBI non-monophyly are the main suspects. Examples:
  - Mortierellomycota among Mucoromycota.
  - A Rhodophyta genome among Bacillariophyta.
  - `Ashace_FM2008_1` (Dothideomycetes) sitting in Saccharomycetes.
- **No indexing bug.** The duplicates were in the stale Sep-29 `merged_input` (copy numbers up to 23).
- **Multicopy code:** the intent is to keep the single top-scoring copy per BUSCO, which the code does. Real problems were fixed: the flag sat at the wrong config path and was never passed to the tree checkpoint, and QC counted rows rather than distinct genes.
- **Top-copy heuristic** matches the most peer-like paralog only 73-77% of the time.
- **Old "good" tree** couldn't be located in git or the logs. The only non-derep JGI+p10k artifact is `backup/final_metapackage.smpkg` (April, 4,910 genomes).

**Changes (commits `b4e9eb9`, `7e9818a`, `c048552`, `155c552`, `e78f085`, `94d4764`, `6aa47b2`, `d1195ef`, `43ee46e`):**
- **Four trees:** `{derep: all|derep} x {copies: single|multi}`. `single` uses Complete BUSCOs only. `multi` adds the top Duplicated copy. Both come from one aggregation through a `<gene>_status.tsv` sidecar.
- **Genome filters:** merge fails on duplicate IDs, stop codons are cleaned, and genomes with fewer than 20 genes are dropped.
- **Subset mode:** `--config subset=family5|jgi_p10k|jgi_p10k_family5` writes to `output_subset_<name>/`. `family5` is 5 random genomes per family (2,209) plus 50 external.
- **Notebooks:**
  - Two global notebooks, `tree_diagnostics.ipynb` and `tree_accuracy.ipynb`, in `output/notebooks/`.
  - Executed via nbconvert to `.ipynb`, mostly figures.
  - `tree_accuracy` includes the old comparison notebook.
  - A rooting-independent, sister-aware monophyly metric.
- **`publish_trees` rule:** copies the final trees, metadata and the two notebooks to a timestamped `~/trees/<date>/` with `MANIFEST.tsv` and `RUN_INFO.txt`.
- **`reroot_tree` bug:** it was a `run:` rule that ignored `conda:` and failed 16 times with `No module named dendropy`. It is now a `script:` rule.
- **`ancient()` on the two upstream inputs of `run_external_busco`:** code or env edits no longer trigger 640 BUSCO reruns.
- **`family5` smoke test (25 steps, all passed):**
  - `all` trees have 2,259 tips and `derep` trees have 936.
  - RF distance between single and multi is ~0.05, and between `all` and `derep` ~0.10.
  - Order-level extra pieces: ~155 for `all` vs ~60-65 for `derep`.
  - `build_tree` ran out of memory at 7 GB on the first attempt, so raise the base memory for the full run (suggested 24 GB).

---

## 8. Introns and SingleM windows (Oct 26-27; commits `54eac49`, `21cc22c`)

**Why:** SingleM windows (20 aa / 60 nt) are matched against metagenomic DNA reads, which cannot span an intron.
A window containing an intron in a genome can never be recovered from that genome.

**Findings (8,931 genomes, 129 genes, ~92% of gene models mappable):**
- Existing windows (from `singlem seqs`) cross an intron in >5% of covering genomes for 49 of 129 genes (13 are >25%; worst `344220at2759` 74%).
- Best-covered window ignoring introns: 51 of 129 over the limit. The window the workflow now picks: 1 of 129 (`248751at2759`, 14.8%, no cleaner window exists).
- Avoiding introns costs little: new windows score a mean 99.8% of the best coverage, never below 97.4%.
- Median marker gene has 2 introns; Basidiomycota and Chlorophyta are intron-rich, Ciliophora nearly intronless.
- **Coordinate bug:** `singlem seqs` ran on GraftM's reference alignment, which has extra insert columns (386 vs 373 HMM columns for one gene; 117 of 129 genes mismatched). `singlem create/pipe` number windows in HMM columns, so for 105 genes the stored window pointed to a different place than scored (up to 25 columns). Root cause: BUSCO proteins contain `*`, which `graftm create --hmm` keeps as nearly empty columns. Recovery effect on existing packages not tested.

**Implemented (`metapackage_workflow/`):**
- `intron_utils.py` (GFF exons -> introns -> HMM columns), `intron_window` rule (per gene), `window_intron_summary` (gates `create_draft_metapackage`; `introns.strict: true` stops on a gene with no clean window).
- `{introns}` wildcard: `corrected` always built; `naive` (negative control) only with `build_naive: true`. Paths of packages, metapackages and transcripts carry it; old `initial/final/regenerated` folders are orphaned.
- `extract_transcripts` rewritten: one record per intron-free exon block (`<gene>|exon<i>of<n>`), output in `outputs/transcripts/exon_blocks/`.
- `aggregate_busco_sequences.py` strips `*`; `graftm_create` fails if alignment width != HMM length.
- Subset mode (`--config subset=family5|jgi_p10k|jgi_p10k_family5`) -> `outputs/subset_<name>/`.
- Notebook `notebooks/intron_survey.ipynb` is generated by `notebooks/build_intron_survey.py`; rendered copies go to `notebooks/rendered/`.
- Snapshots: `scripts/snapshot_intron_survey.py create|compare|render`; baseline `baseline_2026-10-26` (47 MB, gitignored) and `subset_family5_2026-10-27`.
- Subset results (2,209 genomes): genes over the 5% limit 66 (naive) vs 1 (corrected); mean span fraction 12.7% vs 1.5%.
- Boolean `--config` overrides (`"false"` string) are now parsed correctly. Deleted `backup/final_metapackage.smpkg` (as asked) and the old `outputs/singlem/introns/`.
- Fixed `MissingInputException` in `extract_transcripts`: the stale QC manifest still listed 3 deleted genomes; `get_genomes_with_transcripts_factory` now intersects with `genomes_to_process.paths`.

**Transcriptomes (TSA): planned, not started.** Plan in TODO.md. Only discovery code is in `metapackage_workflow/sandbox.ipynb`.

## 9. Decisions and tree work (Oct 27; commits `419b826` ... `6e40118`)

User decisions: MIXED nodes take the majority name (1B); keep polyphyletic splitting shelved (2); replace top-score copy with peer-closeness (3B); unresolved p10k genomes are tree targets only (4); print a summary of undownloaded genomes on exit (5).

- `419b826` download summary at exit (`scripts/download_summary.py`, onsuccess/onerror in `download_data.smk`); manifest warns about absent genomes.
- `5967538` 263 unresolved p10k genomes pass the manifest gate (tree targets); metapackage still skips them via `tax_status`.
- `dcaece0` MIXED nodes above the rank's RED cutoff take the majority reference name; audit in `mixed_nodes.tsv`. On `family5` no real database label changed (32 placeholders replaced).
- `a522d8e`, `b521d17` peer-similarity copy selection for Duplicated BUSCOs (`peer_copy_selection.py`, `copy_selection_lib.py`, `multi_copy_selection: peer`). On `family5` 31% of choices switched away from top score.
- `48e9a66` negative RED branches: cause is PhyloRank taking medians over different rooting sets for parent and child. `monotone_red_tree.py` makes RED non-decreasing root-to-tip in place; raw trees kept as `*.phylorank_raw.tree`. Full tree: 3 negative branches, 8 nodes raised, max 0.040; no rank decision can change.
- `609c85e`, `06f4fc5` tree notebooks to `notebooks/rendered/`, with generated key findings, per-figure takeaways, thresholds cell and what-to-do-next table (`notebook_narrative.py`).
- `75802b9` safeguards: raised nodes must stay 0.05 below the lowest rank cutoff; warning if `enforce_monotone_red` fires; assertion that reference database names never change.
- `26797cd` `multi_copy_min_gain` (0.02): switch rate fell 31% -> 11%.
- `549f791` `annotation_workflow/`: OrthoDB download as Snakemake rules (file names unverified: OrthoDB returned 403 from the sandbox).
- `6e40118` `label_masking` rule: hides 15% of reference labels and re-assigns. `family5` all/single recovery: phylum 99.4%, class 95.1%, order 88.6%, family 74.0%, genus 64.3%; wrong-name rates 0.6-3.9%.

## 10. Context recovery and clean-up (Oct 28)

- Downloads finished: `DOWNLOAD SUMMARY` shows 39,244 expected, 7 outstanding (all p10k zero-byte markers, accepted). On disk: ncbi 33,178; p10k 2,751 (7 empty); algae 260; fungi 3,861.
- Committed the intron work (`54eac49`, `21cc22c`) and removed tracked PBS logs (`7043e78`). Unit tests: 92 pass in the metapackage default env.
- Deleted stale: `jgi_not_ready.*` / `jgi_restore_requested.*` reports, `tree_workflow/output_subset_family5/comparison/`. Refreshed the `phase1_download.done` flag (it was dated Mar 23 and would have made the orchestrator skip phase 1). `LOGGING_HANDOFF_env_predictors.md` no longer exists in this repo.
- Wrote `MEMORY.md`, `HISTORY.md`, `TODO.md`, `ARCHITECTURE.md`, `AGENTS.md` (this file replaces `SESSIONS_SUMMARY.md`).
