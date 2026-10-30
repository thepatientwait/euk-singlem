# TODO

Last updated: 2026-10-28. Mark items done (move to `HISTORY.md` with a line) rather than leaving them here.

## Runs only the user can start (batch queue)

1. **Metapackage phase 2.** ~29,426 `run_busco` jobs (+ QC), then regenerate `dereplicated_genome_representatives.txt`
   (August list has 3,089 genomes). Downloads are finished, so this is unblocked.
   `cd metapackage_workflow && pixi run snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda -s subworkflows/generate_metapackage.smk`
2. **Rerun the full intron survey** after step 1 (`--forcerun intron_window`), create a snapshot and `compare` it with `baseline_2026-10-26`.
3. **Full tree run** after step 1: `cd tree_workflow && pixi run snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda --keep-going`.
   Raise `build_tree` memory to ~24 GB first. Decide a cap for NCBI genomes (the `all` tree may exceed 30k tips).
4. Check the new tree machinery on full data: peer copy selection (switch rate, `copy_selection` notebook), MIXED naming (`mixed_nodes.tsv`), label masking.

## Next to build

5. **TSA transcriptomes** (branch `transcripts`; code done, runs pending).
   - Run phase 1 (`download_data.smk`, no `--config` change needed; TSA is `tsa.enabled: true`): downloads 377 projects (~7 GB), regenerates the unified manifest, taxonomy and genome manifest. Then phase 2 runs BUSCO transcriptome mode on them.
   - Check after the first TSA BUSCO runs: completeness distribution, contamination/isoform redundancy (QC does not filter yet), and how many TSA taxa get `partial` or `unresolved` taxonomy.
   - TSA tips are flagged (`source_group`) and broken out in `tree_accuracy`/`tree_diagnostics` (done). Check the thresholds `TSA_EXCESS` and `TSA_GENES` against the first real run.
   - Intron status for TSA taxa is unknown by design (excluded from `intron_window`).
6. **Derep with galah on marker CDS** (done, branch `transcripts`; full run pending). Rerun the validation on the full set after phase 2 (compare with `genome_clusters.tsv`, tune `ani` 97 and `min_aligned_fraction` 15). `enabled: false` still means the metapackage is not built from representatives; the tree's `derep` mode uses them.
7. **Annotation workflow** (branch `annotation`; REMINDER: the user asked to finish the OrthoDB annotation change): confirm OrthoDB file names/URLs and sizes (403 from sandbox; user must test access), then gene calling (metaeuk vs EukProt, or eggNOG v7 OGs), `family5`-sized pilot (~200 genomes) checking the 129 BUSCO genes map to consistent OGs, then scale. Evaluation of eggNOG v7 euk coverage was not finished (`eggnog.db` timed out).
8. Per-genome parsing speed-up for the 39k-genome intron survey (not done).
9. Test whether the window-coordinate bug (`HISTORY.md` section 8) affected recovery of existing packages (not tested).

## Decisions pending

10. Cap on NCBI genomes in the `all` tree.
11. Whether peer copy selection beats `single` on full data (otherwise drop duplicated copies).
12. Polyphyletic clade splitting: shelved; revisit only after labels/contamination are fixed.

## Housekeeping

D-list reminder (ask the user to pick from items 8, 9, 13-19 once the A/B/C work below is finished).

0. Merge `transcripts` and `annotation` branches into `feature/v2-ncbi-hybrid` when done.

13. Regenerate or retire `taxonomy_mapping.tsv`-dependent tree outputs (tree side still on the old schema until rerun).
14. Simplify the regex re-parse of assign output strings in `tree_workflow/scripts/generate_tree_metadata.py`.
15. Explain the 656 vs 655 job count (tree dry run) and 640 vs 748 BUSCO result directories (unexplained).
16. File the GraftM `*` bug upstream (graftm create --hmm keeps `*` columns; draft in the "Consolidate sessions" opencode session).
17. `benchmarking_workflow`: commit `plan.md` in its nested repo; it needs `final.smpkg`.
18. Update `metapackage_workflow/OUTPUT_LAYOUT.md` for the `{introns}` paths and `exon_blocks/`.
19. `metapackage_workflow/sandbox.ipynb` and `tree_workflow/sandbox.ipynb` are tracked scratch notebooks; delete the metapackage one once TSA discovery is ported.
20. `outputs/orchestration/phase1_download.done` was refreshed on Oct 28; delete it whenever phase 1 must rerun.
