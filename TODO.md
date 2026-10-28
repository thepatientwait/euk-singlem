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

5. **TSA transcriptomes** (not started; asked Oct 26, plan below).
   - `scripts/discover_tsa.py` (checkpoint): non-Metazoa/Viridiplantae `tsa-master` records from `nuccore` (1,162 found in the sandbox notebook); option `tsa.only_novel_taxa` (default true, anti-join against `taxonomy_mapping`); one assembly per TaxID, size/contig filter. Code to reuse: `metapackage_workflow/sandbox.ipynb` cells 3-4.
   - Rule `download_tsa`: bulk `https://ftp.ncbi.nlm.nih.gov/genbank/tsa/<X>/tsa.<PREFIX>.<n>.fsa_nt.gz`, concatenate, md5-check; output `outputs/genomes/tsa/{acc}.fasta.gz`; rate-limit resource like `jgi_api`.
   - Add `tsa` as a source in `build_unified_manifest.py` so `fetch_ncbi_taxonomy.py` and the `tax_status` gate apply.
   - Add `tsa` to `BUSCO_GROUPS`; `run_busco` uses transcriptome mode (MetaEuk) for it. Test on `metapackage_workflow/test/Coolia_malayensis_transcripts.fasta` (gitignored, 55 MB) first.
   - `aggregate_busco_sequences.py` and `busco_results_quality_control.py` must cope with missing `Gene Start/End/Strand`; `extract_busco_transcripts.py` needs a transcript path (emit contig unchanged; no GFF).
   - Intron status for TSA taxa: unknown (spliced); do not let them influence window choice.
6. **Derep with galah on marker CDS** (asked Oct 26): new `build_marker_cds` rule (spliced CDS of Complete + chosen Duplicated BUSCOs per genome/transcriptome), then `galah cluster` replacing `dereplicate_genomes.py`. A galah 0.4.2/skani 0.2.2 install exists under `/mnt/weka/pkg/cmr/sw/conda/envs`. Whole-genome ANI is poor between transcriptome and genome.
7. **Annotation workflow** (to define our own markers): confirm OrthoDB file names/URLs and sizes (403 from sandbox; user must test access), then gene calling (metaeuk vs EukProt, or eggNOG v7 OGs), `family5`-sized pilot (~200 genomes) checking the 129 BUSCO genes map to consistent OGs, then scale. Evaluation of eggNOG v7 euk coverage was not finished (`eggnog.db` timed out).
8. Per-genome parsing speed-up for the 39k-genome intron survey (not done).
9. Test whether the window-coordinate bug (`HISTORY.md` section 8) affected recovery of existing packages (not tested).

## Decisions pending

10. Cap on NCBI genomes in the `all` tree.
11. Whether peer copy selection beats `single` on full data (otherwise drop duplicated copies).
12. Polyphyletic clade splitting: shelved; revisit only after labels/contamination are fixed.

## Housekeeping

13. Regenerate or retire `taxonomy_mapping.tsv`-dependent tree outputs (tree side still on the old schema until rerun).
14. Simplify the regex re-parse of assign output strings in `tree_workflow/scripts/generate_tree_metadata.py`.
15. Explain the 656 vs 655 job count (tree dry run) and 640 vs 748 BUSCO result directories (unexplained).
16. File the GraftM `*` bug upstream (graftm create --hmm keeps `*` columns; draft in the "Consolidate sessions" opencode session).
17. `benchmarking_workflow`: commit `plan.md` in its nested repo; it needs `final.smpkg`.
18. Update `metapackage_workflow/OUTPUT_LAYOUT.md` for the `{introns}` paths and `exon_blocks/`.
19. `metapackage_workflow/sandbox.ipynb` and `tree_workflow/sandbox.ipynb` are tracked scratch notebooks; delete the metapackage one once TSA discovery is ported.
20. `outputs/orchestration/phase1_download.done` was refreshed on Oct 28; delete it whenever phase 1 must rerun.
