# Output layout

The workflow keeps generated products under `outputs/` and separates metadata from analysis products. Paths are defined in
`config/paths.yaml`. `{introns}` is `corrected` (always) or `naive` (only with `build_naive`); `{group}` is a download source
(`fungi`, `algae`, `p10k`, `ncbi`, `tsa`).

## Metadata

- `outputs/metadata/source/`: downloaded source metadata: `jgi/`, `p10k/`, `ncbi/`, and `tsa/tsa_projects.tsv` (NCBI TSA projects chosen by `discover_tsa`).
- `outputs/metadata/manifests/`: `unified_download_manifest.tsv` (all sources, with tax IDs) and `genomes_to_process.paths` (genomes on disk with a taxonomy row).
- `outputs/metadata/taxonomy/`: `taxonomy_mapping.tsv` (ranks, `tax_status`, lineage string) and lineage records.
- `outputs/metadata/qc/`: BUSCO QC manifest (`busco_quality_control_results.tsv`), the CheckM2-style quality report used by dereplication, `genome_clusters.tsv`, `dereplicated_genome_representatives.txt`.
- `outputs/metadata/reports/`: cleanup and download reports (`download_summary.tsv`), `window_intron_summary.tsv`.

## Genomes, BUSCO, dereplication

- `outputs/genomes/{group}/{acc}.fasta.gz`: downloaded genomes; `tsa` holds transcriptomes (one file per project, named by the master accession, e.g. `GLPS00000000.1`).
- `outputs/busco/{group}/{acc}/`: BUSCO result, condensed after the run (`scripts/busco_run.py`): `run_<lineage>/full_table.tsv`,
  `proteins.faa.gz`, `models.gff.gz` (and `cds.fna.gz` for transcriptomes), `short_summary*`, and `logs/busco.log`. Results from before the
  condensing keep the original per-gene files until `scripts/condense_busco_runs.py --apply` is run; both layouts are read.
- `outputs/derep/marker_cds/{group}/{acc}.fasta.gz`: spliced marker-gene CDS of each genome, the sequence galah dereplicates on.

## SingleM / GraftM

- `outputs/singlem/input/`: per-gene protein FASTAs (+ `<gene>_status.tsv`, `multicopy_selection.tsv`); `outputs/singlem/genome_input/`: genome-centric inputs.
- `outputs/singlem/introns/`: per-gene intron survey tables (`<gene>.introns.tsv`, `.genomes.tsv`, `.window_scan.tsv`, `.window_choice.tsv`) and the window files `<gene>.window_position.{naive,corrected}.txt`.
- `outputs/singlem/packages/{introns}/{initial,final,regenerated}/<gene>.spkg`: GraftM/SingleM packages per variant.
- `outputs/singlem/off_target/`: off-target search and preparation products.
- `outputs/singlem/metapackages/{introns}/`: draft and final metapackages.

## Transcripts (nucleotide database)

- `outputs/transcripts/{introns}/extracted/{group}/{acc}.fasta`: `corrected` = one record per intron-free exon block (`<gene>|exon<i>of<n>`), `naive` = whole gene span; TSA = the MetaEuk CDS.
- `outputs/transcripts/{introns}/{singlem,taxonomy}/{group}/{acc}.otu_table.tsv`: per-genome SingleM hits and their taxonomy.
- `outputs/transcripts/{introns}/aggregate/`: combined OTU table and nucleotide database.

## Other

- `outputs/subset_<name>/`: everything above for a quick-check subset (`--config subset=family5|jgi_p10k|jgi_p10k_family5`).
- `outputs/orchestration/`: phase done-flags (delete one to rerun that phase). `outputs/archive/legacy/`: old layouts kept during migration.
- Outside `outputs/`: `notebooks/rendered/` (executed notebooks, gitignored; `intron_survey.ipynb` and `tables/`), `snapshots/intron_survey/` (frozen survey results, gitignored), `logs/<rule>/`, `busco_downloads/` (BUSCO reference cache).

The workflow does not require a persistent `ncbi_dataset/` directory; NCBI archive extraction is temporary.
