# Output layout

The workflow keeps generated products under `outputs/` and separates metadata from analysis products.

- `outputs/metadata/source/`: downloaded JGI, P10K, and NCBI source metadata.
- `outputs/metadata/manifests/`: genome and unified download manifests.
- `outputs/metadata/taxonomy/`: taxonomy mappings and lineage records.
- `outputs/metadata/qc/`: BUSCO, CheckM2, and dereplication records.
- `outputs/metadata/reports/`: cleanup, validation, and diagnostic reports.
- `outputs/genomes/{fungi,algae,p10k,ncbi}/`: downloaded genome FASTA files.
- `outputs/busco/{group}/`: BUSCO results by source group.
- `outputs/singlem/input/`: gene-centric SingleM inputs.
- `outputs/singlem/genome_input/`: genome-centric SingleM inputs.
- `outputs/singlem/packages/`: GraftM, SingleM, window, and regenerated package intermediates.
- `outputs/singlem/off_target/`: off-target search and preparation products.
- `outputs/singlem/metapackages/`: draft and final metapackages.
- `outputs/transcripts/{extracted,singlem,taxonomy}/`: per-genome transcript products.
- `outputs/transcripts/aggregate/`: combined transcript OTU table and nucleotide database.
- `outputs/archive/legacy/`: previous output layouts retained during migration.

`busco_downloads/` remains outside `outputs/` because it is a BUSCO reference cache. The workflow does not require a persistent `ncbi_dataset/` directory; NCBI archive extraction is temporary.
