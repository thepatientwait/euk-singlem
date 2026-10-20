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

# phase 1: metadata + genome download (JGI / NCBI / P10K -> unified manifest)
snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda \
    -s subworkflows/download_data.smk

# phase 2: BUSCO/CheckM2 QC, SingleM packages, off-target screen, metapackage
snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda \
    -s subworkflows/generate_metapackage.smk
```

Notes:

- `config/config.yaml`, `config/paths.yaml`, `config/dryrun.yaml` are versioned;
  `config/secrets.yaml` (JGI API token) is **gitignored** — create it locally:
  ```yaml
  jgi_token: "Bearer <token>"
  ```
- `metapackage_workflow/snakefile` is a legacy top-level entrypoint that still
  points at the old root-level `config.yaml`/`paths.yaml`. Run the subworkflows
  above instead.
- The NCBI/P10K/JGI sources are consolidated by `scripts/build_unified_manifest.py`;
  `scripts/ncbi_taxid_resolver.py` and `scripts/fetch_ncbi_taxonomy_v2.py` resolve
  taxonomy (Datasets API with Entrez fallback for replaced/suppressed assemblies).

### tree_workflow

```bash
cd tree_workflow
snakemake --profile ~/git/hpc_scripts/snakemake_configs/aqua --use-conda -n
```

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
