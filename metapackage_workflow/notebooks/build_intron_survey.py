"""Build notebooks/intron_survey.ipynb (the template is generated; edit this file, then run it with the
`notebook` pixi env: .pixi/envs/notebook/bin/python notebooks/build_intron_survey.py)."""
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

cells = []
def md(t): cells.append(new_markdown_cell(t.strip("\n")))
def code(t): cells.append(new_code_cell(t.strip("\n")))

md("""
# Introns in the BUSCO marker genes and the SingleM windows

SingleM identifies a marker by one short **window** (20 amino acids, 60 nt) of its alignment. For a
*genomic* metagenome the reads come from DNA, and a read cannot cross an intron. If the window of a
gene contains an intron in some genome, reads of that genome can never match the window, so the
genome is invisible at that marker.

This notebook asks, using the BUSCO gene models (miniprot alignments) of every genome that went into
the metapackage:

1. **Presence** - how common are introns in these markers, and how does that differ between lineages?
2. **Conservation** - are intron positions shared across lineages (ancestral) or lineage-specific?
3. **Windows** - how often does the window cross an intron, for the **naive** window (best coverage, introns ignored) and for the
   **corrected** window the workflow chooses (`intron_window` rule)?
4. **Alternatives** - what would other windows of the same gene cost or gain?

**How an intron is called.** An intron is the genomic gap between consecutive exons of a gene model
(gap >= `min_intron`, default 30 nt; shorter gaps are frameshift/indel artefacts). It is placed on the
spliced coding sequence (`offset`, `phase`) and on the HMM alignment (`hmm_col`). A window spans an
intron if the intron lies strictly inside the nucleotides of the window.

**Caveats.** Gene models come from aligning BUSCO reference proteins to the genome with miniprot, so a
very small exon can be missed and divergent genes are modelled less well. Genes with a frameshift or in-frame
stop in the alignment, or whose CDS does not match the protein length, are *unmapped*: they are counted for
coverage but left out of every intron statistic. Prevalences below use the mapped genomes of the group as
denominator (not only those aligned at the column), so they are slight underestimates.
""")

code("""
%matplotlib inline
import gzip, json, os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors
import matplotlib.ticker
import numpy as np
import polars as pl

DEFAULTS = {
    "intron_dir": "outputs/singlem/introns",
    "taxonomy": "outputs/metadata/taxonomy/taxonomy_mapping.tsv",
    "max_span_fraction": 0.05,
    "build_naive": False,
    "subset": "none",
    "out_dir": "outputs/singlem/introns/notebooks/tables",
    "gene_names": "busco_downloads/lineages/eukaryota_odb12/links_to_ODB12.txt",   # BUSCO id <tab> description <tab> url
    "feature_genes": None,
    "feature_rows": 1200,
}
PARAMS = {**DEFAULTS, **(json.load(open(os.environ["NB_PARAMS"])) if "NB_PARAMS" in os.environ else {})}
INTRON_DIR = Path(PARAMS["intron_dir"])
OUT_DIR = Path(PARAMS["out_dir"])
OUT_DIR.mkdir(parents=True, exist_ok=True)
LIMIT = float(PARAMS["max_span_fraction"])
MIN_GENOMES_PER_CLADE = 20     # clades with fewer mapped genomes are not shown
TOP_PHYLA = 12                 # phyla shown in per-clade figures
plt.rcParams.update({"figure.dpi": 100, "axes.spines.top": False, "axes.spines.right": False})
# the same three windows have the same colour in every figure
# magenta and red are kept out of the phylum palette so a window is never confused with a phylum
WINDOW_COLORS = {"naive": "tab:red", "corrected": "magenta"}
WINDOW_LABELS = {"naive": "naive window (best coverage, ignores introns)", "corrected": "corrected window (avoids introns)"}
FEATURE_GENES = PARAMS.get("feature_genes")        # [(label, gene), ...]; default: most covered, cleanest, worst
FEATURE_ROWS = int(PARAMS.get("feature_rows", 1200))   # genomes drawn per alignment figure
SPANNING_ROW_SHARE = 0.5                           # up to this share of each colour group's rows are genomes whose window spans an intron
MIN_ROWS_PER_GROUP = 25                            # every colour group keeps at least this many rows in the subsample
N_PHYLUM_COLORS = 12                               # largest phyla get their own colour; the rest are grouped as "other phyla"


def show_all(df, rows=60):
    \"\"\"Print a polars DataFrame without truncating it.

    Parameters
    ----------
    df : pl.DataFrame
        The table to print.
    rows : int, optional
        Maximum number of rows to print.
    \"\"\"
    with pl.Config(tbl_rows=rows, tbl_cols=30, fmt_str_lengths=60, tbl_width_chars=220):
        print(df)


def save_fig(name):
    \"\"\"Save the current figure as a PNG next to the tables, then show it.\"\"\"
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"{name}.png", dpi=150)
    plt.show()


def save_table(df, name):
    \"\"\"Write a table as TSV next to the figures.\"\"\"
    df.write_csv(OUT_DIR / f"{name}.tsv", separator="\\t")
""")

md("""
## Load the per-gene tables

The `intron_window` rule writes four tables per BUSCO gene: `introns` (one row per intron), `genomes`
(one row per genome), `window_scan` (one row per window start) and `window_choice` (the decision).
Genomes are labelled with their phylum from the taxonomy mapping.
""")

code("""
INTRON_SCHEMA = {"genome_acc": pl.Utf8, "group": pl.Utf8, "busco_id": pl.Utf8, "contig": pl.Utf8, "strand": pl.Utf8,
                 "offset": pl.Int64, "phase": pl.Int64, "residue": pl.Int64, "hmm_col": pl.Int64,
                 "intron_len": pl.Int64, "n_introns": pl.Int64}


def read_all(suffix, **kwargs):
    \"\"\"Concatenate ``<gene>.<suffix>.tsv[.gz]`` over all genes, adding a ``gene`` column.\"\"\"
    frames = []
    for path in sorted(INTRON_DIR.glob(f"*.{suffix}.tsv*")):   # plain or gzipped (snapshots)
        frames.append(pl.read_csv(path, separator="\\t", **kwargs).with_columns(pl.lit(path.name.split(".")[0]).alias("gene")))
    return pl.concat(frames)


taxonomy = (pl.read_csv(PARAMS["taxonomy"], separator="\\t", infer_schema_length=0)
            .select(pl.col("acc").alias("genome_acc"),
                    pl.when(pl.col("phylum").is_null() | (pl.col("phylum") == "")).then(pl.lit("unclassified")).otherwise(pl.col("phylum")).alias("phylum"))
            .unique("genome_acc"))
GENE_NAMES = {}
if Path(PARAMS["gene_names"]).exists():
    opener = gzip.open if str(PARAMS["gene_names"]).endswith(".gz") else open
    with opener(PARAMS["gene_names"], "rt") as handle:
        for line in handle:
            parts = line.rstrip("\\n").split("\\t")
            if len(parts) >= 2:
                GENE_NAMES[parts[0]] = parts[1].strip()
else:
    print("WARNING: no gene name file, plots show BUSCO ids only:", PARAMS["gene_names"])


def gene_name(gene, width=None):
    \"\"\"Protein name of a BUSCO id (empty if unknown), optionally cut to ``width`` characters.\"\"\"
    name = GENE_NAMES.get(gene, "")
    return name if width is None or len(name) <= width else name[:width - 1] + "…"


def gene_label(gene, width=None, sep="\\n"):
    \"\"\"BUSCO id with its protein name underneath (``sep`` = newline) or after it (``sep`` = a space).\"\"\"
    name = gene_name(gene, width)
    return f"{gene}{sep}{name}" if name else gene


def with_name(frame, column="busco_id"):
    \"\"\"Add a ``name`` column (protein name) right after the BUSCO id column.\"\"\"
    named = frame.with_columns(pl.col(column).replace_strict(GENE_NAMES, default="").alias("name"))
    columns = [c for c in named.columns if c != "name"]
    columns.insert(columns.index(column) + 1, "name")
    return named.select(columns)


choices = with_name(read_all("window_choice"))
genomes = read_all("genomes").join(taxonomy, on="genome_acc", how="left").with_columns(pl.col("phylum").fill_null("unclassified"))
introns = (read_all("introns", schema_overrides=INTRON_SCHEMA)
           .join(taxonomy, on="genome_acc", how="left").with_columns(pl.col("phylum").fill_null("unclassified")))
scans = read_all("window_scan")

mapped = genomes.filter(pl.col("mapped"))
# The introns table also lists introns of gene models that failed the frameshift / length check; their positions on the
# alignment are not reliable, so every intron statistic and figure uses mapped genomes only.
n_introns_all = introns.height
introns = introns.join(mapped.select("gene", "genome_acc"), on=["gene", "genome_acc"], how="semi")
top_phyla = mapped.group_by("phylum").agg(pl.col("genome_acc").n_unique().alias("n")).sort(["n", "phylum"], descending=[True, False])["phylum"].head(TOP_PHYLA).to_list()
print(f"{choices.height} genes, {genomes['genome_acc'].n_unique()} genomes, {genomes.height} gene x genome records, "
      f"{introns.height} introns in mapped records ({n_introns_all - introns.height} more in unmapped records are left out)")
""")

md("""
## 1. How much of the data could be analysed?

Each genome x gene record is *mapped* when its gene model could be placed on the protein. Frameshifts and length mismatches
are excluded from intron statistics. A large unmapped share in one source would bias the results towards the other sources.
""")

code("""
status = (genomes.group_by(["group", "status"]).agg(pl.len().alias("n"))
          .with_columns((pl.col("n") / pl.col("n").sum().over("group")).alias("fraction")).sort(["group", "status"]))
show_all(status)
pivot = status.pivot(on="status", index="group", values="fraction").fill_null(0)
fig, ax = plt.subplots(figsize=(6, 3))
bottom = np.zeros(pivot.height)
for column in [c for c in pivot.columns if c != "group"]:
    ax.barh(pivot["group"], pivot[column], left=bottom, label=column)
    bottom += pivot[column].to_numpy()
ax.set_xlabel("fraction of gene x genome records"); ax.set_ylabel("genome source")
ax.set_title("Which gene models could be used for intron statistics?")
ax.legend(title="gene model status", frameon=False, loc="center left", bbox_to_anchor=(1.01, .5))
save_fig("mapping_status")
""")

md("""
## 2. Presence of introns

**Question.** How intron-rich are these markers, and does it depend on the lineage?

**How to read it.** Each point is a genome; its value is the mean number of introns per marker gene among its
mapped genes. Lineages near 0 are essentially intronless in these genes and need no window care; lineages with several
introns per gene are where windows are at risk.
""")

code("""
per_genome = (mapped.group_by(["genome_acc", "group", "phylum"])
              .agg(pl.len().alias("genes"), (pl.col("n_introns") > 0).mean().alias("frac_genes_with_intron"),
                   pl.col("n_introns").mean().alias("mean_introns"))
              .filter(pl.col("genes") >= min(20, max(1, choices.height // 4))).sort("genome_acc"))
save_table(per_genome, "per_genome_intron_summary")
order = (per_genome.filter(pl.col("phylum").is_in(top_phyla)).group_by("phylum")
         .agg(pl.col("mean_introns").median().alias("m")).sort(["m", "phylum"])["phylum"].to_list())
fig, ax = plt.subplots(figsize=(9.5, 4.8))
if not order:
    raise SystemExit("no phylum has enough genomes with mapped genes")
ax.boxplot([per_genome.filter(pl.col("phylum") == p)["mean_introns"].to_numpy() for p in order],
           tick_labels=[f"{p} (n={per_genome.filter(pl.col('phylum') == p).height})" for p in order],
           orientation="horizontal", showfliers=False, medianprops=dict(color="tab:orange", lw=2))
ax.set_xlabel("mean introns per marker gene, one value per genome")
ax.set_ylabel("phylum (n = genomes)")
ax.set_title("How intron-rich are the marker genes, by phylum?")
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc="white", ec="k", label="box: middle half of genomes"),
                   plt.Line2D([], [], color="tab:orange", lw=2, label="median genome"),
                   plt.Line2D([], [], color="k", lw=1, label="whiskers: full range (outliers hidden)")],
          frameon=False, loc="lower right")
save_fig("introns_per_gene_by_phylum")
phylum_table = (per_genome.group_by("phylum")
                .agg(pl.len().alias("genomes"), pl.col("mean_introns").median().alias("median_introns_per_gene"),
                     pl.col("frac_genes_with_intron").median().alias("median_frac_genes_with_intron"))
                .filter(pl.col("genomes") >= MIN_GENOMES_PER_CLADE).sort(["median_introns_per_gene", "phylum"], descending=[True, False]))
save_table(phylum_table, "intron_presence_by_phylum")
show_all(phylum_table)
""")

code("""
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
for group in sorted(per_genome["group"].unique()):
    axes[0].hist(per_genome.filter(pl.col("group") == group)["frac_genes_with_intron"].to_numpy(), bins=40, histtype="step", lw=1.5,
                 label=f"{group} (n={per_genome.filter(pl.col('group') == group).height})", density=True)
axes[0].set_xlabel("fraction of a genome's marker genes with at least one intron"); axes[0].set_ylabel("density of genomes")
axes[0].set_title("Intron-bearing genes per genome, by source"); axes[0].legend(title="genome source", frameon=False)
lengths = introns["intron_len"].to_numpy()
axes[1].hist(lengths, bins=np.logspace(np.log10(max(1, lengths.min())), np.log10(lengths.max()), 60), color="grey")
axes[1].set_xscale("log"); axes[1].axvline(np.median(lengths), color="k", ls="--", lw=1, label=f"median {int(np.median(lengths))} nt")
axes[1].axvline(50, color="tab:red", ls=":", lw=1, label=f"50 nt ({(lengths < 50).mean():.0%} of introns are shorter)")
axes[1].set_xlabel("intron length (nt, log scale)"); axes[1].set_ylabel("introns"); axes[1].set_title("Intron length")
axes[1].legend(frameon=False)
save_fig("intron_presence_and_length")
print("median intron length:", int(np.median(lengths)), "nt;", f"{(lengths < 50).mean():.1%} shorter than 50 nt")
""")

md("""
## 3. Conservation of intron positions

**Question.** Do introns sit at the same place across lineages (inherited from a common ancestor), or are they lineage-specific?

An intron *site* is a (marker gene, HMM column, phase). The prevalence of a site in a phylum is the fraction of that phylum's mapped
genomes with an intron there. Sites present in many phyla are ancient and are the ones that will recur no matter which genomes are
sampled; sites seen in one phylum only are specific to it.

**How to read it.** A histogram piled at low prevalence means most introns are rare. A long tail of sites with prevalence > 0.5 in several
phyla indicates ancestral introns. Windows that avoid these sites protect all lineages at once.
""")

code("""
sites = introns.filter(pl.col("hmm_col") >= 0)
print(f"{sites.height / introns.height:.1%} of introns fall on an HMM match column (the rest are in inserts)")

per_clade = (mapped.group_by(["gene", "phylum"]).agg(pl.len().alias("n_mapped")).filter(pl.col("n_mapped") >= MIN_GENOMES_PER_CLADE))
site_clade = (sites.group_by(["busco_id", "hmm_col", "phase", "phylum"]).agg(pl.len().alias("n"))
              .join(per_clade, left_on=["busco_id", "phylum"], right_on=["gene", "phylum"])
              .with_columns((pl.col("n") / pl.col("n_mapped")).alias("prevalence")))
site_summary = with_name(site_clade.group_by(["busco_id", "hmm_col", "phase"])
                .agg(pl.col("prevalence").max().alias("max_prevalence"),
                     (pl.col("prevalence") >= 0.25).sum().alias("phyla_ge_25pct"),
                     pl.col("phylum").n_unique().alias("phyla_seen")))
save_table(site_summary.sort(["phyla_ge_25pct", "max_prevalence", "busco_id", "hmm_col", "phase"], descending=[True, True, False, False, False]), "intron_sites")

fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
axes[0].hist(site_summary["max_prevalence"].to_numpy(), bins=40, color="grey"); axes[0].set_yscale("log")
axes[0].set_xlabel("highest prevalence of the site in any phylum"); axes[0].set_ylabel("intron sites (log scale)")
axes[0].set_title("Most intron sites are rare")
axes[1].hist(site_summary["phyla_ge_25pct"].to_numpy(), bins=np.arange(0, max(1, site_summary["phyla_ge_25pct"].max()) + 2) - 0.5, color="grey"); axes[1].set_yscale("log")
axes[1].set_xlabel("number of phyla in which the site is present in >= 25% of genomes"); axes[1].set_ylabel("intron sites (log scale)")
axes[1].set_title("How many phyla share a site")
conserved = (site_summary.filter(pl.col("phyla_ge_25pct") >= 3).group_by("busco_id").len().sort(["len", "busco_id"], descending=[True, False]).head(25))
axes[2].barh([gene_label(g, 28, sep="  ") for g in conserved["busco_id"][::-1]], conserved["len"][::-1], color="tab:red")
axes[2].tick_params(axis="y", labelsize=6)
axes[2].set_xlabel("conserved sites per marker\\n(>= 25% of genomes in >= 3 phyla)"); axes[2].set_ylabel("marker gene (top 25)")
axes[2].set_title("Markers with the most conserved sites\\n(BUSCO id, protein name)", fontsize=10)
save_fig("intron_site_conservation")
print(f"{site_summary.height} intron sites; {site_summary.filter(pl.col('phyla_ge_25pct') >= 3).height} are present (>= 25%) in 3 or more phyla")
print("Most widely shared intron sites (prevalence >= 25% in this many phyla):")
show_all(site_summary.sort(["phyla_ge_25pct", "max_prevalence", "busco_id", "hmm_col", "phase"], descending=[True, True, False, False, False]).head(25))
""")

md("""
## 4. Do the windows cross introns?

Two windows per marker are compared, always in **HMM match columns** (the numbering `singlem create --hmm-position` and `singlem pipe` use):

* **naive** - the window covered by the most aligned bases, introns ignored. This is the negative control: a metapackage built with it
  (`singlem.introns.build_naive`) shows what ignoring introns costs.
* **corrected** - the window the workflow uses: the best-covered window that crosses an intron in at most `max_span_fraction` of the genomes
  covering it (dashed line in the figures).

**Question.** For each marker, in what fraction of the genomes a window covers does it span an intron?

**How to read it.** Each panel is a histogram over the markers of one window type; the x-axis is the fraction of covering genomes
in which that window spans an intron. Markers to the right of the dashed limit line cannot be seen in more than the allowed share of
genomes. The naive panel has a long right tail; the corrected panel should sit almost entirely left of the line.
""")

code("""
order = choices.sort(["naive_span_fraction", "busco_id"], descending=[True, False])
columns = [("naive", "naive_span_fraction"), ("corrected", "corrected_span_fraction")]
x_max = max(float(choices[c].max()) for _, c in columns)
bins = np.linspace(0, x_max * 1.02, 34)
fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), sharex=True, sharey=True)
for ax, (key, column) in zip(axes, columns):
    values = choices[column].to_numpy()
    ax.hist(values, bins=bins, color=WINDOW_COLORS[key])
    ax.axvline(LIMIT, color="k", ls="--", lw=1.2, label=f"limit ({LIMIT:.0%})")
    over = int((values > LIMIT).sum())
    ax.set_title(f"{WINDOW_LABELS[key]}\\n{over} of {len(values)} genes over the limit", fontsize=10)
    ax.set_xlabel("fraction of covering genomes whose window spans an intron")
    ax.legend(frameon=False, loc="upper right")
axes[0].set_ylabel("marker genes")
save_fig("window_span_by_gene")

over_naive = choices.filter(pl.col("naive_span_fraction") > LIMIT)
moved = choices.filter(pl.col("corrected_start") != pl.col("naive_start"))
print(f"subset: {PARAMS['subset']}   naive metapackage built: {PARAMS['build_naive']}   limit: {LIMIT:.0%}")
print(f"naive windows over the limit:       {over_naive.height} of {choices.height}")
print(f"corrected window differs from naive: {moved.height}")
print(f"corrected windows still over the limit: {choices.filter(pl.col('violation')).height}")
show_all(order.select("busco_id", "name", "naive_start", "naive_span_fraction", "corrected_start", "corrected_span_fraction",
                      "corrected_score_ratio", "reason", "n_mapped", "n_genomes").head(20))
""")

md("""
### Cost of avoiding introns

`corrected_score_ratio` is the coverage score of the corrected window relative to the best-covered (naive) window of the gene (1 = as well covered as
the naive window). Avoiding an intron usually costs a little coverage; a large drop means the gene has no good intron-free region.
""")

code("""
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
changed = (choices["corrected_start"] != choices["naive_start"]).to_numpy()
for flag, color, label in ((~changed, "grey", "window unchanged"), (changed, WINDOW_COLORS["corrected"], "window moved off the naive position")):
    axes[0].scatter(choices["naive_span_fraction"].to_numpy()[flag], choices["corrected_span_fraction"].to_numpy()[flag], s=14, c=color, label=f"{label} ({int(flag.sum())})")
axes[0].plot([0, 1], [0, 1], "k:", lw=.8, label="no change"); axes[0].axhline(LIMIT, color="k", ls="--", lw=1, label=f"limit ({LIMIT:.0%})")
axes[0].set_xlabel("span fraction of the naive window"); axes[0].set_ylabel("span fraction of the corrected window")
axes[0].set_ylim(-0.01, max(3 * LIMIT, 1.3 * float(choices["corrected_span_fraction"].max())))
axes[0].set_title("What moving the window buys\\n(each dot is a marker gene)", fontsize=10); axes[0].legend(frameon=False, fontsize=8)
zoom = np.linspace(0.95, 1.0, 26)
values = choices["corrected_score_ratio"].to_numpy()
axes[1].hist(np.clip(values, zoom[0], 1.0), bins=zoom, color=WINDOW_COLORS["corrected"])
axes[1].set_xlim(zoom[0], 1.0)
axes[1].text(.03, .95, f"{int((values < zoom[0]).sum())} genes below {zoom[0]} (clipped into the first bar)\\nminimum {values.min():.3f}, median {np.median(values):.4f}",
             transform=axes[1].transAxes, va="top", fontsize=8)
axes[1].set_xlabel("coverage score of the corrected window / score of the naive window"); axes[1].set_ylabel("marker genes")
axes[1].set_title("Coverage of the corrected window\\n(1.0 = as well covered as possible)", fontsize=10)
save_fig("coverage_cost")
cost = pl.DataFrame({
    "window": [WINDOW_LABELS["naive"], WINDOW_LABELS["corrected"]],
    "genes over limit": [over_naive.height, choices.filter(pl.col("corrected_span_fraction") > LIMIT).height],
    "mean span fraction": [choices[c].mean() for c in ("naive_span_fraction", "corrected_span_fraction")],
    "median span fraction": [choices[c].median() for c in ("naive_span_fraction", "corrected_span_fraction")],
    "mean score / naive": [1.0, choices["corrected_score_ratio"].mean()],
    "min score / naive": [1.0, choices["corrected_score_ratio"].min()],
})
show_all(cost)
print("Genes still over the limit with the corrected window:")
show_all(choices.filter(pl.col("violation")).select("busco_id", "name", "corrected_start", "corrected_span_fraction", "corrected_score_ratio", "n_mapped"))
""")

md("""
### Which lineages lose out

For the most affected markers, the fraction of genomes of each phylum in which the **naive** window spans an intron. A marker
that is bad in one phylum only may be acceptable; one that is bad everywhere needs a different window.
""")

code("""
worst = choices.sort(["naive_span_fraction", "busco_id"], descending=[True, False])["busco_id"].head(10).to_list()
by_clade = (genomes.filter(pl.col("gene").is_in(worst) & pl.col("mapped") & pl.col("covers_naive") & pl.col("phylum").is_in(top_phyla))
            .group_by(["gene", "phylum"]).agg(pl.col("spans_naive").mean().alias("span"), pl.len().alias("n"))
            .filter(pl.col("n") >= MIN_GENOMES_PER_CLADE).sort(["gene", "phylum"]))
save_table(by_clade, "naive_window_span_by_phylum")
matrix = by_clade.pivot(on="phylum", index="gene", values="span", sort_columns=True).sort("gene")
fig, ax = plt.subplots(figsize=(11, 4.8))
values = matrix.drop("gene").to_numpy().astype(float)
image = ax.imshow(values, aspect="auto", cmap="Reds", vmin=0, vmax=max(0.3, np.nanmax(values) if values.size else 0.3))
ax.set_yticks(range(matrix.height), [gene_label(g, 30, sep="  ") for g in matrix["gene"]], fontsize=7); ax.set_xticks(range(values.shape[1]), matrix.columns[1:], rotation=60, ha="right", fontsize=7)
ax.set_xlabel("phylum"); ax.set_ylabel("marker gene (10 worst naive windows)")
ax.set_title("Where the naive window crosses introns, by phylum")
plt.colorbar(image, label="fraction of genomes whose\\nnaive window spans an intron")
save_fig("naive_window_span_by_phylum")
print("Fraction of genomes whose naive window spans an intron, worst 10 markers x phylum:")
show_all(with_name(matrix, "gene"), rows=15)

corrected_rate = (genomes.filter(pl.col("mapped") & pl.col("covers_corrected")).group_by("group")
                  .agg(pl.col("spans_corrected").mean().alias("corrected"), pl.len().alias("n")))
naive_rate = (genomes.filter(pl.col("mapped") & pl.col("covers_naive")).group_by("group")
              .agg(pl.col("spans_naive").mean().alias("naive")))
show_all(corrected_rate.join(naive_rate, on="group").select("group", "n", "naive", "corrected").sort("group"))
""")

md("""
## 5. Alternative windows

**Question.** Is a bad window avoidable? For the six markers whose naive window crosses introns most often, one panel per marker shows
what a window would look like at **every** possible start along the HMM.

**How to read a panel.**

* x-axis: the HMM column where a 20-column window starts.
* **Red curve, left axis:** the fraction of genomes covering that window whose window spans an intron. High = invisible to reads for many genomes.
* **Grey area, right axis:** the coverage of that window relative to the best window of the gene (1 = best covered).
* **Dashed line:** the intron limit. Windows whose red curve is below it are allowed.
* Markers sit on the red curve: the naive window (open red circle) and the corrected window (magenta square).

Where the red curve drops below the dashed line while the grey area stays high, a clean and well-covered window exists. If the red curve stays high
everywhere, the gene is intron-rich all along and the marker is a candidate to drop.
""")

code("""
show = choices.sort(["naive_span_fraction", "busco_id"], descending=[True, False])["busco_id"].head(6).to_list()
fig, axes = plt.subplots(2, 3, figsize=(16, 8))
for number, (ax, gene) in enumerate(zip(axes.ravel(), show)):
    s = scans.filter(pl.col("gene") == gene)
    ax2 = ax.twinx()
    ax2.fill_between(s["start"], s["aligned_bases"] / s["aligned_bases"].max(), color="lightgrey")
    ax2.set_ylim(0, 1.05)
    ax2.set_ylabel("relative coverage (grey)" if number % 3 == 2 else "")
    ax.plot(s["start"], s["span_fraction"], color="tab:red", lw=1.5)
    ax.axhline(LIMIT, color="k", ls="--", lw=.9)
    for flag, key, style in (("is_naive", "naive", dict(marker="o", mfc="none", mec="tab:red", mew=2, ms=10)),
                             ("is_corrected", "corrected", dict(marker="s", color=WINDOW_COLORS["corrected"], ms=7))):
        point = s.filter(pl.col(flag))
        ax.plot(point["start"], point["span_fraction"], ls="", **style)
    ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
    ax.set_title(gene_label(gene, 45), fontsize=10); ax.set_xlabel("window start (HMM column)")
    ax.set_ylabel("genomes spanning an intron (red)" if number % 3 == 0 else "")
fig.legend(handles=[plt.Line2D([], [], color="tab:red", label="fraction of genomes whose window spans an intron (left axis)"),
                    plt.Rectangle((0, 0), 1, 1, fc="lightgrey", label="relative coverage of the window (right axis)"),
                    plt.Line2D([], [], color="k", ls="--", label=f"intron limit ({LIMIT:.0%})"),
                    plt.Line2D([], [], ls="", marker="o", mfc="none", mec="tab:red", mew=2, label=WINDOW_LABELS["naive"]),
                    plt.Line2D([], [], ls="", marker="s", color=WINDOW_COLORS["corrected"], label=WINDOW_LABELS["corrected"])],
           loc="lower center", ncol=3, frameon=False)
fig.suptitle("Every possible window along the HMM, six worst markers", fontsize=12)
fig.tight_layout(rect=(0, 0.09, 1, 0.97))
plt.savefig(OUT_DIR / "window_scan_worst_genes.png", dpi=150)
plt.show()
""")

md("""
## 6. The alignment, genome by genome

**Question.** What do the windows and introns look like on the actual alignment? Three genes are drawn: the **most covered** gene (the most
genomes aligned), the gene whose naive window is **cleanest** among well-covered genes, and the gene with the **worst** naive window.

**How to read a figure.**

* x-axis: HMM match column, the same coordinate as the windows.
* **Window track** (top, under the bar strip): each of the three windows is one horizontal bar covering its 20 columns, labelled with its first and last
  column and the share of covering genomes in which it spans an intron. **Red** = naive window (best coverage, introns ignored),
  **magenta** = corrected window the workflow uses. The two edges of each window are repeated as vertical lines across the alignment.
* Each thin horizontal row is one genome, coloured by **phylum** (the 12 largest phyla for this gene get their own colour; the rest are pooled as
  "other phyla"). The coloured bar runs from the first to the last column that genome aligns to; a short bar is a partial gene.
* **Grey-to-black squares** are introns. Their **shade encodes the intron length** (light = short, black = long; log-scale bar at the bottom right).
  **Width is always one column**: the alignment has no coordinates for genomic length, so a long intron is darker, not wider. An intron is a point in the
  spliced sequence (the introns are already removed from the protein), so one column is the honest width. Most of a protein is not part of the HMM
  ("insert" residues between match columns), so many introns fall between two columns; they are drawn on the match column just before them (the strip on
  top shows these separately, in light grey). A window "spans" an intron when the intron lies between the window's first and last residue. Pale rows are
  genes whose model had a frameshift or length mismatch, so their introns are unknown and none is drawn.
* **Flag strip** (right of the alignment): a red mark means the naive window spans an intron in that genome, a magenta mark the corrected window. A window
  is invisible to reads of every flagged genome.
* The strip on top gives, for each column, the share of mapped genomes with an intron there (dark = on the column, light = in the insert after it).

**Which rows are drawn.** The 12 largest phyla, "other phyla" and "unclassified" each keep at least 25 rows. Within each group, genomes whose naive or
corrected window spans an intron are **deliberately oversampled** (up to half of that group's rows, or all of them if there are fewer), so the flagged
genomes are visible. Row counts therefore do not show how common spanning is; the percentages in the labels and titles use all genomes.
""")

code("""
def pick_feature_genes():
    \"\"\"The genes to draw: most covered, cleanest naive window among well-covered genes, worst naive window.

    Returns
    -------
    list of (label, gene)
        Distinct genes; a gene that qualifies twice is drawn once under its first label.
    \"\"\"
    if FEATURE_GENES:
        return [tuple(item) for item in FEATURE_GENES]
    well_covered = choices.filter(pl.col("n_mapped") >= choices["n_mapped"].median())
    picks = [("most covered", choices.sort(["n_mapped", "busco_id"], descending=[True, False])["busco_id"][0]),
             ("cleanest naive window", well_covered.sort(["naive_span_fraction", "n_mapped"], descending=[False, True])["busco_id"][0]),
             ("worst naive window", choices.sort(["naive_span_fraction", "busco_id"], descending=[True, False])["busco_id"][0])]
    seen, out = set(), []
    for label, gene in picks:
        if gene not in seen:
            seen.add(gene); out.append((label, gene))
    return out


def phylum_palette(phyla_by_size, n_colors=N_PHYLUM_COLORS):
    \"\"\"Colours for the largest phyla; the rest are pooled as "other phyla", unclassified is separate.

    The palette leaves out the greys (introns), the reds and greens (naive window) and the pinks (corrected window).
    \"\"\"
    tab = plt.get_cmap("tab20")
    usable = [tab(i) for i in (0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19)]
    colors = {}
    for phylum in [p for p in phyla_by_size if p != "unclassified"][:n_colors]:
        colors[phylum] = usable[len(colors) % len(usable)]
    colors["other phyla"] = (0.55, 0.43, 0.19, 1.0)
    colors["unclassified"] = (0.90, 0.85, 0.72, 1.0)
    return colors


INTRON_LEN_RANGE = (30, 10000)    # nt; shorter / longer introns get the lightest / darkest shade
INTRON_CMAP = matplotlib.colors.LinearSegmentedColormap.from_list("intron_length", [(0.8, 0.8, 0.8), (0.0, 0.0, 0.0)])


def intron_shade(length):
    \"\"\"Grey level (0 = black = long intron, 0.8 = light = short) for an intron length in nt.\"\"\"
    lo, hi = np.log10(INTRON_LEN_RANGE[0]), np.log10(INTRON_LEN_RANGE[1])
    return 0.8 * (1 - (np.log10(np.clip(length, *INTRON_LEN_RANGE)) - lo) / (hi - lo))


def sample_rows(rows, groups, seed):
    \"\"\"Stratified subsample of genome rows that oversamples genomes whose window spans an intron.

    Parameters
    ----------
    rows : pl.DataFrame
        Aligned genomes of one gene with ``color_group``, ``spans_naive`` and ``spans_corrected``.
    groups : pl.DataFrame
        Colour groups with their sizes (columns ``color_group``, ``len``).
    seed : int
        Seed of the random choice.

    Returns
    -------
    pl.DataFrame
        The rows to draw. Every group keeps at least MIN_ROWS_PER_GROUP rows (or all of them); inside a group, up to SPANNING_ROW_SHARE
        of the rows are spanning genomes (all of them if there are fewer), the rest are chosen at random.
    \"\"\"
    rng = np.random.default_rng(seed)
    total = rows.height
    take = []

    def pick(frame, k):
        return frame[rng.choice(frame.height, size=k, replace=False).tolist()] if k < frame.height else frame

    for group, n in zip(groups["color_group"], groups["len"]):
        part = rows.filter(pl.col("color_group") == group).sort("genome_acc")
        quota = min(n, max(MIN_ROWS_PER_GROUP, round(FEATURE_ROWS * n / total)))
        is_flagged = pl.col("spans_naive") | pl.col("spans_corrected")
        flagged, other = part.filter(is_flagged), part.filter(~is_flagged)
        n_flag = min(flagged.height, int(np.ceil(SPANNING_ROW_SHARE * quota)))
        n_other = min(other.height, quota - n_flag)
        n_flag = min(flagged.height, quota - n_other)
        take += [pick(flagged, n_flag), pick(other, n_other)]
    return pl.concat(take)


def draw_alignment(gene, label, seed=1):
    \"\"\"Draw one gene as an alignment of genome rows with introns, a window track and a spanning-flag strip.

    Parameters
    ----------
    gene : str
        BUSCO id.
    label : str
        Why this gene is shown (goes in the title).
    seed : int, optional
        Seed of the stratified row subsample.
    \"\"\"
    if "aln_first" not in genomes.columns or "hmm_col_prev" not in introns.columns:
        raise SystemExit("tables lack aln_first / hmm_col_prev: re-run rule intron_window")
    info = choices.filter(pl.col("busco_id") == gene).row(0, named=True)
    ncol = info["hmm_columns"]
    rows = genomes.filter((pl.col("gene") == gene) & (pl.col("aln_first") > 0))
    sizes = rows.group_by("phylum").len().sort(["len", "phylum"], descending=[True, False])
    palette = phylum_palette(sizes["phylum"].to_list())
    rows = rows.with_columns(pl.when(pl.col("phylum").is_in(list(palette))).then(pl.col("phylum"))
                             .otherwise(pl.lit("other phyla")).alias("color_group"))
    groups = rows.group_by("color_group").len().with_columns(
        pl.when(pl.col("color_group") == "unclassified").then(2).when(pl.col("color_group") == "other phyla").then(1).otherwise(0).alias("_tail")
    ).sort(["_tail", "len", "color_group"], descending=[False, True, False])

    total = rows.height
    order_of = {g: i for i, g in enumerate(groups["color_group"].to_list())}
    drawn = (sample_rows(rows, groups, seed).with_columns(pl.col("color_group").replace_strict(order_of).alias("_rank"))
             .sort(["_rank", "aln_first", "aln_last", "genome_acc"]))
    row_of = {acc: i for i, acc in enumerate(drawn["genome_acc"].to_list())}
    n_flag_all = rows.filter(pl.col("spans_naive") | pl.col("spans_corrected")).height
    n_flag_drawn = drawn.filter(pl.col("spans_naive") | pl.col("spans_corrected")).height

    image = np.ones((drawn.height, ncol, 4))
    for i, (group, first, last, ok) in enumerate(zip(drawn["color_group"], drawn["aln_first"], drawn["aln_last"], drawn["mapped"])):
        color = np.array(palette[group], dtype=float)
        if not ok:
            color = 0.45 * color + 0.55   # pale: no usable gene model
        image[i, first - 1:last, :] = color
    gene_introns = introns.filter(pl.col("busco_id") == gene)
    on_column = gene_introns.filter(pl.col("hmm_col") >= 0)
    in_insert = gene_introns.filter((pl.col("hmm_col") < 0) & (pl.col("hmm_col_prev") >= 0))
    before_first = gene_introns.height - on_column.height - in_insert.height
    n_drawn = {}
    for name, table in (("insert", in_insert), ("column", on_column)):   # column introns on top
        column_name = "hmm_col" if name == "column" else "hmm_col_prev"
        sub = table.filter(pl.col("genome_acc").is_in(list(row_of)))
        n_drawn[name] = sub.height
        for acc, col, length in zip(sub["genome_acc"], sub[column_name], sub["intron_len"]):
            shade = intron_shade(length)
            image[row_of[acc], col, :] = (shade, shade, shade, 1.0)

    # which drawn genomes have a window spanning an intron
    flags = np.ones((drawn.height, 2, 4))
    flags[drawn["spans_naive"].to_numpy(), 0] = matplotlib.colors.to_rgba(WINDOW_COLORS["naive"])
    flags[drawn["spans_corrected"].to_numpy(), 1] = matplotlib.colors.to_rgba(WINDOW_COLORS["corrected"])

    # share of mapped genomes with an intron at each column
    n_mapped = max(1, rows.filter(pl.col("mapped")).height)
    per_column = np.bincount(on_column["hmm_col"].to_numpy(), minlength=ncol)[:ncol] / n_mapped
    per_insert = np.bincount(in_insert["hmm_col_prev"].to_numpy(), minlength=ncol)[:ncol] / n_mapped

    fig = plt.figure(figsize=(13.5, 10.5))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.1, 0.8, 7], width_ratios=[60, 1.6], hspace=0.06, wspace=0.02)
    top = fig.add_subplot(grid[0, 0])
    track = fig.add_subplot(grid[1, 0], sharex=top)
    ax = fig.add_subplot(grid[2, 0], sharex=top)
    flag_ax = fig.add_subplot(grid[2, 1], sharey=ax)
    top.bar(np.arange(ncol), per_column, width=1.0, color="0.2")
    top.bar(np.arange(ncol), per_insert, bottom=per_column, width=1.0, color="0.62")
    top.set_ylabel("genomes with\\nan intron here"); top.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    top.tick_params(labelbottom=False)
    ax.imshow(image, aspect="auto", interpolation="nearest", extent=(-0.5, ncol - 0.5, drawn.height - 0.5, -0.5))
    flag_ax.imshow(flags, aspect="auto", interpolation="nearest", extent=(-0.5, 1.5, drawn.height - 0.5, -0.5))
    flag_ax.set_xticks([0, 1], ["naive", "corrected"], rotation=90, fontsize=8)
    flag_ax.tick_params(labelleft=False, left=False)
    flag_ax.set_xlabel("window spans\\nan intron", fontsize=8)
    ticks, names, start = [], [], 0
    for group, n_rows in drawn.group_by("color_group", maintain_order=True).len().iter_rows():
        ax.axhline(start - 0.5, color="white", lw=1)
        n_group = groups.filter(pl.col("color_group") == group)["len"][0]
        if n_rows >= 0.04 * drawn.height:   # small groups are only named in the legend
            ticks.append(start + n_rows / 2); names.append(f"{group}\\n{n_rows} of {n_group} drawn")
        start += n_rows
    ax.set_yticks(ticks, names, fontsize=8)
    ax.set_ylabel("phylum (groups too small to label are named in the legend)"); ax.set_xlabel("HMM match column")

    # one bar per window in the track, and its two edges as lines across the alignment
    handles = []
    track.set_ylim(-0.6, 1.6); track.set_yticks([1, 0], ["naive", "corrected"], fontsize=8)
    track.tick_params(labelbottom=False)
    for spine in ("left", "bottom"):
        track.spines[spine].set_visible(False)
    for y, (key, start_col, linestyle) in zip((1, 0), (("naive", info["naive_start"], "--"), ("corrected", info["corrected_start"], "-"))):
        lo, hi = start_col - 1 - 0.5, start_col - 1 + 20 - 0.5
        span = info[f"{key}_span_fraction"]
        track.add_patch(plt.Rectangle((lo, y - 0.35), 20, 0.7, fc=WINDOW_COLORS[key], ec="none"))
        right_side = hi < 0.62 * ncol
        track.text(hi + 1 if right_side else lo - 1, y, f"{start_col}-{start_col + 19}: spans an intron in {span:.1%}", color=WINDOW_COLORS[key],
                   fontsize=8, va="center", ha="left" if right_side else "right")
        for edge in (lo, hi):
            ax.axvline(edge, color=WINDOW_COLORS[key], lw=1.6, ls=linestyle)
        handles.append(plt.Rectangle((0, 0), 1, 1, fc=WINDOW_COLORS[key],
                                     label=f"{WINDOW_LABELS[key]}: columns {start_col}-{start_col + 19}, spans an intron in {span:.1%} of genomes"))
    track.set_xlim(-0.5, ncol - 0.5)
    handles += [plt.Rectangle((0, 0), 1, 1, fc="0.85", ec="0.6", label="pale row: no usable gene model, introns unknown")]
    handles += [plt.Rectangle((0, 0), 1, 1, fc=palette[g], label=g) for g in groups["color_group"]]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.94, 0.93), frameon=False, fontsize=8, title="windows and phyla (row colour)")
    cax = fig.add_axes([0.945, 0.10, 0.17, 0.012])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=matplotlib.colors.LogNorm(*INTRON_LEN_RANGE), cmap=INTRON_CMAP), cax=cax, orientation="horizontal")
    colorbar.set_label("intron length (nt, log scale; shade of the intron squares)", fontsize=8)
    colorbar.ax.tick_params(labelsize=8)
    top.set_title(f"{gene_label(gene, sep='  ')}  ({label})\\n{drawn.height} of {total} aligned genomes drawn\\n"
                  f"{n_flag_drawn} of the {n_flag_all} genomes whose naive or corrected window spans an intron are drawn\\n"
                  f"introns drawn: {n_drawn['column']} on a column, {n_drawn['insert']} in inserts; {before_first} lie before the first match column (not drawn)",
                  fontsize=9, loc="left")
    plt.savefig(OUT_DIR / f"alignment_{gene}.png", dpi=150, bbox_inches="tight")
    plt.show()


FEATURE = pick_feature_genes()
print("genes drawn:", FEATURE)
show_all(choices.filter(pl.col("busco_id").is_in([g for _, g in FEATURE]))
         .select("busco_id", "name", "n_mapped", "n_with_intron", "naive_start", "naive_span_fraction", "corrected_start", "corrected_span_fraction"))
""")

code("""
draw_alignment(FEATURE[0][1], FEATURE[0][0])
""")

code("""
if len(FEATURE) > 1:
    draw_alignment(FEATURE[1][1], FEATURE[1][0])
""")

code("""
if len(FEATURE) > 2:
    draw_alignment(FEATURE[2][1], FEATURE[2][0])
""")

md("""
## Per-gene decisions

One row per marker: the two windows, how often each crosses an intron, and what the corrected window costs in coverage. Also saved as `window_choices.tsv`.
""")

code("""
save_table(choices, "window_choices")
show_all(choices.sort(["naive_span_fraction", "busco_id"], descending=[True, False])
         .select("busco_id", "name", "n_mapped", "n_with_intron", "naive_start", "naive_span_fraction", "corrected_start",
                 "corrected_span_fraction", "corrected_score_ratio", "reason"),
         rows=140)
""")

md("""
## Summary
""")

code("""
summary = {
    "genes": choices.height,
    "genomes": genomes["genome_acc"].n_unique(),
    "records mapped": f"{mapped.height / genomes.height:.1%}",
    "introns called": introns.height,
    "intron sites": site_summary.height,
    "sites shared by >= 3 phyla": site_summary.filter(pl.col("phyla_ge_25pct") >= 3).height,
    "median introns per mapped gene": float(mapped["n_introns"].median()),
    "subset": PARAMS["subset"], "naive metapackage built": PARAMS["build_naive"], "limit": LIMIT,
    "naive windows over the limit": over_naive.height,
    "corrected window differs from naive": moved.height,
    "corrected windows over the limit": choices.filter(pl.col("violation")).height,
    "mean naive span fraction": round(float(choices["naive_span_fraction"].mean()), 4),
    "mean corrected span fraction": round(float(choices["corrected_span_fraction"].mean()), 4),
    "mean corrected score ratio": round(float(choices["corrected_score_ratio"].mean()), 4),
}
for key, value in summary.items():
    print(f"{key:36s} {value}")
save_table(pl.DataFrame({"metric": list(summary), "value": [str(v) for v in summary.values()]}), "summary")
""")

nb = new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}})
for i, c in enumerate(nb.cells):
    c["id"] = f"intron-{i:02d}"
nbf.write(nb, str(Path(__file__).resolve().parent / "intron_survey.ipynb"))
print(len(cells), "cells")
