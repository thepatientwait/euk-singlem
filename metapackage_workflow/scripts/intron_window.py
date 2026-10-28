#!/usr/bin/env python3
"""Intron survey and intron-aware window choice for one BUSCO gene (rule ``intron_window``).

For every genome that contributes a sequence to the gene this finds the BUSCO gene model
(miniprot GFF), calls its introns, aligns the protein to the BUSCO HMM and checks which
20-column (60 nt) windows would cross an intron. Windows are numbered in HMM match columns, the
coordinate ``singlem create --hmm-position`` and ``singlem pipe`` use. It writes one window per
metapackage variant (wildcard ``introns``):

* ``naive``     - the best-covered window, introns ignored (ties: the earliest);
* ``corrected`` - the best-covered window that crosses an intron in at most ``max_span_fraction``
  of the genomes covering it;

as ``window_position.<variant>.txt`` (1-based), plus

* ``introns``  - one row per intron (HMM column, phase, length) for conservation analysis. It also lists the
  introns of genomes whose gene model failed the frameshift/length check (``mapped`` false in ``genomes``);
  their HMM columns are unreliable, so consumers must keep mapped genomes only;
* ``genomes``  - one row per genome (exon/intron counts, mapping status, aligned extent, and whether the naive
  and corrected windows span an intron in it);
* ``scan``     - for every window start: coverage score and intron-spanning fraction;
* ``choice``   - the decision for the gene.

Runs inside Snakemake, or standalone: ``intron_window.py --help``.
"""
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

sys.path.append(str(Path(__file__).parent))
import intron_utils as iu  # noqa: E402

try:
    snakemake  # type: ignore[name-defined]
except NameError:
    snakemake = None

if snakemake is not None:
    from log_config import configure_logger
    logger = configure_logger(__file__)
else:
    import logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logger = logging.getLogger("intron_window")


def read_fasta(path: Path) -> dict:
    records, name, chunks = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def read_busco_results(manifest: Path) -> dict:
    """genome accession -> BUSCO full_table path (later manifest entries win, like the aggregator)."""
    table = pl.read_csv(manifest, separator="\t")
    return {Path(p).parts[-3]: Path(p) for p in sorted(table["busco_result"].to_list())}


def top_hit(full_table: Path, busco_id: str, use_multi_copy: bool):
    """The row the aggregator uses for this gene: best Score, then longest, then position."""
    header, rows = None, []
    with open(full_table) as handle:
        for line in handle:
            if line.startswith("# Busco id"):
                header = line.lstrip("# ").rstrip("\n").split("\t")
            elif not line.startswith("#"):
                parts = line.rstrip("\n").split("\t")
                if parts[0] == busco_id and header:
                    rows.append(dict(zip(header, parts)))
    statuses = ("Complete", "Duplicated") if use_multi_copy else ("Complete",)
    rows = [r for r in rows if r.get("Status") in statuses]
    if not rows:
        return None

    def key(r):
        return (-float(r["Score"]), -int(r.get("Length") or 0), r["Sequence"], int(r["Gene Start"]))

    return sorted(rows, key=key)[0]


def gene_model(full_table: Path, busco_id: str, use_multi_copy: bool):
    """(status, model, row) for one genome and gene; ``model`` is None when it can't be found."""
    if not full_table.exists():
        return "no_full_table", None
    row = top_hit(full_table, busco_id, use_multi_copy)
    if row is None:
        return "no_hit", None
    subdir = "single_copy_busco_sequences" if row["Status"] == "Complete" else "multi_copy_busco_sequences"
    gff = full_table.parent / "busco_sequences" / subdir / f"{busco_id}.gff"
    if not gff.exists():
        return "no_gff", None
    models = iu.parse_gff(gff)
    if not models:
        return "no_model", None
    model = iu.select_model(models, row["Sequence"], int(row["Gene Start"]), int(row["Gene End"]), row.get("Strand"))
    if model is None and len(models) == 1:
        model = models[0]
    return ("ok", model) if model is not None else ("no_model", None)


def run(a) -> dict:
    busco_id = a.busco_id
    proteins_raw = read_fasta(a.sequences)
    results = read_busco_results(a.busco_manifest)
    names = list(proteins_raw)
    genome_of = {n: n.split("|")[0] for n in names}
    proteins = {n: iu.clean_protein(s) for n, s in proteins_raw.items()}

    def analyse(name):
        genome = genome_of[name]
        table = results.get(genome)
        if table is None:
            return name, "not_in_manifest", None, None
        status, model = gene_model(table, busco_id, a.use_multi_copy)
        return name, status, model, table

    with ThreadPoolExecutor(max(1, a.threads)) as pool:
        analysed = list(pool.map(analyse, names))

    n = len(names)
    status = [""] * n
    group = [""] * n
    n_exons = np.zeros(n, dtype=int)
    cds_len = np.zeros(n, dtype=int)
    introns = [[] for _ in range(n)]
    gaps = [[] for _ in range(n)]
    contig = [""] * n
    strand = [""] * n
    valid = np.zeros(n, dtype=bool)
    for i, (name, st, model, table) in enumerate(analysed):
        status[i] = st
        group[i] = table.parts[-4] if table is not None else ""
        if model is None:
            continue
        pairs, cds = iu.intron_gaps(model, a.min_intron)
        contig[i], strand[i] = model.contig, model.strand
        n_exons[i], cds_len[i] = len(model.exons), cds
        introns[i] = [o for o, _ in pairs]
        gaps[i] = [g for _, g in pairs]
        if iu.cds_matches_protein(model, cds, len(proteins[name])):
            valid[i] = True
        else:
            status[i] = "frameshift" if model.irregular else "length_mismatch"

    hmm_names, matrix = iu.hmmalign_matrix(a.hmm, proteins, a.hmmalign)
    assert hmm_names == names
    ncol = matrix.shape[1]
    window_cols = a.window_nt // 3
    scan = iu.window_scan(matrix, introns, valid, window_cols)

    corrected, reason, violation, naive = iu.choose_window(
        scan, a.max_span_fraction, a.score_tolerance, a.min_valid)

    def frac(start):
        v = scan["n_cover_valid"][start]
        return float(scan["n_span"][start] / v) if v else float("nan")

    corrected_frac = frac(corrected)
    corrected_violation = bool(scan["n_cover_valid"][corrected] >= a.min_valid and corrected_frac > a.max_span_fraction)
    best_score = float(scan["aligned_bases"].max())

    # --- outputs ---------------------------------------------------------------------
    for path in (a.out_position_naive, a.out_position_corrected, a.out_introns, a.out_genomes, a.out_scan, a.out_choice):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out_position_naive).write_text(f"{naive + 1}\n")
    Path(a.out_position_corrected).write_text(f"{corrected + 1}\n")

    intron_rows = []
    for i in range(n):
        for offset, gap in zip(introns[i], gaps[i]):
            residue = offset // 3
            res_cols = iu.residue_columns(matrix[i], len(proteins[names[i]]))
            column = int(res_cols[residue]) if residue < len(res_cols) else -1
            # an intron inside an insert has no match column of its own: also record the last match
            # column before it (-1 if it lies before the first one), so it can be placed on the alignment
            before = res_cols[:residue + 1]
            before = before[before >= 0]
            intron_rows.append({
                "genome_acc": genome_of[names[i]], "group": group[i], "busco_id": busco_id,
                "contig": contig[i], "strand": strand[i], "offset": offset, "phase": offset % 3,
                "residue": residue,
                "hmm_col": column, "hmm_col_prev": column if column >= 0 else (int(before[-1]) if len(before) else -1),
                "intron_len": gap, "n_introns": len(introns[i]),
            })
    pl.DataFrame(intron_rows, schema={
        "genome_acc": pl.Utf8, "group": pl.Utf8, "busco_id": pl.Utf8, "contig": pl.Utf8, "strand": pl.Utf8,
        "offset": pl.Int64, "phase": pl.Int64, "residue": pl.Int64, "hmm_col": pl.Int64, "hmm_col_prev": pl.Int64,
        "intron_len": pl.Int64, "n_introns": pl.Int64}).write_csv(a.out_introns, separator="\t")

    def window_flags(start):
        return scan["covers"][:, start], scan["spans"][:, start]

    aln_first, aln_last, aln_cols = iu.aligned_extent(matrix)
    cov_n, span_n = window_flags(naive)
    cov_c, span_c = window_flags(corrected)
    pl.DataFrame({
        "genome_acc": [genome_of[x] for x in names], "group": group, "status": status,
        "n_exons": n_exons, "n_introns": [len(x) for x in introns],
        "protein_len": [len(proteins[x]) for x in names], "cds_len": cds_len, "mapped": valid,
        # 1-based HMM match columns the protein spans (-1 if none) and how many of them it occupies
        "aln_first": np.where(aln_first >= 0, aln_first + 1, -1), "aln_last": np.where(aln_last >= 0, aln_last + 1, -1),
        "aln_cols": aln_cols,
        "covers_naive": cov_n, "spans_naive": span_n & valid,
        "covers_corrected": cov_c, "spans_corrected": span_c & valid,
    }).write_csv(a.out_genomes, separator="\t")

    starts = np.arange(len(scan["n_cover"]))
    pl.DataFrame({
        "start": starts + 1,  # 1-based
        "n_cover": scan["n_cover"], "aligned_bases": scan["aligned_bases"],
        "n_cover_valid": scan["n_cover_valid"], "n_span": scan["n_span"],
        "span_fraction": np.where(scan["n_cover_valid"] > 0, scan["n_span"] / np.maximum(scan["n_cover_valid"], 1), np.nan),
        "is_naive": starts == naive, "is_corrected": starts == corrected,
    }).write_csv(a.out_scan, separator="\t")

    def ratio(start):
        return float(scan["aligned_bases"][start]) / best_score if best_score else float("nan")

    choice = {
        "busco_id": busco_id, "hmm_columns": ncol, "n_genomes": n, "n_mapped": int(valid.sum()),
        "n_with_intron": int(sum(1 for x in introns if x)),
        "naive_start": naive + 1, "naive_span_fraction": frac(naive), "naive_n_cover": int(scan["n_cover_valid"][naive]),
        "corrected_start": corrected + 1, "corrected_span_fraction": corrected_frac,
        "corrected_score_ratio": ratio(corrected), "reason": reason,
        "max_span_fraction": a.max_span_fraction, "violation": corrected_violation,
    }
    pl.DataFrame([choice]).write_csv(a.out_choice, separator="\t")
    logger.info("%s: naive window %d spans an intron in %.1f%% of %d genomes; corrected window %d (%.1f%%, %.1f%% of best score): %s",
                busco_id, naive + 1, 100 * choice["naive_span_fraction"], choice["naive_n_cover"],
                corrected + 1, 100 * corrected_frac, 100 * choice["corrected_score_ratio"], reason)
    if corrected_violation:
        logger.warning("%s: window %d still spans an intron in %.1f%% of genomes (limit %.1f%%)",
                       busco_id, corrected + 1, 100 * corrected_frac, 100 * a.max_span_fraction)
    return choice


def args_from_snakemake(smk):
    p = smk.params
    return argparse.Namespace(
        busco_id=smk.wildcards.busco_gene, sequences=smk.input.sequences, hmm=smk.input.hmm,
        busco_manifest=smk.input.busco_manifest,
        out_position_naive=smk.output.position_naive, out_position_corrected=smk.output.position_corrected,
        out_introns=smk.output.introns, out_genomes=smk.output.genomes,
        out_scan=smk.output.scan, out_choice=smk.output.choice,
        min_intron=int(p.min_intron), max_span_fraction=float(p.max_span_fraction),
        score_tolerance=float(p.score_tolerance), min_valid=int(p.min_valid), window_nt=int(p.window_nt),
        use_multi_copy=bool(p.use_multi_copy), hmmalign="hmmalign", threads=int(smk.threads))


def args_from_cli():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--busco-id", required=True)
    ap.add_argument("--sequences", required=True, help="aggregated gene FASTA (outputs/singlem/input/<gene>.faa)")
    ap.add_argument("--hmm", required=True)
    ap.add_argument("--busco-manifest", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--min-intron", type=int, default=iu.DEFAULT_MIN_INTRON)
    ap.add_argument("--max-span-fraction", type=float, default=0.05)
    ap.add_argument("--score-tolerance", type=float, default=0.005)
    ap.add_argument("--min-valid", type=int, default=20)
    ap.add_argument("--window-nt", type=int, default=iu.WINDOW_NT)
    ap.add_argument("--no-multi-copy", dest="use_multi_copy", action="store_false")
    ap.add_argument("--hmmalign", default="hmmalign")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    d = Path(a.out_dir)
    a.out_position_naive, a.out_position_corrected = (d / f"{a.busco_id}.window_position.{v}.txt" for v in ("naive", "corrected"))
    a.out_introns, a.out_genomes = d / f"{a.busco_id}.introns.tsv", d / f"{a.busco_id}.genomes.tsv"
    a.out_scan, a.out_choice = d / f"{a.busco_id}.window_scan.tsv", d / f"{a.busco_id}.window_choice.tsv"
    return a


if __name__ == "__main__":
    run(args_from_snakemake(snakemake) if snakemake is not None else args_from_cli())
