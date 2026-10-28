#!/usr/bin/env python3
"""Discover NCBI TSA (transcriptome shotgun assembly) projects to add as transcriptome inputs.

Steps
-----
1. Entrez: every ``tsa-master`` nuccore record of Eukaryota, minus Metazoa and Viridiplantae.
2. Drop projects whose taxon is already covered by a genome (``tsa.only_novel_taxa``): the
   genomes are trusted more. Coverage is judged on species name (``tsa.novel_level: species``)
   or on the exact tax ID (``taxid``) against the taxonomy mapping.
3. Look up the bulk files NCBI publishes for each project
   (``genbank/tsa/<first letter>/tsa.<PREFIX>.<part>.fsa_nt.gz``) in the FTP directory listing.
4. Keep one project per tax ID (the largest), above a minimum size.

Output (snakemake.output.metadata): TSV with the columns of ``ncbi_genomes.tsv``
(``genome_acc, tax_id, original_name, download_source, group``) plus ``prefix``,
``version``, ``parts`` (semicolon-separated file names), ``size_bytes``, ``title``.
``genome_acc`` is the master accession with its version, e.g. ``GLPS00000000.1``.

Snakemake params: ``email``, ``only_novel_taxa``, ``novel_level``, ``min_size_mb``,
``max_projects`` (0 = all).  ``taxonomy`` (path of taxonomy_mapping.tsv, read if it exists); inputs: the genome tables
(ncbi_genomes.tsv, p10k_meta.csv) that define which taxa a genome covers.
"""
import csv
import os
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

from log_config import configure_logger

logger = configure_logger(os.path.splitext(os.path.basename(__file__))[0])

FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/genbank/tsa"
SEARCH_TERM = ('("Eukaryota"[Organism] NOT "Metazoa"[Organism] NOT "Viridiplantae"[Organism]) '
               "AND tsa-master[prop]")
COLUMNS = ["genome_acc", "tax_id", "original_name", "download_source", "group",
           "prefix", "version", "parts", "size_bytes", "title"]

_TITLE_RE = re.compile(r"^TSA:\s*(?P<name>.+?)(?:\s+(?:strain|isolate|cultivar|voucher|clone|"
                       r"ecotype|breed|genotype)\b.*|,.*|\s+transcriptome shotgun assembly.*)?$",
                       re.IGNORECASE)
_SIZE_RE = re.compile(r"^(?P<n>[\d.]+)(?P<u>[KMGT]?)$")
_UNITS = {"": 1, "K": 1024, "M": 1024 ** 2, "G": 1024 ** 3, "T": 1024 ** 4}


def organism_from_title(title: str) -> str:
    """'TSA: Coolia malayensis strain DF364, transcriptome shotgun assembly' -> 'Coolia malayensis'."""
    match = _TITLE_RE.match((title or "").strip())
    return match.group("name").strip() if match else ""


def split_accession(accession: str):
    """'GLPS00000000.1' -> ('GLPS', 1): leading letters and the assembly version."""
    base, _, ver = accession.partition(".")
    match = re.match(r"^([A-Z]+)\d*$", base)
    if not match:
        return "", 0
    return match.group(1), int(ver) if ver.isdigit() else 1


def parse_size(text: str) -> int:
    match = _SIZE_RE.match(text.strip())
    if not match:
        return 0
    return int(float(match.group("n")) * _UNITS[match.group("u")])


def parse_listing(html: str):
    """{prefix: [(part_file_name, size_bytes), ...]} from an Apache index of genbank/tsa/<letter>/."""
    parts = defaultdict(list)
    pattern = re.compile(r'href="(tsa\.([A-Z]+)\.(\d+)\.fsa_nt\.gz)"[^\n]*?\s(\d+(?:\.\d+)?[KMGT]?)\s*$',
                         re.MULTILINE)
    for name, prefix, part, size in pattern.findall(html):
        parts[prefix].append((name, parse_size(size), int(part)))
    return {p: [(n, s) for n, s, _ in sorted(v, key=lambda x: x[2])] for p, v in parts.items()}


def select_projects(records, covered_taxids, covered_species, listings, only_novel=True,
                    level="species", min_size_bytes=0):
    """Pure selection step.

    records: dicts with ``accession``, ``tax_id``, ``title``.  listings: ``{prefix: [(file, size)]}``.
    Returns (selected rows, dropped counts by reason).
    """
    dropped = defaultdict(int)
    best = {}
    for rec in records:
        prefix, version = split_accession(rec["accession"])
        name = organism_from_title(rec["title"])
        tax_id = str(rec["tax_id"])
        if not tax_id or tax_id == "0":
            dropped["no_taxid"] += 1
            continue
        if only_novel:
            if level == "taxid" and tax_id in covered_taxids:
                dropped["covered_by_genome"] += 1
                continue
            if level == "species" and (tax_id in covered_taxids or _species_key(name) in covered_species):
                dropped["covered_by_genome"] += 1
                continue
        files = listings.get(prefix)
        if not files:
            dropped["no_bulk_file"] += 1
            continue
        size = sum(s for _, s in files)
        if size < min_size_bytes:
            dropped["too_small"] += 1
            continue
        row = {
            "genome_acc": rec["accession"], "tax_id": tax_id, "original_name": name,
            "download_source": "TSA", "group": "tsa", "prefix": prefix, "version": version,
            "parts": ";".join(n for n, _ in files), "size_bytes": size, "title": rec["title"],
        }
        current = best.get(tax_id)
        if current is None or size > current["size_bytes"]:
            if current is not None:
                dropped["smaller_same_taxid"] += 1
            best[tax_id] = row
        else:
            dropped["smaller_same_taxid"] += 1
    return sorted(best.values(), key=lambda r: r["genome_acc"]), dict(dropped)


def _species_key(name: str) -> str:
    """'Ochromonas species_1' / 'Homo sapiens strain X' -> 'ochromonas species' / 'homo sapiens' (two words)."""
    name = re.sub(r"_\d+$", "", (name or "").strip())
    words = name.split()
    return " ".join(words[:2]).lower() if len(words) >= 2 else ""


def read_covered(taxonomy_path=None, genome_tables=(), exclude_accs=()):
    """Tax IDs and species names that a genome already covers.

    Sources: the taxonomy mapping when it exists (columns ``acc, ncbi_id, species``), plus the
    pre-taxonomy genome tables (``tax_id``/``original_name`` or ``Scientific name``) so a first run
    that has no mapping yet is still correct.  ``exclude_accs`` (earlier TSA accessions, which a
    previous run may have put into the mapping) never count as genomes.
    """
    exclude = set(exclude_accs)
    taxids, species = set(), set()

    def add(tax_id, name):
        if tax_id and str(tax_id) not in ("", "nan", "None"):
            taxids.add(str(tax_id))
        key = _species_key(name)
        if key:
            species.add(key)

    if taxonomy_path and Path(taxonomy_path).exists():
        with open(taxonomy_path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("acc") not in exclude:
                    add(row.get("ncbi_id"), row.get("species"))
    for table in genome_tables:
        if not table or not Path(table).exists():
            continue
        delim = "," if str(table).endswith(".csv") else "\t"
        with open(table, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter=delim):
                if (row.get("genome_acc") or row.get("Accession")) in exclude:
                    continue
                add(row.get("tax_id"), row.get("original_name") or row.get("Scientific name"))
    return taxids, species


def search_tsa_masters(email: str, batch: int = 200):
    from Bio import Entrez
    Entrez.email = email
    handle = Entrez.esearch(db="nuccore", term=SEARCH_TERM, usehistory="y", retmax=0)
    found = Entrez.read(handle)
    handle.close()
    count = int(found["Count"])
    logger.info("Entrez: %d TSA master records", count)
    records = []
    for start in range(0, count, batch):
        for attempt in range(5):
            try:
                handle = Entrez.esummary(db="nuccore", query_key=found["QueryKey"], WebEnv=found["WebEnv"],
                                         retstart=start, retmax=batch)
                summaries = Entrez.read(handle)
                handle.close()
                break
            except Exception as err:  # transient NCBI errors
                logger.warning("esummary %d failed (%s), retry %d", start, err, attempt + 1)
                time.sleep(2 ** attempt)
        else:
            raise RuntimeError(f"esummary failed repeatedly at {start}")
        for s in summaries:
            records.append({"accession": s.get("AccessionVersion", ""), "tax_id": str(int(s.get("TaxId", 0) or 0)),
                            "title": s.get("Title", "")})
        time.sleep(0.4)
    return records


EXCLUDED_CLADES = ("Metazoa", "Viridiplantae")


def excluded_by_lineage(lineages: dict, excluded=EXCLUDED_CLADES):
    """Tax IDs whose lineage string (``'Eukaryota; Opisthokonta; Metazoa; ...'``) contains an excluded clade."""
    out = set()
    for tax_id, lineage in lineages.items():
        names = {n.strip() for n in lineage.split(";")}
        if names & set(excluded):
            out.add(tax_id)
    return out


def fetch_lineages(tax_ids, email: str, batch: int = 200) -> dict:
    """{tax_id: lineage string} from the Entrez taxonomy database."""
    from Bio import Entrez
    Entrez.email = email
    ids = sorted(set(tax_ids))
    lineages = {}
    for start in range(0, len(ids), batch):
        chunk = ids[start:start + batch]
        for attempt in range(5):
            try:
                handle = Entrez.efetch(db="taxonomy", id=",".join(chunk), retmode="xml")
                taxa = Entrez.read(handle)
                handle.close()
                break
            except Exception as err:  # transient NCBI errors
                logger.warning("taxonomy efetch failed (%s), retry %d", err, attempt + 1)
                time.sleep(2 ** attempt)
        else:
            raise RuntimeError("taxonomy efetch failed repeatedly")
        for taxon in taxa:
            lineages[str(taxon["TaxId"])] = taxon.get("Lineage", "")
            for alias in taxon.get("AkaTaxIds", []):
                lineages[str(alias)] = taxon.get("Lineage", "")
        time.sleep(0.4)
    return lineages


def fetch_listings(letters):
    import requests
    listings = {}
    for letter in sorted(letters):
        url = f"{FTP_ROOT}/{letter}/"
        for attempt in range(4):
            try:
                response = requests.get(url, timeout=120)
                if response.status_code == 404:  # no TSA files under this letter: those projects are dropped
                    logger.warning("no listing at %s", url)
                    response = None
                    break
                response.raise_for_status()
                break
            except requests.RequestException as err:
                logger.warning("listing %s failed (%s)", url, err)
                time.sleep(2 ** attempt)
        else:
            raise RuntimeError(f"could not list {url}")
        if response is not None:
            listings.update(parse_listing(response.text))
        logger.info("listing %s: %d prefixes so far", url, len(listings))
    return listings


def main() -> None:
    p = snakemake.params  # type: ignore[name-defined]
    out = Path(snakemake.output.metadata)  # type: ignore[name-defined]
    previous = set()
    if out.exists():
        with open(out, newline="", encoding="utf-8") as handle:
            previous = {r["genome_acc"] for r in csv.DictReader(handle, delimiter="\t")}
    records = [r for r in search_tsa_masters(p.email) if r["accession"]]
    listings = fetch_listings({split_accession(r["accession"])[0][:1] for r in records})
    # The taxonomy mapping is read by path, not declared as an input: it depends on this rule's output
    # (TSA taxa are added to it), which would be a cycle. Rerun with --forcerun discover_tsa after it changes.
    taxids, species = read_covered(p.taxonomy, list(snakemake.input), previous)  # type: ignore[name-defined]
    rows, dropped = select_projects(
        records, taxids, species, listings, only_novel=bool(p.only_novel_taxa), level=p.novel_level,
        min_size_bytes=int(float(p.min_size_mb) * 1024 ** 2))
    # the Entrez organism query is not reliable at excluding animals and plants: check the lineages
    lineages = fetch_lineages([r["tax_id"] for r in rows], p.email)
    bad = excluded_by_lineage(lineages)
    unknown = [r["tax_id"] for r in rows if r["tax_id"] not in lineages]
    if unknown:
        logger.warning("no lineage for tax IDs %s; kept", unknown[:10])
    if bad:
        logger.info("dropping %d projects in %s", sum(r["tax_id"] in bad for r in rows), "/".join(EXCLUDED_CLADES))
        dropped["metazoa_or_plant"] = sum(r["tax_id"] in bad for r in rows)
        rows = [r for r in rows if r["tax_id"] not in bad]
    if int(p.max_projects) > 0:
        rows = rows[:int(p.max_projects)]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    logger.info("%d records -> %d TSA projects selected; dropped: %s", len(records), len(rows), dropped)


if __name__ == "__main__":
    try:
        snakemake  # type: ignore[name-defined]
    except NameError:
        sys.exit("run through snakemake")
    main()
