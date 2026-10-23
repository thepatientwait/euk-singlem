#!/usr/bin/env python3
"""Query NCBI for the taxonomy of every genome in a manifest (single source of truth).

This is the ONLY place where NCBI lineages are fetched, cleaned and formatted.
Both the metapackage and the tree workflow use it; downstream scripts read the
mapping it writes and do not re-parse or re-clean names.

Input  (snakemake.input.unified_manifest): TSV with columns ``genome_acc``, ``tax_id``
Output (snakemake.output.taxonomy): TSV with columns

    acc, ncbi_id, domain, phylum, class, order, family, genus, species,
    tax_status, parsed_taxon_string

* ``parsed_taxon_string`` is ``d__X;p__X;c__X;o__X;f__X;g__X;s__X`` (``;`` with no
  spaces, always seven ranks, blank ranks left empty).
* ``tax_status`` is the one definition of "usable taxonomy":
    ``complete``    all seven ranks populated
    ``partial``     domain plus at least one other rank
    ``unresolved``  nothing below domain
  Genomes with no NCBI record are absent from the mapping (i.e. "missing").

Cleaning (applied once, here):
* parenthetical text is removed; ``sp. X`` -> ``X``; a bare ``sp.`` -> ``species``;
  ``:`` and ``;`` -> ``_``
* placeholder names (``uncultured``, ``Eukaryota`` used as a rank name) are blanked
* the domain is always ``Eukarya``
* different taxa that end up with the same species name are indexed
  (``Ochromonas species_1``, ``Ochromonas species_2``, ordered by tax ID).
  Genomes of the same taxon share one name and are not indexed.

Snakemake params: ``email`` (required); ``strict`` (optional, default False) makes the
script fail if any requested tax ID could not be resolved.
"""
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, NamedTuple, Optional

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from log_config import configure_logger  # noqa: E402

logger = configure_logger(os.path.splitext(os.path.basename(__file__))[0])

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = {rank: f"{rank[0]}__" for rank in RANKS}
DOMAIN_NAME = "Eukarya"
USABLE_STATUSES = ("complete", "partial")

# NCBI rank names that map onto our ranks (kingdom and everything else is ignored)
RANK_ALIASES = {"superkingdom": "domain", "domain": "domain", **{r: r for r in RANKS[1:]}}

DATASETS_BATCH_SIZE = 200
ENTREZ_BATCH_SIZE = 200


class Lineage(NamedTuple):
    taxon_id: int                        # primary tax ID of the queried taxon
    ranks: Dict[str, str]                # rank -> cleaned name (blank ranks omitted)
    species_id: Optional[int] = None     # tax ID of the species-level taxon, if any

    @property
    def species_key(self) -> int:
        """Identity of the species: strains of one species share it."""
        return self.species_id or self.taxon_id


# --------------------------------------------------------------------------
# Cleaning / formatting (pure functions)
# --------------------------------------------------------------------------
def clean_name(name: Optional[str]) -> str:
    """Normalise one NCBI rank name; returns '' if it is not a usable name."""
    if not name:
        return ""
    name = re.sub(r"\s*\([^)]*\)", "", str(name))
    name = name.replace("sp. ", "").replace("sp.", "species")
    name = name.replace(":", "_").replace(";", "_")
    name = " ".join(name.split())
    lowered = name.lower()
    if "eukaryota" in lowered or "uncultured" in lowered:
        return ""
    return name


def clean_ranks(raw: Dict[str, str]) -> Dict[str, str]:
    """raw rank -> NCBI name  =>  cleaned dict (domain forced, blanks dropped)."""
    cleaned = {}
    for rank in RANKS[1:]:
        value = clean_name(raw.get(rank))
        if value:
            cleaned[rank] = value
    cleaned["domain"] = DOMAIN_NAME
    return cleaned


def tax_status(ranks: Dict[str, str]) -> str:
    below_domain = [bool(ranks.get(rank)) for rank in RANKS[1:]]
    if all(below_domain):
        return "complete"
    if any(below_domain):
        return "partial"
    return "unresolved"


def format_lineage(ranks: Dict[str, str]) -> str:
    return ";".join(f"{PREFIXES[rank]}{ranks.get(rank, '')}" for rank in RANKS)


def index_duplicate_species(lineages: Dict[int, Lineage]) -> Dict[int, Lineage]:
    """Give taxa that share a species name a numeric suffix (ordered by tax ID).

    ``lineages`` maps requested tax ID -> Lineage. Taxa are compared by species tax ID, so
    strains of one species and secondary (merged) IDs keep the same single name.
    """
    species_by_name = defaultdict(set)
    for lineage in lineages.values():
        species = lineage.ranks.get("species")
        if species:
            species_by_name[species].add(lineage.species_key)

    new_names = {}
    for name, species_keys in species_by_name.items():
        if len(species_keys) > 1:
            for index, key in enumerate(sorted(species_keys), start=1):
                new_names[(name, key)] = f"{name}_{index}"

    result = {}
    for tax_id, lineage in lineages.items():
        new_name = new_names.get((lineage.ranks.get("species"), lineage.species_key))
        if new_name:
            lineage = lineage._replace(ranks=dict(lineage.ranks, species=new_name))
        result[tax_id] = lineage
    return result


# --------------------------------------------------------------------------
# Record parsing (pure functions, easy to test)
# --------------------------------------------------------------------------
def lineage_from_datasets_record(record: dict) -> Optional[tuple]:
    """NCBI Datasets JSON-lines record -> (Lineage, secondary_ids)."""
    taxonomy = record.get("taxonomy", {}) or {}
    tax_id = int(taxonomy.get("tax_id", 0) or 0)
    if not tax_id:
        return None
    raw = {}
    species_id = None
    for rank, item in (taxonomy.get("classification", {}) or {}).items():
        mapped = RANK_ALIASES.get(rank.lower())
        if mapped and isinstance(item, dict):
            raw[mapped] = item.get("name", "")
            if mapped == "species" and item.get("id"):
                species_id = int(item["id"])
    secondary = [int(x) for x in taxonomy.get("secondary_tax_ids", []) or []]
    return Lineage(tax_id, clean_ranks(raw), species_id), secondary


def lineage_from_entrez_record(record: dict) -> Optional[tuple]:
    """Entrez Taxonomy record -> (Lineage, secondary_ids)."""
    tax_id = int(record.get("TaxId", 0) or 0)
    if not tax_id:
        return None
    raw = {}
    species_id = None
    for item in record.get("LineageEx", []) or []:
        mapped = RANK_ALIASES.get(str(item.get("Rank", "")).lower())
        if mapped:
            raw[mapped] = item.get("ScientificName", "")
            if mapped == "species" and item.get("TaxId"):
                species_id = int(item["TaxId"])
    # LineageEx excludes the record itself; species-level records carry the species name
    own_rank = RANK_ALIASES.get(str(record.get("Rank", "")).lower())
    if own_rank:
        raw[own_rank] = record.get("ScientificName", "")
        if own_rank == "species":
            species_id = tax_id
    secondary = [int(x) for x in record.get("AkaTaxIds", []) or []]
    return Lineage(tax_id, clean_ranks(raw), species_id), secondary


# --------------------------------------------------------------------------
# Querying NCBI
# --------------------------------------------------------------------------
def _batches(values: List[int], size: int) -> Iterable[List[int]]:
    for start in range(0, len(values), size):
        yield values[start:start + size]


def _register(found: Dict[int, Lineage], parsed: Optional[tuple]) -> None:
    if parsed is None:
        return
    lineage, secondary = parsed
    for tax_id in [lineage.taxon_id, *secondary]:
        found[tax_id] = lineage


def query_datasets(tax_ids: List[int]) -> Dict[int, Lineage]:
    """Query the NCBI Datasets CLI. Returns {tax_id: Lineage}."""
    found: Dict[int, Lineage] = {}
    for number, batch in enumerate(_batches(tax_ids, DATASETS_BATCH_SIZE), start=1):
        cmd = ["datasets", "summary", "taxonomy", "taxon", ",".join(map(str, batch)), "--as-json-lines"]
        logger.info("NCBI Datasets batch %d (%d IDs)", number, len(batch))
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.splitlines():
            if not line.strip():
                continue
            try:
                _register(found, lineage_from_datasets_record(json.loads(line)))
            except (json.JSONDecodeError, ValueError, KeyError, AttributeError) as exc:
                logger.debug("Skipping unparsable Datasets record: %s", exc)
    return found


def query_entrez(tax_ids: List[int], email: str) -> Dict[int, Lineage]:
    """Query NCBI Entrez Taxonomy. Returns {tax_id: Lineage}."""
    from Bio import Entrez

    Entrez.email = email
    found: Dict[int, Lineage] = {}
    for number, batch in enumerate(_batches(tax_ids, ENTREZ_BATCH_SIZE), start=1):
        logger.info("NCBI Entrez batch %d (%d IDs)", number, len(batch))
        handle = Entrez.efetch(db="Taxonomy", id=",".join(map(str, batch)), retmode="xml")
        try:
            records = Entrez.read(handle)
        finally:
            handle.close()
        for record in records:
            _register(found, lineage_from_entrez_record(record))
    return found


def resolve_lineages(tax_ids: List[int], email: str) -> Dict[int, Lineage]:
    """Datasets first; anything it did not return (or all, if the CLI fails) via Entrez."""
    found: Dict[int, Lineage] = {}
    try:
        found.update(query_datasets(tax_ids))
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        detail = getattr(exc, "stderr", None) or exc
        logger.error("NCBI Datasets query failed (%s); using Entrez for all tax IDs", str(detail)[:500])

    remaining = sorted(set(tax_ids) - set(found))
    if remaining:
        logger.info("Querying Entrez for %d tax IDs not resolved by Datasets", len(remaining))
        try:
            found.update(query_entrez(remaining, email))
        except Exception as exc:  # network/XML errors: report, let strict mode decide
            logger.error("NCBI Entrez query failed: %s", exc)

    wanted = set(tax_ids)
    return {tax_id: lineage for tax_id, lineage in found.items() if tax_id in wanted}


# --------------------------------------------------------------------------
# Mapping
# --------------------------------------------------------------------------
MAPPING_COLUMNS = ["acc", "ncbi_id", *RANKS, "tax_status", "parsed_taxon_string"]


def build_mapping(manifest: pl.DataFrame, lineages: Dict[int, Lineage]) -> pl.DataFrame:
    rows = []
    missing = []
    for genome_acc, tax_id in manifest.select("genome_acc", "tax_id").iter_rows():
        lineage = lineages.get(tax_id)
        if lineage is None:
            missing.append((genome_acc, tax_id))
            continue
        ranks = lineage.ranks
        row = {"acc": genome_acc, "ncbi_id": tax_id}
        row.update({rank: ranks.get(rank, "") for rank in RANKS})
        row["tax_status"] = tax_status(ranks)
        row["parsed_taxon_string"] = format_lineage(ranks)
        rows.append(row)
    if missing:
        logger.warning("No taxonomy found for %d genomes (e.g. %s)", len(missing), missing[:5])
    schema = {"acc": pl.Utf8, "ncbi_id": pl.Int64, **{r: pl.Utf8 for r in RANKS},
              "tax_status": pl.Utf8, "parsed_taxon_string": pl.Utf8}
    return pl.DataFrame(rows, schema=schema)


def load_manifest(path: str) -> pl.DataFrame:
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)  # all strings
    missing = {"genome_acc", "tax_id"} - set(frame.columns)
    if missing:
        raise ValueError(f"Manifest {path} is missing columns: {sorted(missing)}")
    frame = frame.select(
        pl.col("genome_acc").str.strip_chars(),
        pl.col("tax_id").str.strip_chars().cast(pl.Int64, strict=False),
    )
    dropped = frame.filter(pl.col("tax_id").is_null()).height
    if dropped:
        logger.warning("Ignoring %d manifest rows without a numeric tax ID", dropped)
    return frame.filter(pl.col("tax_id").is_not_null()).unique(subset=["genome_acc"], keep="first")


def main() -> None:
    manifest = load_manifest(str(snakemake.input["unified_manifest"]))
    email = str(snakemake.params["email"])
    strict = bool(snakemake.params.get("strict", False))

    tax_ids = sorted(manifest["tax_id"].unique().to_list())
    logger.info("Found %d genomes and %d unique tax IDs", manifest.height, len(tax_ids))

    lineages = resolve_lineages(tax_ids, email)
    unresolved = sorted(set(tax_ids) - set(lineages))
    if unresolved:
        message = f"{len(unresolved)} tax IDs could not be resolved (e.g. {unresolved[:10]})"
        if strict:
            raise ValueError(message)
        logger.warning(message)
    if not lineages:
        raise ValueError("No taxonomy records retrieved; cannot proceed")

    mapping = build_mapping(manifest, index_duplicate_species(lineages))
    if mapping.height == 0:
        raise ValueError("No taxonomy records matched the manifest; cannot proceed")

    output = Path(str(snakemake.output["taxonomy"]))
    output.parent.mkdir(parents=True, exist_ok=True)
    mapping.write_csv(output, separator="\t")
    counts = dict(mapping["tax_status"].value_counts().iter_rows())
    logger.info("Wrote %d rows to %s (tax_status: %s)", mapping.height, output, counts)


if __name__ == "__main__":
    main()
