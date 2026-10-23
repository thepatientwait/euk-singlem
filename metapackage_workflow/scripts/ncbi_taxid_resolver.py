#!/usr/bin/env python3
"""
Resolve assembly accessions to NCBI taxonomy using the NCBI Datasets REST API,
with an Entrez fallback for replaced/suppressed accessions.

Why this exists:
- The Datasets API omits replaced/suppressed assemblies from genome dataset
  reports; Entrez esummary still returns them (often under the paired RefSeq
  accession), with a valid tax_id.
- Assembly accessions in the metadata may reference older versions than the
  current record; results are indexed by both full and version-stripped
  accession so lookups are version-insensitive.
"""

import json
import ssl
import time
import urllib.request
from typing import Dict, List, Optional, Tuple

from log_config import configure_logger

logger = configure_logger(__file__)

DATASETS_BASE = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
EUTILS_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
# The genome dataset_report POST endpoint returns at most 20 reports per call.
DATASETS_POST_BATCH = 20
RATE_LIMIT_SLEEP = 0.4


def base_accession(acc: str) -> str:
    """Strip a trailing numeric .version from an assembly accession."""
    if '.' in acc:
        base, _, version = acc.rpartition('.')
        if version.isdigit():
            return base
    return acc


def accession_core(acc: str) -> str:
    """
    Assembly identity core shared by GenBank/RefSeq pairs.

    GCA_000002825.1 and its paired GCF_000002825.2 are the same assembly
    under different databases; both reduce to '000002825'.
    """
    base = base_accession(acc)
    if base.startswith(('GCA_', 'GCF_')):
        return base.split('_', 1)[1]
    return base
    return acc


def ssl_context() -> ssl.SSLContext:
    """Build an SSL context that works on hosts without a system CA bundle."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        pass
    for path in (
        '/etc/pki/tls/certs/ca-bundle.crt',      # RHEL-family
        '/etc/ssl/certs/ca-certificates.crt',    # Debian-family
        '/etc/ssl/cert.pem',
    ):
        try:
            return ssl.create_default_context(cafile=path)
        except (OSError, ssl.SSLError):
            continue
    return ssl.create_default_context()


_CTX = None


def get_context() -> ssl.SSLContext:
    global _CTX
    if _CTX is None:
        _CTX = ssl_context()
    return _CTX


def _http_json(url: str, payload: Optional[dict] = None, timeout: int = 60):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(
        url,
        data=data,
        headers={'Content-Type': 'application/json'} if data else {},
        method='POST' if data else 'GET',
    )
    with urllib.request.urlopen(req, context=get_context(), timeout=timeout) as response:
        return json.loads(response.read().decode())


def datasets_accs_to_taxids(accessions: List[str]) -> Dict[str, Tuple[str, str, str]]:
    """
    Resolve assembly accessions via the Datasets genome dataset_report API.

    Returns dict: accession -> (tax_id, organism_name, current_accession).
    Only accessions with a report are included.
    """
    resolved: Dict[str, Tuple[str, str, str]] = {}

    for start in range(0, len(accessions), DATASETS_POST_BATCH):
        chunk = accessions[start:start + DATASETS_POST_BATCH]
        try:
            data = _http_json(
                f'{DATASETS_BASE}/genome/dataset_report',
                payload={'accessions': chunk},
            )
        except Exception as exc:
            logger.warning('Datasets genome query failed for chunk starting %s: %s',
                           chunk[0] if chunk else '?', exc)
            time.sleep(RATE_LIMIT_SLEEP * 2)
            continue

        for report in data.get('reports', []):
            acc = report.get('accession')
            organism = report.get('organism', {}) or {}
            tax_id = organism.get('tax_id')
            if acc and tax_id:
                resolved[acc] = (
                    str(tax_id),
                    organism.get('organism_name') or '',
                    report.get('current_accession') or acc,
                )

        time.sleep(RATE_LIMIT_SLEEP)

    return resolved


def entrez_accs_to_taxids(accessions: List[str], email: str,
                          chunk_size: int = 20) -> Dict[str, Tuple[str, str, str]]:
    """
    Resolve assembly accessions via Entrez esearch + esummary.

    Needed for replaced/suppressed assemblies that the Datasets API omits.
    esummary may report the paired RefSeq accession instead of the requested
    GenBank one, so each hit is attributed to the requested accession that
    produced it.

    Returns dict: accession -> (tax_id, organism_name, current_accession).
    """
    from Bio import Entrez

    resolved: Dict[str, Tuple[str, str, str]] = {}
    Entrez.email = email

    for start in range(0, len(accessions), chunk_size):
        chunk = accessions[start:start + chunk_size]

        for acc in chunk:
            try:
                handle = Entrez.esearch(
                    db='assembly', term=f'{acc}[Assembly Accession]', retmax=5)
                uids = Entrez.read(handle).get('IdList', [])
                handle.close()
            except Exception as exc:
                logger.warning('Entrez esearch failed for %s: %s', acc, exc)
                time.sleep(RATE_LIMIT_SLEEP)
                continue

            if not uids:
                continue

            try:
                time.sleep(RATE_LIMIT_SLEEP)
                handle = Entrez.esummary(db='assembly', id=','.join(uids))
                docs = Entrez.read(handle).get('DocumentSummarySet', {}).get('DocumentSummary', [])
                handle.close()
            except Exception as exc:
                logger.warning('Entrez esummary failed for %s: %s', acc, exc)
                continue

            core = accession_core(acc)
            best = None
            for doc in docs:
                reported = str(doc.get('AssemblyAccession', '') or '')
                # Accept the exact accession, another version of it, or its
                # paired GenBank/RefSeq record (replaced assemblies often only
                # surface under the paired accession).
                if (reported == acc
                        or base_accession(reported) == base_accession(acc)
                        or accession_core(reported) == core):
                    best = doc
                    break
            if best is None:
                logger.info('Entrez: no matching assembly summary for %s (hits: %s)',
                            acc, [d.get('AssemblyAccession') for d in docs])
                continue

            tax_id = str(best.get('Taxid', '') or '')
            if tax_id and tax_id != '0':
                resolved[acc] = (
                    tax_id,
                    str(best.get('Organism', '') or ''),
                    str(best.get('AssemblyAccession', '') or acc),
                )

    return resolved


def resolve_accessions_to_taxids(accessions: List[str], email: str) -> Dict[str, Tuple[str, str, str]]:
    """
    Resolve accessions with the Datasets API, then Entrez for the remainder.

    Returns dict: requested accession -> (tax_id, organism_name, current_accession).
    """
    accessions = sorted({a for a in accessions if a})
    logger.info('Resolving %d accessions to tax IDs (datasets API)...', len(accessions))

    resolved = datasets_accs_to_taxids(accessions)
    logger.info('Datasets API resolved %d/%d accessions', len(resolved), len(accessions))

    missing = [acc for acc in accessions if acc not in resolved]
    if missing:
        logger.info('Falling back to Entrez for %d unresolved (replaced/suppressed?) accessions...',
                    len(missing))
        entrez_hits = entrez_accs_to_taxids(missing, email)
        resolved.update(entrez_hits)
        logger.info('Entrez resolved %d of the remaining %d; total %d/%d',
                    len(entrez_hits), len(missing), len(resolved), len(accessions))

    # Index by version-stripped accession so callers can look up older/newer
    # versions of the same assembly.
    indexed: Dict[str, Tuple[str, str, str]] = {}
    for acc, value in resolved.items():
        indexed[acc] = value
        indexed.setdefault(base_accession(acc), value)

    return indexed


def lookup_taxid(acc_to_taxid: Dict[str, Tuple[str, str, str]], acc: str):
    """Version-insensitive lookup; returns (tax_id, name, current_acc) or None."""
    return acc_to_taxid.get(acc) or acc_to_taxid.get(base_accession(acc))
