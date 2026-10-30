"""Which source a genome came from (fungi, algae, p10k, ncbi, tsa, external, ...).

The BUSCO quality-control manifests (metapackage and external) list ``group`` and ``genome_acc`` for every genome;
``group`` is the download source (``tsa`` = NCBI transcriptome shotgun assembly, ``external`` = the tree workflow's
own genomes). Used for the ``source_group`` column of the tree metadata and for the by-source breakdowns of the
tree notebooks, so transcriptome tips can be judged separately from genome tips.
"""
import csv
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence

UNKNOWN = "unknown"


def read_source_groups(manifests: Iterable) -> Dict[str, str]:
    """``{genome_acc: group}`` from QC manifests (tab-separated, columns ``group`` and ``genome_acc``).

    Missing files are skipped; for a genome in several manifests the first one wins.
    """
    groups: Dict[str, str] = {}
    for manifest in manifests:
        path = Path(str(manifest))
        if not path.is_file():
            continue
        with open(path, newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                acc, group = row.get("genome_acc"), row.get("group")
                if acc and group and acc not in groups:
                    groups[acc] = group
    return groups


def group_of(accession: str, groups: Mapping[str, str]) -> str:
    return groups.get(accession, UNKNOWN)


def by_source(tips: Sequence[str], groups: Mapping[str, str], flagged: Iterable[str] = ()) -> list:
    """Per source group: tips in the tree, how many are in ``flagged`` and the flagged share.

    Rows are sorted by group name; the share is flagged / tips.
    """
    flagged = set(flagged)
    total, bad = Counter(), Counter()
    for tip in tips:
        group = group_of(tip, groups)
        total[group] += 1
        bad[group] += tip in flagged
    return [{"group": g, "tips": total[g], "flagged": bad[g], "share": bad[g] / total[g]} for g in sorted(total)]
