"""Read (and shrink) the result directory of one BUSCO run.

A BUSCO run directory (``outputs/busco/<group>/<acc>/run_<lineage>/``) written by BUSCO itself holds about 540 files and
170 MB, of which a genome's own miniprot index (``miniprot_output/ref.mpi``, 167 MB) is the bulk. Downstream steps use only:

* ``full_table.tsv`` (and the ``short_summary*`` files);
* per hit, its protein (``.faa``), its gene model (``.gff``, one ``CDS`` line per exon) and, in transcriptome mode,
  its coding sequence (``.fna``).

``condense_run`` keeps exactly that, in three files per run, and deletes the rest (also the genome folder's ``tmp/`` and
all logs but ``busco.log``):

* ``proteins.faa.gz``  every ``busco_sequences/{single,multi}_copy_busco_sequences/<id>.faa``, concatenated;
* ``models.gff.gz``    the same for ``.gff``;
* ``cds.fna.gz``       the same for ``.fna`` (only when there are any).

Nothing is lost by concatenating: each FASTA header is ``><busco id>_<n>:<hash>|<contig>:<start>-<end>|<strand>`` and
each GFF line carries ``Target=<busco id>_<n>:<hash> ...``, so a record is assigned to its BUSCO gene by the text before
the first ``_``. Whether a gene was ``Complete`` or ``Duplicated`` is in ``full_table.tsv``.

``BuscoRun`` reads either layout (the original per-gene files or the condensed ones), so code written against it works
on old and new results alike.
"""
import gzip
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SINGLE, MULTI = "single_copy_busco_sequences", "multi_copy_busco_sequences"
FILES = {"proteins": "proteins.faa.gz", "models": "models.gff.gz", "cds": "cds.fna.gz"}
EXT = {"proteins": ".faa", "models": ".gff", "cds": ".fna"}
# directories BUSCO writes that nothing downstream reads (ref.mpi alone is ~167 MB per genome)
DISPOSABLE = ("busco_sequences", "miniprot_output", "hmmer_output", "metaeuk_output", "blast_output",
              "augustus_output", "prodigal_output", "tmp")
_TARGET = re.compile(r"Target=([^_;\s]+)_")


class Rec:
    """A FASTA record: ``id`` (first word of the header), ``description`` (whole header) and ``seq`` (str)."""
    __slots__ = ("id", "description", "seq")

    def __init__(self, header: str, seq: str):
        self.description = header
        self.id = header.split(None, 1)[0] if header.strip() else ""
        self.seq = seq

    def __repr__(self):
        return f"Rec({self.description!r}, {len(self.seq)} aa/nt)"


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_fasta_lines(lines) -> List[Rec]:
    records, header, chunks = [], None, []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                records.append(Rec(header, "".join(chunks)))
            header, chunks = line[1:], []
        elif header is not None:
            chunks.append(line.strip())
    if header is not None:
        records.append(Rec(header, "".join(chunks)))
    return records


def busco_id_of_header(header: str) -> str:
    return header.split("_", 1)[0]


class BuscoRun:
    """The sequences and gene models of one BUSCO run, whichever layout it has (``layout``: condensed/legacy/none)."""

    def __init__(self, run_dir):
        self.dir = Path(run_dir)
        if any((self.dir / name).exists() for name in FILES.values()):
            self.layout = "condensed"
        elif (self.dir / "busco_sequences").is_dir():
            self.layout = "legacy"
        else:
            self.layout = "none"
        self._cache: Dict[str, Dict[str, List]] = {}

    # -- loading ------------------------------------------------------------------------------------------------
    def _legacy_files(self, kind: str, busco_id: Optional[str] = None):
        for sub in (SINGLE, MULTI):
            folder = self.dir / "busco_sequences" / sub
            if not folder.is_dir():
                continue
            if busco_id is not None:
                path = folder / f"{busco_id}{EXT[kind]}"
                if path.exists():
                    yield path
            else:
                yield from sorted(folder.glob(f"*{EXT[kind]}"))

    def _load(self, kind: str) -> Dict[str, List]:
        """``{busco id: [records or gff lines]}`` for ``kind`` (cached for this instance)."""
        if kind in self._cache:
            return self._cache[kind]
        grouped: Dict[str, List] = defaultdict(list)
        if self.layout == "condensed":
            path = self.dir / FILES[kind]
            sources = [path] if path.exists() else []
        elif self.layout == "legacy":
            sources = list(self._legacy_files(kind))
        else:
            sources = []
        for path in sources:
            with _open_text(path) as handle:
                if kind == "models":
                    for line in handle:
                        if line.startswith("#") or not line.strip():
                            continue
                        match = _TARGET.search(line)
                        if self.layout == "legacy":
                            grouped[path.name[: -len(EXT[kind])]].append(line.rstrip("\n"))   # the file name is the gene
                        elif match:
                            grouped[match.group(1)].append(line.rstrip("\n"))
                else:
                    for record in parse_fasta_lines(handle):
                        key = path.name[: -len(EXT[kind])] if self.layout == "legacy" else busco_id_of_header(record.description)
                        grouped[key].append(record)
        self._cache[kind] = dict(grouped)
        return self._cache[kind]

    # -- queries ------------------------------------------------------------------------------------------------
    def protein_records(self, busco_id: str) -> List[Rec]:
        """Every protein BUSCO kept for the gene (one per copy of a Duplicated gene)."""
        return list(self._load("proteins").get(busco_id, []))

    def cds_records(self, busco_id: str) -> List[Rec]:
        """Coding sequences (transcriptome mode only)."""
        return list(self._load("cds").get(busco_id, []))

    def gff_lines(self, busco_id: str) -> List[str]:
        return list(self._load("models").get(busco_id, []))

    def gff_models(self, busco_id: str):
        """Gene models (``intron_utils.Model``) of the gene; empty when there is no model."""
        from intron_utils import parse_gff_text
        return parse_gff_text(self.gff_lines(busco_id))

    def has_models(self) -> bool:
        return bool(self._load("models"))


def _tag_targets(gff_text: str, busco_id: str) -> str:
    """Make sure every GFF line names its gene (``Target=<busco id>_...``): BUSCO's lines do, so this is a no-op for them."""
    out = []
    for line in gff_text.splitlines():
        if line.strip() and not line.startswith("#") and not _TARGET.search(line):
            line = f"{line};Target={busco_id}_0" if len(line.split("\t")) >= 9 else line
        out.append(line)
    return "\n".join(out) + ("\n" if out else "")


def _tidy_genome_dir(genome_dir: Path, stats: Dict[str, int]) -> None:
    """The folder above ``run_*``: drop ``tmp/`` (a copy of the lineage proteins) and every log but ``busco.log``.

    The tool logs (miniprot alignment output alone is ~120 MB for a large genome) are only useful while debugging a run.
    """
    def drop(path: Path):
        try:
            stats["removed_bytes"] += path.stat().st_size
            stats["removed_files"] += 1
            path.unlink()
        except OSError:
            pass

    tmp = genome_dir / "tmp"
    if tmp.is_dir():
        for sub, _, files in os.walk(tmp):
            for f in files:
                drop(Path(sub) / f)
        shutil.rmtree(tmp, ignore_errors=True)
    logs = genome_dir / "logs"
    if logs.is_dir():
        for path in logs.iterdir():
            if path.is_file() and path.name != "busco.log":
                drop(path)


def condense_run(run_dir, remove: bool = True) -> Dict[str, int]:
    """Concatenate the per-gene files into the three condensed files and (``remove``) delete everything else big.

    Safe to repeat: an already condensed run is left alone. The condensed files are written under a temporary name and
    renamed first, so an interrupted call never leaves a run without its sequences. Returns file counts.
    """
    run_dir = Path(run_dir)
    run = BuscoRun(run_dir)
    stats = {"genes": 0, "removed_files": 0, "removed_bytes": 0}
    if run.layout == "legacy":
        for kind, name in FILES.items():
            files = list(run._legacy_files(kind))
            if not files:
                continue
            if kind == "proteins":
                stats["genes"] = len(files)
            tmp = run_dir / (name + ".tmp")
            with gzip.open(tmp, "wt", compresslevel=6) as out:
                for path in files:
                    text = path.read_text()
                    if kind == "models":
                        text = _tag_targets(text, path.name[: -len(EXT[kind])])
                    out.write(text if text.endswith("\n") or not text else text + "\n")
            os.replace(tmp, run_dir / name)
    if remove:
        _tidy_genome_dir(run_dir.parent, stats)
        for name in DISPOSABLE:
            path = run_dir / name
            if path.is_dir():
                for sub, _, files in os.walk(path):
                    for f in files:
                        stats["removed_files"] += 1
                        try:
                            stats["removed_bytes"] += os.path.getsize(os.path.join(sub, f))
                        except OSError:
                            pass
                shutil.rmtree(path, ignore_errors=True)
    return stats
