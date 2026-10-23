"""Multicopy BUSCO handling: one top-scoring sequence per gene, QC counts genes (no network)."""
import builtins
import sys
import tempfile
import types
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
builtins.snakemake = types.SimpleNamespace()
import aggregate_busco_sequences as ab  # noqa: E402
import busco_results_quality_control as qc  # noqa: E402

HEADER = "# BUSCO version is: 5.8.2\n# The lineage dataset is: eukaryota_odb12\n"
COLUMNS = "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
ROWS = [
    ("g1", "Complete", "c1", 1, 30, "+", 100.0, 10),
    ("g2", "Duplicated", "c1", 100, 160, "+", 50.0, 20),
    ("g2", "Duplicated", "c2", 5, 95, "+", 80.0, 30),   # best score
    ("g2", "Duplicated", "c3", 5, 65, "-", 80.0, 20),   # same score, shorter: loses the tie
    ("g3", "Fragmented", "c1", 200, 230, "+", 10.0, 10),
    ("g4", "Missing", "", "", "", "", "", ""),
]


def write_table(root):
    path = Path(root) / "busco" / "grp" / "GCA_1" / "run_x" / "full_table.tsv"
    path.parent.mkdir(parents=True)
    path.write_text(HEADER + COLUMNS + "".join("\t".join(str(v) for v in r) + "\n" for r in ROWS))
    return path


class TestTopHit(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.table = write_table(self.tmp.name)

    def selected(self, multi):
        _, acc, df = ab.parse_busco_file(self.table, multi)
        self.assertEqual(acc, "GCA_1")
        return {r["Busco id"]: r for r in df.iter_rows(named=True)}

    def test_flag_off_drops_duplicated_genes(self):
        self.assertEqual(set(self.selected(False)), {"g1"})

    def test_flag_on_keeps_one_top_copy_per_gene(self):
        sel = self.selected(True)
        self.assertEqual(set(sel), {"g1", "g2"})
        self.assertEqual((sel["g2"]["Sequence"], sel["g2"]["Gene Start"]), ("c2", 5))

    def test_qc_counts_genes_not_rows(self):
        self.assertEqual(qc._count_complete_markers(self.table, False), 1)
        self.assertEqual(qc._count_complete_markers(self.table, True), 2)  # g2 has 3 rows but counts once

    def test_qc_agrees_with_aggregation(self):
        for multi in (False, True):
            self.assertEqual(qc._count_complete_markers(self.table, multi), len(self.selected(multi)))


class TestStatusSidecar(unittest.TestCase):
    def test_sidecar_lists_status_per_sequence(self):
        seqs = {
            "A|g|1": {"sequence": ">x\nMK", "sequence_id": "A|g|1", "genome_acc": "A", "busco_id": "g",
                      "status": "Complete", "score": 100.0},
            "B|g|1": {"sequence": ">x\nMK", "sequence_id": "B|g|1", "genome_acc": "B", "busco_id": "g",
                      "status": "Duplicated", "score": 80.0},
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = ab.write_status_file("g", seqs, Path(tmp))
            lines = path.read_text().splitlines()
        self.assertEqual(lines[0], "sequence_id\tstatus\tscore")
        self.assertEqual(lines[1:], ["A|g|1\tComplete\t100.0", "B|g|1\tDuplicated\t80.0"])


if __name__ == "__main__":
    unittest.main()
