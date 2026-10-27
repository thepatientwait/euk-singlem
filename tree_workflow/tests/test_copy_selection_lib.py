"""copy_selection_lib on a small synthetic multicopy_selection.tsv."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import copy_selection_lib as lib  # noqa: E402

HEADER = "busco_id\tgenome_acc\tn_candidates\tchosen_copy\tchosen_identity\ttop_score_copy\ttop_score_identity\tagrees_with_top_score\tpeer_level\tn_peers\n"
ROWS = [
    "g1\tA\t2\t2\t0.80\t1\t0.50\t0\t6\t8",   # switched, gain .30
    "g1\tB\t2\t1\t0.70\t1\t0.70\t1\t6\t8",   # agrees
    "g2\tA\t3\t1\t0.30\t1\t0.30\t1\t5\t2",   # weak: 2 peers, low identity
    "g2\tC\t2\t1\t\t1\t\t1\t0\t0",           # fallback, no peers
]


class TestLib(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = Path(tmp.name) / "multicopy_selection.tsv"
        self.path.write_text(HEADER + "\n".join(ROWS) + "\n")
        self.frame = lib.load_selection(self.path, "t")

    def test_load_flags(self):
        self.assertEqual(self.frame["decided_by_peers"].to_list(), [True, True, True, False])
        self.assertEqual(self.frame["switched"].to_list(), [True, False, False, False])

    def test_overview(self):
        row = lib.overview(self.frame).row(0, named=True)
        self.assertEqual((row["decisions"], row["by_peers"], row["switched"]), (4, 3, 1))
        self.assertAlmostEqual(row["share_switched"], 1 / 3)

    def test_gains_only_switched(self):
        g = lib.gains(self.frame)
        self.assertEqual(g.height, 1)
        self.assertAlmostEqual(g["gain"][0], 0.30)

    def test_peer_level_table_names_levels(self):
        t = lib.by_peer_level(self.frame)
        self.assertIn("genus", t["level"].to_list())
        self.assertIn("no peers (top score)", t["level"].to_list())

    def test_weak_choices(self):
        w = lib.weak_choices(self.frame)
        self.assertEqual(w["decided by peers"], 3)
        self.assertEqual(w["fewer than 3 peers"], 1)
        self.assertEqual(w["best copy identity < 0.4"], 1)

    def test_load_all_reports_missing(self):
        frame, missing = lib.load_all([{"label": "x", "path": str(self.path)}, {"label": "y", "path": "/nonexistent"}])
        self.assertEqual(frame.height, 4)
        self.assertEqual(missing, ["y"])
        self.assertEqual(lib.load_all([{"label": "y", "path": "/nonexistent"}]), (None, ["y"]))

    def test_phylum(self):
        f = lib.attach_phylum(self.frame, {"A": "d__E;p__Ascomycota;c__x", "B": "d__E"})
        self.assertEqual(dict(zip(f["genome_acc"], f["phylum"]))["A"], "Ascomycota")
        self.assertEqual(dict(zip(f["genome_acc"], f["phylum"]))["B"], "unknown")

    def test_gene_and_genome_tables_respect_minimum(self):
        self.assertEqual(lib.gene_table(self.frame, min_decisions=2).height, 1)   # only g1 has 2 peer decisions
        self.assertEqual(lib.genome_table(self.frame, min_decisions=2)["genome_acc"].to_list(), ["A"])


if __name__ == "__main__":
    unittest.main()
