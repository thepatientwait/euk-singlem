"""Tests for scripts/tree_diagnostics_lib.py (synthetic data, no network)."""
import random
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import tree_diagnostics_lib as lib  # noqa: E402


class TestSeparation(unittest.TestCase):
    def test_extremes_and_ties(self):
        self.assertEqual(lib.separation([0.1, 0.2], [0.5, 0.6]), 1.0)
        self.assertEqual(lib.separation([0.5, 0.6], [0.1, 0.2]), 0.0)
        self.assertEqual(lib.separation([0.3, 0.3], [0.3, 0.3]), 0.5)
        self.assertTrue(np.isnan(lib.separation([], [1.0])))

    def test_matches_brute_force(self):
        rng = random.Random(1)
        w = [rng.random() for _ in range(40)]
        b = [rng.random() + 0.2 for _ in range(50)]
        brute = np.mean([(x < y) + 0.5 * (x == y) for x in w for y in b])
        self.assertAlmostEqual(lib.separation(w, b), brute)


class TestPairs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = Path(self.tmp.name)
        # two clades of 3 genomes: clade X sequences are identical-ish, clade Y differ from X
        x = "ACDEFGHIKLMNPQRSTVWY" * 3
        y = "YWVTSRQPNMLKIHGFEDCA" * 3
        def mut(s, pos):
            s = list(s); s[pos] = "A" if s[pos] != "A" else "C"; return "".join(s)
        rows = {"x1": x, "x2": mut(x, 3), "x3": mut(x, 9), "y1": y, "y2": mut(y, 5), "y3": mut(y, 20)}
        rows["x3"] = rows["x3"][:30] + "-" * 30                      # x3 lacks gene 2
        (d / "aln.fasta").write_text("".join(f">{k}\n{v}\n" for k, v in rows.items()))
        (d / "partitions.tsv").write_text("gene\tstart\tend\tgenomes\ng1\t1\t30\t6\ng2\t31\t60\t5\n")
        self.ids, self.matrix = lib.load_alignment(d / "aln.fasta")
        self.parts = lib.read_partitions(d / "partitions.tsv")
        self.tax = {k: {"phylum": "p__X" if k[0] == "x" else "p__Y"} for k in rows}

    def test_load_and_presence(self):
        self.assertEqual(self.matrix.shape, (6, 60))
        occ, present = lib.occupancy(self.ids, self.matrix, self.parts)
        self.assertEqual(dict(zip(occ["genome"], occ["genes"]))["x3"], 1)
        self.assertEqual(present.shape, (6, 2))
        self.assertAlmostEqual(dict(zip(occ["genome"], occ["gap_fraction"]))["x3"], 0.5)

    def test_pdist(self):
        self.assertAlmostEqual(lib.pdist(self.matrix[0], self.matrix[1], 5), 1 / 60)
        self.assertTrue(np.isnan(lib.pdist(self.matrix[0], self.matrix[1], 100)))

    def test_clades_separate_and_gene_overlap(self):
        _, present = lib.occupancy(self.ids, self.matrix, self.parts)
        report = lib.rank_report(self.ids, self.matrix, present, self.parts, self.tax, ranks=["phylum"],
                                 n_pairs=60, seed=3)
        s = report["phylum"]["summary"]
        self.assertEqual(s["separation_pdist"], 1.0)
        self.assertLess(s["median_pdist_within"], s["median_pdist_between"])
        genes = lib.gene_separation(report, self.parts, "phylum")
        self.assertEqual([g for g, _, _ in genes], ["g1", "g2"])
        self.assertEqual(genes[0][1], 1.0)

    def test_sample_pairs_respects_labels(self):
        labels = ["a", "a", "b", "", "b"]
        within, between = lib.sample_pairs(labels, 30, random.Random(0))
        self.assertTrue(all(labels[i] == labels[j] != "" for i, j in within))
        self.assertTrue(all(labels[i] != labels[j] and labels[i] and labels[j] for i, j in between))


class TestMulticopy(unittest.TestCase):
    def test_prefers_peer_like_copy_when_score_disagrees(self):
        peer = "MKVLAAGIVGLLLAQPSTEEKD" * 3
        close = peer[:-3] + "AAA"
        far = "WWPPHHYYCCGGNNQQRRSSTT" * 3
        outcome = lib.evaluate_copy_choice([far, close], [90.0, 50.0], [peer, peer, close])
        self.assertEqual((outcome["chosen"], outcome["closest"], outcome["agree"]), (0, 1, False))
        self.assertGreater(outcome["similarity_gap"], 0.3)
        outcome = lib.evaluate_copy_choice([close, far], [90.0, 50.0], [peer, peer])
        self.assertTrue(outcome["agree"])
        self.assertIsNone(lib.evaluate_copy_choice([close], [1.0], [peer]))
        self.assertIsNone(lib.evaluate_copy_choice([close, far], [1.0, 0.5], []))

    def test_parse_duplicated_rows_and_sequences(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "run"
            (run / "busco_sequences" / "multi_copy_busco_sequences").mkdir(parents=True)
            table = run / "full_table.tsv"
            table.write_text(
                "# BUSCO\n# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
                "g1\tDuplicated\tc1\t10\t40\t+\t50.0\t10\n"
                "g1\tDuplicated\tc2\t5\t35\t-\t90.0\t10\n"
                "g2\tComplete\tc1\t1\t9\t+\t70.0\t3\n")
            (run / "busco_sequences" / "multi_copy_busco_sequences" / "g1.faa").write_text(
                ">g1|c1:10-40|+ x\nMMMM*\n>g1|c2:5-35|- y\nKKKK*\n")
            cases = lib.duplicated_cases(table)
            self.assertEqual(list(cases), ["g1"])
            self.assertEqual([r["Sequence"] for r in cases["g1"]], ["c2", "c1"])   # best score first
            self.assertEqual(lib.copy_sequences(table, "g1", cases["g1"]), ["KKKK", "MMMM"])


if __name__ == "__main__":
    unittest.main()
