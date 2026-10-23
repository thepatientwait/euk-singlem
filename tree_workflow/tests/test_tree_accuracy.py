"""Tests for scripts/tree_accuracy_lib.py (small synthetic trees)."""
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import Phylo

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import tree_accuracy_lib as lib  # noqa: E402


def tree(nwk):
    return Phylo.read(io.StringIO(nwk), "newick")


def tax(**kw):
    return {k: {"order": v[0], "family": v[1]} for k, v in kw.items()}


class TestFragmentation(unittest.TestCase):
    NWK = "(((a1:1,a2:1):1,(b1:1,b2:1):1):1,(a3:1,(c1:1,(c2:1,c3:1):1):1):1);"
    TAX = tax(a1=("A", "f"), a2=("A", "f"), a3=("A", "f"), b1=("B", "f"), b2=("B", "f"),
              c1=("C", "f"), c2=("C", "f"), c3=("C", "f"))

    def test_extra_piece_found(self):
        pieces = lib.pieces_unrooted(tree(self.NWK), self.TAX, "order")
        self.assertEqual({n: len(p) for n, p in pieces.items()}, {"A": 2, "B": 1, "C": 1})
        self.assertEqual([len(t) for _, t in pieces["A"]], [2, 1])

    def test_independent_of_root(self):
        # same unrooted topology written with a different root
        other = "((a1:1,a2:1):1,(b1:1,b2:1):1,(a3:1,(c1:1,(c2:1,c3:1):1):1):1);"
        a = lib.pieces_unrooted(tree(self.NWK), self.TAX, "order")
        b = lib.pieces_unrooted(tree(other), self.TAX, "order")
        self.assertEqual({n: sorted(len(t) for _, t in p) for n, p in a.items()},
                         {n: sorted(len(t) for _, t in p) for n, p in b.items()})

    def test_midpoint_style_root_does_not_create_fragments(self):
        # A is monophyletic in the unrooted sense but a rooting inside it would split it
        nwk = "((a1:1,a2:1):1,(a3:1,(b1:1,b2:1):1):1);"
        taxonomy = tax(a1=("A", "f"), a2=("A", "f"), a3=("A", "f"), b1=("B", "f"), b2=("B", "f"))
        pieces = lib.pieces_unrooted(tree(nwk), taxonomy, "order")
        self.assertEqual(len(pieces["A"]), 1)
        self.assertEqual(len(pieces["B"]), 1)

    def test_summary(self):
        summary, names = lib.fragmentation(tree(self.NWK), self.TAX, ranks=["order"])
        self.assertEqual(summary[0]["extra_pieces"], 1)
        self.assertEqual(summary[0]["fragmented_names"], 1)
        self.assertEqual(summary[0]["tips_in_fragmented"], 3)
        self.assertEqual({n["name"] for n in names}, {"A", "B", "C"})

    def test_extra_piece_support(self):
        t = tree("(((a1:1,a2:1)0.9:1,(b1:1,b2:1)0.8:1)0.7:1,(a3:1,(c1:1,(c2:1,c3:1)0.95:1)0.5:1)0.6:1);")
        pieces = lib.pieces_unrooted(t, self.TAX, "order")
        rows = lib.extra_piece_supports(pieces)
        self.assertEqual([(n, k) for n, k, _ in rows], [("A", 1)])


class TestMisplaced(unittest.TestCase):
    def test_flags_non_sister_but_not_sister_swaps(self):
        # clade X (order OX, class K1) surrounds a tip labelled order OY; OY in class K2 => flagged;
        # tip labelled OS (class K1, sister of OX) inside OX => not flagged
        inner = ",".join(f"x{i}:1" for i in range(10))
        nwk = f"(({inner},y1:1,s1:1):1,(z1:1,z2:1):1);"
        taxonomy = {f"x{i}": {"order": "OX", "class": "K1"} for i in range(10)}
        taxonomy["y1"] = {"order": "OY", "class": "K2"}
        taxonomy["s1"] = {"order": "OS", "class": "K1"}
        taxonomy["z1"] = {"order": "OZ", "class": "K3"}
        taxonomy["z2"] = {"order": "OZ", "class": "K3"}
        rows = {r["tip"]: r for r in lib.misplaced_tips(tree(nwk), taxonomy, "order", min_neighbours=8, frac=0.8)}
        self.assertFalse(rows["y1"]["sister"])
        self.assertTrue(rows["s1"]["sister"])
        self.assertEqual(rows["y1"]["neighbours_say"], "OX")


class TestRF(unittest.TestCase):
    def test_identical_and_different(self):
        a = tree("(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
        same = tree("((a:1,b:1):1,(c:1,d:1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
        diff = tree("(((a:1,c:1):1,(b:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
        self.assertEqual(lib.rf_distance(a, same)[0], 0)
        rf, norm, n = lib.rf_distance(a, diff)
        self.assertEqual((rf, n), (4, 8))
        self.assertGreater(norm, 0)

    def test_uses_common_tips_only(self):
        a = tree("(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);")
        b = tree("(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,x:1):1):1);")
        self.assertEqual(lib.rf_distance(a, b)[2], 7)


class TestRed(unittest.TestCase):
    def test_monotonicity_and_negative_lengths(self):
        t = tree("(((a|RED=1.000:0.5,b|RED=1.000:0.5)0.9|RED=0.30:-0.1,c|RED=1.000:0.4)0.8|RED=0.50:0.2,d|RED=1.000:1)0.7|RED=0.20:0;")
        n, lower, neg, worst = lib.red_monotonicity(t)
        self.assertEqual((n, lower, neg), (3, 1, 1))
        self.assertAlmostEqual(worst, -0.1)

    def test_red_by_rank_and_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tsv = Path(tmp) / "r.tsv"
            tsv.write_text("Taxa\tGTDB taxonomy\tMedian distance\tx\ty\tz\n"
                           "p__A\td__E\t0.3\t0\tphylum\tOK\np__B\td__E\t0.4\t0\tphylum\tOK\n"
                           "c__C\tp__A\t0.6\t0\tclass\tOK\nc__D\tp__A\t0.8\t0\tclass\tOK\n")
            cut = Path(tmp) / "c.dict"
            cut.write_text(json.dumps({"phylum": 0.35, "class": 0.7}))
            values = lib.red_by_rank(tsv, cut)
            self.assertEqual(values["phylum"], (0.35, [0.3, 0.4]))
            overlap = lib.rank_overlap(values)
            self.assertEqual(overlap["class"][0], 0.0)    # none below the phylum median (0.35)
            self.assertEqual(overlap["phylum"][1], 0.0)   # none above the class median (0.7)


class TestBusco(unittest.TestCase):
    def test_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            f = Path(tmp) / "full_table.tsv"
            f.write_text("# x\n# Busco id\tStatus\n1\tComplete\n2\tDuplicated\n3\tDuplicated\n4\tMissing\n")
            self.assertEqual(dict(lib.busco_status_counts(f)), {"Complete": 1, "Duplicated": 2, "Missing": 1})


if __name__ == "__main__":
    unittest.main()
