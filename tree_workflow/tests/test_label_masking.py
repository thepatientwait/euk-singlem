import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import label_masking as lm  # noqa: E402

CUTOFFS = {"phylum": 0.5, "class": 0.6, "order": 0.7, "family": 0.8, "genus": 0.9, "species": 0.95}


def full(p, c, o, f, g, s):
    return f"d__Eukarya;p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}"


def clade(prefix, p, n):
    """n leaves of one genus, all sharing every rank (a pure clade)."""
    return {f"{prefix}{i}": full(p, p + "c", p + "o", p + "f", p + "g", f"{p}g sp{i}") for i in range(n)}


def newick(groups):
    parts = []
    for prefix, n in groups:
        leaves = ",".join(f"{prefix}{i}|RED=1:0.1" for i in range(n))
        parts.append(f"({leaves})|RED=0.97:0.1")
    return "(" + ",".join(parts) + ")|RED=0.2;"


class TestMasking(unittest.TestCase):
    def test_pure_clades_recover_every_rank(self):
        known = {**clade("a", "A", 20), **clade("b", "B", 20)}
        table = lm.run(newick([("a", 20), ("b", 20)]), known, CUTOFFS, fraction=0.2, replicates=2)
        for row in table.iter_rows(named=True):
            self.assertEqual(row["correct"], row["tested"], row)
            self.assertEqual((row["novel"], row["wrong"]), (0, 0))

    def test_mislabelled_neighbours_cause_wrong_names(self):
        # the A clade is half labelled B at phylum: hidden A genomes get either name, never a placeholder
        known = {**clade("a", "A", 20), **clade("b", "B", 20)}
        for i in range(10):
            known[f"a{i}"] = known[f"b{i}"]
        table = lm.run(newick([("a", 20), ("b", 20)]), known, CUTOFFS, fraction=0.3, replicates=3)
        phylum = table.filter(table["rank"] == "phylum")
        self.assertGreater(phylum["wrong"].sum() + phylum["correct"].sum(), 0)
        self.assertEqual(phylum["tested"].sum(), phylum["correct"].sum() + phylum["wrong"].sum() + phylum["novel"].sum())

    def test_counts_add_up_and_are_seeded(self):
        known = {**clade("a", "A", 15), **clade("b", "B", 15)}
        nwk = newick([("a", 15), ("b", 15)])
        first, second = lm.run(nwk, known, CUTOFFS, seed=5), lm.run(nwk, known, CUTOFFS, seed=5)
        self.assertTrue(first.equals(second))
        for r in first.iter_rows(named=True):
            self.assertEqual(r["tested"], r["correct"] + r["novel"] + r["wrong"])


if __name__ == "__main__":
    unittest.main()
