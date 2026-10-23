"""Tests for scripts/subset_genomes.py."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import subset_genomes as sg  # noqa: E402

ROWS = [("fungi", f"f{i}") for i in range(12)] + [("ncbi", f"n{i}") for i in range(12)] + [("p10k", f"p{i}") for i in range(4)]
TAX = {}
for i in range(12):
    TAX[f"f{i}"] = ("Fam%d" % (i % 2), "complete")
    TAX[f"n{i}"] = ("Fam%d" % (i % 2), "complete")
for i in range(4):
    TAX[f"p{i}"] = ("Fam2", "partial")
TAX["f0"] = ("Fam0", "unresolved")      # unusable
TAX["n11"] = ("", "complete")            # no family


class TestSubset(unittest.TestCase):
    def test_family5(self):
        base, ext = sg.select_subset("family5", ROWS, TAX, ["e1", "e2", "e3"], per_family=5, n_external=2)
        fams = {}
        for a in base:
            fams.setdefault(TAX[a][0], []).append(a)
        self.assertEqual({f: len(v) for f, v in fams.items()}, {"Fam0": 5, "Fam1": 5, "Fam2": 4})
        self.assertNotIn("f0", base)
        self.assertNotIn("n11", base)
        self.assertEqual(len(ext), 2)

    def test_deterministic_and_order_independent(self):
        a = sg.select_subset("family5", ROWS, TAX, ["e1", "e2"], 3, 1, seed=7)
        b = sg.select_subset("family5", list(reversed(ROWS)), TAX, ["e2", "e1"], 3, 1, seed=7)
        self.assertEqual(a, b)
        self.assertNotEqual(a[0], sg.select_subset("family5", ROWS, TAX, ["e1", "e2"], 3, 1, seed=8)[0])

    def test_jgi_p10k(self):
        base, _ = sg.select_subset("jgi_p10k", ROWS, TAX, [])
        self.assertTrue(all(a[0] in "fp" for a in base))
        self.assertEqual(len(base), 11 + 4)   # f0 unusable
        base, _ = sg.select_subset("jgi_p10k_family5", ROWS, TAX, [], per_family=2)
        self.assertEqual(len(base), 6)
        self.assertTrue(all(a[0] in "fp" for a in base))

    def test_bad_mode(self):
        with self.assertRaises(ValueError):
            sg.select_subset("nope", ROWS, TAX, [])


if __name__ == "__main__":
    unittest.main()
