import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import notebook_narrative as nn  # noqa: E402


class TestRule(unittest.TestCase):
    def test_higher_is_worse(self):
        r = nn.Rule(look=0.1, bad=0.3)
        self.assertEqual([r.status(x) for x in (0.0, 0.1, 0.2, 0.3, 0.9)], [nn.OK, nn.LOOK, nn.LOOK, nn.BAD, nn.BAD])

    def test_lower_is_worse(self):
        r = nn.Rule(look=0.8, bad=0.5, higher_is_worse=False)
        self.assertEqual([r.status(x) for x in (0.95, 0.8, 0.6, 0.5, 0.1)], [nn.OK, nn.LOOK, nn.LOOK, nn.BAD, nn.BAD])

    def test_missing_value_is_a_look(self):
        self.assertEqual(nn.Rule(1, 2).status(None), nn.LOOK)
        self.assertEqual(nn.Rule(1, 2).status(float("nan")), nn.LOOK)


class TestFindings(unittest.TestCase):
    def test_overall_is_worst(self):
        f = nn.Findings()
        f.add(nn.OK, "fine")
        f.add(nn.BAD, "broken", "3 of 4", "fix it")
        self.assertEqual(f.overall(), nn.BAD)
        self.assertIn("broken", f.actions())
        self.assertNotIn("fine", f.actions())

    def test_actions_when_all_ok(self):
        f = nn.Findings()
        f.add(nn.OK, "fine")
        self.assertEqual(f.actions(), "Nothing needs attention.")
        self.assertEqual(nn.Findings().overall(), nn.OK)

    def test_actions_sorted_worst_first_and_pipes_escaped(self):
        f = nn.Findings()
        f.add(nn.LOOK, "minor")
        f.add(nn.BAD, "major | bad")
        lines = f.actions().splitlines()
        self.assertIn("major / bad", lines[2])
        self.assertIn("minor", lines[3])

    def test_formatting(self):
        self.assertEqual(nn.pct(0.314), "31%")
        self.assertEqual(nn.pct(None), "n/a")
        self.assertEqual(nn.count(1, "genome"), "1 genome")
        self.assertEqual(nn.count(2, "genome"), "2 genomes")
        self.assertEqual(nn.count(3, "family", "families"), "3 families")

    def test_status_table(self):
        t = nn.status_table([{"status": nn.OK, "v": "a"}], ["v"])
        self.assertTrue(t.splitlines()[2].endswith("| a |"))


if __name__ == "__main__":
    unittest.main()


class TestMoveToTop(unittest.TestCase):
    def test_tagged_cells_move_below_title(self):
        import json, tempfile
        import run_notebook as rn
        cells = [{"id": n, "metadata": ({"tags": ["findings-top"]} if n in "de" else {})} for n in "abcde"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "n.ipynb"
            path.write_text(json.dumps({"cells": cells}))
            self.assertEqual(rn.move_tagged_to_top(path), 2)
            self.assertEqual([c["id"] for c in json.loads(path.read_text())["cells"]], list("adebc"))
            self.assertEqual(rn.move_tagged_to_top(path, tag="none"), 0)


class TestAcrossAndVerdict(unittest.TestCase):
    def test_across_names_the_ends(self):
        status, text = nn.across(nn.Rule(0.1, 0.3), {"a": 0.05, "b": 0.2, "c": 0.12}, nn.pct)
        self.assertEqual(status, nn.LOOK)
        self.assertEqual(text, "5% (`a`) to 20% (`b`)")

    def test_across_single_or_equal(self):
        self.assertEqual(nn.across(nn.Rule(1, 2), {"a": 0.5})[1], "0.50")
        self.assertEqual(nn.across(nn.Rule(1, 2), {"a": 0.5, "b": 0.501})[1], "0.50")

    def test_verdict_counts(self):
        f = nn.Findings()
        f.add(nn.OK, "a"); f.add(nn.LOOK, "b"); f.add(nn.BAD, "c")
        v = f.verdict()
        self.assertIn("3 checks: 1 fine, 1 worth a look, 1 problems", v)
        self.assertIn("problems to fix", v)
        self.assertIn("No problems found", nn.Findings().verdict())
