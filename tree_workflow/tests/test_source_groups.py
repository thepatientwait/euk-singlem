"""Source group lookup and the by-source breakdown used by the tree metadata and notebooks."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import source_groups as sg  # noqa: E402

HEADER = "group\tgenome_acc\tbusco_result\tcomplete_count\n"


class TestReadGroups(unittest.TestCase):
    def write(self, tmp, name, rows):
        path = Path(tmp) / name
        path.write_text(HEADER + "".join(f"{g}\t{a}\tx/full_table.tsv\t50\n" for g, a in rows))
        return path

    def test_reads_several_manifests_first_wins_and_skips_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = self.write(tmp, "a.tsv", [("fungi", "F1"), ("tsa", "GLPS00000000.1"), ("ncbi", "dup")])
            b = self.write(tmp, "b.tsv", [("external", "E1"), ("external", "dup")])
            groups = sg.read_source_groups([a, Path(tmp) / "missing.tsv", b])
        self.assertEqual(groups, {"F1": "fungi", "GLPS00000000.1": "tsa", "dup": "ncbi", "E1": "external"})
        self.assertEqual(sg.group_of("nope", groups), "unknown")


class TestBySource(unittest.TestCase):
    def test_counts_and_shares(self):
        groups = {"a": "fungi", "b": "fungi", "c": "tsa", "d": "tsa", "e": "tsa"}
        rows = sg.by_source(["a", "b", "c", "d", "e", "x"], groups, flagged={"b", "c", "d", "x"})
        self.assertEqual([r["group"] for r in rows], ["fungi", "tsa", "unknown"])
        by = {r["group"]: r for r in rows}
        self.assertEqual((by["fungi"]["tips"], by["fungi"]["flagged"], by["fungi"]["share"]), (2, 1, 0.5))
        self.assertAlmostEqual(by["tsa"]["share"], 2 / 3)
        self.assertEqual(by["unknown"]["flagged"], 1)

    def test_no_flagged(self):
        self.assertEqual(sg.by_source(["a"], {"a": "p10k"})[0]["share"], 0.0)


if __name__ == "__main__":
    unittest.main()
