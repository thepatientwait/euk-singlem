import gzip
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import orthodb_tables as ot  # noqa: E402

ROWS = [
    "100at2759\t9606_0:001\tMKVL*\n",
    "100at2759\t4932_0:002\tMKIL\n",
    "200at2\t562_0:001\tMSSS\n",          # bacterial level: must be dropped
    "300at33154\t9606_0:009\tMAAA\n",      # opisthokonta level: must be dropped
    "1000at2759\t7227_0:003\t\n",          # empty sequence: skipped in FASTA
    "12at27590\t1_0:1\tMKV\n",             # 'at27590' is not 'at2759'
]


class TestTables(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.src = self.dir / "t.tab.gz"
        with gzip.open(self.src, "wt") as handle:
            handle.writelines(ROWS)

    def test_subset_keeps_only_the_level(self):
        out = self.dir / "o.tab.gz"
        kept, total = ot.subset(self.src, out, 0, 2759)
        self.assertEqual((kept, total), (3, 6))
        with gzip.open(out, "rt") as handle:
            self.assertEqual([l.split("\t")[0] for l in handle], ["100at2759", "100at2759", "1000at2759"])

    def test_fasta_headers_and_cleaning(self):
        out = self.dir / "o.faa"
        records, _ = ot.to_fasta(self.src, out, 0, 1, 2, 2759)
        self.assertEqual(records, 2)
        self.assertEqual(out.read_text(), ">9606_0:001|100at2759\nMKVL\n>4932_0:002|100at2759\nMKIL\n")

    def test_wrong_column_is_an_error_not_an_empty_file(self):
        with self.assertRaisesRegex(ValueError, "column 1"):
            ot.subset(self.src, self.dir / "x.gz", 1, 2759)

    def test_other_level(self):
        kept, _ = ot.subset(self.src, self.dir / "x.gz", 0, 2)
        self.assertEqual(kept, 1)


if __name__ == "__main__":
    unittest.main()
