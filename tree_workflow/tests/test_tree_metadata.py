"""generate_tree_metadata: rank columns from the lineage string, and novelty as written by assign_taxonomy."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import generate_tree_metadata as gm  # noqa: E402
import polars as pl  # noqa: E402


class TestMetadata(unittest.TestCase):
    def test_rank_columns(self):
        frame = pl.DataFrame({"genome_id": ["a", "b"], "taxonomy": [
            "d__Eukarya;p__Asco;c__Sord;o__Hyp;f__Nec;g__Fus;s__Fus sp",
            "d__Eukarya;p__Asco;c__;o__;f__;g__;s__"]})
        out = gm.expand_taxonomy(frame)
        self.assertEqual(out.row(0, named=True)["family"], "f__Nec")
        self.assertEqual(out.row(0, named=True)["species"], "s__Fus sp")
        self.assertEqual(out.row(1, named=True)["class"], "c__")
        self.assertEqual(out.row(1, named=True)["phylum"], "p__Asco")

    def test_missing_rank_item_is_just_the_prefix(self):
        out = gm.expand_taxonomy(pl.DataFrame({"genome_id": ["a"], "taxonomy": ["d__Eukarya;p__Asco"]}))
        self.assertEqual(out.row(0, named=True)["genus"], "g__")

    def test_novelty_column_is_used_and_old_files_still_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            new, old = Path(tmp) / "new.tsv", Path(tmp) / "old.tsv"
            new.write_text("genome_id\ttaxonomy\tnovelty\ng1\td__Eukarya;p__A\tgenus\ng2\td__Eukarya;p__A\tknown\n")
            old.write_text("genome_id\ttaxonomy\ng3\td__Eukarya;p__A\n")
            got = gm.read_taxonomy(new, source="new", novelty="novel")
            self.assertEqual(dict(zip(got["genome_id"], got["novelty"])), {"g1": "genus", "g2": "known"})
            self.assertEqual(gm.read_taxonomy(old, source="new", novelty="novel")["novelty"].to_list(), ["novel"])
            self.assertEqual(gm.read_taxonomy(old, source="reference")["novelty"].to_list(), ["known"])


if __name__ == "__main__":
    unittest.main()
