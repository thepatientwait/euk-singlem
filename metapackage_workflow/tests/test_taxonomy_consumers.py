"""Tests for how downstream scripts use the tax_status column (no network)."""
import builtins
import sys
import tempfile
import types
import unittest
from pathlib import Path

import polars as pl

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
builtins.snakemake = types.SimpleNamespace()  # the scripts reference it only inside main()
import aggregate_busco_sequences as ab  # noqa: E402

MAPPING = pl.DataFrame({
    "acc": ["g_complete", "g_partial", "g_unres"],
    "parsed_taxon_string": ["d__Eukarya;p__P;c__C;o__O;f__F;g__G;s__G x", "d__Eukarya;p__P;c__;o__;f__;g__;s__",
                            "d__Eukarya;p__;c__;o__;f__;g__;s__"],
    "tax_status": ["complete", "partial", "unresolved"],
})


class TestUsableTaxonomy(unittest.TestCase):
    def test_unresolved_and_missing_are_unusable(self):
        usable = ab.load_usable_taxonomy(MAPPING)
        self.assertEqual(set(usable), {"g_complete", "g_partial"})

    def test_old_schema_fails_clearly(self):
        with self.assertRaisesRegex(ValueError, "tax_status"):
            ab.load_usable_taxonomy(MAPPING.drop("tax_status"))

    def test_filter_resolved_sequences(self):
        usable = ab.load_usable_taxonomy(MAPPING)
        seqs = {"a": {"original_genome_acc": "g_complete"}, "b": {"original_genome_acc": "g_unres"},
                "c": {"original_genome_acc": "nope"}}
        self.assertEqual(set(ab.filter_resolved_sequences(seqs, usable, "BUSCO1")), {"a"})
        self.assertEqual(set(ab.filter_resolved_sequences(seqs, usable, "BUSCO1", allow_missing_taxonomy=True)),
                         {"a", "b", "c"})

    def test_taxonomy_file_only_lists_usable(self):
        usable = ab.load_usable_taxonomy(MAPPING)
        seqs = {"a": {"original_genome_acc": "g_complete"}, "b": {"original_genome_acc": "g_unres"}}
        with tempfile.TemporaryDirectory() as tmp:
            path = ab.write_taxonomy_file("BUSCO1", seqs, usable, Path(tmp))
            lines = path.read_text().splitlines()
        self.assertEqual(lines[0], "sequence_id\ttaxonomy")
        self.assertEqual([l.split("\t")[0] for l in lines[1:]], ["a"])


if __name__ == "__main__":
    unittest.main()
