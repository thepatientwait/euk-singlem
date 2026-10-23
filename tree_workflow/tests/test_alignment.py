"""Tests for merge_gene_fastas.py and generate_alignment.py (no network)."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import generate_alignment as ga  # noqa: E402
import merge_gene_fastas as mg  # noqa: E402


def fasta(path, records):
    Path(path).write_text("".join(f">{k}\n{v}\n" for k, v in records.items()))
    return Path(path)


class TestMerge(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def test_cleans_stops_and_merges(self):
        base = fasta(self.dir / "b.faa", {"A|g|1": "MKV*", "B|g|1": "MK*V*"})
        ext = fasta(self.dir / "e.faa", {"E|g|1": "MMM"})
        out = self.dir / "m.faa"
        stats = mg.merge(base, ext, out, True)
        self.assertEqual(dict(mg.read_records(out)), {"A|g|1": "MKV", "B|g|1": "MKXV", "E|g|1": "MMM"})
        self.assertEqual((stats["base"], stats["external"]), (2, 1))

    def test_duplicate_ids_raise(self):
        base = fasta(self.dir / "b.faa", {"A|g|1": "MK"})
        ext = fasta(self.dir / "e.faa", {"A|g|1": "MK"})
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            mg.merge(base, ext, self.dir / "m.faa", True)

    def test_external_ignored_when_disabled_and_keep_filter(self):
        base = fasta(self.dir / "b.faa", {"A|g|1": "MK", "B|g|1": "MV"})
        ext = fasta(self.dir / "e.faa", {"E|g|1": "MM"})
        mg.merge(base, ext, self.dir / "m.faa", False)
        self.assertEqual(set(dict(mg.read_records(self.dir / "m.faa"))), {"A|g|1", "B|g|1"})
        stats = mg.merge(base, ext, self.dir / "k.faa", True, keep={"A"})
        self.assertEqual(set(dict(mg.read_records(self.dir / "k.faa"))), {"A|g|1", "E|g|1"})
        self.assertEqual(stats["dropped_not_kept"], 1)


class TestCopyAndSubset(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)
        self.base = fasta(self.dir / "g.faa", {"A|g|1": "MK", "B|g|1": "MV", "C|g|1": "MM"})
        (self.dir / "g_status.tsv").write_text("sequence_id\tstatus\tscore\nA|g|1\tComplete\t9\nB|g|1\tDuplicated\t8\nC|g|1\tComplete\t7\n")
        (self.dir / "e").mkdir()
        self.ext = fasta(self.dir / "e" / "g.faa", {"E|g|1": "MQ", "F|g|1": "MW"})
        (self.dir / "e" / "g_status.tsv").write_text("sequence_id\tstatus\tscore\nE|g|1\tDuplicated\t5\nF|g|1\tComplete\t5\n")

    def ids(self, **kw):
        out = self.dir / "m.faa"
        stats = mg.merge(self.base, self.ext, out, True, **kw)
        return set(dict(mg.read_records(out))), stats

    def test_multi_keeps_everything(self):
        ids, _ = self.ids(copy="multi")
        self.assertEqual(ids, {"A|g|1", "B|g|1", "C|g|1", "E|g|1", "F|g|1"})

    def test_single_drops_duplicated_in_both_inputs(self):
        ids, stats = self.ids(copy="single")
        self.assertEqual(ids, {"A|g|1", "C|g|1", "F|g|1"})
        self.assertEqual(stats["dropped_duplicated"], 2)

    def test_subset_applies_to_base_and_external_keep_to_base_only(self):
        ids, _ = self.ids(copy="multi", subset={"A", "E"})
        self.assertEqual(ids, {"A|g|1", "E|g|1"})
        ids, _ = self.ids(copy="multi", keep={"A"})
        self.assertEqual(ids, {"A|g|1", "E|g|1", "F|g|1"})
        ids, _ = self.ids(copy="single", keep={"A", "B"}, subset={"A", "B", "F"})
        self.assertEqual(ids, {"A|g|1", "F|g|1"})

    def test_missing_sidecar_counts_everything_as_complete(self):
        (self.dir / "g_status.tsv").unlink()
        ids, stats = self.ids(copy="single")
        self.assertIn("B|g|1", ids)
        self.assertEqual(stats["no_sidecar"], 1)

    def test_bad_copy(self):
        with self.assertRaises(ValueError):
            self.ids(copy="both")


class TestAlignment(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = Path(self.tmp.name)
        self.g1 = fasta(d / "g1_bmge.fasta", {"A|g1|1": "AAA", "B|g1|1": "CCC", "C|g1|1": "DDD"})
        self.g2 = fasta(d / "g2_bmge.fasta", {"A|g2|1": "EE", "B|g2|1": "FF"})
        # the best copy (lowest number) wins even though it comes later in the file
        self.g3 = fasta(d / "g3_bmge.fasta", {"A|g3|2": "XX", "A|g3|1": "YY", "B|g3|1": "ZZ"})

    def test_gap_fill_and_best_copy(self):
        seqs, stats = ga.concatenate([self.g1, self.g2, self.g3])
        self.assertEqual(seqs["A"], "AAA" + "EE" + "YY")
        self.assertEqual(seqs["C"], "DDD" + "--" + "--")
        self.assertEqual(stats["genes_with_repeated_genomes"], 1)
        self.assertEqual(stats["columns"], 7)
        self.assertEqual([(q["gene"], q["start"], q["end"], q["genomes"]) for q in stats["partitions"]],
                         [("g1", 1, 3, 3), ("g2", 4, 5, 2), ("g3", 6, 7, 2)])

    def test_min_genes_and_exclude(self):
        seqs, stats = ga.concatenate([self.g1, self.g2, self.g3], min_genes_per_genome=2)
        self.assertEqual(set(seqs), {"A", "B"})
        self.assertEqual(stats["genomes_dropped_few_genes"], 1)
        seqs, _ = ga.concatenate([self.g1, self.g2, self.g3], exclude={"B"})
        self.assertEqual(set(seqs), {"A", "C"})

    def test_ragged_alignment_is_an_error(self):
        bad = fasta(Path(self.tmp.name) / "g9_bmge.fasta", {"A|g9|1": "AA", "B|g9|1": "A"})
        with self.assertRaisesRegex(ValueError, "different lengths"):
            ga.concatenate([bad])


if __name__ == "__main__":
    unittest.main()
