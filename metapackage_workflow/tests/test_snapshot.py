"""Snapshot create / compare on a tiny synthetic survey (no workflow outputs needed)."""
import contextlib
import shutil
import io
import json
import sys
import tempfile
import types
import unittest
import unittest.mock
from pathlib import Path

import polars as pl

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
import snapshot_intron_survey as snap  # noqa: E402


def choice(gene, corrected_start, corrected_span, n_mapped=100):
    return {"busco_id": gene, "hmm_columns": 50, "n_genomes": n_mapped, "n_mapped": n_mapped, "n_with_intron": 10,
            "naive_start": 5, "naive_span_fraction": 0.2, "naive_n_cover": n_mapped,
            "corrected_start": corrected_start, "corrected_span_fraction": corrected_span, "corrected_score_ratio": 0.99,
            "reason": "moved_for_introns", "max_span_fraction": 0.05, "violation": False}


class TestSnapshot(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        for name in ("config", "survey", "input", "hmm", "qc"):
            (self.root / name).mkdir()
        (self.root / "config" / "config.yaml").write_text("singlem:\n  introns:\n    max_span_fraction: 0.05\n")
        (self.root / "config" / "paths.yaml").write_text("x: 1\n")
        (self.root / "qc" / "manifest.tsv").write_text("group\tgenome_acc\tbusco_result\tcomplete_count\nfungi\tG1\tp\t5\nfungi\tG2\tp\t5\n")
        (self.root / "qc" / "taxonomy.tsv").write_text("acc\tphylum\nG1\tAscomycota\nG2\tAscomycota\n")
        for attr, value in {
            "ROOT": self.root, "SNAPSHOT_ROOT": self.root / "snapshots", "SOURCE_DIR": self.root / "survey",
"GENE_INPUT_DIR": self.root / "input", "HMM_DIR": self.root / "hmm",
            "BUSCO_MANIFEST": self.root / "qc" / "manifest.tsv", "TAXONOMY": self.root / "qc" / "taxonomy.tsv",
            "SUMMARY": self.root / "qc" / "summary.tsv", "GENE_NAMES": self.root / "qc" / "names.tsv", "NOTEBOOK_DIR": self.root / "nonexistent"}.items():
            patcher = unittest.mock.patch.object(snap, attr, value)
            patcher.start()
            self.addCleanup(patcher.stop)

    def survey(self, final_start, final_span, genes=("g1", "g2")):
        for path in (self.root / "survey").glob("*"):
            path.unlink()
        for gene in genes:
            pl.DataFrame([choice(gene, final_start, final_span)]).write_csv(self.root / "survey" / f"{gene}.window_choice.tsv", separator="\t")
            for suffix in ("introns", "genomes", "window_scan"):
                (self.root / "survey" / f"{gene}.{suffix}.tsv").write_text("a\tb\n1\t2\n")
            (self.root / "survey" / f"{gene}.window_position.naive.txt").write_text("5\n")
            (self.root / "survey" / f"{gene}.window_position.corrected.txt").write_text(f"{final_start}\n")
            (self.root / "input" / f"{gene}.faa").write_text(f">{gene}\nMK\n")
            (self.root / "hmm" / f"{gene}.hmm").write_text("HMMER3\nLENG 50\n")

    def create(self, name, **survey_args):
        self.survey(**survey_args)
        with contextlib.redirect_stdout(io.StringIO()):
            return snap.create(types.SimpleNamespace(name=name, note="test"))

    def test_create_is_complete_and_self_describing(self):
        out = self.create("one", final_start=10, final_span=0.01)
        info = json.loads((out / "RUN_INFO.json").read_text())
        self.assertEqual((info["n_genes"], info["n_genomes_in_busco_manifest"]), (2, 2))
        self.assertEqual(info["intron_settings"], {"max_span_fraction": 0.05})
        self.assertTrue((out / "data" / "g1.window_choice.tsv.gz").exists())
        self.assertTrue((out / "data" / "g1.window_position.naive.txt").exists())
        self.assertTrue((out / "data" / "g1.window_position.corrected.txt").exists())
        manifest = pl.read_csv(out / "MANIFEST.tsv", separator="\t")
        self.assertIn("RUN_INFO.json", manifest["path"].to_list())
        inputs = pl.read_csv(out / "inputs.tsv", separator="\t")
        self.assertEqual(sorted(inputs["kind"].unique().to_list()), ["busco_manifest", "gene_fasta", "hmm", "taxonomy_mapping"])
        self.assertEqual((self.root / "snapshots" / "latest").resolve(), out.resolve())

    def test_existing_snapshot_is_never_overwritten(self):
        self.create("one", final_start=10, final_span=0.01)
        with self.assertRaises(SystemExit), contextlib.redirect_stdout(io.StringIO()):
            snap.create(types.SimpleNamespace(name="one", note=""))

    def test_gzip_is_deterministic(self):
        source = self.root / "x.tsv"
        source.write_text("a\tb\n" * 1000)
        snap.gzip_copy(source, self.root / "a.gz")
        snap.gzip_copy(source, self.root / "b.gz")
        self.assertEqual((self.root / "a.gz").read_bytes(), (self.root / "b.gz").read_bytes())

    def test_compare_reports_moved_windows_and_changed_inputs(self):
        self.create("one", final_start=10, final_span=0.01)
        self.survey(final_start=12, final_span=0.03)
        (self.root / "hmm" / "g1.hmm").write_text("HMMER3\nLENG 51\n")          # an input changed
        with contextlib.redirect_stdout(io.StringIO()):
            snap.create(types.SimpleNamespace(name="two", note=""))
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            snap.compare(types.SimpleNamespace(old="one", new="two", out=None))
        text = buffer.getvalue()
        self.assertIn("corrected window moved", text)
        self.assertIn("2 in both", text)
        self.assertIn("'hmm'", text)
        diff = pl.read_csv(self.root / "snapshots" / "two" / "comparison_with_one.tsv", separator="\t")
        self.assertEqual(diff["corrected_start_new"].to_list(), [12, 12])
        self.assertAlmostEqual(diff["corrected_span_change"][0], 0.02)

    def test_old_snapshots_with_legacy_column_names_still_compare(self):
        self.create("one", final_start=10, final_span=0.01)
        old = self.root / "snapshots" / "old"
        shutil.copytree(self.root / "snapshots" / "one", old)
        for path in (old / "data").glob("*.window_choice.tsv.gz"):
            frame = pl.read_csv(path, separator="\t").rename({"naive_start": "best_coverage_start", "naive_span_fraction": "best_coverage_span_fraction",
                                                              "corrected_start": "final_start", "corrected_span_fraction": "final_span_fraction",
                                                              "corrected_score_ratio": "final_score_ratio"})
            frame.write_csv(path.with_suffix(""), separator="\t")
            path.unlink()
            snap.gzip_copy(path.with_suffix(""), path)
            path.with_suffix("").unlink()
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            snap.compare(types.SimpleNamespace(old="old", new="one", out=str(self.root / "cmp.tsv")))
        self.assertIn("corrected window moved", buffer.getvalue())


if __name__ == "__main__":
    import unittest.mock  # noqa: F401
    unittest.main()
