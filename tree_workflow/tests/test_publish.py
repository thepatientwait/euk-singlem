"""Tests for scripts/publish_trees.py."""
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import publish_trees as pt  # noqa: E402


class TestPublish(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        d = Path(self.tmp.name)
        (d / "a.nwk").write_text("(a,b);")
        (d / "m.tsv").write_text("x\ty\n")
        self.files = [(d / "a.nwk", "all/single/tree.rerooted.nwk"), (d / "m.tsv", "all/single/tree_metadata.tsv")]
        self.dest = d / "trees"

    def test_copies_files_and_writes_manifest(self):
        snap = pt.publish(self.files, self.dest, info={"modes": ["all"]}, repo=self.tmp.name)
        self.assertEqual((snap / "all/single/tree.rerooted.nwk").read_text(), "(a,b);")
        manifest = (snap / "MANIFEST.tsv").read_text().splitlines()
        self.assertEqual(manifest[0], "file\tbytes\tmd5")
        self.assertEqual(len(manifest), 3)
        info = (snap / "RUN_INFO.txt").read_text()
        self.assertIn("files       : 2", info)
        self.assertIn('modes       : ["all"]', info)

    def test_each_call_gets_its_own_directory_and_latest_points_to_newest(self):
        now = time.localtime(1_700_000_000)
        a = pt.publish(self.files, self.dest, now=now, repo=self.tmp.name)
        b = pt.publish(self.files, self.dest, now=now, repo=self.tmp.name)   # same second
        self.assertNotEqual(a, b)
        self.assertTrue(b.name.endswith("_2"))
        self.assertEqual((self.dest / "latest").resolve(), b.resolve())

    def test_subset_in_name(self):
        self.assertTrue(pt.snapshot_name("family5", time.localtime(0)).endswith("_subset_family5"))
        self.assertNotIn("subset", pt.snapshot_name("none", time.localtime(0)))

    def test_missing_source_publishes_nothing(self):
        with self.assertRaisesRegex(FileNotFoundError, "nope.nwk"):
            pt.publish(self.files + [(Path(self.tmp.name) / "nope.nwk", "x")], self.dest)
        self.assertFalse(self.dest.exists())


if __name__ == "__main__":
    unittest.main()
