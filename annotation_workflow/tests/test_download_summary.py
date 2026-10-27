import io
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import download_summary as ds  # noqa: E402


class TestSummary(unittest.TestCase):
    def test_statuses(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = Path(tmp) / "raw"
            raw.mkdir()
            (raw / "a.gz").write_bytes(b"x" * 500)
            (raw / "b.gz").write_bytes(b"x" * 5)
            (raw / "c.gz.partial").write_bytes(b"x" * 50)
            config = {"min_download_bytes": 100, "orthodb": {"files": {k: {"name": f"{k}.gz"} for k in "abcd"}},
                      "paths": {"orthodb_raw": str(raw), "summary": str(Path(tmp) / "r" / "s.tsv")}}
            out = io.StringIO()
            rows = ds.summarise(config, out)
            self.assertEqual([r[2] for r in rows], ["present", "too_small", "partial", "missing"])
            self.assertIn("1 of 4 files present", out.getvalue())
            self.assertTrue((Path(tmp) / "r" / "s.tsv").exists())


if __name__ == "__main__":
    unittest.main()
