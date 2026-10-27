import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import fetch_file as ff  # noqa: E402


class TestFetch(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.src = self.dir / "src.txt"
        self.src.write_bytes(b"orthodb-test\n" * 100)
        self.url = f"file://{self.src}"
        self.out = self.dir / "out" / "dst.txt"

    def test_downloads_and_leaves_no_partial(self):
        self.assertEqual(ff.download(self.url, self.out, sleep=lambda s: None), 1)
        self.assertEqual(self.out.read_bytes(), self.src.read_bytes())
        self.assertFalse(self.out.with_name("dst.txt.partial").exists())

    def test_md5_checked(self):
        good = hashlib.md5(self.src.read_bytes()).hexdigest()
        ff.download(self.url, self.out, md5=good, sleep=lambda s: None)
        self.out.unlink()
        with self.assertRaisesRegex(RuntimeError, "md5 mismatch"):
            ff.download(self.url, self.out, md5="0" * 32, retries=2, sleep=lambda s: None)
        self.assertFalse(self.out.exists())
        self.assertFalse(self.out.with_name("dst.txt.partial").exists())

    def test_empty_file_refused(self):
        empty = self.dir / "empty.txt"
        empty.write_bytes(b"")
        with self.assertRaisesRegex(RuntimeError, "only 0 bytes"):
            ff.download(f"file://{empty}", self.out, retries=2, sleep=lambda s: None)
        self.assertFalse(self.out.exists())

    def test_missing_source_fails_after_retries(self):
        waits = []
        with self.assertRaisesRegex(RuntimeError, "after 3 attempts"):
            ff.download(f"file://{self.dir}/nope", self.out, retries=3, sleep=waits.append)
        self.assertEqual(len(waits), 2)       # backs off between attempts, not after the last
        self.assertFalse(self.out.exists())

    def test_resumes_a_partial_file(self):
        partial = self.out.with_name("dst.txt.partial")
        partial.parent.mkdir(parents=True)
        partial.write_bytes(self.src.read_bytes()[:50])
        ff.download(self.url, self.out, md5=hashlib.md5(self.src.read_bytes()).hexdigest(), sleep=lambda s: None)
        self.assertEqual(self.out.read_bytes(), self.src.read_bytes())


if __name__ == "__main__":
    unittest.main()
