"""Tests for the download robustness changes (no network)."""
import builtins
import gzip
import io
import json
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
builtins.snakemake = types.SimpleNamespace()  # the scripts reference it only inside main()

import download_genome as dg  # noqa: E402
import helpers  # noqa: E402
import jgi_common  # noqa: E402
import refresh_jgi_metadata as rjm  # noqa: E402


def org(acc, status, size=10):
    return {"portal_detail_id": acc, "files": [{"_id": acc, "file_name": f"{acc}.fasta.gz",
                                                "file_size": size, "file_status": status}]}


class FakeLogger:
    def __getattr__(self, name):
        return lambda *a, **k: None


class TestBackoff(unittest.TestCase):
    def test_jitter_and_growth(self):
        for attempt, top in [(1, 15), (2, 30), (3, 60), (10, 450)]:
            for _ in range(50):
                delay = jgi_common.backoff_delay(attempt)
                self.assertGreaterEqual(delay, top / 3 - 1e-9)
                self.assertLessEqual(delay, top + 1e-9)

    def test_retry_after_wins(self):
        response = types.SimpleNamespace(headers={"Retry-After": "120"})
        self.assertTrue(120 <= jgi_common.backoff_delay(1, response) <= 125)
        huge = types.SimpleNamespace(headers={"Retry-After": "99999"})
        self.assertLessEqual(jgi_common.backoff_delay(1, huge), 605)


class TestFetchAllPages(unittest.TestCase):
    @staticmethod
    def fake_pages(total, size=50, drop_page=None):
        def request(session, method, url, logger, **kwargs):
            page = int(url.rsplit("&p=", 1)[1])
            start = (page - 1) * size
            orgs = [{"portal_detail_id": f"o{i}"} for i in range(start, min(start + size, total))]
            if page == drop_page:
                orgs = orgs[:-1]
            body = {"organisms": orgs, "total": total, "next_page": page + 1 if start + size < total else False}
            return types.SimpleNamespace(raise_for_status=lambda: None, json=lambda: body)
        return request

    def test_merges_pages_in_order(self):
        with mock.patch.object(jgi_common, "request_with_retries", self.fake_pages(257)):
            data = jgi_common.fetch_all_pages("http://x?a=1", FakeLogger())
        self.assertEqual([o["portal_detail_id"] for o in data["organisms"]], [f"o{i}" for i in range(257)])
        self.assertFalse(data["next_page"])

    def test_single_page(self):
        with mock.patch.object(jgi_common, "request_with_retries", self.fake_pages(7)):
            self.assertEqual(len(jgi_common.fetch_all_pages("http://x?a=1", FakeLogger())["organisms"]), 7)

    def test_short_page_is_an_error(self):
        with mock.patch.object(jgi_common, "request_with_retries", self.fake_pages(257, drop_page=3)):
            with self.assertRaisesRegex(RuntimeError, "Expected 257"):
                jgi_common.fetch_all_pages("http://x?a=1", FakeLogger())


class TestRefresh(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "fungi.json"
        self.old = {"organisms": [org("a", "PURGED"), org("b", "RESTORE_REGISTERED")]}
        self.path.write_text(json.dumps(self.old))
        os.utime(self.path, (1000, 1000))

    def tearDown(self):
        self.tmp.cleanup()

    def refresh(self, new=None, error=None):
        with mock.patch.object(rjm, "fetch_all_pages", side_effect=error, return_value=new):
            return rjm.refresh_group("fungi", "http://x", self.path, FakeLogger())

    def test_unchanged_keeps_mtime(self):
        self.assertFalse(self.refresh({"organisms": list(self.old["organisms"]), "next_page": False}))
        self.assertEqual(self.path.stat().st_mtime, 1000)

    def test_changed_status_replaces(self):
        new = {"organisms": [org("a", "RESTORE_REGISTERED"), org("b", "RESTORED")]}
        self.assertTrue(self.refresh(new))
        saved = json.loads(self.path.read_text())["organisms"]
        self.assertEqual(saved[1]["files"][0]["file_status"], "RESTORED")
        self.assertGreater(self.path.stat().st_mtime, 1000)
        self.assertFalse(list(Path(self.tmp.name).glob("*.tmp")))

    def test_failure_keeps_existing(self):
        self.assertFalse(self.refresh(error=RuntimeError("HTTP 500")))
        self.assertEqual(json.loads(self.path.read_text()), self.old)

    def test_truncated_answer_refused(self):
        self.assertFalse(self.refresh({"organisms": [org("a", "RESTORED")]}))
        self.assertEqual(json.loads(self.path.read_text()), self.old)

    def test_missing_file_is_left_to_the_checkpoint(self):
        self.path.unlink()
        self.assertFalse(self.refresh({"organisms": []}))
        self.assertFalse(self.path.exists())


class TestWhitelist(unittest.TestCase):
    def test_only_ready_files_get_jobs(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = os.getcwd()
            os.chdir(tmp)
            try:
                Path("outputs/metadata/source/jgi").mkdir(parents=True)
                Path("outputs/genomes/fungi").mkdir(parents=True)
                orgs = [org("tape", "PURGED"), org("pending", "RESTORE_REGISTERED"), org("ready", "RESTORED"),
                        org("bk", "BACKUP_COMPLETE"), org("have", "PURGED")]
                Path("outputs/metadata/source/jgi/fungi.json").write_text(json.dumps({"organisms": orgs}))
                Path("outputs/genomes/fungi/have.fasta.gz").write_bytes(b"x")

                not_ready = helpers.get_jgi_not_ready_paths(["fungi"])
                self.assertEqual({Path(p).stem.split(".")[0]: s for p, s in not_ready.items()},
                                 {"tape": "PURGED", "pending": "RESTORE_REGISTERED", "have": "PURGED"})

                with mock.patch.object(helpers, "get_expected_genome_paths", return_value=[
                        f"outputs/genomes/fungi/{a}.fasta.gz" for a in ["tape", "pending", "ready", "bk", "have"]]):
                    todo = helpers.build_all_genomes_to_download(["fungi"], None, None, None)(None)
                    self.assertEqual([Path(p).name for p in todo], ["ready.fasta.gz", "bk.fasta.gz"])
            finally:
                os.chdir(old_cwd)


class TestResources(unittest.TestCase):
    def test_threads_double_and_cap(self):
        fn = helpers.get_threads(1, max_threads=4)
        self.assertEqual([fn(None, attempt=a) for a in (1, 2, 3, 4)], [1, 2, 4, 4])

    def test_fixed_memory(self):
        profile = helpers.create_resource_profile_factory(4000, "4h", 3)
        self.assertEqual(profile(mem_scale=0.5, fixed_mem=True)["mem_mb"], 2000)
        doubling = profile()["mem_mb"]
        self.assertEqual([doubling(None, attempt=a) for a in (1, 2)], [4000, 8000])


class TestNcbi(unittest.TestCase):
    def test_ftp_parent_url(self):
        self.assertEqual(dg.ftp_parent_url("GCF_019279795.1"),
                         "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/279/795/")
        self.assertIsNone(dg.ftp_parent_url("not-an-accession"))

    def test_find_dir(self):
        html = '<a href="GCA_000000001.1_asm/">x</a> <a href="GCA_000000001.2_other/">y</a>'
        session = mock.Mock()
        session.get.return_value = types.SimpleNamespace(status_code=200, text=html, raise_for_status=lambda: None)
        url = dg.find_ftp_assembly_dir("GCA_000000001.2", session)
        self.assertTrue(url.endswith("/GCA/000/000/001/GCA_000000001.2_other/"))
        session.get.return_value = types.SimpleNamespace(status_code=404, text="", raise_for_status=lambda: None)
        self.assertIsNone(dg.find_ftp_assembly_dir("GCA_000000001.1", session))

    def test_gzip_roundtrip_with_and_without_pigz(self):
        data = b">s\nACGT\n" * 1000
        for pigz in (None, "real"):
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / "g.gz"
                which = (lambda _: None) if pigz is None else dg.shutil.which
                if pigz == "real" and not dg.shutil.which("pigz"):
                    continue
                with mock.patch.object(dg.shutil, "which", which):
                    dg.gzip_stream(io.BytesIO(data), out, 2)
                self.assertEqual(gzip.open(out).read(), data)


if __name__ == "__main__":
    unittest.main()
