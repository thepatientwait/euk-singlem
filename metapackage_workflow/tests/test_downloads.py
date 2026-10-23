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


def fake_jwt(exp_offset=7200):
    import base64, time
    enc = lambda d: base64.urlsafe_b64encode(json.dumps(d).encode()).rstrip(b"=").decode()
    return f"{enc({'alg': 'none'})}.{enc({'exp': time.time() + exp_offset})}.sig"


class TestToken(unittest.TestCase):
    def test_double_bearer_is_normalised(self):
        jwt = fake_jwt()
        for raw in (jwt, f"Bearer {jwt}", f"Bearer Bearer {jwt}", f"  bearer   {jwt}\n"):
            self.assertEqual(jgi_common.normalize_token(raw), jwt)
            self.assertGreater(jgi_common.check_token(raw), 3600)
            self.assertEqual(jgi_common.make_session(raw).headers["Authorization"], f"Bearer {jwt}")

    def test_garbage_and_expired_fail(self):
        with self.assertRaises(jgi_common.JgiAuthError):
            jgi_common.check_token("Bearer not-a-jwt")
        with self.assertRaises(jgi_common.JgiAuthError):
            jgi_common.check_token(fake_jwt(-10))
        with self.assertRaises(jgi_common.JgiAuthError):
            jgi_common.check_token("")


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


class TestWriteIfChanged(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "fungi.json"
        self.old = {"organisms": [org("a", "PURGED"), org("b", "RESTORE_REGISTERED")]}
        self.path.write_text(json.dumps(self.old))
        os.utime(self.path, (1000, 1000))

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, new):
        return rjm.write_json_if_changed(self.path, new, FakeLogger(), label="fungi")

    def test_unchanged_keeps_mtime(self):
        self.assertFalse(self.write({"organisms": list(self.old["organisms"]), "next_page": False}))
        self.assertEqual(self.path.stat().st_mtime, 1000)

    def test_changed_status_replaces(self):
        new = {"organisms": [org("a", "RESTORE_REGISTERED"), org("b", "RESTORED")]}
        self.assertTrue(self.write(new))
        self.assertEqual(json.loads(self.path.read_text())["organisms"][1]["files"][0]["file_status"], "RESTORED")
        self.assertGreater(self.path.stat().st_mtime, 1000)
        self.assertFalse(list(Path(self.tmp.name).glob("*.tmp")))

    def test_truncated_answer_refused(self):
        with self.assertRaisesRegex(RuntimeError, "refusing"):
            self.write({"organisms": [org("a", "RESTORED")]})
        self.assertEqual(json.loads(self.path.read_text()), self.old)

    def test_missing_file_is_created(self):
        self.path.unlink()
        self.assertTrue(self.write({"organisms": [org("a", "RESTORED")]}))
        self.assertTrue(self.path.exists())

    def test_bytes(self):
        csv = Path(self.tmp.name) / "m.csv"
        self.assertTrue(rjm.write_bytes_if_changed(csv, b"x"))
        os.utime(csv, (1000, 1000))
        self.assertFalse(rjm.write_bytes_if_changed(csv, b"x"))
        self.assertEqual(csv.stat().st_mtime, 1000)
        self.assertTrue(rjm.write_bytes_if_changed(csv, b"y"))


class TestWaitForRestore(unittest.TestCase):
    """download_jgi_genome.main polls until the file is restored, then downloads."""

    def run_main(self, statuses, wait=True, report=None, request_restore=False):
        import download_jgi_genome as dj
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        out = Path(tmp.name) / "g.fasta.gz"
        sleeps, posts = [], []
        queue = list(statuses)

        def fake_query(acc, name, anonymous, authed):
            status = queue.pop(0) if len(queue) > 1 else queue[0]
            return {"_id": "o"}, {"_id": "f", "file_status": status, "file_name": name}

        class Download:
            status_code = 200
            def __enter__(self_inner):
                return self_inner
            def __exit__(self_inner, *a):
                return False
            def iter_content(self_inner, chunk_size):
                import zipfile
                buf = io.BytesIO()
                with zipfile.ZipFile(buf, "w") as zf:
                    zf.writestr("g.fa.gz", gzip.compress(b">s\nACGT\n"))
                yield buf.getvalue()

        def fake_request(session, method, url, logger, **kwargs):
            posts.append(url)
            return Download() if url == dj.JDP_DOWNLOAD_URL else types.SimpleNamespace(status_code=200, text="{}")

        secrets = Path(tmp.name) / "secrets.yaml"
        secrets.write_text("jgi_token: tok\n")
        builtins.snakemake = types.SimpleNamespace(
            wildcards=types.SimpleNamespace(genome_acc="g"),
            output=types.SimpleNamespace(genome=str(out)),
            input=types.SimpleNamespace(get=lambda key, default=None: report),
            params=types.SimpleNamespace(name="g.fa.gz", token="t", secrets=str(secrets), wait=wait,
                                         poll_seconds=900, request_restore=request_restore))
        with mock.patch.object(dj, "query_file", fake_query), \
                mock.patch.object(dj, "request_with_retries", fake_request), \
                mock.patch.object(dj, "check_token", lambda *a, **k: 1e6), \
                mock.patch.object(dj.time, "sleep", sleeps.append):
            dj.main()
        return out, sleeps, posts

    def test_polls_until_restored_then_downloads(self):
        out, sleeps, posts = self.run_main(["PURGED", "RESTORE_REGISTERED", "RESTORE_REGISTERED", "RESTORED"])
        self.assertEqual(len(sleeps), 3)
        self.assertTrue(all(810 <= x <= 990 for x in sleeps))  # 15 min +-10 %
        self.assertTrue(out.exists() and out.stat().st_size > 0)
        self.assertEqual(posts, [__import__("download_jgi_genome").JDP_DOWNLOAD_URL])

    def test_ready_file_does_not_wait(self):
        out, sleeps, _ = self.run_main(["RESTORED"])
        self.assertEqual(sleeps, [])
        self.assertTrue(out.exists())

    def test_without_wait_fails_at_once(self):
        import jgi_common
        with self.assertRaises(jgi_common.JgiNotRestoredError):
            self.run_main(["RESTORE_REGISTERED"], wait=False)

    def test_refused_restore_fails_at_once(self):
        import jgi_common
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        report = Path(tmp.name) / "nr.tsv"
        report.write_text("genome_acc\tfile_status\tstate\ng\tPURGED\trestore_refused\n")
        with self.assertRaisesRegex(jgi_common.JgiNotRestoredError, "refused"):
            self.run_main(["PURGED"], report=str(report))

    def test_requests_restore_itself_when_asked(self):
        out, sleeps, posts = self.run_main(["PURGED", "RESTORED"], request_restore=True)
        dj = __import__("download_jgi_genome")
        self.assertEqual(posts, [dj.JDP_RESTORE_URL, dj.JDP_DOWNLOAD_URL])


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
                    todo = helpers.build_all_genomes_to_download(["fungi"], None, None, None, skip_jgi_unrestored=True)(None)
                    self.assertEqual([Path(p).name for p in todo], ["ready.fasta.gz", "bk.fasta.gz"])
            finally:
                os.chdir(old_cwd)


class TestWaitingStatus(unittest.TestCase):
    def test_status_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = os.getcwd()
            os.chdir(tmp)
            try:
                Path("outputs/metadata/source/jgi").mkdir(parents=True)
                Path("outputs/metadata/source/jgi/fungi.json").write_text(json.dumps({"organisms": [
                    org("tape", "PURGED"), org("ready", "RESTORED")]}))
                self.assertTrue(helpers.jgi_is_waiting("fungi", "tape"))
                self.assertFalse(helpers.jgi_is_waiting("fungi", "ready"))
                self.assertEqual(helpers.get_jgi_file_status("fungi", "tape"), "PURGED")
            finally:
                os.chdir(old_cwd)

    def test_default_includes_waiting_genomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = os.getcwd()
            os.chdir(tmp)
            try:
                Path("outputs/metadata/source/jgi").mkdir(parents=True)
                Path("outputs/metadata/source/jgi/fungi.json").write_text(json.dumps({"organisms": [
                    org("tape", "PURGED"), org("ready", "RESTORED")]}))
                with mock.patch.object(helpers, "get_expected_genome_paths", return_value=[
                        "outputs/genomes/fungi/tape.fasta.gz", "outputs/genomes/fungi/ready.fasta.gz"]):
                    todo = helpers.build_all_genomes_to_download(["fungi"], None, None, None)(None)
                self.assertEqual(len(todo), 2)
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
