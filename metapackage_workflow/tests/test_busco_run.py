"""BUSCO run directories: the original per-gene layout and the condensed one read identically."""
import builtins
import gzip
import sys
import tempfile
import types
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
builtins.snakemake = types.SimpleNamespace()
import aggregate_busco_sequences as agg  # noqa: E402
import busco_run as br  # noqa: E402
import extract_busco_transcripts as ex  # noqa: E402
import run_busco  # noqa: E402
from Bio import SeqIO  # noqa: E402

TABLE = ("# BUSCO version is: 6.1.0\n# The lineage dataset is: eukaryota_odb12\n"
         "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
         "1at2759\tComplete\tc1\t11\t100\t+\t100.0\t30\n"
         "2at2759\tDuplicated\tc1\t201\t290\t-\t90.0\t30\n"
         "2at2759\tDuplicated\tc2\t5\t94\t+\t80.0\t30\n")
GFF1 = ("c1\tminiprot\tmRNA\t11\t100\t100\t+\t.\tID=M1;Rank=1;Target=1at2759_7_0:aa 1 30\n"
        "c1\tminiprot\tCDS\t11\t55\t50\t+\t0\tParent=M1;Rank=1;Target=1at2759_7_0:aa 1 15\n"
        "c1\tminiprot\tCDS\t76\t100\t50\t+\t0\tParent=M1;Rank=1;Target=1at2759_7_0:aa 16 30\n")
GFF2 = ("c1\tminiprot\tmRNA\t201\t290\t100\t-\t.\tID=M2;Rank=1;Target=2at2759_8_0:bb 1 30\n"
        "c1\tminiprot\tCDS\t201\t290\t50\t-\t0\tParent=M2;Rank=1;Target=2at2759_8_0:bb 1 30\n"
        "c2\tminiprot\tmRNA\t5\t94\t100\t+\t.\tID=M3;Rank=1;Target=2at2759_8_0:cc 1 30\n"
        "c2\tminiprot\tCDS\t5\t94\t50\t+\t0\tParent=M3;Rank=1;Target=2at2759_8_0:cc 1 30\n")
FAA1 = ">1at2759_7_0:aa|c1:11-100|+\nMKKLLAAGGMKKLLAAGGMKKLLAAGGMKK*\n"
FAA2 = ">2at2759_8_0:bb|c1:201-290|-\nMAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n>2at2759_8_0:cc|c2:5-94|+\nMCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"


def make_legacy(root: Path) -> Path:
    run = root / "fungi" / "G1" / "run_eukaryota_odb12"
    single = run / "busco_sequences" / "single_copy_busco_sequences"
    multi = run / "busco_sequences" / "multi_copy_busco_sequences"
    frag = run / "busco_sequences" / "fragmented_busco_sequences"
    for d in (single, multi, frag, run / "miniprot_output", run / "hmmer_output"):
        d.mkdir(parents=True)
    (run / "full_table.tsv").write_text(TABLE)
    (run / "short_summary.txt").write_text("summary\n")
    (single / "1at2759.faa").write_text(FAA1)
    (single / "1at2759.gff").write_text(GFF1)
    (multi / "2at2759.faa").write_text(FAA2)
    (multi / "2at2759.gff").write_text(GFF2)
    (run / "miniprot_output" / "ref.mpi").write_bytes(b"x" * 100000)
    (run.parent / "tmp").mkdir()
    (run.parent / "tmp" / "refseq_db.faa").write_text(">x\nMK\n")
    (run.parent / "logs").mkdir()
    (run.parent / "logs" / "busco.log").write_text("log\n")
    (run.parent / "logs" / "miniprot_align_out.log").write_text("big\n" * 1000)
    (run / "hmmer_output" / "a.out").write_text("hmm\n")
    return run


class TestBuscoRun(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.run = make_legacy(self.root)

    def test_legacy_reads(self):
        r = br.BuscoRun(self.run)
        self.assertEqual(r.layout, "legacy")
        self.assertEqual([x.description for x in r.protein_records("2at2759")],
                         ["2at2759_8_0:bb|c1:201-290|-", "2at2759_8_0:cc|c2:5-94|+"])
        self.assertEqual(r.protein_records("1at2759")[0].id, "1at2759_7_0:aa|c1:11-100|+")
        models = r.gff_models("2at2759")
        self.assertEqual([(m.contig, m.strand, m.exons) for m in models], [("c1", "-", [(201, 290)]), ("c2", "+", [(5, 94)])])
        self.assertEqual(len(r.gff_models("1at2759")[0].exons), 2)
        self.assertEqual(r.protein_records("9at2759"), [])
        self.assertEqual(r.gff_models("9at2759"), [])

    def test_condensed_reads_the_same_and_files_are_gone(self):
        before = br.BuscoRun(self.run)
        want = {g: ([(x.description, x.seq) for x in before.protein_records(g)], before.gff_lines(g))
                for g in ("1at2759", "2at2759")}
        stats = br.condense_run(self.run)
        self.assertEqual(stats["genes"], 2)
        self.assertGreaterEqual(stats["removed_bytes"], 100000)
        names = sorted(p.name for p in self.run.iterdir())
        self.assertEqual(names, ["full_table.tsv", "models.gff.gz", "proteins.faa.gz", "short_summary.txt"])
        after = br.BuscoRun(self.run)
        self.assertEqual(after.layout, "condensed")
        got = {g: ([(x.description, x.seq) for x in after.protein_records(g)], after.gff_lines(g))
               for g in ("1at2759", "2at2759")}
        self.assertEqual(got, want)

    def test_gff_lines_without_a_target_attribute_survive_condensing(self):
        gff = self.run / "busco_sequences" / "single_copy_busco_sequences" / "1at2759.gff"
        gff.write_text("c1\tminiprot\tmRNA\t11\t100\t100\t+\t.\tID=M1\nc1\tminiprot\tCDS\t11\t100\t50\t+\t0\tParent=M1\n")
        br.condense_run(self.run)
        models = br.BuscoRun(self.run).gff_models("1at2759")
        self.assertEqual([m.exons for m in models], [[(11, 100)]])

    def test_condense_is_idempotent(self):
        br.condense_run(self.run)
        again = br.condense_run(self.run)
        self.assertEqual(again["removed_files"], 0)
        self.assertEqual(br.BuscoRun(self.run).layout, "condensed")
        self.assertEqual(len(br.BuscoRun(self.run).protein_records("2at2759")), 2)

    def test_genome_folder_is_tidied(self):
        br.condense_run(self.run)
        self.assertFalse((self.run.parent / "tmp").exists())
        self.assertEqual(sorted(p.name for p in (self.run.parent / "logs").iterdir()), ["busco.log"])

    def test_keep_originals_when_not_removing(self):
        br.condense_run(self.run, remove=False)
        self.assertTrue((self.run / "miniprot_output" / "ref.mpi").exists())
        self.assertEqual(br.BuscoRun(self.run).layout, "condensed")   # condensed files win once present

    def test_aggregation_gives_identical_sequences_for_both_layouts(self):
        table = self.run / "full_table.tsv"
        legacy, err = agg.process_single_busco_file(str(table), use_multi_copy=True)
        self.assertIsNone(err)
        br.condense_run(self.run)
        condensed, err = agg.process_single_busco_file(str(table), use_multi_copy=True)
        self.assertIsNone(err)
        self.assertEqual(legacy, condensed)
        self.assertEqual(set(legacy), {"1at2759", "2at2759"})

    def test_extraction_of_exon_blocks_matches_for_both_layouts(self):
        genome = self.root / "g.fasta"
        genome.write_text(">c1\n" + "".join("ACGT"[i % 4] for i in range(400)) + "\n>c2\n" + "ACGT" * 30 + "\n")
        out1, out2 = self.root / "o1.fa", self.root / "o2.fa"
        ex.extract_busco_transcripts(genome, self.run / "full_table.tsv", out1, min_intron=10, min_exon_nt=10)
        br.condense_run(self.run)
        ex.extract_busco_transcripts(genome, self.run / "full_table.tsv", out2, min_intron=10, min_exon_nt=10)
        a = [(r.id, str(r.seq)) for r in SeqIO.parse(out1, "fasta")]
        b = [(r.id, str(r.seq)) for r in SeqIO.parse(out2, "fasta")]
        self.assertEqual(a, b)
        self.assertGreaterEqual(len(a), 2)


class TestFinalize(unittest.TestCase):
    def test_result_is_condensed_and_moved(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            produced = tmp / "scratch" / "G1"
            make_legacy(tmp / "mk")
            (tmp / "scratch").mkdir()
            (tmp / "mk" / "fungi" / "G1").rename(produced)
            (produced / "short_summary.specific.eukaryota_odb12.G1.json").write_text("{}")
            final = tmp / "outputs" / "busco" / "fungi" / "G1"
            final.mkdir(parents=True)
            (final / "stale.txt").write_text("old")
            stats = run_busco.finalize(produced, final)
            self.assertFalse(produced.exists())
            self.assertFalse((final / "stale.txt").exists())
            self.assertTrue((final / "run_eukaryota_odb12" / "full_table.tsv").exists())
            self.assertTrue((final / "short_summary.specific.eukaryota_odb12.G1.json").exists())
            self.assertFalse((final / "run_eukaryota_odb12" / "miniprot_output").exists())
            self.assertGreaterEqual(stats["removed_bytes"], 100000)

    def test_can_keep_everything(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            make_legacy(tmp / "mk")
            final = tmp / "out" / "G1"
            run_busco.finalize(tmp / "mk" / "fungi" / "G1", final, condense=False)
            self.assertTrue((final / "run_eukaryota_odb12" / "miniprot_output" / "ref.mpi").exists())


class TestCondenseScript(unittest.TestCase):
    def test_report_then_apply(self):
        import condense_busco_runs as cb
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            make_legacy(tmp / "outputs" / "busco")
            roots = [tmp / "outputs" / "busco"]
            self.assertEqual(len(cb.find_runs(roots)), 1)
            self.assertEqual(cb.main([str(roots[0])]), 0)    # report only
            run = cb.find_runs(roots)[0]
            self.assertTrue((run / "miniprot_output" / "ref.mpi").exists())
            self.assertEqual(cb.main([str(roots[0]), "--apply"]), 0)
            self.assertFalse((run / "miniprot_output").exists())
            self.assertEqual(br.BuscoRun(run).layout, "condensed")
            self.assertEqual(cb.process(run, True)[1], "condensed")   # a second pass changes nothing


if __name__ == "__main__":
    unittest.main()
