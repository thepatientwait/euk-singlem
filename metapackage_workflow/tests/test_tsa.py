"""TSA transcriptomes: discovery selection, download filtering, BUSCO transcriptome-mode tables (no network)."""
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
import busco_table as bt  # noqa: E402
import discover_tsa as dt  # noqa: E402
import download_tsa as dl  # noqa: E402
import extract_busco_transcripts as ex  # noqa: E402
from Bio import SeqIO  # noqa: E402

LISTING = """
<a href="tsa.GLPS.1.fsa_nt.gz">tsa.GLPS.1.fsa_nt.gz</a>    2026-09-14 15:47   16M  
<a href="tsa.GLPS.1.gbff.gz">tsa.GLPS.1.gbff.gz</a>    2026-09-14 15:47   30M  
<a href="tsa.GAAA.1.fsa_nt.gz">tsa.GAAA.1.fsa_nt.gz</a>    2026-09-14 15:47   700M  
<a href="tsa.GAAA.2.fsa_nt.gz">tsa.GAAA.2.fsa_nt.gz</a>    2026-09-14 15:47   1.5G  
<a href="tsa.GAAB.1.fsa_nt.gz">tsa.GAAB.1.fsa_nt.gz</a>    2026-09-14 15:47   512K  
"""


class TestDiscover(unittest.TestCase):
    def test_title_and_accession(self):
        self.assertEqual(dt.organism_from_title("TSA: Coolia malayensis strain DF364, transcriptome shotgun assembly"),
                         "Coolia malayensis")
        self.assertEqual(dt.organism_from_title("TSA: Fragilariopsis cylindrus, transcriptome shotgun assembly"),
                         "Fragilariopsis cylindrus")
        self.assertEqual(dt.split_accession("GLPS00000000.1"), ("GLPS", 1))
        self.assertEqual(dt.split_accession("GAAA00000000.2"), ("GAAA", 2))

    def test_listing_sizes_and_parts(self):
        parts = dt.parse_listing(LISTING)
        self.assertEqual(parts["GLPS"], [("tsa.GLPS.1.fsa_nt.gz", 16 * 1024 ** 2)])
        self.assertEqual([n for n, _ in parts["GAAA"]], ["tsa.GAAA.1.fsa_nt.gz", "tsa.GAAA.2.fsa_nt.gz"])
        self.assertEqual(parts["GAAB"][0][1], 512 * 1024)
        self.assertNotIn("gbff", str(parts))

    def records(self):
        return [
            {"accession": "GLPS00000000.1", "tax_id": "119754", "title": "TSA: Coolia malayensis strain DF364, transcriptome shotgun assembly"},
            {"accession": "GAAA00000000.1", "tax_id": "119754", "title": "TSA: Coolia malayensis strain X, transcriptome shotgun assembly"},
            {"accession": "GAAB00000000.1", "tax_id": "5", "title": "TSA: Tiny alga, transcriptome shotgun assembly"},
            {"accession": "GXXX00000000.1", "tax_id": "7", "title": "TSA: Covered species strain A, transcriptome shotgun assembly"},
            {"accession": "GYYY00000000.1", "tax_id": "8", "title": "TSA: No file here, transcriptome shotgun assembly"},
        ]

    def test_selection(self):
        listings = dt.parse_listing(LISTING)
        rows, dropped = dt.select_projects(self.records(), set(), {"covered species"}, listings,
                                           min_size_bytes=1024 ** 2)
        self.assertEqual([r["genome_acc"] for r in rows], ["GAAA00000000.1"])  # largest per tax ID
        self.assertEqual(rows[0]["parts"], "tsa.GAAA.1.fsa_nt.gz;tsa.GAAA.2.fsa_nt.gz")
        self.assertEqual(rows[0]["group"], "tsa")
        self.assertEqual(dropped["covered_by_genome"], 1)
        self.assertEqual(dropped["too_small"], 1)
        self.assertEqual(dropped["no_bulk_file"], 1)
        self.assertEqual(dropped["smaller_same_taxid"], 1)

    def test_covered_by_taxid_and_not_novel_only(self):
        listings = dt.parse_listing(LISTING)
        rows, _ = dt.select_projects(self.records(), {"119754"}, set(), listings, level="taxid")
        self.assertNotIn("GAAA00000000.1", [r["genome_acc"] for r in rows])
        rows, dropped = dt.select_projects(self.records(), {"119754"}, {"covered species"}, listings, only_novel=False)
        self.assertIn("GAAA00000000.1", [r["genome_acc"] for r in rows])
        self.assertNotIn("covered_by_genome", dropped)

    def test_lineage_exclusion(self):
        out = dt.excluded_by_lineage({"1": "Eukaryota; Opisthokonta; Metazoa; Annelida", "2": "Eukaryota; Sar; Alveolata",
                                      "3": "Eukaryota; Viridiplantae; Chlorophyta"})
        self.assertEqual(out, {"1", "3"})

    def test_read_covered_ignores_earlier_tsa_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            tax = Path(tmp) / "tax.tsv"
            tax.write_text("acc\tncbi_id\tspecies\nGCA_1.1\t10\tOchromonas species_2\nGLPS00000000.1\t119754\tComp x\n")
            ncbi = Path(tmp) / "ncbi.tsv"
            ncbi.write_text("genome_acc\ttax_id\toriginal_name\nGCA_2.1\t11\tHomo sapiens strain Z\n")
            taxids, species = dt.read_covered(tax, [ncbi], exclude_accs={"GLPS00000000.1"})
            self.assertEqual(taxids, {"10", "11"})
            self.assertEqual(species, {"ochromonas species", "homo sapiens"})


class TestDownload(unittest.TestCase):
    def test_filter_keeps_only_requested_version(self):
        with tempfile.TemporaryDirectory() as tmp:
            part1, part2, out = (Path(tmp) / n for n in ("a.gz", "b.gz", "out.gz"))
            with gzip.open(part1, "wt") as h:
                h.write(">GLPS01000001.1 x\nACGT\n>GLPS02000001.1 y\nTTTT\n")
            with gzip.open(part2, "wt") as h:
                h.write(">GLPS01000002.1 z\nGGGG\n")
            kept = dl.filter_version([part1, part2], out, "GLPS", 1)
            self.assertEqual(kept, 2)
            with gzip.open(out, "rt") as h:
                text = h.read()
            self.assertIn("GLPS01000001.1", text)
            self.assertIn("GLPS01000002.1", text)
            self.assertNotIn("GLPS02", text)

    def test_old_gi_headers(self):
        self.assertEqual(dl.header_accession(">gi|760317201|gb|GBYM01000001.1| TSA: Coriolopsis gallica cotig_1"), "GBYM01000001.1")
        self.assertEqual(dl.header_accession(">GLPS01000001.1 TSA: x"), "GLPS01000001.1")
        with tempfile.TemporaryDirectory() as tmp:
            part, out = Path(tmp) / "a.gz", Path(tmp) / "out.gz"
            with gzip.open(part, "wt") as h:
                h.write(">gi|1|gb|GBYM01000001.1| TSA: a\nACGT\n>gi|2|gb|GBYM02000001.1| TSA: b\nTT\n")
            self.assertEqual(dl.filter_version([part], out, "GBYM", 1), 1)


TABLE = (
    "# BUSCO version is: 6.1.0\n# The lineage dataset is: eukaryota_odb12\n"
    "# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription\n"
    "1at2759\tComplete\tc1:11-40\t100.0\t10\turl\td\n"
    "2at2759\tMissing\n"
    "3at2759\tDuplicated\tc2:5-34\t90.0\t10\turl\td\n"
    "3at2759\tDuplicated\tc3:7-36\t80.0\t10\turl\td\n"
)


class TestTranscriptomeMode(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        run = Path(tmp.name) / "outputs" / "busco" / "tsa" / "GLPS00000000.1" / "run_eukaryota_odb12"
        for sub in ("single_copy_busco_sequences", "multi_copy_busco_sequences"):
            (run / "busco_sequences" / sub).mkdir(parents=True)
        (run / "full_table.tsv").write_text(TABLE)
        sc, mc = run / "busco_sequences" / "single_copy_busco_sequences", run / "busco_sequences" / "multi_copy_busco_sequences"
        (sc / "1at2759.faa").write_text(">1at2759_1_0:aa|c1:11-40|+\nMKKLLAAGG*\n")
        (sc / "1at2759.fna").write_text(">1at2759_1_0:aa|c1:11-40|+\n" + "ATG" * 10 + "\n")
        (mc / "3at2759.faa").write_text(">3at2759_1_0:bb|c2:5-34|-\nMAAAAAAAAA\n>3at2759_1_0:cc|c3:7-36|+\nMCCCCCCCCC\n")
        self.table = run / "full_table.tsv"
        self.out = Path(tmp.name) / "out.fa"

    def test_table_normalised(self):
        df = bt.read_full_table(self.table)
        row = df.filter(df["Busco id"] == "1at2759").row(0, named=True)
        self.assertEqual((row["Sequence"], row["Gene Start"], row["Gene End"], row["Strand"], row["Transcriptome"]),
                         ("c1", 11, 40, "", True))
        self.assertIsNone(df.filter(df["Busco id"] == "2at2759")["Gene Start"][0])

    def test_genome_table_unchanged(self):
        import polars as pl
        df = pl.DataFrame({"Busco id": ["a"], "Status": ["Complete"], "Sequence": ["ctg"], "Gene Start": [1],
                           "Gene End": [9], "Strand": ["+"]})
        out = bt.normalise_full_table(df)
        self.assertEqual(out["Sequence"][0], "ctg")
        self.assertFalse(out["Transcriptome"][0])

    def test_aggregation_matches_records_without_strand(self):
        genes, err = agg.process_single_busco_file(str(self.table), use_multi_copy=True)
        self.assertIsNone(err)
        self.assertEqual(set(genes), {"1at2759", "3at2759"})
        one = next(iter(genes["1at2759"].values()))
        self.assertEqual(one["length"], 9)  # '*' stripped
        three = genes["3at2759"]
        self.assertEqual(len(three), 1)  # top-scoring copy of the Duplicated gene
        self.assertIn("MAAAAAAAAA", next(iter(three.values()))["sequence"])

    def test_extract_writes_transcript_cds(self):
        n = ex.extract_busco_transcripts("/does/not/exist.fasta.gz", self.table, self.out)
        self.assertEqual(n, 1)  # only the Complete gene
        rec = next(SeqIO.parse(self.out, "fasta"))
        self.assertEqual(rec.id, "1at2759|exon1of1")
        self.assertEqual(str(rec.seq), "ATG" * 10)


class TestQualityPenalty(unittest.TestCase):
    def test_transcriptome_penalty(self):
        import build_checkm2_quality_report as q
        groups = {"tsa"}
        self.assertEqual(q.adjusted_completeness("tsa", 90.0, 50, groups), 40.0)
        self.assertEqual(q.adjusted_completeness("tsa", 30.0, 50, groups), 0.0)
        self.assertEqual(q.adjusted_completeness("fungi", 90.0, 50, groups), 90.0)
        self.assertEqual(q.adjusted_completeness("tsa", 90.0, 0, groups), 90.0)


class TestMarkerCds(unittest.TestCase):
    def test_spliced_cds_orientation(self):
        import extract_marker_cds as mc
        import intron_utils as iu
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        contig = SeqRecord(Seq("AAACCCGGGTTTAAACCCGGGTTT"), id="c1")  # exons 1-6 and 13-18
        plus = iu.Model("c1", 1, 18, "+", [(1, 6), (13, 18)])
        self.assertEqual(str(mc.spliced_cds(plus, contig)), "AAACCC" + "AAACCC")
        minus = iu.Model("c1", 1, 18, "-", [(1, 6), (13, 18)])
        self.assertEqual(str(mc.spliced_cds(minus, contig)), "GGGTTT" + "GGGTTT")
        self.assertIsNone(mc.spliced_cds(iu.Model("c1", 1, 99, "+", [(1, 99)]), contig))


if __name__ == "__main__":
    unittest.main()
