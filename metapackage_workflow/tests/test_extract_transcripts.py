"""Transcript extraction cuts BUSCO genes at introns (synthetic genome, no BUSCO needed)."""
import builtins
import sys
import tempfile
import types
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
builtins.snakemake = types.SimpleNamespace()
import extract_busco_transcripts as ex  # noqa: E402
from Bio import SeqIO  # noqa: E402

HEADER = "# BUSCO version is: 6.1.0\n# The lineage dataset is: eukaryota_odb12\n"
COLUMNS = "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"


def gff_lines(contig, strand, exons):
    rows = [f"{contig}\tminiprot\tmRNA\t{min(s for s, _ in exons)}\t{max(e for _, e in exons)}\t100\t{strand}\t.\tID=M1;Rank=1"]
    for s, e in exons:
        rows.append(f"{contig}\tminiprot\tCDS\t{s}\t{e}\t50\t{strand}\t0\tParent=M1;Rank=1")
    return "\n".join(rows) + "\n"


class TestExtract(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        # contig of 600 nt: exon A 1-90, intron, exon B 201-290 (+ tiny 5 nt indel gap) 296-320, intron, exon C 401-430
        self.contig = "".join("ACGT"[i % 4] for i in range(600))
        self.genome = root / "g.fasta"
        self.genome.write_text(f">c1\n{self.contig}\n")
        run = root / "busco" / "fungi" / "G1" / "run_x"
        (run / "busco_sequences" / "single_copy_busco_sequences").mkdir(parents=True)
        self.table = run / "full_table.tsv"
        self.out = root / "out.fasta"

    def write(self, exons, strand="+"):
        s, e = min(a for a, _ in exons), max(b for _, b in exons)
        self.table.write_text(HEADER + COLUMNS + f"g1\tComplete\tc1\t{s - 1}\t{e - 3}\t{strand}\t100.0\t50\n")
        (self.table.parent / "busco_sequences" / "single_copy_busco_sequences" / "g1.gff").write_text(gff_lines("c1", strand, exons))

    def records(self, **kw):
        ex.extract_busco_transcripts(self.genome, self.table, self.out, **kw)
        return list(SeqIO.parse(self.out, "fasta"))

    def test_cut_at_each_intron_and_merge_short_gap(self):
        self.write([(1, 90), (201, 290), (296, 320), (401, 430)])
        recs = self.records(min_exon_nt=20)
        self.assertEqual([len(r.seq) for r in recs], [90, 120, 30])         # 201-320 merged across the 5 nt gap
        self.assertEqual(str(recs[0].seq), self.contig[0:90])
        self.assertEqual(str(recs[1].seq), self.contig[200:320])
        self.assertTrue(recs[0].id.endswith("exon1of3"))

    def test_minus_strand_is_reverse_complemented_in_transcript_order(self):
        self.write([(1, 90), (201, 290)], strand="-")
        recs = self.records(min_exon_nt=20)
        first_exon = self.contig[200:290]    # transcript starts at the highest coordinates
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        self.assertEqual(str(recs[0].seq), "".join(comp[b] for b in reversed(first_exon)))
        self.assertEqual(len(recs), 2)

    def test_short_blocks_are_dropped(self):
        self.write([(1, 90), (201, 230), (401, 490)])
        self.assertEqual([len(r.seq) for r in self.records(min_exon_nt=60)], [90, 90])

    def test_no_record_contains_an_intron(self):
        exons = [(1, 90), (201, 290), (401, 430)]
        self.write(exons)
        joined = "".join(str(r.seq) for r in self.records(min_exon_nt=1))
        self.assertEqual(joined, "".join(self.contig[s - 1:e] for s, e in exons))

    def test_gene_span_mode_keeps_the_introns(self):
        self.write([(11, 100), (201, 290), (401, 430)])
        recs = self.records(mode="gene_span")
        self.assertEqual(len(recs), 1)
        start, end = 10, 427                      # Gene Start / Gene End written by write()
        self.assertEqual(str(recs[0].seq), self.contig[start - 1:end])
        self.assertEqual(recs[0].id, "g1")
        self.assertGreater(len(recs[0].seq), 90 + 90 + 30)    # longer than the exons alone

    def test_gene_span_minus_strand_is_reverse_complemented(self):
        self.write([(11, 100), (201, 290)], strand="-")
        rec = self.records(mode="gene_span")[0]
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        self.assertEqual(str(rec.seq)[:5], "".join(comp[b] for b in reversed(self.contig[282:287])))

    def test_missing_gff_skips_the_gene(self):
        self.write([(1, 90)])
        (self.table.parent / "busco_sequences" / "single_copy_busco_sequences" / "g1.gff").unlink()
        self.assertEqual(self.records(), [])


if __name__ == "__main__":
    unittest.main()
