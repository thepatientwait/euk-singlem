"""Intron detection, window spanning and window choice (no network, no BUSCO)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
import intron_utils as iu  # noqa: E402


def gff(rows):
    lines = []
    for feature, start, end, strand, attrs in rows:
        lines.append("\t".join(["c1", "miniprot", feature, str(start), str(end), "10", strand, ".", attrs]))
    return "\n".join(lines) + "\n"


def model(exons, strand="+", irregular=False):
    return iu.Model("c1", min(s for s, _ in exons), max(e for _, e in exons), strand, list(exons), irregular)


class TestGff(unittest.TestCase):
    def test_parse_models_and_frameshift_flag(self):
        text = gff([
            ("mRNA", 100, 400, "+", "ID=M1"), ("CDS", 100, 150, "+", "Parent=M1"), ("CDS", 300, 400, "+", "Parent=M1;Frameshift=1"),
            ("mRNA", 900, 1000, "-", "ID=M2"), ("CDS", 900, 1000, "-", "Parent=M2"),
        ])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.gff"
            path.write_text(text)
            models = iu.parse_gff(path)
        self.assertEqual(len(models), 2)
        self.assertEqual(models[0].exons, [(100, 150), (300, 400)])
        self.assertTrue(models[0].irregular)
        self.assertFalse(models[1].irregular)

    def test_select_model_tolerates_stop_codon_shift(self):
        a, b = model([(100, 200)]), model([(500, 600)])
        self.assertIs(iu.select_model([a, b], "c1", 501, 597, "+"), b)   # BUSCO span ~3 nt shorter
        self.assertIsNone(iu.select_model([a, b], "c1", 1000, 1100, "+"))


class TestIntrons(unittest.TestCase):
    def test_plus_strand_offsets_and_gap(self):
        introns, cds = iu.intron_gaps(model([(100, 149), (300, 359), (500, 539)]))
        self.assertEqual(introns, [(50, 150), (110, 140)])     # offsets are exon-length sums
        self.assertEqual(cds, 150)

    def test_minus_strand_is_read_in_transcript_order(self):
        # same gene on the other strand: exons listed by genomic position, transcript runs high -> low
        introns, cds = iu.intron_gaps(model([(500, 539), (300, 359), (100, 149)], strand="-"))
        self.assertEqual(introns, [(40, 140), (100, 150)])
        self.assertEqual(cds, 150)

    def test_short_gap_is_not_an_intron_and_segments_merge(self):
        m = model([(100, 149), (155, 200), (400, 450)])        # 5 nt gap = indel, 199 nt gap = intron
        introns, _ = iu.intron_gaps(m)
        self.assertEqual([o for o, _ in introns], [96])
        self.assertEqual(iu.exon_segments(m), [(100, 200), (400, 450)])

    def test_segments_minus_strand(self):
        m = model([(400, 450), (155, 200), (100, 149)], strand="-")
        self.assertEqual(iu.exon_segments(m), [(400, 450), (100, 200)])

    def test_cds_matches_protein(self):
        m = model([(1, 30)])
        self.assertTrue(iu.cds_matches_protein(m, 30, 10))      # exact
        self.assertTrue(iu.cds_matches_protein(m, 33, 10))      # + stop codon
        self.assertTrue(iu.cds_matches_protein(m, 27, 10))      # one residue past the last CDS line
        self.assertFalse(iu.cds_matches_protein(m, 31, 10))     # not whole codons
        self.assertFalse(iu.cds_matches_protein(m, 60, 10))
        self.assertFalse(iu.cds_matches_protein(model([(1, 30)], irregular=True), 30, 10))


class TestSpan(unittest.TestCase):
    def test_boundaries(self):
        # window covers residues 10..29 => nt [30, 90)
        self.assertFalse(iu.spans_intron(10, 29, [30]))   # exactly at the window start
        self.assertTrue(iu.spans_intron(10, 29, [31]))    # inside the first codon
        self.assertTrue(iu.spans_intron(10, 29, [60]))
        self.assertTrue(iu.spans_intron(10, 29, [89]))
        self.assertFalse(iu.spans_intron(10, 29, [90]))   # exactly at the window end


class TestScan(unittest.TestCase):
    def make(self, introns):
        # 3 genomes, 40 columns, residue i at column i (no deletions), window of 20 columns
        n, ncol = len(introns), 40
        matrix = np.tile(np.arange(ncol, dtype=np.int32), (n, 1))
        return matrix, np.ones(n, dtype=bool)

    def test_scan_matches_scalar_definition(self):
        introns = [[], [90], [31], [200]]
        matrix, valid = self.make(introns)
        scan = iu.window_scan(matrix, introns, valid, 20)
        self.assertEqual(len(scan["n_cover"]), 21)
        for p in range(21):
            expected = sum(iu.spans_intron(p, p + 19, offs) for offs in introns)
            self.assertEqual(scan["n_span"][p], expected, f"start {p}")
        self.assertEqual(scan["n_cover"].tolist(), [4] * 21)

    def test_invalid_genomes_count_for_coverage_only(self):
        introns = [[31], [31]]
        matrix, _ = self.make(introns)
        scan = iu.window_scan(matrix, introns, np.array([True, False]), 20)
        self.assertEqual(scan["n_cover"][0], 2)
        self.assertEqual(scan["n_cover_valid"][0], 1)
        self.assertEqual(scan["n_span"][0], 1)

    def test_deletion_at_window_end_means_not_covered(self):
        matrix = np.tile(np.arange(40, dtype=np.int32), (1, 1))
        matrix[0, 19] = -1
        scan = iu.window_scan(matrix, [[]], np.ones(1, dtype=bool), 20)
        self.assertEqual(scan["n_cover"][0], 0)
        self.assertEqual(scan["n_cover"][1], 1)


class TestExtent(unittest.TestCase):
    def test_extent_with_internal_deletions_and_empty_rows(self):
        m = np.array([[-1, 0, 1, -1, 2, -1],       # spans columns 1..4, occupies 3
                      [0, 1, 2, 3, 4, 5],          # full
                      [-1, -1, -1, -1, -1, -1]], dtype=np.int32)
        first, last, n = iu.aligned_extent(m)
        self.assertEqual(first.tolist(), [1, 0, -1])
        self.assertEqual(last.tolist(), [4, 5, -1])
        self.assertEqual(n.tolist(), [3, 6, 0])


class TestChoose(unittest.TestCase):
    def scan(self, spans, scores, valid=100):
        n = len(spans)
        return {"n_span": np.array(spans), "aligned_bases": np.array(scores), "n_cover_valid": np.full(n, valid),
                "n_cover": np.full(n, valid)}

    def test_best_coverage_wins_when_clean(self):
        start, reason, violation, best = iu.choose_window(self.scan([0, 0, 0], [100, 120, 98]), 0.05)
        self.assertEqual((start, reason, violation, best), (1, "best_coverage", False, 1))

    def test_intron_crossing_best_window_is_skipped(self):
        scan = self.scan([20, 30, 1, 0], [100, 120, 90, 80])
        start, reason, violation, best = iu.choose_window(scan, 0.05)
        self.assertEqual((start, reason, violation, best), (2, "moved_for_introns", False, 1))

    def test_equally_covered_windows_prefer_fewer_introns(self):
        scan = self.scan([3, 1, 2], [100, 100, 100])
        start, reason, _, _ = iu.choose_window(scan, 0.05)
        self.assertEqual((start, reason), (1, "best_coverage_fewer_introns"))

    def test_near_best_coverage_counts_as_equal(self):
        scan = self.scan([4, 0], [1000, 996])           # 0.4% lower coverage, no introns
        self.assertEqual(iu.choose_window(scan, 0.05, score_tolerance=0.005)[0], 1)
        self.assertEqual(iu.choose_window(scan, 0.05, score_tolerance=0.001)[0], 0)

    def test_no_window_under_threshold_takes_least_spanning(self):
        scan = self.scan([40, 25, 30], [100, 50, 90])
        start, reason, violation, _ = iu.choose_window(scan, 0.05)
        self.assertEqual((start, violation, reason), (1, True, "no_window_under_threshold"))

    def test_too_little_data_falls_back_to_best_coverage(self):
        scan = self.scan([0, 0], [10, 20], valid=5)
        self.assertEqual(iu.choose_window(scan, 0.05, min_valid=20)[:3], (1, "insufficient_intron_data", False))


if __name__ == "__main__":
    unittest.main()
