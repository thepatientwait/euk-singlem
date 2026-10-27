"""Peer-similarity choice among duplicated BUSCO copies (no network)."""
import builtins
import random
import sys
import tempfile
import types
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS))
builtins.snakemake = types.SimpleNamespace()
import aggregate_busco_sequences as ab  # noqa: E402
import peer_copy_selection as pcs  # noqa: E402

AA = "ACDEFGHIKLMNPQRSTVWY"


def random_protein(rng, n=120):
    return "".join(rng.choice(AA) for _ in range(n))


def mutate(rng, seq, fraction):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < fraction:
            out[i] = rng.choice(AA)
    return "".join(out)


def info(genome, busco, copy, protein, status="Complete", score=100.0):
    return {"sequence": f">{genome}\n{protein}", "sequence_id": f"{genome}|{busco}|{copy}",
            "genome_acc": genome, "busco_id": busco, "copy_number": copy, "score": score,
            "status": status, "length": len(protein)}


def tax(family, genus):
    return f"d__Eukarya;p__P;c__C;o__O;f__{family};g__{genus};s__{genus} x"


class TestPeerSelection(unittest.TestCase):
    def setUp(self):
        rng = random.Random(1)
        self.ortholog = random_protein(rng)
        self.paralog = random_protein(rng)
        self.taxonomy = {}
        seqs = {}
        for i in range(5):  # five Complete peers, all orthologs
            genome = f"peer{i}"
            self.taxonomy[genome] = tax("F1", "G1")
            seqs[f"{genome}|g1|1"] = info(genome, "g1", 1, mutate(rng, self.ortholog, 0.1))
        # candidate genome: copy 1 is the paralog with the better BUSCO score, copy 2 the ortholog
        self.taxonomy["cand"] = tax("F1", "G1")
        for copy, (protein, score) in enumerate([(self.paralog, 200.0), (mutate(rng, self.ortholog, 0.05), 150.0)], 1):
            seqs[f"cand|g1|{copy}"] = info("cand", "g1", copy, protein, "Duplicated", score)
        self.genes = {"g1": seqs}

    def test_ortholog_chosen_over_better_score(self):
        selected, report = pcs.select_peer_copies(self.genes, self.taxonomy)
        chosen = selected["g1"]["cand|g1|1"]
        self.assertEqual(chosen["score"], 150.0)  # the ortholog, not the 200-score paralog
        self.assertEqual(report[0]["chosen_copy"], 2)
        self.assertEqual(report[0]["top_score_copy"], 1)
        self.assertEqual(report[0]["agrees_with_top_score"], 0)
        self.assertEqual(len(selected["g1"]), 6)  # 5 peers + one copy for cand
        self.assertEqual(len([k for k in selected["g1"] if k.startswith("cand|")]), 1)

    def test_no_peers_falls_back_to_top_score(self):
        del self.taxonomy["cand"]
        selected, report = pcs.select_peer_copies(self.genes, self.taxonomy)
        self.assertEqual(selected["g1"]["cand|g1|1"]["score"], 200.0)
        self.assertEqual(report[0]["n_peers"], 0)
        self.assertEqual(report[0]["chosen_identity"], "")

    def test_uses_deepest_level_with_enough_peers(self):
        # only two genus peers (< min_peers=3): the family level supplies the peers
        self.taxonomy["cand"] = tax("F1", "G1")
        for genome in ("peer2", "peer3", "peer4"):
            self.taxonomy[genome] = tax("F1", "Gother")
        _, report = pcs.select_peer_copies(self.genes, self.taxonomy)
        self.assertEqual(report[0]["peer_level"], 5)  # family
        self.assertEqual(report[0]["n_peers"], 5)

    def test_single_copy_genomes_untouched(self):
        selected, report = pcs.select_peer_copies(self.genes, self.taxonomy)
        self.assertEqual(selected["g1"]["peer0|g1|1"], self.genes["g1"]["peer0|g1|1"])
        self.assertEqual(len(report), 1)

    def test_deterministic(self):
        first = pcs.select_peer_copies(self.genes, self.taxonomy)[1]
        second = pcs.select_peer_copies(self.genes, self.taxonomy)[1]
        self.assertEqual(first, second)

    def test_identity_basics(self):
        self.assertAlmostEqual(pcs.identity("ACDEFGHIK", "ACDEFGHIK"), 1.0)
        self.assertEqual(pcs.identity("", "ACD"), 0.0)
        self.assertGreater(pcs.identity("ACDEFGHIK", "ACDEFGHIK"), pcs.identity("ACDEFGHIK", "WWWWWWWWW"))
        self.assertEqual(pcs.clean("acd*UOJ"), "ACDXXX")  # stops removed, unknown letters become X


class TestKeepAllCopies(unittest.TestCase):
    HEADER = "# BUSCO version is: 5.8.2\n# The lineage dataset is: eukaryota_odb12\n"
    COLUMNS = "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"

    def table(self, rows):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "busco" / "grp" / "GCA_1" / "run_x" / "full_table.tsv"
        path.parent.mkdir(parents=True)
        path.write_text(self.HEADER + self.COLUMNS + "".join("\t".join(str(v) for v in r) + "\n" for r in rows))
        return path

    def test_all_copies_kept_best_first_and_capped(self):
        rows = [("g1", "Complete", "c1", 1, 30, "+", 100.0, 10)]
        rows += [("g2", "Duplicated", f"c{i}", 1, 30, "+", float(i), 10) for i in range(1, 10)]
        _, _, df = ab.parse_busco_file(self.table(rows), True, keep_all_copies=True)
        g2 = df.filter(df["Busco id"] == "g2")
        self.assertEqual(g2.height, ab.MAX_CANDIDATE_COPIES)
        self.assertEqual(g2["Score"].to_list(), sorted(g2["Score"].to_list(), reverse=True))
        self.assertEqual(df.filter(df["Busco id"] == "g1").height, 1)

    def test_default_still_one_per_gene(self):
        rows = [("g2", "Duplicated", f"c{i}", 1, 30, "+", float(i), 10) for i in range(1, 4)]
        _, _, df = ab.parse_busco_file(self.table(rows), True)
        self.assertEqual(df.height, 1)


if __name__ == "__main__":
    unittest.main()


class TestMinGain(unittest.TestCase):
    def build(self, rng_seed=3, mutate_fraction=0.05):
        rng = random.Random(rng_seed)
        ortholog = random_protein(rng)
        genes, taxonomy, seqs = {}, {}, {}
        for i in range(5):
            genome = f"peer{i}"
            taxonomy[genome] = tax("F1", "G1")
            seqs[f"{genome}|g1|1"] = info(genome, "g1", 1, mutate(rng, ortholog, 0.1))
        taxonomy["cand"] = tax("F1", "G1")
        # both candidates are nearly equally similar to the peers; copy 1 has the better score
        seqs["cand|g1|1"] = info("cand", "g1", 1, mutate(rng, ortholog, mutate_fraction), "Duplicated", 200.0)
        seqs["cand|g1|2"] = info("cand", "g1", 2, mutate(rng, ortholog, mutate_fraction), "Duplicated", 150.0)
        return {"g1": seqs}, taxonomy

    def test_small_difference_keeps_top_score(self):
        genes, taxonomy = self.build()
        _, report = pcs.select_peer_copies(genes, taxonomy, min_gain=1.0)   # impossible margin
        self.assertEqual(report[0]["chosen_copy"], report[0]["top_score_copy"])
        self.assertEqual(report[0]["agrees_with_top_score"], 1)

    def test_margin_zero_follows_identity(self):
        genes, taxonomy = self.build()
        _, report = pcs.select_peer_copies(genes, taxonomy, min_gain=0.0)
        best = max(report[0]["chosen_identity"], report[0]["top_score_identity"])
        self.assertEqual(report[0]["chosen_identity"], best)
