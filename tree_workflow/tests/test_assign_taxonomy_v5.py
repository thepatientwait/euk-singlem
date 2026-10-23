"""Edge-case tests for scripts/assign_taxonomy_v5.py (stdlib unittest).

Run with the taxonomy conda env (dendropy + polars):
    python -m unittest discover -s tests -v
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import assign_taxonomy_v5 as at  # noqa: E402

CUTOFFS = {"phylum": 0.5, "class": 0.6, "order": 0.7, "family": 0.8, "genus": 0.9, "species": 0.95}


def full(p="", c="", o="", f="", g="", s=""):
    return f"d__Eukarya;p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}"


def run(newick, known, cutoffs=CUTOFFS):
    tree = at.read_red_tree(newick)
    ref, tgt, stats = at.assign_taxonomy(tree, known, cutoffs)
    parsed = {r["genome_id"]: at.parse_tax(r["taxonomy"]) for r in ref + tgt}
    return parsed, {r["genome_id"] for r in tgt}, stats


class TestMonotoneRed(unittest.TestCase):
    def test_child_raised_to_parent_only_when_lower(self):
        tree = at.read_red_tree("(((a:1,b:1)n3|RED=0.7:0.1,c:1)n2|RED=0.4:0.1,d:1)n1|RED=0.5:0.1;")
        for n in tree.postorder_node_iter():
            n.red = None if n.is_leaf() else at.parse_node_label(n)[1]
        raised = at.enforce_monotone_red(tree)
        reds = {n.label.split("|")[0]: n.red for n in tree.internal_nodes()}
        self.assertEqual(raised, 1)
        self.assertEqual(reds, {"n1": 0.5, "n2": 0.5, "n3": 0.7})

    def test_assign_reports_the_count(self):
        tree = at.read_red_tree("((a|RED_1:1,b|RED_1:1)x|RED_0.2:1,c|RED_1:1)root|RED_0.3:1;")
        _, _, stats = at.assign_taxonomy(tree, {"a": "d__Eukarya;p__P;c__;o__;f__;g__;s__"}, CUTOFFS)
        self.assertEqual(stats["red_raised_to_parent"], 1)


class TestNegativeBranches(unittest.TestCase):
    def test_counts_negative_edges(self):
        tree = at.read_red_tree("((a:0.1,b:-0.2)x|RED=0.5:-0.05,c:0.3)root;")
        self.assertEqual(at.report_negative_branches(tree), (2, -0.2))
        tree = at.read_red_tree("((a:0.1,b:0.2)x|RED=0.5:0.05,c:0.3)root;")
        self.assertEqual(at.report_negative_branches(tree), (0, 0.0))


class TestParsing(unittest.TestCase):
    def test_parse_and_format_roundtrip(self):
        t = at.parse_tax("d__Eukarya;p__A;c__;o__O;f__;g__G;s__G x")
        self.assertEqual(t, {"domain": "Eukarya", "phylum": "A", "order": "O", "genus": "G", "species": "G x"})
        self.assertEqual(at.format_tax(t), "d__Eukarya;p__A;c__;o__O;f__;g__G;s__G x")

    def test_empty(self):
        self.assertEqual(at.parse_tax(""), {})
        self.assertEqual(at.parse_tax(None), {})

    def test_red_label_parsing(self):
        tree = at.read_red_tree("((a|RED=1.000:0.1,b|RED=1.000:0.1)0.9|RED=0.95:0.1,c|RED=1.000:0.2)d__Eukarya|RED=0.000;")
        reds = sorted(at.parse_node_label(n)[1] for n in tree.internal_nodes())
        self.assertEqual(reds, [0.0, 0.95])
        self.assertEqual(sorted(at.parse_node_label(n)[0] for n in tree.leaf_node_iter()), ["a", "b", "c"])


class TestTargets(unittest.TestCase):
    def test_blank_taxonomy_row_is_target(self):
        """Bug 4: a row with nothing below domain is a target, not a reference."""
        known = {"a": full("P", "C", "O", "F", "G", "G a"), "b": full(), "c": ""}
        parsed, targets, _ = run("((a|RED=1:0.1,b|RED=1:0.1)|RED=0.96:0.1,c|RED=1:0.1)|RED=0.0;", known)
        self.assertEqual(targets, {"b", "c"})
        self.assertEqual(parsed["a"]["phylum"], "P")

    def test_every_leaf_fully_resolved(self):
        known = {"a": full("P", "C", "O", "F", "G", "G a")}
        parsed, _, _ = run("((a|RED=1:0.1,x|RED=1:0.1)|RED=0.96:0.1,y|RED=1:0.1)|RED=0.0;", known)
        for taxa in parsed.values():
            self.assertEqual(set(taxa), set(at.RANKS))


class TestConflictPropagation(unittest.TestCase):
    def test_mixed_is_sticky_upward(self):
        """Bug 2: N=(X, (A:P1, B:P2)) must not be labelled P1 with a target X beside the mixed clade."""
        known = {"A": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "B": full("P2", "C2", "O2", "F2", "G2", "G2 b")}
        # X is a target sister to the mixed clade (A,B); RED values below phylum cutoff at the top
        nwk = "(X|RED=1:0.1,(A|RED=1:0.1,B|RED=1:0.1)|RED=0.3:0.1)|RED=0.2;"
        parsed, targets, _ = run(nwk, known)
        self.assertEqual(targets, {"X"})
        self.assertNotIn(parsed["X"]["phylum"], {"P1", "P2"})

    def test_single_value_clade_labels_target_and_numbering_has_no_gaps(self):
        known = {"A": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "B": full("P1", "C1", "O1", "F1", "G1", "G1 b")}
        # (A,B) is species-mixed but deeper than the species cutoff, so it holds an
        # unused species placeholder; X must still get sp_1, not sp_2.
        nwk = "((A|RED=1:0.1,B|RED=1:0.1)|RED=0.97:0.1,X|RED=1:0.1)|RED=0.92:0.1;"
        parsed, targets, _ = run(nwk, known)
        self.assertEqual(parsed["X"]["genus"], "G1")
        self.assertEqual(parsed["X"]["species"], "G1 sp_1")


class TestInheritanceLeak(unittest.TestCase):
    def test_changed_rank_clears_lower_inherited_ranks(self):
        """Bug 3: a leaf with its own phylum must not keep the parent's class."""
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "b": full("P1", "C1", "O1", "F1", "G1", "G1 b"),
                 "z": "d__Eukarya;p__P2;c__;o__;f__;g__;s__"}
        nwk = "((a|RED=1:0.1,b|RED=1:0.1)|RED=0.99:0.1,z|RED=1:0.1)|RED=0.3;"
        parsed, _, _ = run(nwk, known)
        self.assertEqual(parsed["z"]["phylum"], "P2")
        self.assertNotEqual(parsed["z"]["class"], "C1")
        self.assertNotIn("C1", parsed["z"].values())

    def test_apply_state_unit(self):
        assigned = {"phylum": "P1", "class": "C1", "order": "O1"}
        at.apply_state(assigned, {"phylum": "P2"})
        self.assertEqual(assigned, {"phylum": "P2"})


class TestRedTruncation(unittest.TestCase):
    def test_below_phylum_cutoff_keeps_labels(self):
        """Decision 2: nodes older than the phylum cutoff keep database labels."""
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "b": full("P1", "C1", "O1", "F1", "G1", "G1 b")}
        parsed, _, _ = run("(a|RED=1:0.1,b|RED=1:0.1)|RED=0.1;", known)
        self.assertEqual(parsed["a"]["phylum"], "P1")
        self.assertEqual(parsed["b"]["genus"], "G1")

    def test_truncation_creates_new_family_for_target(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "b": full("P1", "C1", "O1", "F1", "G1", "G1 b")}
        # node RED 0.85 is above order(0.7) / family(0.8) cutoffs, below genus(0.9)
        nwk = "((a|RED=1:0.1,b|RED=1:0.1)|RED=0.97:0.1,t|RED=1:0.1)|RED=0.85;"
        parsed, targets, _ = run(nwk, known)
        self.assertEqual(parsed["t"]["family"], "F1")
        self.assertEqual(parsed["t"]["genus"], "F1_g1")
        self.assertEqual(parsed["t"]["species"], "F1_g1 sp_1")
        self.assertEqual(parsed["a"]["genus"], "G1")

    def test_missing_red_does_not_crash(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a")}
        parsed, _, stats = run("(a|RED=1:0.1,t|RED=1:0.1);", known)
        self.assertEqual(stats["internal_nodes_without_red"], 1)
        self.assertEqual(set(parsed["t"]), set(at.RANKS))


class TestHarmonise(unittest.TestCase):
    def test_majority_parent_deterministic_with_ties(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G", "G a"),
                 "b": full("P2", "C2", "O2", "F2", "G", "G b")}
        nwk = "(a|RED=1:0.1,b|RED=1:0.1)|RED=0.99;"
        first, _, _ = run(nwk, known)
        second, _, _ = run("(b|RED=1:0.1,a|RED=1:0.1)|RED=0.99;", known)
        self.assertEqual(first["a"]["family"], first["b"]["family"])
        self.assertEqual(first["a"]["family"], second["a"]["family"])  # tie broken alphabetically, order-independent
        self.assertEqual(first["a"]["family"], "F1")

    def test_blank_higher_rank_filled(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G", "G a"),
                 "b": "d__Eukarya;p__;c__;o__;f__;g__G;s__G b"}
        parsed, _, _ = run("(a|RED=1:0.1,b|RED=1:0.1)|RED=0.99;", known)
        self.assertEqual(parsed["b"]["family"], "F1")
        self.assertEqual(parsed["b"]["phylum"], "P1")


class TestPlaceholders(unittest.TestCase):
    def test_no_collision_with_database_names(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "b": full("P1", "C1", "O1", "F1", "G1", "G1 b"),
                 "c": full("P1", "C1", "O1", "F1", "G1", "G1 sp_1")}
        nwk = "((a|RED=1:0.1,b|RED=1:0.1,t|RED=1:0.1)|RED=0.5:0.1,c|RED=1:0.1)|RED=0.3;"
        parsed, _, _ = run(nwk, known)
        self.assertEqual(parsed["t"]["genus"], "G1")
        self.assertEqual(parsed["t"]["species"], "G1 sp_2")

    def test_unique_species_per_target(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "b": full("P1", "C1", "O1", "F1", "G1", "G1 b")}
        nwk = "(((a|RED=1:0.1,b|RED=1:0.1)|RED=0.5:0.1,t1|RED=1:0.1)|RED=0.4:0.1,t2|RED=1:0.1)|RED=0.3;"
        parsed, targets, _ = run(nwk, known)
        self.assertEqual({parsed["t1"]["species"], parsed["t2"]["species"]}, {"G1 sp_1", "G1 sp_2"})

    def test_nested_placeholders_resolve_parent_first(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a"),
                 "c": full("P2", "C2", "O2", "F2", "G2", "G2 c")}
        # root is phylum-mixed; the (t1,t2) clade is deeper than every cutoff -> a whole new lineage
        nwk = "((a|RED=1:0.1,c|RED=1:0.1)|RED=0.1:0.1,(t1|RED=1:0.1,t2|RED=1:0.1)|RED=0.99:0.1)|RED=0.0;"
        parsed, _, _ = run(nwk, known)
        t1, t2 = parsed["t1"], parsed["t2"]
        for rank in at.RANKS[2:]:
            self.assertTrue(t1[rank].startswith(t1[at.RANKS[at.RANKS.index(rank) - 1]]), (rank, t1))
        self.assertEqual(t1["genus"], t2["genus"])
        self.assertEqual(t1["species"], t2["species"])  # RED 0.99 > species cutoff: same species

    def test_deep_ladder_tree_no_recursion_error(self):
        """Processing must not rely on recursion (dendropy's newick parser does, so build via API)."""
        import dendropy
        n = 3000
        tree = dendropy.Tree()
        node = tree.seed_node
        for i in range(n):
            leaf = node.new_child(label=f"t{i}|RED_1.000")
            leaf.taxon = None
            nxt = node.new_child(label="x|RED_0.5")
            node = nxt
        node.label = "tlast|RED_1.000"
        known = {"t0": full("P", "C", "O", "F", "G", "G a")}
        ref, tgt, _ = at.assign_taxonomy(tree, known, CUTOFFS)
        self.assertEqual(len(ref) + len(tgt), n + 1)


class TestDeterminism(unittest.TestCase):
    def test_repeatable(self):
        known = {"a": full("P1", "C1", "O1", "F1", "G1", "G1 a")}
        nwk = "((a|RED=1:0.1,t1|RED=1:0.1)|RED=0.96:0.1,t2|RED=1:0.1)|RED=0.5;"
        self.assertEqual(run(nwk, known)[0], run(nwk, known)[0])


if __name__ == "__main__":
    unittest.main()
