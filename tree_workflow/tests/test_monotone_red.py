"""monotone_red_tree: raise RED to the parent's, no negative branches (no PhyloRank needed)."""
import re
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import monotone_red_tree as mr  # noqa: E402

# root(0) -> A(RED .5) -> B(RED .4, below its parent) -> leaves x,y (RED 1); and leaf z under root
BAD = "(((x|RED=1.000:0.6,y|RED=1.000:0.6)0.9|RED=0.400:-0.1)0.9|RED=0.500:0.5,z|RED=1.000:1.0)|RED=0.000;"
GOOD = "(((x|RED=1.000:0.5,y|RED=1.000:0.5)0.9|RED=0.500:0.1)0.9|RED=0.400:0.4,z|RED=1.000:1.0)|RED=0.000;"


def run(text):
    tree = mr.read_tree(text)
    stats = mr.make_monotone(tree)
    return tree, stats


def lengths(tree):
    return [e.length for e in tree.preorder_edge_iter() if e.length is not None]


class TestMonotone(unittest.TestCase):
    def test_negative_branch_fixed(self):
        tree, stats = run(BAD)
        self.assertEqual(stats["negative_before"], 1)
        self.assertEqual(stats["raised"], 1)
        self.assertAlmostEqual(stats["max_raise"], 0.1)
        self.assertTrue(all(x >= 0 for x in lengths(tree)))
        labels = [n.label for n in tree.internal_nodes() if n.label]
        self.assertIn("0.9|RED_0.500", labels)           # B raised to its parent's .5
        self.assertNotIn("0.9|RED_0.400", labels)

    def test_child_red_never_below_parent(self):
        tree, _ = run(BAD)
        for node in tree.preorder_node_iter():
            if node.parent_node is not None:
                self.assertGreaterEqual(node.red_fixed, node.parent_node.red_fixed - 1e-12)

    def test_clean_tree_unchanged(self):
        tree, stats = run(GOOD)
        self.assertEqual((stats["negative_before"], stats["raised"]), (0, 0))
        self.assertAlmostEqual(sum(lengths(tree)), 0.4 + 0.1 + 0.5 + 0.5 + 1.0)

    def test_root_to_tip_distances_preserved_where_not_raised(self):
        tree, _ = run(BAD)
        x = [l for l in tree.leaf_node_iter() if l.taxon.label.startswith("x|")][0]
        self.assertAlmostEqual(sum(n.edge_length or 0 for n in [x] + list(x.ancestor_iter())), 1.0)

    def test_file_roundtrip_is_idempotent_and_keeps_raw(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "t.red_decorated.tree"
            path.write_text(BAD)
            mr.fix_file(path)
            first = path.read_text()
            self.assertEqual(mr.raw_path(path).read_text(), BAD)
            mr.fix_file(path)
            self.assertEqual(path.read_text(), first)
            self.assertNotIn("RED_", first)
            self.assertTrue(re.search(r"x\|RED=1\.000", first))
            self.assertFalse(any(float(v) < 0 for v in re.findall(r":(-?[\d.]+(?:e-?\d+)?)", first)))

    def test_scaled_tree_without_labels(self):
        text = "(((x:0.6,y:0.6):-0.1,z:0.0):0.5,w:1.0);"
        tree, stats = run(text)
        self.assertEqual(stats["negative_before"], 1)
        self.assertTrue(all(x >= 0 for x in lengths(tree)))


if __name__ == "__main__":
    unittest.main()


class TestCutoffCheck(unittest.TestCase):
    def _cut(self, tmp, low):
        p = Path(tmp) / "c.dict"
        p.write_text('{"phylum": %s, "class": 0.8}' % low)
        return p

    def test_passes_when_raised_nodes_are_far_below_cutoffs(self):
        _, stats = run(BAD)           # raised node ends at 0.5 of RED scale
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(mr.check_cutoffs(stats, self._cut(tmp, 0.9)), 0.8)

    def test_fails_when_a_raised_node_reaches_the_cutoffs(self):
        _, stats = run(BAD)
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit):
                mr.check_cutoffs(stats, self._cut(tmp, 0.52))

    def test_nothing_raised_nothing_checked(self):
        _, stats = run(GOOD)
        self.assertIsNone(mr.check_cutoffs(stats, "/nonexistent"))
