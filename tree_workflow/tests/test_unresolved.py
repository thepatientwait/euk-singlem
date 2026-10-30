"""Unresolved metapackage genomes are routed through the external genomes of the tree workflow."""
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import helpers  # noqa: E402

MAPPING = ("acc\tncbi_id\ttax_status\tparsed_taxon_string\n"
           "GCA_1.1\t1\tunresolved\td__Eukarya\n"
           "GCA_2.1\t2\tcomplete\td__Eukarya;p__A\n"
           "GCA_3.1\t3\tunresolved\td__Eukarya\n"
           "GLPS00000000.1\t4\tunresolved\td__Eukarya\n"
           "Fung1\t5\tunresolved\td__Eukarya\n"
           "Gone1\t6\tunresolved\td__Eukarya\n")
MANIFEST = ("outputs/genomes/p10k/GCA_1.1.fasta.gz\n"
            "outputs/genomes/p10k/GCA_2.1.fasta.gz\n"
            "outputs/genomes/ncbi/GCA_3.1.fasta.gz\n"
            "outputs/genomes/tsa/GLPS00000000.1.fasta.gz\n"
            "# comment\n"
            "outputs/genomes/fungi/Fung1.fasta.gz\n")


class TestUnresolved(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.mapping = self.root / "taxonomy_mapping.tsv"
        self.mapping.write_text(MAPPING)
        self.manifest = self.root / "genomes_to_process.paths"
        self.manifest.write_text(MANIFEST)

    def test_selects_unresolved_genomes_in_the_manifest_without_transcriptomes(self):
        found = helpers.unresolved_genome_paths(self.mapping, self.manifest, Path("/meta"))
        self.assertEqual(set(found), {"GCA_1.1", "GCA_3.1", "Fung1"})   # not complete, not tsa, not absent from the manifest
        self.assertEqual(found["GCA_1.1"], "/meta/outputs/genomes/p10k/GCA_1.1.fasta.gz")
        self.assertEqual(list(found), sorted(found))

    def test_transcriptome_groups_can_be_included(self):
        found = helpers.unresolved_genome_paths(self.mapping, self.manifest, Path("/meta"), exclude_groups=())
        self.assertIn("GLPS00000000.1", found)

    def test_missing_inputs_give_nothing(self):
        self.assertEqual(helpers.unresolved_genome_paths(self.root / "no.tsv", self.manifest, Path("/m")), {})
        self.assertEqual(helpers.unresolved_genome_paths(self.mapping, self.root / "no.paths", Path("/m")), {})

    def test_absolute_manifest_paths_are_kept(self):
        self.manifest.write_text("/abs/p10k/GCA_1.1.fasta.gz\n")
        self.assertEqual(helpers.unresolved_genome_paths(self.mapping, self.manifest, Path("/m")),
                         {"GCA_1.1": "/abs/p10k/GCA_1.1.fasta.gz"})


if __name__ == "__main__":
    unittest.main()
