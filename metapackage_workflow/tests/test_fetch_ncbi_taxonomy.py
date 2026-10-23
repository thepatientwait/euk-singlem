"""Tests for scripts/fetch_ncbi_taxonomy.py (stdlib unittest, no network).

Run from metapackage_workflow/ with an env that has polars:
    python -m unittest discover -s tests -v
"""
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import fetch_ncbi_taxonomy as ft  # noqa: E402
import polars as pl  # noqa: E402


def ds_record(tax_id, secondary=(), species_id=None, **ranks):
    classification = {k: {"name": v} for k, v in ranks.items()}
    if species_id and "species" in classification:
        classification["species"]["id"] = species_id
    return {"taxonomy": {"tax_id": tax_id, "secondary_tax_ids": list(secondary), "classification": classification}}


FULL = dict(superkingdom="Eukaryota", phylum="Ascomycota", **{"class": "Sordariomycetes"},
            order="Sordariales", family="Chaetomiaceae", genus="Chaetomium", species="Chaetomium globosum")


class TestCleaning(unittest.TestCase):
    def test_sp_rules(self):
        self.assertEqual(ft.clean_name("Chaetomium sp."), "Chaetomium species")
        self.assertEqual(ft.clean_name("Umbelopsis sp. PMI_123"), "Umbelopsis PMI_123")
        self.assertEqual(ft.clean_name("Foo bar (strain X1)"), "Foo bar")
        self.assertEqual(ft.clean_name("a:b;c"), "a_b_c")

    def test_placeholders_blanked(self):
        self.assertEqual(ft.clean_name("uncultured fungus"), "")
        self.assertEqual(ft.clean_name("Eukaryota incertae sedis"), "")
        self.assertEqual(ft.clean_name(None), "")

    def test_domain_forced(self):
        ranks = ft.clean_ranks({"domain": "Eukaryota", "phylum": "P"})
        self.assertEqual(ranks, {"phylum": "P", "domain": "Eukarya"})


class TestStatusAndFormat(unittest.TestCase):
    def test_status(self):
        full = ft.clean_ranks({r: r.upper() for r in ft.RANKS})
        self.assertEqual(ft.tax_status(full), "complete")
        self.assertEqual(ft.tax_status(ft.clean_ranks({"phylum": "P"})), "partial")
        self.assertEqual(ft.tax_status(ft.clean_ranks({})), "unresolved")

    def test_format_has_seven_ranks_no_spaces(self):
        s = ft.format_lineage(ft.clean_ranks({"phylum": "P", "species": "S x"}))
        self.assertEqual(s, "d__Eukarya;p__P;c__;o__;f__;g__;s__S x")


class TestRecordParsing(unittest.TestCase):
    def test_datasets(self):
        lineage, sec = ft.lineage_from_datasets_record(ds_record(5, [6, 7], species_id=4, **FULL))
        self.assertEqual((lineage.taxon_id, lineage.species_id, sec), (5, 4, [6, 7]))
        self.assertEqual(lineage.ranks["domain"], "Eukarya")
        self.assertEqual(lineage.ranks["species"], "Chaetomium globosum")
        self.assertEqual(ft.tax_status(lineage.ranks), "complete")

    def test_datasets_bad_record(self):
        self.assertIsNone(ft.lineage_from_datasets_record({"taxonomy": {}}))

    def test_entrez_uses_domain_rank_and_own_species(self):
        record = {"TaxId": "5", "Rank": "species", "ScientificName": "Chaetomium globosum", "AkaTaxIds": ["9"],
                  "LineageEx": [{"Rank": "domain", "ScientificName": "Eukaryota"},
                                {"Rank": "kingdom", "ScientificName": "Fungi"},
                                {"Rank": "phylum", "ScientificName": "Ascomycota"},
                                {"Rank": "genus", "ScientificName": "Chaetomium"}]}
        lineage, sec = ft.lineage_from_entrez_record(record)
        self.assertEqual((lineage.taxon_id, lineage.species_id, sec), (5, 5, [9]))
        self.assertEqual(lineage.ranks, {"domain": "Eukarya", "phylum": "Ascomycota", "genus": "Chaetomium",
                                         "species": "Chaetomium globosum"})

    def test_entrez_strain_record_takes_species_id_from_lineage(self):
        record = {"TaxId": "50", "Rank": "strain", "ScientificName": "Chaetomium globosum CBS 1",
                  "LineageEx": [{"Rank": "species", "TaxId": "5", "ScientificName": "Chaetomium globosum"}]}
        lineage, _ = ft.lineage_from_entrez_record(record)
        self.assertEqual((lineage.taxon_id, lineage.species_id), (50, 5))
        self.assertEqual(lineage.ranks["species"], "Chaetomium globosum")

    def test_entrez_non_species_record_does_not_invent_species(self):
        record = {"TaxId": "3", "Rank": "genus", "ScientificName": "Chaetomium",
                  "LineageEx": [{"Rank": "phylum", "ScientificName": "Ascomycota"}]}
        lineage, _ = ft.lineage_from_entrez_record(record)
        self.assertNotIn("species", lineage.ranks)
        self.assertEqual(lineage.ranks["genus"], "Chaetomium")


class TestSpeciesIndexing(unittest.TestCase):
    L = ft.Lineage

    def test_distinct_taxa_same_name_are_indexed(self):
        L = self.L
        lineages = {10: L(10, {"genus": "Och", "species": "Och species"}, 10),
                    30: L(30, {"genus": "Och", "species": "Och species"}, 30),
                    20: L(20, {"genus": "Och", "species": "Och real"}, 20)}
        out = ft.index_duplicate_species(lineages)
        self.assertEqual(out[10].ranks["species"], "Och species_1")
        self.assertEqual(out[30].ranks["species"], "Och species_2")
        self.assertEqual(out[20].ranks["species"], "Och real")

    def test_strains_of_one_species_not_indexed(self):
        L = self.L
        lineages = {101: L(101, {"species": "S cerevisiae"}, 4932), 102: L(102, {"species": "S cerevisiae"}, 4932)}
        out = ft.index_duplicate_species(lineages)
        self.assertEqual({v.ranks["species"] for v in out.values()}, {"S cerevisiae"})

    def test_secondary_id_of_same_taxon_not_indexed(self):
        L = self.L
        lin = L(10, {"species": "Foo bar"}, 10)
        out = ft.index_duplicate_species({10: lin, 11: lin})
        self.assertEqual({v.ranks["species"] for v in out.values()}, {"Foo bar"})

    def test_deterministic_regardless_of_input_order(self):
        L = self.L
        a = {2: L(2, {"species": "X y"}, 2), 1: L(1, {"species": "X y"}, 1)}
        b = {1: L(1, {"species": "X y"}, 1), 2: L(2, {"species": "X y"}, 2)}
        self.assertEqual(ft.index_duplicate_species(a)[1].ranks, ft.index_duplicate_species(b)[1].ranks)


class TestMapping(unittest.TestCase):
    def test_build_mapping_columns_and_missing(self):
        manifest = pl.DataFrame({"genome_acc": ["a", "b", "c"], "tax_id": [1, 1, 99]})
        lineages = ft.index_duplicate_species(
            {1: ft.Lineage(1, ft.clean_ranks({"phylum": "P", "species": "S x"}))})
        out = ft.build_mapping(manifest, lineages)
        self.assertEqual(out.columns, ft.MAPPING_COLUMNS)
        self.assertEqual(out["acc"].to_list(), ["a", "b"])  # 99 absent == missing
        self.assertEqual(out["tax_status"].to_list(), ["partial", "partial"])
        self.assertEqual(out["species"].to_list(), ["S x", "S x"])  # same taxon, not indexed
        self.assertEqual(out["parsed_taxon_string"][0], "d__Eukarya;p__P;c__;o__;f__;g__;s__S x")

    def test_resolve_falls_back_to_entrez_for_unreturned_ids(self):
        got_ds = {1: ft.Lineage(1, {"domain": "Eukarya", "phylum": "A"})}
        got_en = {2: ft.Lineage(2, {"domain": "Eukarya", "phylum": "B"})}
        with mock.patch.object(ft, "query_datasets", return_value=got_ds), \
                mock.patch.object(ft, "query_entrez", return_value=got_en) as en:
            out = ft.resolve_lineages([1, 2], "x@y")
        en.assert_called_once_with([2], "x@y")
        self.assertEqual(set(out), {1, 2})

    def test_resolve_survives_missing_datasets_cli(self):
        got_en = {1: ft.Lineage(1, {"domain": "Eukarya"})}
        with mock.patch.object(ft, "query_datasets", side_effect=FileNotFoundError("datasets")), \
                mock.patch.object(ft, "query_entrez", return_value=got_en):
            self.assertEqual(set(ft.resolve_lineages([1], "x@y")), {1})


if __name__ == "__main__":
    unittest.main()
