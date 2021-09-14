import unittest
from typing import List
from pthr_db_caller.models.panther import RefProtPantherMapping
from pthr_db_caller.models import paint, metadata, orthoxml
from pthr_db_caller.models.refprot_file import RefProtGeneAccFile, RefProtIdmappingFile, RefProtFastaFile
from pthr_db_caller.panther_tree_graph import PantherTreeGraph


class TestRefProtMapping(unittest.TestCase):

    def test_swissprot_status(self):
        current_mapping_path = "resources/test/refProteomePANTHERmapping_swissprot_status_test_current"
        mapping_w_status_path = "resources/test/refProteomePANTHERmapping_swissprot_status_test_data"
        current_mapping = RefProtPantherMapping.parse(current_mapping_path)
        mapping_w_status = RefProtPantherMapping.parse(mapping_w_status_path)

        for entry in current_mapping:
            status_entry = mapping_w_status.find_uniprot(entry.uniprot_id)
            entry.extras = status_entry.extras  # Copy status info

        self.assertEqual(current_mapping.entries[0].extras, ['tr'])


class TestRefProtPantherIdMapping(unittest.TestCase):
    def test_mgi(self):
        gene_2_acc_file = RefProtGeneAccFile.parse("resources/test/UP000000589_10090_MOUSE.gene2acc")
        idmapping_file = RefProtIdmappingFile.parse("resources/test/UP000000589_10090_MOUSE.idmapping")

        self.assertEqual(len(gene_2_acc_file.entries), 1000)
        self.assertEqual(len(idmapping_file.entries), 28236)

    def test_tair(self):
        pass


class TestRefProtFastaFile(unittest.TestCase):
    def test_reviewed_status(self):
        banana_fasta = "resources/test/UP000012960_214687.fasta"
        fasta_file = RefProtFastaFile.parse(banana_fasta)
        seq_to_status = {}
        for entry in fasta_file:
            seq_to_status[entry.uniprot_id] = entry.reviewed_status
        self.assertEqual(100, len(seq_to_status))
        self.assertEqual("tr", seq_to_status["M0RE52"])


class TestPantherTreeGraph(unittest.TestCase):
    def test_species_tree(self):
        tree = PantherTreeGraph.parse(tree_file="resources/test/species_pthr16_annot.nhx")
        n = tree.node("HUMAN")
        self.assertIsNotNone(n, "HUMAN leaf node not found")
        n = tree.node("Viridiplantae")
        self.assertIsNotNone(n, "Viridiplantae internal node not found")

    def test_family_tree(self):
        tree = PantherTreeGraph.parse(tree_file="resources/test/PTHR10000.tree")
        leaf_node = tree.node("AN96")
        self.assertEqual("LISMO|Gene=CAC98245|UniProtKB=Q8YAT3", leaf_node.get("long_id"),
                         msg="ID Q8YAT3 not found for node AN96")
        n = tree.node("AN98")
        self.assertEqual("Bacillus", n.get("species"), msg="Species Bacillus not found for node AN98")

    def test_pruning(self):
        tree = PantherTreeGraph.parse(tree_file="resources/test/PTHR10013.divided.tree.00")
        tree.prune_species(taxon_list=["STAA8", "BACCR"])  # The only two species in this tree
        self.assertEqual(len(tree), 5)  # 3 leaves + 2 internal = 5
        tree.prune_species(taxon_list=["HUMAN", "MOUSE"])
        self.assertEqual(len(tree), 0)  # Should prune all nodes in tree


class TestXmlToGaf(unittest.TestCase):
    ASPECT_FILE = "resources/test/go_aspects.tsv"
    COMPLEX_FILE = "resources/test/complex_terms.tsv"
    WRITER = paint.PaintIbaWriter(go_aspect=ASPECT_FILE, complex_termlist=COMPLEX_FILE)

    def test_not_w_contributes_to(self):
        xml_file = "resources/test/PTHR12548.xml"
        annotated_node_collection = paint.PaintIbaXmlParser.parse(xml_file)
        anode = annotated_node_collection.find_persistent_id("PTN002649017")
        # Find annotation to GO:0000977
        annot = anode.annotations.find_term("GO:0000977")[0]
        self.assertEqual(annot.qualifiers, ["NOT", "contributes_to"])

    def run_term_and_qualifiers_test(self, term: str, qualifiers: List, expected: List):
        annot = paint.Annotation(evidence_code="IBA", term=term, qualifiers=qualifiers, evidence_list=[])
        self.assertEqual(self.WRITER.get_qualifiers(annot.qualifiers, annot.term), expected)

    def test_default_output_relations(self):
        # regulation of transcription by RNA polymerase II (P)
        self.run_term_and_qualifiers_test(term="GO:0006357", qualifiers=["NOT"], expected=["NOT", "involved_in"])
        # kinase activity (F)
        self.run_term_and_qualifiers_test(term="GO:0016301", qualifiers=[], expected=["enables"])
        self.run_term_and_qualifiers_test(term="GO:0016301", qualifiers=["contributes_to"], expected=["contributes_to"])
        # GINS complex (complex)
        self.run_term_and_qualifiers_test(term="GO:0000811", qualifiers=[], expected=["part_of"])
        self.run_term_and_qualifiers_test(term="GO:0000811", qualifiers=["NOT", "colocalizes_with"], expected=["NOT", "colocalizes_with"])

    def test_iba_metadata_file_parse(self):
        iba_files = metadata.parse_iba_metadata_file("resources/test/paint_iba_files.tsv")
        self.assertEqual(len(iba_files), 13)


class TestOrthoXml(unittest.TestCase):
    def test_parse_pthr_orthoxml(self):
        # Parse preliminary-orthoXML format produced by divideHTtrees
        xml_file = "resources/test/orthoxml_pthr/PTHR21234.xml"
        groups = orthoxml.PthrOrthoXmlParser.parse(xml_file)
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups.genes), 7)

        xml_dir = "resources/test/orthoxml_pthr/"
        groups = orthoxml.PthrOrthoXmlParser.parse(xml_dir)
        self.assertEqual(len(groups), 3)  # Should equal number of input files?
        self.assertEqual(len(groups.genes), 133)

    def test_parse_orthoxml(self):
        xml_file = "resources/test/orthoxml/PTHR21234.divided.tree.00.nhx.xml"
        groups = orthoxml.PthrOrthoXmlParser.parse(xml_file)
        self.assertEqual(len(groups), 2)
        self.assertEqual(len(groups.genes), 7)

        xml_dir = "resources/test/orthoxml/"
        groups = orthoxml.PthrOrthoXmlParser.parse(xml_dir)
        self.assertEqual(len(groups), 6)
        self.assertEqual(len(groups.genes), 215)


if __name__ == "__main__":
    unittest.main()