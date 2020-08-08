import unittest
from pthr_db_caller.models.panther import RefProtPantherMapping
from ref_prot_id_mapping import RefProtGeneAccFile, RefProtIdmappingFile


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


if __name__ == "__main__":
    unittest.main()