import unittest
from pthr_db_caller.models import RefProtPantherMapping, RefProtPantherMappingEntry
import csv

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


if __name__ == "__main__":
    unittest.main()