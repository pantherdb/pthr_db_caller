import csv


class TaxonomyRecord:
    def __init__(self, taxon_id: int, organism: str = None, species_code: str = None, pthr_version=None,
                 mod_prefix: str = None, up_id: str = None):
        self.taxon_id = taxon_id
        # $taxId, $org, $symbol, $panther_version, $MOD, $uid
        self.organism = organism
        self.species_code = species_code
        self.pthr_version = pthr_version
        self.mod_prefix = mod_prefix
        self.up_id = up_id


class TaxonomyFile:
    # A TaxonomyFile contains TaxonomyRecords

    def __init__(self):
        self.records = []
        self.up_id_to_taxon_id = {}
        self.taxon_id_to_organism = {}
        self.taxon_id_to_species_code = {}
        self.taxon_id_to_mod_type = {}
        self.species_code_to_taxon_id = {}

    def add_record(self, record: TaxonomyRecord):
        self.records.append(record)
        self.up_id_to_taxon_id[record.up_id] = record.taxon_id
        self.taxon_id_to_organism[record.taxon_id] = record.organism
        self.taxon_id_to_species_code[record.taxon_id] = record.species_code
        if record.mod_prefix and record.mod_prefix != "":
            self.taxon_id_to_mod_type[record.taxon_id] = record.mod_prefix
        self.species_code_to_taxon_id[record.species_code] = record.taxon_id

    @staticmethod
    def parse_file(file_path):
        taxonomy_file = TaxonomyFile()
        with open(file_path) as tf:
            reader = csv.reader(tf, delimiter="\t")
            for r in reader:
                record = TaxonomyRecord(taxon_id=r[0],
                                        organism=r[1],
                                        species_code=r[2],
                                        pthr_version=r[3],
                                        mod_prefix=r[4],
                                        up_id=r[5])
                taxonomy_file.add_record(record)
        return taxonomy_file

    def __iter__(self):
        return iter(self.records)