import csv
from typing import List
from pthr_db_caller.config import BuildConfig


class PthrSequence:
    def __init__(self, long_id):
        self.long_id = long_id

        species_abbr, gene_id, uniprot = self.long_id.split("|")
        self.species_abbr = species_abbr
        self.gene_id = gene_id
        self.uniprot = uniprot

        self.uniprot_id = self.uniprot.split("=")[1]

    def __str__(self):
        return self.long_id

    def __eq__(self, other):
        if isinstance(other, PthrSequence):
            return self.long_id == other.long_id
        return False


class RefProtPantherMappingEntry:
    def __init__(self, uniprot_id, long_id: PthrSequence, symbol, family, description, mapping_method, old_long_id: PthrSequence = None):
        self.uniprot_id = uniprot_id
        self.long_id = long_id
        self.symbol = symbol
        self.family = family
        self.description = description
        self.old_long_id = old_long_id
        self.mapping_method = mapping_method

    @classmethod
    def from_row(cls, csv_row: List[str]):
        vals = [el.strip() for el in csv_row]

        # Mapping file is expected to be 7 cols long
        if len(vals) > 7:
            vals = vals[:7]
            
        uniprot_id, long_id, symbol, family, description, old_long_id, mapping_method = vals
        long_id = PthrSequence(long_id)
        if len(old_long_id) > 0:
            old_long_id = PthrSequence(old_long_id)
        entry = cls(uniprot_id, long_id, symbol, family, description, mapping_method, old_long_id=old_long_id)
        return entry

    def __str__(self):
        # Recreate mapping file format
        if not isinstance(self.old_long_id, PthrSequence):
            old_long_id = ""
        else:
            old_long_id = self.old_long_id.long_id
        line_elements = [
            self.uniprot_id,
            self.long_id.long_id,
            self.symbol,
            self.family,
            self.description,
            old_long_id,
            self.mapping_method
        ]
        return "\t".join(line_elements)


class RefProtPantherMapping:
    def __init__(self):
        self.entries = []  # [RefProtPantherMappingEntry]

    @classmethod
    def parse(cls, mapping_file):
        mapping_obj = cls()
        with open(mapping_file) as mf:
            reader = csv.reader(mf, delimiter="\t")
            for r in reader:
                mapping_obj.entries.append(RefProtPantherMappingEntry.from_row(r))
        return mapping_obj

    def find_uniprot(self, uniprot_id):
        for entry in self.entries:
            if entry.uniprot_id == uniprot_id:
                return entry

    def find_long_id(self, long_id: PthrSequence):
        for entry in self.entries:
            if entry.long_id == long_id:
                return entry

    def find_entries_by_base_family(self, family):
        entries = []
        for entry in self.entries:
            base_family = entry.family.split(":")[0]
            if base_family == family:
                entries.append(entry)
        return entries

    def add_entry(self, entry: RefProtPantherMappingEntry):
        self.entries.append(entry)

    def write(self, outfile):
        with open(outfile, "w+") as out_f:
            for entry in self.entries:
                out_f.write(str(entry) + "\n")

    def __iter__(self):
        return iter(self.entries)


class TaxonomyReadmeDetail:
    def __init__(self, proteome_id, tax_id, oscode, species_name):
        # Headers - Proteome_ID Tax_ID  OSCODE     #(1)    #(2)    #(3)  Species Name
        self.proteome_id = proteome_id
        self.tax_id = tax_id
        self.oscode = oscode
        if self.oscode == "None":
            self.oscode = None
        self.species_name = species_name


class TaxonomyDetails:
    # What data should this class handle?
    #  README_Reference_Proteome.txt
    #  README_QfO_release.txt
    #  RP_taxonomy_organism_lib##.txt

    def __init__(self, ref_prot_readme=None):
        if not ref_prot_readme:
            cfg = BuildConfig()
            ref_prot_readme = cfg.properties['README_Reference_Proteome']
        self.readme_tax_details = self.parse_readme(ref_prot_readme)

    @staticmethod
    def parse_readme(readme_file):
        details = []
        with open(readme_file) as rf:
            for l in rf.readlines():
                if l.startswith("UP"):
                    row = l.split()
                    detail = TaxonomyReadmeDetail(
                        proteome_id=row[0],
                        tax_id=row[1],
                        oscode=row[2],
                        species_name=row[6]
                    )
                    details.append(detail)

        return details

    @staticmethod
    def find_readme_detail(details: List[TaxonomyReadmeDetail], search_term, field):
        for d in details:
            if getattr(d, field) == search_term:
                return d

    def find_ref_prot_detail(self, search_term, field):
        return self.find_readme_detail(self.readme_tax_details, search_term, field)

    def find_by_oscode(self, oscode):
        return self.find_ref_prot_detail(oscode, 'oscode')

    def find_by_tax_id(self, tax_id):
        return self.find_ref_prot_detail(tax_id, 'tax_id')