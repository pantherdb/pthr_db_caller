import os
import re
import csv
from typing import List
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.SeqRecord import SeqRecord
from pthr_db_caller.config import BuildConfig


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


class RefProtEntry:
    def __int__(self, uniprot_id: str, taxon_id: str):
        self.uniprot_id = uniprot_id
        self.taxon_id = taxon_id

    @classmethod
    def parse_row(cls, row: List, taxon_id: str):
        pass


class RefProtFile:
    ENTRY_TYPE = RefProtEntry

    def __init__(self, filename):
        self.filename = filename
        self.up_id, self.taxon_id, self.file_type = self.parse_filename(self.filename)
        self.entries: List[self.ENTRY_TYPE] = []
        self.by_uniprot = {}

    @staticmethod
    def parse_filename(filename):
        # e.g. UP000009136_9913.fasta
        up_id, taxon_id = os.path.basename(filename).split("_", maxsplit=1)
        taxon_id = re.match("[0-9]+", taxon_id).group(0)
        file_type = os.path.splitext(filename)[1].lstrip(".")
        return up_id, taxon_id, file_type

    def add_entry(self, entry: RefProtEntry):
        self.entries.append(entry)
        if entry.uniprot_id not in self.by_uniprot:
            self.by_uniprot[entry.uniprot_id] = []
        self.by_uniprot[entry.uniprot_id].append(entry)

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)


class RefProtFastaEntry(RefProtEntry):
    def __init__(self, uniprot_id: str, sequence: str, taxon_id: str, gene_name: str, definition: str, reviewed_status: str, seq_version: str):
        RefProtEntry.__int__(self, uniprot_id, taxon_id)
        self.sequence = sequence
        self.gene_name = gene_name
        self.definition = definition
        self.reviewed_status = reviewed_status
        self.seq_version = seq_version

    @staticmethod
    def parse_header(header_description):
        # tr|A0A0A0MP85|A0A0A0MP85_BOVIN Isoform of Q2YDN1, G protein-coupled receptor 161 OS=Bos taurus OX=9913 GN=GPR161 PE=3 SV=1
        ids_and_def, attributes = header_description.split("OS=")
        ids, definition = ids_and_def.split(" ", maxsplit=1)
        definition = definition.rstrip()
        reviewed_status = ids.split("|")[0]
        uniprot = ids.split("|")[1]

        # Pull out the key-value pairs of interest
        gene_name = None
        seq_version = None
        for rp_attr in attributes.split(" "):
            if "=" in rp_attr:
                attr_key, attr_val = rp_attr.split("=")
                if attr_key == "GN":
                    gene_name = attr_val
                elif attr_key == "SV":
                    seq_version = attr_val
        return uniprot, definition, gene_name, reviewed_status, seq_version

    @staticmethod
    def parse_fasta_record(record: SeqRecord, taxon_id: str):
        uniprot, definition, gene_name, reviewed_status, seq_version = RefProtFastaEntry.parse_header(record.description)
        return RefProtFastaEntry(uniprot, record.seq, taxon_id, gene_name, definition, reviewed_status, seq_version)


class RefProtFastaFile(RefProtFile):
    def __init__(self, filename: str):
        RefProtFile.__init__(self, filename)
        self.entries: List[RefProtFastaEntry] = []

    @staticmethod
    def parse(fasta_file):
        ref_prot_fasta_file = RefProtFastaFile(fasta_file)
        with open(ref_prot_fasta_file.filename) as ff:
            for record in FastaIterator(ff):
                entry = RefProtFastaEntry.parse_fasta_record(record, ref_prot_fasta_file.taxon_id)
                ref_prot_fasta_file.add_entry(entry)
        return ref_prot_fasta_file


class RefProtMappingEntry:
    def __init__(self, uniprot_id: str, taxon_id: str):
        self.uniprot_id = uniprot_id
        self.taxon_id = taxon_id

    @classmethod
    def parse_row(cls, row: List, taxon_id: str):
        pass


class RefProtMappingFile(RefProtFile):
    ENTRY_TYPE = RefProtMappingEntry

    @classmethod
    def parse(cls, filename: str):
        ref_prot_file = cls(filename)
        with open(ref_prot_file.filename) as gf:
            reader = csv.reader(gf, delimiter="\t")
            for r in reader:
                entry = cls.ENTRY_TYPE.parse_row(row=r, taxon_id=ref_prot_file.taxon_id)
                if entry:
                    ref_prot_file.add_entry(entry)
        return ref_prot_file


class IdMappingEntry(RefProtMappingEntry):
    def __init__(self, uniprot_id: str, source_type: str, gene_id: str, taxon_id: str):
        RefProtMappingEntry.__init__(self, uniprot_id, taxon_id)
        self.source_type = source_type
        self.gene_id = gene_id

    @classmethod
    def parse_row(cls, row: List, taxon_id: str):
        uniprot, source_type, gene_id = row
        if gene_id is None or gene_id == "" or source_type == "Araport":
            return None
        entry = cls(uniprot, source_type, gene_id, taxon_id)
        return entry

    def __str__(self):
        return "\t".join([self.uniprot_id, self.source_type, self.gene_id])


class RefProtIdmappingFile(RefProtMappingFile):
    ENTRY_TYPE = IdMappingEntry

    def __init__(self, filename):
        RefProtMappingFile.__init__(self, filename)
        self.by_source_type = {}

    def add_entry(self, entry: IdMappingEntry):
        RefProtMappingFile.add_entry(self, entry)
        if entry.source_type not in self.by_source_type:
            self.by_source_type[entry.source_type] = []
        self.by_source_type[entry.source_type].append(entry)


class GeneAccEntry(RefProtMappingEntry):
    TAXON_ID_TO_MOD_TYPE = {}

    def __init__(self, gene_id_1: str, uniprot_id: str, gene_id_2: str, taxon_id: str):
        RefProtMappingEntry.__init__(self, uniprot_id, taxon_id)
        # gene_id_1 - a unique gene symbol that is chosen with the following order of
        # preference from the annotation found in:
        # 1) Model Organism Database (MOD)
        # 2) Ensembl or Ensembl Genomes database
        # 3) UniProt Ordered Locus Name (OLN)
        # 4) UniProt Open Reading Frame (ORF)
        # 5) UniProt Gene Name
        self.gene_id_1 = gene_id_1

        # gene_id_2 - the gene symbol of the canonical accession used to represent the
        # respective gene group
        self.gene_id_2 = gene_id_2

    @classmethod
    def parse_row(cls, row: List, taxon_id):
        gene_id_1, uniprot, gene_id_2 = row
        if gene_id_1 in ["-", ""] or gene_id_1.isspace()\
                or gene_id_2 in ["-", ""] or gene_id_2.isspace():
            # These are blank so skip
            return None
        gene_id_1 = gene_id_1.lstrip().rstrip()
        gene_id_2 = gene_id_2.lstrip().rstrip()
        entry = GeneAccEntry(gene_id_1, uniprot, gene_id_2, taxon_id)
        return entry

    def __str__(self):
        return "\t".join([self.gene_id_1, self.uniprot_id, self.gene_id_2])


class RefProtGeneAccFile(RefProtMappingFile):
    ENTRY_TYPE = GeneAccEntry

    def __init__(self, filename):
        RefProtMappingFile.__init__(self, filename)
        self.by_gene_id = {}

    def add_entry(self, entry: GeneAccEntry):
        RefProtMappingFile.add_entry(self, entry)
        for gid in [entry.gene_id_1, entry.gene_id_2]:
            if gid not in self.by_gene_id:
                self.by_gene_id[gid] = []
            self.by_gene_id[gid].append(entry)


class RefProtFileSet:
    def __init__(self, taxon_id: str):
        self.taxon_id = taxon_id
        self.fasta : List[RefProtFastaFile] = []
        self.idmapping : RefProtIdmappingFile = None
        self.gene2acc : RefProtGeneAccFile = None

    def add_file(self, ref_prot_file: RefProtFile):
        if isinstance(ref_prot_file, RefProtFastaFile):
            self.fasta.append(ref_prot_file)
        elif isinstance(ref_prot_file, RefProtIdmappingFile):
            self.idmapping = ref_prot_file
        elif isinstance(ref_prot_file, RefProtGeneAccFile):
            self.gene2acc = ref_prot_file
