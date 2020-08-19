import csv
from typing import List, Dict
import logging

logger = logging.getLogger(__name__)


class PthrSequence:
    def __init__(self, long_id):
        self.long_id: str = long_id

        species_abbr, gene_id, uniprot = self.long_id.split("|")
        self.species_abbr: str = species_abbr
        self.gene_id: str = gene_id
        self.uniprot: str = uniprot

        self.uniprot_id: str = self.uniprot.split("=")[1]

    def __str__(self):
        return self.long_id

    def __eq__(self, other):
        if isinstance(other, PthrSequence):
            return self.long_id == other.long_id
        return False


class RefProtPantherMappingEntry:
    def __init__(self, uniprot_id, long_id: PthrSequence, symbol, family, description, mapping_method, old_long_id: PthrSequence = None, extras = None):
        self.uniprot_id = uniprot_id
        self.long_id = long_id
        self.symbol = symbol
        self.family = family
        self.description = description
        self.old_long_id = old_long_id
        self.mapping_method = mapping_method
        self.extras = extras  # any trailing fields - consumer scripts should know what to do with this

    @classmethod
    def from_row(cls, csv_row: List[str]):
        vals = [el.strip() for el in csv_row]

        # Mapping file is expected to be 7 cols long - but do save the 'extras'
        extras = None
        if len(vals) > 7:
            extras = vals[7:]
            vals = vals[:7]
            
        uniprot_id, long_id, symbol, family, description, old_long_id, mapping_method = vals
        long_id = PthrSequence(long_id)
        if len(old_long_id) > 0:
            old_long_id = PthrSequence(old_long_id)
        entry = cls(uniprot_id, long_id, symbol, family, description, mapping_method, old_long_id=old_long_id, extras=extras)
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
        if self.extras:
            line_elements = line_elements + self.extras
        return "\t".join(line_elements)


class RefProtPantherMapping:
    def __init__(self):
        self.entries: List[RefProtPantherMappingEntry] = []
        self.by_uniprot: Dict[str, List[RefProtPantherMappingEntry]] = {}

    @classmethod
    def parse(cls, mapping_file):
        mapping_obj = cls()
        with open(mapping_file) as mf:
            reader = csv.reader(mf, delimiter="\t")
            for r in reader:
                entry = RefProtPantherMappingEntry.from_row(r)
                if entry:
                    mapping_obj.add_entry(entry)
        return mapping_obj

    def find_uniprot(self, uniprot_id):
        if uniprot_id in self.by_uniprot:
            return self.by_uniprot[uniprot_id][0]
        return None
        # for entry in self.entries:
        #     if entry.uniprot_id == uniprot_id:
        #         return entry

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
        if entry.uniprot_id not in self.by_uniprot:
            self.by_uniprot[entry.uniprot_id] = []
        self.by_uniprot[entry.uniprot_id].append(entry)
        self.entries.append(entry)

    def write(self, outfile):
        with open(outfile, "w+") as out_f:
            for entry in self.entries:
                out_f.write(str(entry) + "\n")

    def __iter__(self):
        return iter(self.entries)


class DatEntry:
    def __init__(self, *args):
        pass

    @classmethod
    def parse_row(cls, row: List[str]):
        pass


class DatFile:
    ENTRY_TYPE = DatEntry

    def __init__(self, filename: str):
        self.filename = filename
        self.entries: List[self.ENTRY_TYPE] = []

    def add_entry(self, entry: DatEntry):
        self.entries.append(entry)

    @classmethod
    def parse(cls, filename: str):
        dat_file = cls(filename)
        with open(filename) as f:
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                entry = cls.ENTRY_TYPE.parse_row(r)
                if entry:
                    dat_file.add_entry(entry)
        return dat_file

    def __iter__(self):
        return iter(self.entries)

    def __len__(self):
        return len(self.entries)


class GeneDatEntry(DatEntry):
    def __init__(self, long_id: PthrSequence, description: str, synonym: str, mod_id: str, *args):
        super().__init__(self, long_id, description, synonym, mod_id)
        self.long_id = long_id
        self.description = description
        self.synonym = synonym
        self.mod_id = mod_id

    @classmethod
    def parse_row(cls, row: List[str]):
        pthr_sequence = PthrSequence(row[0])
        return cls(pthr_sequence, *row[1:4])

    def __str__(self):
        return "\t".join([str(self.long_id), self.description, self.synonym, self.mod_id])


class GeneDatFile(DatFile):
    ENTRY_TYPE = GeneDatEntry


class NodeDatEntry(DatEntry):
    def __init__(self, an_id: str, ptn: str, node_type: str, event_type: str, branch_length: float):
        super().__init__(self, an_id, ptn, node_type, event_type, branch_length)
        self.an_id = an_id
        self.ptn = ptn
        self.node_type = node_type
        self.event_type = event_type
        self.branch_length = branch_length

    @classmethod
    def parse_row(cls, row: List[str]):
        return cls(*row[0:4], float(row[4]))

    def __str__(self):
        event_type = self.event_type
        if event_type is None:
            # Will be blank if node_type == LEAF
            event_type = ""
        return "\t".join([self.an_id, self.ptn, self.node_type, event_type, self.branch_length])


class NodeDatFile(DatFile):
    ENTRY_TYPE = NodeDatEntry
