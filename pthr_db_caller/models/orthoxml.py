import os
import io
from typing import List, Dict
from dataclasses import dataclass
from lxml import etree
from ete3 import orthoxml
from pthr_db_caller.models.panther import OrganismDatFile, PthrSequence


@dataclass
class Gene:
    orthoxml_id: str
    gene_id: str
    oscode: str

    @classmethod
    def from_element(Gene, element: etree.Element):
        orthoxml_id = element.attrib["id"]
        prot_id = element.attrib["protId"]
        oscode, gene_id = prot_id.split("_", maxsplit=1)
        return Gene(orthoxml_id=orthoxml_id, gene_id=gene_id, oscode=oscode)

    @classmethod
    def from_pthr_sequence(Gene, pthr_sequence: PthrSequence, orthoxml_id: str):
        gene_id = pthr_sequence.uniprot_id
        oscode = pthr_sequence.species_abbr
        return Gene(orthoxml_id=orthoxml_id, gene_id=gene_id, oscode=oscode)

@dataclass
class GeneCollection:
    genes: Dict = None
    species: Dict = None

    def add_gene(self, gene: Gene):
        # Track by OrthoXML ID
        if self.genes is None:
            self.genes = {}
        self.genes[gene.orthoxml_id] = gene
        # Track by species/oscode
        if self.species is None:
            self.species = {}
        if gene.oscode not in self.species:
            self.species[gene.oscode] = []
        self.species[gene.oscode].append(gene)

    def add_genes_from_species_element(self, species_element: etree.Element):
        db_element = species_element.find("database")
        genes_element = db_element.find("genes")
        for c in genes_element.getchildren():
            self.add_gene(Gene.from_element(c))

    def new_gene_from_long_id(self, pthr_long_id: str):
        pthr_seq = PthrSequence(pthr_long_id)
        next_orthoxml_id = self.max_orthoxml_id() + 1
        new_gene = Gene.from_pthr_sequence(pthr_seq, orthoxml_id=str(next_orthoxml_id))
        self.add_gene(new_gene)
        return new_gene

    def get_gene_by_orthoxml_id(self, orthoxml_id):
        return self.genes.get(orthoxml_id)

    def __len__(self):
        return len(self.genes)

    def max_orthoxml_id(self):
        if self.genes is None:
            return 0
        return int(sorted(self.genes.keys(), key=lambda x: int(x), reverse=True)[0])


@dataclass
class OrthoXmlGroup:
    genes: List[Gene] = None
    groups: List = None

    def add_gene(self, gene: Gene):
        if self.genes is None:
            self.genes = []
        self.genes.append(gene)

    def add_group(self, group):
        if self.groups is None:
            self.groups = []
        self.groups.append(group)

    def to_orthoxml_obj(self, parent_group=None):
        if parent_group is None:
            parent_group = orthoxml.group()
        if self.genes:
            for gene in self.genes:
                gene_ref = orthoxml.geneRef(gene.orthoxml_id)
                parent_group.add_geneRef(gene_ref)
        if self.groups:
            for g in self.groups:
                child_group = orthoxml.group()
                if isinstance(g, OrthologGroup):
                    g.to_orthoxml_obj(parent_group=child_group)
                    parent_group.add_orthologGroup(child_group)
                elif isinstance(g, ParalogGroup):
                    g.to_orthoxml_obj(parent_group=child_group)
                    parent_group.add_paralogGroup(child_group)
        return parent_group

    def __len__(self):
        num_genes = 0
        if self.genes:
            num_genes = len(self.genes)
        num_groups = 0
        if self.groups:
            num_groups = len(self.groups)
        return num_genes + num_groups


@dataclass
class OrthologGroup(OrthoXmlGroup):
    groups: List[OrthoXmlGroup] = None


@dataclass
class ParalogGroup(OrthoXmlGroup):
    groups: List[OrthoXmlGroup] = None


@dataclass
class GroupCollection:
    groups: List = None
    genes: GeneCollection = None

    def add_group(self, group: OrthoXmlGroup):
        if self.groups is None:
            self.groups = []
        self.groups.append(group)

    def remove_groups(self, groups: List[OrthoXmlGroup]):
        new_groups_list = []
        for g in self.groups:
            if g not in groups:
                # This group can stay
                new_groups_list.append(g)
        self.groups = new_groups_list

    def group_from_group_element(self, group_element: etree.Element):
        group = None
        if group_element.tag == "orthologGroup":
            group = OrthologGroup()
        elif group_element.tag == "paralogGroup":
            group = ParalogGroup()
        for c in group_element.getchildren():
            if c.tag in ["geneRef", "gene"]:
                if c.tag == "geneRef":
                    gene = self.genes.get_gene_by_orthoxml_id(c.attrib["id"])
                else:
                    gene = self.genes.new_gene_from_long_id(c.text)
                group.add_gene(gene)
            elif c.tag.endswith("Group"):
                group.add_group(self.group_from_group_element(c))
        return group

    def add_groups_from_groups_element(self, groups_element: etree.Element):
        for g in groups_element.getchildren():
            self.add_group(self.group_from_group_element(g))

    def merge_collection(self, collection):
        for group in collection.groups:
            self.add_group(group)

    def __len__(self):
        return len(self.groups)

    def __iter__(self):
        return iter(self.groups)

    def to_orthoxml_str(self, pthr_version: str, database_version: str, organism_dat: str = None):
        # Ready some element data
        oscode_taxid_lkp = {}
        if organism_dat:
            oscode_taxid_lkp = OrganismDatFile.parse_organism_dat(organism_dat)
        ### Write out compiled OrthoXML
        # Write out all genes
        oxml = orthoxml.orthoXML()
        oxml.set_version(0.3)
        oxml.set_origin("PANTHER")
        oxml.set_originVersion(pthr_version)
        for oscode, gene_list in self.genes.species.items():
            taxon_id = oscode_taxid_lkp.get(oscode)
            ortho_species = orthoxml.species(name=oscode, NCBITaxId=taxon_id)
            ortho_db = orthoxml.database(name="UniProt", version=database_version)
            ortho_genes = orthoxml.genes()
            for gene in gene_list:
                ortho_genes.add_gene(orthoxml.gene(protId=gene.gene_id, id=gene.orthoxml_id))
            ortho_db.set_genes(ortho_genes)
            ortho_species.add_database(ortho_db)
            oxml.add_species(ortho_species)

        # Write out all groups
        ortho_groups = orthoxml.groups()
        oxml.set_groups(ortho_groups)
        for g in self:
            ortho_groups.add_orthologGroup(g.to_orthoxml_obj())

        # print this to STDOUT
        out_str = io.StringIO()
        oxml.export(out_str, level=0, namespace_="", namespacedef_="xmlns=\"http://orthoXML.org/2011/\"")
        final_orthoxml_str = sanitize_xml_str(out_str.getvalue())
        return final_orthoxml_str


def sanitize_xml_str(xml_str: str):
    # Remove bytes syntax around strings. Ex: <gene protId=b'"ECOLI_P21829"' id="43"/>
    sanitized = xml_str.replace("b\'", "").replace("\'", "")
    sanitized = sanitized.replace("version=\"0.300000\"", "version=\"0.3\"")  # Very hacky, sorry
    return sanitized


class PthrOrthoXmlParser:
    @staticmethod
    def parse(xml_path: str):
        if os.path.isdir(xml_path):
            xml_files = [os.path.join(xml_path, xf_basename) for xf_basename in os.listdir(xml_path)]
        else:
            # Maybe don't allow single files cuz what's the point?
            xml_files = [xml_path]

        # Parse into Genes+Groups DS
        gene_id_index = 1
        all_genes = GeneCollection()
        all_groups = GroupCollection(genes=all_genes)
        for xf in xml_files:
            # Gotta fix ete3.orthoxml's bytes-encoding quirk (I think it's a python2 thing):
            xml_string = ""
            orthoxml.orthoXML()
            with open(xf) as xml_f:
                for l in xml_f.readlines():
                    xml_string += sanitize_xml_str(l)

            file_genes = GeneCollection()
            file_groups = GroupCollection(genes=file_genes)
            try:
                tree = etree.fromstring(xml_string, parser=etree.XMLParser(recover=True))
            except etree.XMLSyntaxError as e:
                # Some input files generated by divideHTtrees can be empty
                # TODO: Log out exception and xf (filename)
                continue
            if tree.tag != "orthoXML":
                # This is not standard orthoXML, assume it's output from divideHTtrees
                orthoxml_root = etree.Element("orthoXML")
                orthoxml_root.append(tree)
                tree = orthoxml_root
            for node in tree.getchildren():
                if node.tag == "species":
                    file_genes.add_genes_from_species_element(node)
            for node in tree.getchildren():
                if node.tag == "groups":
                    file_groups.add_groups_from_groups_element(node)
                elif node.tag in ["orthologGroup", "paralogGroup"]:
                    # This is coming from divideHTtrees output orthoXML
                    file_groups.add_group(file_groups.group_from_group_element(node))

            # Extra filter to remove singleton groups produced by etree2orthoxml.py
            groups_to_remove = []
            for g in file_groups:
                if len(g) < 2:
                    groups_to_remove.append(g)
            file_groups.remove_groups(groups_to_remove)

            # Remint orthoxml_ids to avoid collisions across files
            for orthoxml_id in sorted(file_genes.genes.keys(), key=int):
                gene = file_genes.genes[orthoxml_id]
                gene.orthoxml_id = str(gene_id_index)
                gene_id_index += 1
                all_genes.add_gene(gene)

            all_groups.merge_collection(file_groups)

        return all_groups
