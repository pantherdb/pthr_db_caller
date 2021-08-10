#!/usr/bin/python3

import argparse
import csv
import os
import io
from typing import List, Dict
from dataclasses import dataclass
from lxml import etree
from ete3 import orthoxml


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--xml_path', help="Path to directory containing OrthoXML files to merge")
parser.add_argument('-p', '--pthr_version', help="PANTHER version from where the input files originate")
parser.add_argument('-d', '--database_version', help="DB where gene IDs were minted")
parser.add_argument('-o', '--organism_dat', help="Oscode-to-taxonID lookup from PANTHER build process")


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

    def __len__(self):
        return len(self.genes)


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

    def to_orthoxml(self, parent_group=None):
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
                    g.to_orthoxml(parent_group=child_group)
                    parent_group.add_orthologGroup(child_group)
                elif isinstance(g, ParalogGroup):
                    g.to_orthoxml(parent_group=child_group)
                    parent_group.add_paralogGroup(child_group)
        return parent_group


@dataclass
class OrthologGroup(OrthoXmlGroup):
    groups: List[OrthoXmlGroup] = None


@dataclass
class ParalogGroup(OrthoXmlGroup):
    groups: List[OrthoXmlGroup] = None


@dataclass
class GroupCollection:
    groups: List = None
    gene_collection: GeneCollection = None

    def add_group(self, group: OrthoXmlGroup):
        if self.groups is None:
            self.groups = []
        self.groups.append(group)

    def group_from_group_element(self, group_element: etree.Element):
        group = None
        if group_element.tag == "orthologGroup":
            group = OrthologGroup()
        elif group_element.tag == "paralogGroup":
            group = ParalogGroup()
        for c in group_element.getchildren():
            if c.tag == "geneRef":
                gene = self.gene_collection.genes.get(c.attrib["id"])
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


def sanitize_xml_str(xml_str: str):
    # Remove bytes syntax around strings. Ex: <gene protId=b'"ECOLI_P21829"' id="43"/>
    sanitized = xml_str.replace("b\'", "").replace("\'", "")
    sanitized = sanitized.replace("version=\"0.300000\"", "version=\"0.3\"")  # Very hacky, sorry
    return sanitized


def parse_organism_dat(organism_dat_path: str):
    oscode_taxon_lkp = {}
    with open(organism_dat_path) as odf:
        reader = csv.reader(odf, delimiter="\t")
        for r in reader:
            oscode_taxon_lkp[r[2]] = r[5]
    return oscode_taxon_lkp


if __name__ == "__main__":
    args = parser.parse_args()

    if os.path.isdir(args.xml_path):
        xml_files = [os.path.join(args.xml_path, xf_basename) for xf_basename in os.listdir(args.xml_path)]
    else:
        # Maybe don't allow single files cuz what's the point?
        xml_files = [args.xml_path]

    # Parse into Genes+Groups DS
    gene_id_index = 1
    all_genes = GeneCollection()
    all_groups = GroupCollection()
    for xf in xml_files:
        # Gotta fix ete3.orthoxml's bytes-encoding quirk (I think it's a python2 thing):
        xml_string = ""
        orthoxml.orthoXML()
        with open(xf) as xml_f:
            for l in xml_f.readlines():
                xml_string += sanitize_xml_str(l)

        file_genes = GeneCollection()
        file_groups = GroupCollection(gene_collection=file_genes)
        tree = etree.fromstring(xml_string, parser=etree.XMLParser(recover=True))
        for node in tree.getchildren():
            if node.tag == "species":
                file_genes.add_genes_from_species_element(node)
        for node in tree.getchildren():
            if node.tag == "groups":
                file_groups.add_groups_from_groups_element(node)

        # Remint orthoxml_ids to avoid collisions across files
        for orthoxml_id in sorted(file_genes.genes.keys(), key=int):
            gene = file_genes.genes[orthoxml_id]
            gene.orthoxml_id = str(gene_id_index)
            gene_id_index += 1
            all_genes.add_gene(gene)

        all_groups.merge_collection(file_groups)

    # Ready some element data
    pthr_version = args.pthr_version
    database_version = args.database_version
    oscode_taxid_lkp = {}
    if args.organism_dat:
        oscode_taxid_lkp = parse_organism_dat(args.organism_dat)
    ### Write out compiled OrthoXML
    # Write out all genes
    oxml = orthoxml.orthoXML()
    oxml.set_version(0.3)
    oxml.set_origin("PANTHER")
    oxml.set_originVersion(pthr_version)
    for oscode, gene_list in all_genes.species.items():
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
    for g in all_groups.groups:
        ortho_groups.add_orthologGroup(g.to_orthoxml())

    # print this to STDOUT
    out_str = io.StringIO()
    oxml.export(out_str, level=0, namespace_="")
    print(sanitize_xml_str(out_str.getvalue()))
