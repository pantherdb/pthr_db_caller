import os
import csv
from typing import List
from lxml import etree
from dataclasses import dataclass
from pthr_db_caller.models import panther

PAINT_PMID = "PMID:21873635"
MOD_LIST = [
    ""
]
MOD_SUFFIXES = {
    'MOUSE': "mgi",
    'HUMAN': "human",
    'RAT': "rgd",
    'DROME': "fb",
    'ARATH': "tair",
    'CAEEL': "wb",
    'CHiCK': "chicken",
    'ECOLI': "ecocyc",
    'YEAST': "sgd",
    'DICDI': "dictyBase",
    'SCHPO': "pombase",
    'DANRE': "zfin",
    'CANAL': "cgd",
    'CHICK': "chicken"
}
DEFAULT_UNIPROT_ID_OUTPUT_FOR_PREFIXES = [
    "Gene",
    "Ensembl",
    "HGNC",
    "EcoGene"
]
PREFIX_TRANSFORMATIONS = {
    "WormBase": "WB",
    "FlyBase": "FB"
}
DEFAULT_RELATIONS = {
    "F": "enables",
    "P": "involved_in",
    "C": "is_active_in",
    "complex": "part_of"
}
ASPECT_LABEL_TO_SYMBOL = {
    "biological_process": "P",
    "molecular_function": "F",
    "cellular_component": "C"
}


def go_appropriate_id(long_id: panther.PthrSequence):
    # Extract from long ID, and perhaps external lookup tables, the ID expected by GO MODs and other consumers.
    #  Could be a MOD gene ID or UniProtKB depending on to-be-coded factors.
    gene_id = long_id.gene_id.replace("=", ":")
    gene_id_prefix, gene_id_suffix = gene_id.split(":", maxsplit=1)
    for prefix in DEFAULT_UNIPROT_ID_OUTPUT_FOR_PREFIXES:
        if gene_id_prefix.startswith(prefix):
            return long_id.uniprot.replace("=", ":")
    for prefix, new_prefix in PREFIX_TRANSFORMATIONS.items():
        if gene_id_prefix.startswith(prefix):
            gene_id = ":".join([new_prefix, gene_id_suffix])
    return gene_id


@dataclass
class WithBase:
    with_ids: List

    def go_appropriate_ids(self):
        return []


@dataclass
class ExperimentalWith(WithBase):
    with_ids: List[panther.PthrSequence]

    def go_appropriate_ids(self):
        return [go_appropriate_id(wid) for wid in self.with_ids]


@dataclass
class AncestralWith(WithBase):
    with_ids: List[str]  # A single PTN
    evidence_code: str

    def go_appropriate_ids(self):
        return ["PANTHER:{}".format(self.with_ids[0])]


def make_with(evidence_code: str, with_element: etree.Element):
    # with_ids could be list of gene_long_id elements (if IBD) or one persistent_id and one evidence_code (if IKR)
    if evidence_code in ["IKR"]:
        persistent_id = with_element.find("persistent_id").text
        with_evidence_code = with_element.find("evidence_code").text
        return AncestralWith([persistent_id], with_evidence_code)
    else:
        with_ids = [panther.PthrSequence(wi.text) for wi in with_element.getchildren()]
        return ExperimentalWith(with_ids)


@dataclass
class WithAnnotation:
    persistent_id: str
    evidence_code: str
    creation_date: str
    with_ids: WithBase

    @classmethod
    def from_element(WithAnnotation, element: etree.Element):
        persistent_id = element.find("persistent_id").text
        evidence_code = element.find("evidence_code").text
        creation_date = element.find("creation_date").text
        with_ids = make_with(evidence_code, element.find("with"))

        return WithAnnotation(persistent_id, evidence_code, creation_date, with_ids)


def determine_qualifiers(qualifier_elements: List[etree.Element]):
    qualifiers = []
    if qualifier_elements:
        negated = False
        for q_ele in qualifier_elements:
            qual = q_ele.text
            if qual == "NOT":
                negated = True
                # Don't append yet, add at the end
            else:
                qual = qual.lower()
                qualifiers.append(qual)
        # 'NOT' needs to be first
        if negated:
            qualifiers.insert(0, "NOT")

    # TODO: Set default qualifiers per GAF 2.2 spec

    return qualifiers


@dataclass
class Annotation:
    evidence_code: str
    term: str
    evidence_list: List[WithAnnotation]
    qualifiers: List[str]

    @classmethod
    def from_element(Annotation, element: etree.Element):
        evidence_code = element.find("evidence_code").text
        term = element.find("term").text
        evidence_list = element.find("evidence_list").getchildren()
        qualifiers = determine_qualifiers(element.findall("qualifier"))

        return Annotation(evidence_code, term, [WithAnnotation.from_element(we) for we in evidence_list], qualifiers)


@dataclass
class AnnotationCollection:
    annotations: List[Annotation]

    @classmethod
    def initial(AnnotationCollection):
        return AnnotationCollection([])

    def __iter__(self):
        return iter(self.annotations)

    def add(self, annotation: Annotation):
        self.annotations.append(annotation)

    def find_term(self, term):
        annotations = []
        for annotation in self:
            if annotation.term == term:
                annotations.append(annotation)
        return annotations


@dataclass()
class AnnotatedNode:
    persistent_id: str
    gene_long_id: panther.PthrSequence
    gene_name: str
    gene_symbol: str
    taxon_id: str
    annotations: AnnotationCollection

    @classmethod
    def from_element(AnnotatedNode, element: etree.Element):
        # Node will have fields like: persistent_id, gene_long_id, gene_name, gene_symbol, taxon_id
        persistent_id = element.find("persistent_id").text
        gene_long_id = panther.PthrSequence(element.find("gene_long_id").text)
        gene_name = element.find("gene_name").text
        gene_symbol = element.find("gene_symbol").text
        taxon_id = element.find("taxon_id").text
        annotations = AnnotationCollection.initial()

        return AnnotatedNode(persistent_id, gene_long_id, gene_name, gene_symbol, taxon_id, annotations)

# ibd_column_vals = [
#                         ibd_ptn,
#                         # "|".join(ibd_quals),
#                         sorted(ibd_quals),
#                         ibd_term,
#                         # "PMID:21873635",
#                         ibd_ev_code,
#                         # "|".join(ev.with_ids.go_appropriate_ids()),
#                         sorted(ev.with_ids),
#                         "taxon:",
#                         ev.creation_date
#                     ]
#
"""
Currently a combination of node and annotation 
"""
@dataclass
class Ibd:
    persistent_id: str
    qualifiers: List[str]
    term: str
    evidence_code: str
    evidence_with_ids: List[str]
    taxon_id: str
    creation_date: str

    @classmethod
    def from_list(Ibd, ibd_data: List):
        return Ibd(
            persistent_id=ibd_data[0],
            qualifiers=ibd_data[1],
            term=ibd_data[2],
            evidence_code=ibd_data[3],
            evidence_with_ids=ibd_data[4],
            taxon_id=ibd_data[5],
            creation_date=ibd_data[6]
        )


@dataclass
class AnnotatedNodeCollection:
    annotated_nodes: List[AnnotatedNode]

    @classmethod
    def initial(AnnotatedNodeCollection):
        return AnnotatedNodeCollection([])

    def add(self, annotated_node: AnnotatedNode):
        self.annotated_nodes.append(annotated_node)

    """
    Effectively "merges" two AnnotatedNodeCollection objects together
    """
    def merge_collection(self, collection):
        for annotated_node in collection:
            self.add(annotated_node)

    """
    Finds the AnnotatedNode by persistent_id (PTN). There should only be one.
    """
    def find_persistent_id(self, persistent_id: str):
        for annotated_node in self:
            if annotated_node.persistent_id == persistent_id:
                return annotated_node

    def __iter__(self):
        return iter(self.annotated_nodes)

    def __len__(self):
        return len(self.annotated_nodes)

    """
    Extracts IBD data from existing annotated_nodes
    """
    def ibd_nodes(self):
        # First, just collect lists of the data, ensuring they're unique
        ibd_data = []
        # Then, create Ibd objects that can be used by the IbaWriter (to access aspect, etc.)
        for anode in self.annotated_nodes:
            for a in anode.annotations:
                ibd_quals = a.qualifiers
                ibd_term = a.term
                for ev in a.evidence_list:
                    ibd_ptn = ev.persistent_id
                    ibd_ev_code = ev.evidence_code
                    ibd_column_vals = [
                        ibd_ptn,
                        sorted(ibd_quals),
                        ibd_term,
                        ibd_ev_code,
                        sorted(ev.with_ids.go_appropriate_ids()),
                        "taxon:",
                        ev.creation_date
                    ]
                    if ibd_column_vals not in ibd_data:
                        ibd_data.append(ibd_column_vals)
        ibd_objs = [Ibd.from_list(n) for n in ibd_data]
        return ibd_objs


class PaintIbaWriter:
    def __init__(self, go_aspect: str, complex_termlist: str, file_format: str = "GAF", obsolete_uniprots: str = None):
        self.go_aspects = self.parse_go_aspect(go_aspect)
        self.complex_terms = self.parse_list_file(complex_termlist)
        self.file_format = file_format
        self.obsolete_uniprot_ids = None
        if obsolete_uniprots:
            self.obsolete_uniprot_ids = self.parse_list_file(obsolete_uniprots)

    def get_aspect(self, term):
        try:
            aspect = self.go_aspects[term]
        except KeyError as ex:
            print("ERROR: Missing GO aspect for term '{}'".format(term))
            raise ex
        return aspect

    def get_qualifiers(self, qualifiers: List[str], term: str):
        # Start with assuming format is GAF 2.2
        aspect = self.get_aspect(term)
        # Determine default relation in case we need it
        if aspect == "C" and term in self.complex_terms:
            default_relation = DEFAULT_RELATIONS["complex"]
        else:
            default_relation = DEFAULT_RELATIONS[aspect]
        # Apply default relation if qualifiers are blank or only NOT
        if len(qualifiers) == 0 or qualifiers == ["NOT"]:
            qualifiers.append(default_relation)
        return qualifiers

    @staticmethod
    def parse_go_aspect(go_aspect):
        go_aspects = {}
        with open(go_aspect) as af:
            reader = csv.reader(af, delimiter="\t")
            for r in reader:
                term = r[0]
                aspect = r[1]
                if term and aspect:
                    go_aspects[term] = ASPECT_LABEL_TO_SYMBOL[aspect]
        return go_aspects

    @staticmethod
    def parse_list_file(list_file):
        # Generic for single-column, non-redundant lists
        list_of_things = set()
        with open(list_file) as lf:
            for l in lf.readlines():
                list_of_things.add(l.rstrip())
        return list_of_things

    def gaf_line(self, annotation: Annotation, annotated_node: AnnotatedNode):
        # GAF 2.2:
        # PomBase	SPAC959.04c	omh6	involved_in	GO:0006493	PMID:21873635	IBA	PANTHER:PTN000779407|SGD:S000002891|CGD:CAL0000188662|SGD:S000005625|SGD:S000000409
        # P	O-glycoside alpha-1,2-mannosyltransferase homolog 6	UniProtKB:Q9P4X2|PTN001258804	protein	taxon:284812	20170228	GO_Central
        # UniProtKB   F7HDM2  TFDP3   NOT|contributes_to      GO:0000977      PMID:21873635   IBA
        # PANTHER:PTN000284512|PANTHER:PTN000284480       F       Transcription factor    UniProtKB:F7HDM2|PTN000284516   protein taxon:9544      20200914        GO_Central
        s = go_appropriate_id(annotated_node.gene_long_id)  # MGI=MGI=12345 -> MGI:MGI:12345
        subject_prefix, subject_id = s.split(":", maxsplit=1)  # MGI:MGI:12345 -> MGI, MGI:12345

        # Use uniprot_id is no gene_symbol  # TODO: Fix PomBase symbols
        gene_symbol = annotated_node.gene_symbol
        if gene_symbol is None:
            gene_symbol = annotated_node.gene_long_id.uniprot_id

        qualifiers = self.get_qualifiers(annotation.qualifiers, annotation.term)
        qualifiers_str = "|".join(qualifiers)

        first_with = annotation.evidence_list[0]
        with_ptn = first_with.persistent_id  # PTN of IBD or IKR
        with_id_list = first_with.with_ids.with_ids
        if isinstance(first_with.with_ids, ExperimentalWith):
            with_ids = [go_appropriate_id(long_id) for long_id in with_id_list]
        else:
            # IKR with - PTN in list will be for IBD
            with_ids = ["PANTHER:{}".format(ptn) for ptn in with_id_list]

        aspect = self.get_aspect(annotation.term)

        return "\t".join([
            subject_prefix,
            subject_id,
            gene_symbol,
            qualifiers_str,
            annotation.term,
            PAINT_PMID,
            annotation.evidence_code,
            "|".join(["PANTHER:{}".format(with_ptn)] + with_ids),
            aspect,
            annotated_node.gene_name,
            "{}|{}".format(annotated_node.gene_long_id.uniprot.replace("=", ":"), annotated_node.persistent_id),
            "protein",
            "taxon:{}".format(annotated_node.taxon_id),
            first_with.creation_date,
            "GO_Central",
            "",  # Annotation Extension placeholder
            "",  # Gene Product Form ID placeholder
        ])

    def annotation_lines(self, annotated_node_collection: AnnotatedNodeCollection):
        lines = []
        for anode in annotated_node_collection:
            uniprot_id = anode.gene_long_id.uniprot_id
            if self.obsolete_uniprot_ids and uniprot_id in self.obsolete_uniprot_ids:
                print("\t".join(["Skipping - obsolete ID missing from latest uniprot_protein.gpi", "taxon:{}".format(anode.taxon_id), uniprot_id]))
                continue
            for annot in anode.annotations:
                if annot.evidence_code == "IBA":
                    if self.file_format.upper() == "GAF":
                        line = self.gaf_line(annot, anode)
                        lines.append(line)
        return lines

    def print(self, annotated_node_collection: AnnotatedNodeCollection):
        lines = self.annotation_lines(annotated_node_collection)
        for l in lines:
            print(l)

    def ibd_line(self, ibd_node: Ibd):
        return "\t".join([
            "PANTHER",
            ibd_node.persistent_id,
            ibd_node.persistent_id,
            "|".join(self.get_qualifiers(ibd_node.qualifiers, ibd_node.term)),
            ibd_node.term,
            "PMID:21873635",
            ibd_node.evidence_code,
            "|".join(ibd_node.evidence_with_ids),
            self.get_aspect(ibd_node.term),
            "",  # DB Object Name placeholder
            "",  # DB Object Synonym placeholder
            "protein",
            ibd_node.taxon_id,
            ibd_node.creation_date,
            "GO_Central",
            "",  # Annotation Extension placeholder
            "",  # Gene Product Form ID placeholder
        ])

    def ibd_lines(self, ibd_nodes: List[Ibd]):
        lines = []
        for node in ibd_nodes:
            line = self.ibd_line(node)
            lines.append(line)
        return lines


class PaintIbaXmlParser:
    PARSER = etree.XMLParser(recover=True)

    def extract_annotations(self, node: etree.Element, annotations: AnnotationCollection = None, is_leaf=None):
        if annotations is None:
            annotations = AnnotationCollection.initial()
        if not is_leaf:
            # Check
            is_leaf = node.tag == 'node' and 'children' not in [c.tag for c in node]
        for c in node:
            if c.tag == "persistent_id":
                pass
            elif c.tag == "annotation" and is_leaf:
                annotations.add(Annotation.from_element(c))
            else:
                self.extract_annotations(c, annotations, is_leaf=is_leaf)
        return annotations

    def parse_xml(self, xml_path: str):
        # This is where we have access to complex_terms and aspects
        annotated_node_collection = AnnotatedNodeCollection.initial()

        try:
            tree = etree.parse(xml_path, PaintIbaXmlParser.PARSER)
            node_list = tree.find("node_list")
            if node_list is not None:
                for node in node_list.getchildren():
                    anode = AnnotatedNode.from_element(node)
                    anode.annotations = self.extract_annotations(node)
                    annotated_node_collection.add(anode)
        except Exception as e:
            print(xml_path, e)

        return annotated_node_collection

    @staticmethod
    def parse(xml_path: str):
        parser = PaintIbaXmlParser()
        annotated_node_collection = AnnotatedNodeCollection.initial()

        # Sort out if xml_path is file or directory
        # TODO: Ensure these are .xml?
        if os.path.isdir(xml_path):
            xml_dir = xml_path
            xml_files = [os.path.join(xml_dir, xf) for xf in os.listdir(xml_path)]
        else:
            xml_files = [xml_path]

        for xf in xml_files:
            annotated_node_collection.merge_collection(parser.parse_xml(xf))

        return annotated_node_collection
