from typing import List
from lxml import etree
from dataclasses import dataclass
from pthr_db_caller.models.panther import PthrSequence


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
DEFAULT_QUALIFIERS = {
    "F": "enables",
    "P": "involved_in",
    "C": "is_active_in",
    "complex": "part_of"
}


@dataclass
class WithBase:
    with_ids: List


@dataclass
class ExperimentalWith(WithBase):
    with_ids: List[PthrSequence]


@dataclass
class AncestralWith(WithBase):
    with_ids: List[str]  # A single PTN
    evidence_code: str


def make_with(evidence_code: str, with_element: etree.Element):
    # with_ids could be list of gene_long_id elements (if IBD) or one persistent_id and one evidence_code (if IKR)
    if evidence_code in ["IKR"]:
        persistent_id = with_element.find("persistent_id").text
        with_evidence_code = with_element.find("evidence_code").text
        return AncestralWith([persistent_id], with_evidence_code)
    else:
        with_ids = [PthrSequence(wi.text) for wi in with_element.getchildren()]
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
    gene_long_id: PthrSequence
    gene_name: str
    gene_symbol: str
    taxon_id: str
    annotations: AnnotationCollection

    @classmethod
    def from_element(AnnotatedNode, element: etree.Element):
        # Node will have fields like: persistent_id, gene_long_id, gene_name, gene_symbol, taxon_id
        persistent_id = element.find("persistent_id").text
        gene_long_id = PthrSequence(element.find("gene_long_id").text)
        gene_name = element.find("gene_name").text
        gene_symbol = element.find("gene_symbol").text
        taxon_id = element.find("taxon_id").text
        annotations = AnnotationCollection.initial()

        return AnnotatedNode(persistent_id, gene_long_id, gene_name, gene_symbol, taxon_id, annotations)


@dataclass
class AnnotatedNodeCollection:
    annotated_nodes: List[AnnotatedNode]

    @classmethod
    def initial(AnnotatedNodeCollection):
        return AnnotatedNodeCollection([])

    def add(self, annotated_node: AnnotatedNode):
        self.annotated_nodes.append(annotated_node)

    """
    Finds the AnnotatedNode by persistent_id (PTN). There should only be one.
    """
    def find_persistent_id(self, persistent_id: str):
        for annotated_node in self:
            if annotated_node.persistent_id == persistent_id:
                return annotated_node

    def __iter__(self):
        return iter(self.annotated_nodes)


def go_appropriate_id(long_id: PthrSequence):
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


def gaf_line(annotation: Annotation, annotated_node: AnnotatedNode):
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

    qualifiers = "|".join(annotation.qualifiers)

    first_with = annotation.evidence_list[0]
    with_ptn = first_with.persistent_id  # PTN of IBD or IKR
    with_id_list = first_with.with_ids.with_ids
    if isinstance(first_with.with_ids, ExperimentalWith):
        with_ids = [go_appropriate_id(long_id) for long_id in with_id_list]
    else:
        # IKR with - PTN in list will be for IBD
        with_ids = ["PANTHER:{}".format(ptn) for ptn in with_id_list]

    aspect = "P"  # TODO: calculate aspect

    return "\t".join([
        subject_prefix,
        subject_id,
        gene_symbol,
        qualifiers,
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


class PaintIbaXmlParser:
    PARSER = etree.XMLParser(recover=True)

    @staticmethod
    def extract_annotations(node: etree.Element, annotations: AnnotationCollection = None, is_leaf=None):
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
                PaintIbaXmlParser.extract_annotations(c, annotations, is_leaf=is_leaf)
        return annotations

    @staticmethod
    def parse(xml_path: str):
        annotated_node_collection = AnnotatedNodeCollection.initial()

        tree = etree.parse(xml_path, PaintIbaXmlParser.PARSER)
        node_list = tree.find("node_list")
        for node in node_list.getchildren():
            anode = AnnotatedNode.from_element(node)
            anode.annotations = PaintIbaXmlParser.extract_annotations(node)
            annotated_node_collection.add(anode)

        return annotated_node_collection