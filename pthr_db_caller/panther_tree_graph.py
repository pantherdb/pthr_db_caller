import networkx
from networkx import MultiDiGraph
from typing import List, Dict
from Bio import Phylo
from io import StringIO
from pthr_db_caller.models.panther import NodeDatFile


# Unfortunately, this only uses AN# node IDs instead of PTNs due to parsing from tree files.
# If we parsed the enormous node.dat we could add a PTN option.

class PantherTreePhylo:

    def __init__(self, tree_file):
        # Load graph from tree file
        with open(tree_file) as tf:
            tree_line = tf.readline()
            tree_string = StringIO(tree_line)
            # tree_phylo = next(PantherNewickIOParser(tree_string).parse())
            tree_phylo = next(Phylo.parse(tree_string, "newick"))
            # Leaves parse clean due to not having species name in 'S:'

        self.tree = tree_phylo


def extract_clade_name(clade_comment):
    if clade_comment is None:
        clade_comment = ""
    # Name-parsing should be more robust
    new_comment = ""
    comment_bits = clade_comment.split(":")
    for b in comment_bits:
        if b.startswith("S="):
            new_comment = b.replace("S=", "")
            break
    # Also grab ID
    an_id = ""
    for b in comment_bits:
        if b.startswith("ID="):
            an_id = b.replace("ID=", "")
            break
    ###
    new_comment = new_comment.replace("&&NHX:S=", "")
    new_comment = new_comment.replace("&&NXH:S=", "")
    if new_comment == "Opisthokonts":
        new_comment = "Opisthokonta"
    return new_comment, an_id


class PantherTreeGraph:
    def __init__(self, tree_name: str = None):
        self.graph = MultiDiGraph()
        self.name: str = tree_name
        self.ptn_to_an: Dict = None

    # Recursive method to fill graph from Phylo clade, parsing out node accession and species name (if present)
    def add_children(self, parent_clade):
        self.add_node_from_clade(parent_clade)
        for child in parent_clade.clades:
            self.add_node_from_clade(child)
            self.graph.add_edge(parent_clade.name, child.name)

            if len(child.clades) > 0:
                self.add_children(child)

    def add_node_from_clade(self, clade):
        # Few cases:
        #  1. Leaf - no comment; name=AN#
        #       Add node child.name to graph
        #  2. Internal - AN# in comment; name not set; no species name in comment
        #       Parse ID from comment; set name=AN#; Add node child.name to graph
        #  3. Internal - AN# in comment; name not set; species name in comment
        #       Parse ID and species from comment; set name=AN#; Add node child.name to graph; Add node
        #       property of species
        species, id = extract_clade_name(clade.comment)
        if clade.name is None:
            clade.name = id
        if clade.name not in self.graph.nodes():
            self.graph.add_node(clade.name)
        if species:
            self.graph.nodes[clade.name]["species"] = species

    def extract_leaf_ids(self, tree_file):
        with open(tree_file) as tf:
            tf.readline()  # ignore first line, it's parsed already
            for l in tf.readlines():
                an_id, long_id = l.split(":", maxsplit=1)
                long_id = long_id.rstrip().rstrip(";")
                if an_id in self.graph:
                    self.graph.node[an_id]["long_id"] = long_id

    def extract_node_properties(self, node_file: str):
        node_dat_file = NodeDatFile.parse(node_file)
        if self.ptn_to_an is None:
            self.ptn_to_an = {}
        for entry in node_dat_file:
            family_name, an_id = entry.an_id.split(":", maxsplit=1)
            if self.name == family_name and an_id in self.graph:
                self.ptn_to_an[entry.ptn] = an_id
                self.graph.node[an_id]["ptn"] = entry.ptn
                self.graph.node[an_id]["node_type"] = entry.node_type
                # event_type?
                # branch_length?

    @staticmethod
    def parse(tree_file: str, tree_name: str = None, node_file: str = None):
        pthr_tree_graph = PantherTreeGraph(tree_name)

        # Parse Newick line
        phylo = PantherTreePhylo(tree_file)
        # Fill networkx graph from Phylo obj
        pthr_tree_graph.add_children(phylo.tree.clade)
        # Fill in long IDs on leaf nodes
        pthr_tree_graph.extract_leaf_ids(tree_file)
        # Fill in PTNs if node_file specified
        if node_file:
            pthr_tree_graph.extract_node_properties(node_file)


        return pthr_tree_graph

    def node(self, node):
        return self.graph.node.get(node)

    def root(self):
        for n in self.graph.nodes():
            # This criteria is asking to be abused
            if self.parents(n) == []:
                return n

    def ancestors(self, node, reflexive=False):
        nodes = list(networkx.ancestors(self.graph, node))
        if reflexive:
            nodes.append(node)
        return nodes

    def descendants(self, node, reflexive=False):
        nodes = list(networkx.descendants(self.graph, node))
        if reflexive:
            nodes.append(node)
        return nodes

    def parents(self, node):
        return list(self.graph.predecessors(node))

    def children(self, node):
        return list(self.graph.successors(node))

    def leaves(self, node=None):
        leaves = []
        if node is None:
            node = self.root()
        for n in self.descendants(node):
            if self.children(n):
                continue  # Has children so not a leaf
            leaves.append(n)
        return leaves

    def subgraph(self, nodes: List):
        return self.graph.subgraph(nodes).copy()

    def subtree(self, split_node):
        pthr_tree_graph = PantherTreeGraph(tree_name=self.name)
        pthr_tree_graph.graph = self.subgraph(
            self.descendants(split_node, reflexive=True)
        )
        return pthr_tree_graph

    def nodes_between(self, ancestor_node, descendant_node):
        descendants_of_anc = self.descendants(ancestor_node)
        ancestors_of_desc = self.ancestors(descendant_node)
        return list(set(descendants_of_anc) & set(ancestors_of_desc))

    def __len__(self):
        return len(self.graph)
