import networkx
from networkx import MultiDiGraph
from Bio import Phylo
from io import StringIO


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


class PantherTreeGraph(MultiDiGraph):

    def __init__(self, tree_file):
        MultiDiGraph.__init__(self)

        self.phylo = PantherTreePhylo(tree_file)

        # Fill graph from Phylo obj
        self.add_children(self.phylo.tree.clade)

    # Recursive method to fill graph from Phylo clade, parsing out node accession and species name (if present)
    def add_children(self, parent_clade):
        self.add_node_from_clade(parent_clade)
        for child in parent_clade.clades:
            self.add_node_from_clade(child)
            self.add_edge(parent_clade.name, child.name)

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
        if clade.name not in self.nodes():
            self.add_node(clade.name)
        if species:
            self.nodes[clade.name]["species"] = species

    def ancestors(self, node):
        return list(networkx.ancestors(self, node))

    def descendants(self, node):
        return list(networkx.descendants(self, node))

    def parents(self, node):
        return list(self.predecessors(node))

    def children(self, node):
        return list(self.successors(node))

    def nodes_between(self, ancestor_node, descendant_node):
        descendants_of_anc = self.descendants(ancestor_node)
        ancestors_of_desc = self.ancestors(descendant_node)
        return list(set(descendants_of_anc) & set(ancestors_of_desc))