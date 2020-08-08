import networkx
import csv
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger(__name__)


class ProteinClassGraph(networkx.MultiDiGraph):
    def __init__(self):
        networkx.MultiDiGraph.__init__(self)

    @staticmethod
    def parse_class_and_rel_files(class_file, rel_file):
        pc_graph = ProteinClassGraph()

        # PC-to-PC
        pc_relationships = ProteinClassRelationship.parse_relationship_file(rel_file)
        for r in pc_relationships:
            pc_graph.add_edge(r.parent.id, r.child.id) # parent, child

        # PC, name
        in_pc = open(class_file)
        pc_reader = csv.reader(in_pc, delimiter=";")
        pcs = {}
        next(pc_reader)  # skip header
        for pc in pc_reader:
            pcs[pc[0]] = pc[1]
            if pc[0] not in pc_graph:
                pc_graph.add_node(pc[0])
        in_pc.close()

        return pc_graph

    @staticmethod
    def show_graph(graph):
        pos = networkx.spring_layout(graph)
        networkx.drawing.nx_pylab.draw(graph, pos=pos)
        networkx.draw_networkx_labels(graph, pos)
        plt.show()


class ProteinClass:
    def __init__(self, id, description=None, is_valid=True):
        self.id = id
        self.description = description
        # is_valid currently denotes whether assignment was commented out (is_valid==False) or not
        self.is_valid = is_valid

    def __str__(self):
        if self.description:
            return "{} - {}".format(self.id, self.description)
        else:
            return self.id

    def __eq__(self, other):
        return self.id == other.id and self.description == other.description and self.is_valid == other.is_valid


class ProteinClassRelationship:
    def __init__(self, parent, child):
        self.parent = parent
        self.child = child

    @staticmethod
    def parse_relationship_file(filename):
        pc_relationships = []
        with open(filename) as rf:
            reader = csv.reader(rf, delimiter="\t")
            for r in reader:
                if len(r) > 0 and r[0].startswith("PC"):
                    parent_class = ProteinClass(r[2], description=r[3])
                    child_class = ProteinClass(r[0], description=r[1])
                    pc_relationships.append(ProteinClassRelationship(parent=parent_class,
                                                                    child=child_class))
        return pc_relationships


class PthrToPc:
    PC_NAMES = {}

    def __init__(self, family, status, fam_name, pc_annots, source_line=None):
        self.family = family
        self.status = None
        self.fam_name = fam_name
        self.new_name = None
        self.pc_annots = self.handle_pc_annots(pc_annots)
        self.source_line = source_line
        self.handle_status(status)

    def handle_status(self, status):
        if len(status) > 0:
            if status in ["REMOVE", "DIVIDE"]:
                self.status = status
            else:
                self.new_name = status

    def get_pc_annots(self, is_valid=None):
        annots_list = []
        for a in self.pc_annots:
            if is_valid:
                # Check validity
                if is_valid == a.is_valid:
                    annots_list.append(a)
            else:
                # Return valid and invalid
                annots_list.append(a)
        return annots_list

    def get_valid_annots(self):
        return self.get_pc_annots(is_valid=True)

    def get_invalid_annots(self):
        return self.get_pc_annots(is_valid=False)

    def handle_pc_annots(self, pc_annots):
        pc_annots_list = []
        for pca in pc_annots:
            if len(pca) > 0:
                pc_bits = pca.split(" - ")  # TODO: Allow for extra whitespace
                pc_id = pc_bits[0].rstrip()
                # pc_description = None
                is_valid = True
                if len(pc_bits) > 1:
                    pc_description = pc_bits[1]
                    PthrToPc.PC_NAMES[pc_id] = pc_description
                else:
                    pc_description = PthrToPc.PC_NAMES.get(pc_id)
                if pc_id.startswith("#") or (pc_description and pc_description.startswith("#")):
                    # pc_annots_ds["invalid_annots"].append(pc_id.replace("#", ""))
                    pc_id = pc_id.replace("#", "")
                    is_valid=False
                    logging.info("{} - {} is commented out".format(self.family, pc_id))
                # else:
                    # pc_annots_ds["valid_annots"].append(pc_id)
                pc_class = ProteinClass(pc_id, pc_description, is_valid=is_valid)
                if pc_class not in pc_annots_list:
                    pc_annots_list.append(pc_class)
        return pc_annots_list

    def __str__(self):
        col_2 = self.status
        if col_2 is None:
            col_2 = self.new_name
        if col_2 is None:
            col_2 = ""
        return "\t".join([self.family, col_2, self.fam_name] + [str(a) for a in self.get_valid_annots()])

    @staticmethod
    def parse_pthr_to_pc(filename):
        with open(filename) as f:
            pthr_to_pcs = []
            reader = csv.reader(f, delimiter="\t")
            for r in reader:
                pthr_to_pcs.append(PthrToPc(family=r[0],
                                    status=r[1],
                                    fam_name=r[2],
                                    pc_annots=r[3:],
                                    source_line="\t".join(r))
                )
            return pthr_to_pcs

    @staticmethod
    def find_by_fam(pthr_to_pcs, family):
        for r in pthr_to_pcs:
            if r.family == family:
                return r
