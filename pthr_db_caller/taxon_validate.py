from Bio import Phylo
import csv
import argparse
import logging

logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output_file")
parser.add_argument("-s", "--slim_terms")
parser.add_argument("-t", "--taxon_term_table")
parser.add_argument("-p", "--panther_tree_nhx")

subparsers = parser.add_subparsers(dest='subcommand', help='sub-command help')
parser_n = subparsers.add_parser('lkp', help='Direct lookup of input taxon-term table')
parser_n.add_argument('-n', '--taxon', type=str)
parser_n.add_argument("-g", "--term", type=str)
parser_n.set_defaults(function="taxon_term_lookup")

THE_REST = []


def extract_clade_name(clade_comment):
    if clade_comment is None:
        clade_comment = ""
    ### Name-parsing should be more robust
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
    return new_comment, an_id


def name_children(parent_clade):
    # print(parent_clade.name)
    for child in parent_clade.clades:
        child.name, child.id = extract_clade_name(child.comment)
        if len(child.clades) > 0:
            name_children(child)


def find_taxon_clade(taxon_name, root_clade):
    # print(root_clade)
    if root_clade.name == taxon_name:
        return root_clade
    elif len(root_clade.clades) > 0:
        for c in root_clade.clades:
            result = find_taxon_clade(taxon_name, c)
            if result is not None:
                return result


class TaxonTermValidator:
    def __init__(self, taxon_term_table, panther_tree_nhx=None, slim_terms=None):
        self.term_constraint_lists = {}
        self.taxon_indexes = {}
        self.slim_terms = []
        self.tree = None

        if slim_terms:
            # Get list of slim terms to filter for
            slim_file = open(slim_terms)
            for t in slim_file.readlines():
                self.slim_terms.append(t.rstrip())
            slim_file.close()
        
        with open(taxon_term_table) as t3f:
            header = t3f.readline().rstrip()
            headers = header.split("\t")
            index_count = 0
            for h in headers[1:len(headers)]:
                self.taxon_indexes[h] = index_count
                index_count += 1

            for l in t3f.readlines():
                cols = l.split("\t")
                go_term = cols[0]
                if len(self.slim_terms) == 0 or go_term in self.slim_terms:
                    self.term_constraint_lists[go_term] = cols[1:len(cols)]

        logger.debug("taxon_indexes: {}".format(len(self.taxon_indexes)))
        logger.debug("term_constraint_lists: {}".format(len(self.term_constraint_lists)))

        # Parse species_tree
        if panther_tree_nhx:
            self.tree = next(Phylo.parse(panther_tree_nhx, "newick"))
            self.tree.clade.name, self.tree.clade.id = extract_clade_name(self.tree.clade.comment)
            name_children(self.tree.clade)

    # This is used to replace NCBITaxon:##### column labels with species codes used in PANTHER and PAINT
    def replace_taxon_header_labels(self, new_label_lookup: dict):
        new_taxon_indexes = {}
        for t in self.taxon_indexes:
            taxon_label = t
            if t in new_label_lookup:
                taxon_label = new_label_lookup[t]
            new_taxon_indexes[taxon_label] = self.taxon_indexes[t]
        self.taxon_indexes = new_taxon_indexes

    def validate_taxon_term(self, taxon, term):
        # node_path[-2] won't work for LUCA. LUCA should equal NCBITaxon:131567 for "cellular organisms".
        # Need to rerun gaferencer to include this taxon, then convert "cellular organisms" header to "LUCA" in
        # taxon_to_oscode.py
        # print(taxon)
        while taxon in THE_REST and taxon != "LUCA":
            # print(taxon)
            # Get parent of taxon - handy BioPython trick
            taxon_clade = find_taxon_clade(taxon, self.tree.clade)
            node_path = self.tree.get_path(taxon_clade)
            if len(node_path) > 1:
                parent_clade = node_path[-2]
            else:
                parent_clade = self.tree.clade
            taxon = parent_clade.name

        # Remove after getting LUCA's real values - assuming most everything is cool with LUCA
        if taxon == "LUCA":
            # virus stuff
            if term in ("GO:0019012", "GO:0039679", "GO:0044423"):
                return False
            else:
                return True

        try:
            result = self.taxon_term_lookup(taxon, term)
        except IndexError:
            print(taxon)
            result = self.taxon_term_lookup(taxon, term)
        if result == '0':
            return False
        return True

    def taxon_term_lookup(self, taxon, term):
        return self.term_constraint_lists[term][self.taxon_indexes[taxon]]


def append_species_to_table(validator : TaxonTermValidator, output_file):
    # List hyphenated species
    # Reconstruct entire table file
    # with open("new_table_file", "w+") as nf:
    with open(output_file, "w+") as nf:
        writer = csv.writer(nf, delimiter="\t")
        header = ["GOterm"]
        out_term_values = {}
        for tk in list(validator.taxon_indexes.keys()) + THE_REST:
            header.append(tk)
            for term in validator.term_constraint_lists:
                if validator.validate_taxon_term(tk, term):
                    result = 1
                else:
                    result = 0
                if term in out_term_values:
                    out_term_values[term].append(result)
                else:
                    out_term_values[term] = [term, result]

        writer.writerow(header)
        for otv in out_term_values:
            if len(validator.slim_terms) == 0 or otv in validator.slim_terms:
                writer.writerow(out_term_values[otv])


def get_all_species_from_tree(validator : TaxonTermValidator):
    for c in validator.tree.find_clades():
        if len(c.name) > 0 and c.name not in list(validator.taxon_indexes.keys()) + THE_REST:
            THE_REST.append(c.name)


if __name__ == "__main__":
    args = parser.parse_args()

    func = args.function
    if func and func == "taxon_term_lookup":
        validator = TaxonTermValidator(args.taxon_term_table)
        result = validator.taxon_term_lookup(args.taxon, args.term)
        print(result)
    else:
        validator = TaxonTermValidator(args.taxon_term_table, args.panther_tree_nhx, args.slim_terms)
        get_all_species_from_tree(validator)
        append_species_to_table(validator, args.output_file)
