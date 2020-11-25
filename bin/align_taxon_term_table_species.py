#!/usr/bin/env python3

import csv
import argparse
from pthr_db_caller.taxon_validate import TaxonTermValidator, get_all_species_from_tree, append_species_to_table


parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxon_term_table", help="Gaferencer output taxon-term table")
parser.add_argument("-o", "--organism_dat")
parser.add_argument("-f", "--output_file")
parser.add_argument("-p", "--panther_tree_nhx")
parser.add_argument("-s", "--slim_terms", help="Optional - filters output table to only include supplied terms")
args = parser.parse_args()

# This script converts NCBITaxon ID column headers to PANTHER/PAINT organism labels from provided organism.dat
#  RefProt OSCODE (e.g. 'HUMAN', 'MOUSE', 'SACS2') for leaf organisms
#  Full taxon name (e.g. 'Eukaryota', 'Metazoa-Choanoflagellida') for ancestor organisms
# This script also fills in species missing a taxon ID and computes taxon constraints using provided panther_tree.nhx

if __name__ == "__main__":
    taxon_names = {}
    # organism.dat will hold all species in PAINT
    # If you don't like the name in output table, then change it in organism.dat
    with open(args.organism_dat) as org_f:
        reader = csv.reader(org_f, delimiter="\t")
        for r in reader:
            org_name = r[2]
            taxon_id = r[5]
            if taxon_id:
                taxon_names[f"NCBITaxon:{taxon_id}"] = org_name

    validator = TaxonTermValidator(args.taxon_term_table, args.panther_tree_nhx, args.slim_terms)
    validator.replace_taxon_header_labels(taxon_names)
    get_all_species_from_tree(validator)
    append_species_to_table(validator, args.output_file)
