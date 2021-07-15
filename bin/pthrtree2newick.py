#!/usr/bin/env python3

import argparse
from pthr_db_caller.panther_tree_graph import PantherTreeGraph


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree_file')
parser.add_argument('-o', '--out_file')
parser.add_argument('-p', '--prune_species')


if __name__ == "__main__":
    args = parser.parse_args()
    tree = PantherTreeGraph.parse(tree_file=args.tree_file)
    if args.prune_species:
        taxon_list = []
        with open(args.prune_species) as spf:
            for l in spf.readlines():
                taxon_list.append(l.rstrip())
        tree.prune_species(taxon_list=taxon_list)
    if len(tree) > 0:
        tree.write(args.out_file)
    else:
        print("ERROR: Empty tree so nothing to write - {}".format(args.tree_file))
