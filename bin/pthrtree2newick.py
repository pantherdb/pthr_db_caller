#!/usr/bin/env python3

import argparse
from pthr_db_caller.panther_tree_graph import PantherTreeGraph


parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree_file')
parser.add_argument('-o', '--out_file')


if __name__ == "__main__":
    args = parser.parse_args()
    tree = PantherTreeGraph.parse(tree_file=args.tree_file)
    tree.write(args.out_file)
