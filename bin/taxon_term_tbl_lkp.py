#!/usr/bin/env python3

import argparse
from pthr_db_caller.taxon_validate import TaxonTermValidator


parser = argparse.ArgumentParser()
parser.add_argument("-t", "--taxon_term_table")
parser.add_argument('-n', '--taxon', type=str)
parser.add_argument("-g", "--term", type=str)


if __name__ == "__main__":
    args = parser.parse_args()

    validator = TaxonTermValidator(args.taxon_term_table)
    result = validator.taxon_term_lookup(args.taxon, args.term)
    print(result)
