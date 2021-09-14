#!/usr/bin/python3

import argparse
from pthr_db_caller.models.orthoxml import PthrOrthoXmlParser


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--xml_path', help="Path to directory containing OrthoXML files to merge")
parser.add_argument('-p', '--pthr_version', help="PANTHER version from where the input files originate")
parser.add_argument('-d', '--database_version', help="DB where gene IDs were minted")
parser.add_argument('-o', '--organism_dat', help="Oscode-to-taxonID lookup from PANTHER build process")


if __name__ == "__main__":
    args = parser.parse_args()

    all_groups = PthrOrthoXmlParser.parse(args.xml_path)

    pthr_version = args.pthr_version
    database_version = args.database_version
    organism_dat = None
    if args.organism_dat:
        organism_dat = args.organism_dat
    orthoxml_out = all_groups.to_orthoxml_str(pthr_version, database_version, organism_dat)
    print(orthoxml_out)
