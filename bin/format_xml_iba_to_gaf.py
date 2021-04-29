import argparse
from pthr_db_caller.models.paint import PaintIbaXmlParser, gaf_line


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file_xml')


if __name__ == "__main__":
    args = parser.parse_args()

    anodes = PaintIbaXmlParser.parse(args.file_xml)
    for anode in anodes:
        for annot in anode.annotations:
            if annot.evidence_code == "IBA":
                print(gaf_line(annot, anode))
