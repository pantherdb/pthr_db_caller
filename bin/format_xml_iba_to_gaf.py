import argparse
from pthr_db_caller.models.paint import PaintIbaXmlParser, PaintIbaWriter


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file_xml')
parser.add_argument('-g', '--go_aspect', help="Filepath to TSV of GO term ID -> GO aspect. "
                                              "E.g. 'GO:0009507\tcellular_component'")
parser.add_argument('-c', '--complex_termlist', help="Filepath to (pre-computed from GO ontology) list of "
                                                     "protein-containing complex (GO:0032991) and all its descendant "
                                                     "terms.")


if __name__ == "__main__":
    args = parser.parse_args()

    anodes = PaintIbaXmlParser.parse(args.file_xml)
    writer = PaintIbaWriter(go_aspect=args.go_aspect, complex_termlist=args.complex_termlist)
    writer.write(anodes)
