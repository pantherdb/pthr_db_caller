#!/usr/bin/python3

import argparse
import os
from pthr_db_caller.models import paint, metadata


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file_xml')
parser.add_argument('-g', '--go_aspect', help="Filepath to TSV of GO term ID -> GO aspect. "
                                              "E.g. 'GO:0009507\tcellular_component'")
parser.add_argument('-c', '--complex_termlist', help="Filepath to (pre-computed from GO ontology) list of "
                                                     "protein-containing complex (GO:0032991) and all its descendant "
                                                     "terms.")
parser.add_argument('-s', '--split_by_species', help="Filepath to 'filename,taxon_id,oscode' TSV. Write to STDOUT"
                                                     "if this option is omitted.")
parser.add_argument('-d', '--out_directory', help="Destination directory for split files. Only used if"
                                                  "--split_by_species is specified.")
parser.add_argument('-a', '--file_format', help="GO annotation format to output. Default is 'GAF' (version 2.2)")
parser.add_argument('-p', '--panther_version', help="PANTHER library version. E.g. 15.0, 16.0")
parser.add_argument('-r', '--go_release_date', help="GO release date in YYYY-MM-DD format")
parser.add_argument('-u', '--obsolete_uniprots', help="Filepath to list of PANTHER UniProt IDs not in UniProt GPI")
parser.add_argument('-b', '--ibd_file_outpath', help="If supplied, filepath to write IBD file to")


if __name__ == "__main__":
    args = parser.parse_args()

    if args.file_format:
        file_format = args.file_format.upper()
    else:
        file_format = "GAF"
    writer = paint.PaintIbaWriter(go_aspect=args.go_aspect,
                                  complex_termlist=args.complex_termlist,
                                  file_format=file_format,
                                  obsolete_uniprots=args.obsolete_uniprots)

    anodes = paint.PaintIbaXmlParser.parse(args.file_xml)
    if args.ibd_file_outpath:
        ibd_nodes = anodes.ibd_nodes()
        ibd_file = metadata.PaintIbdFile(writer=writer,
                                         panther_version=args.panther_version,
                                         go_release_date=args.go_release_date,
                                         ibd_nodes=ibd_nodes)
        ibd_file.write(args.ibd_file_outpath)

    # Split anodes by file to write to; by taxon
    if args.split_by_species:
        iba_file_data = metadata.parse_iba_metadata_file(args.split_by_species)
        # Add the catch-all, fallback file
        iba_file_data.append({"basename": "gene_association.paint_other", "taxon_id": "other", "oscode": None})
        taxon_to_file = {}
        for pif in iba_file_data:
            iba_file = metadata.PaintIbaFile(writer=writer,
                                             panther_version=args.panther_version,
                                             go_release_date=args.go_release_date,
                                             basename=pif["basename"],
                                             taxon_id=pif["taxon_id"],
                                             oscode=pif["oscode"])
            taxon_to_file[iba_file.taxon_id] = iba_file
        for node in anodes:
            iba_file = taxon_to_file.get(node.taxon_id)
            if iba_file is None:
                iba_file = taxon_to_file.get("other")
            iba_file.add_node(node)
        # Now that the iba_files have their annotated_nodes
        for taxon_id, iba_file in taxon_to_file.items():
            # Specify format (gaf) and outdir and
            if iba_file.annotated_nodes:
                full_filepath = os.path.join(args.out_directory, iba_file.basename)
                full_filepath = "{}.{}".format(full_filepath, file_format.lower())
                print(iba_file.basename, len(iba_file.annotated_nodes))
                iba_file.write(full_filepath)
    else:
        writer.print(anodes)
