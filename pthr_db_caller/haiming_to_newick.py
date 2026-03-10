#!/usr/bin/env python3
"""
ToNewick Python Module
Converted (by Claude) from Perl HAIMING::ToNewick module
Handles phylogenetic tree operations and Newick format conversions
"""

import re
import os
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from io import StringIO
import tempfile


def format_branch_length(length):
    """Format branch length to remove unnecessary precision and zeros"""
    if length is None or length == 0:
        return None

    # Round to 3 decimal places and remove trailing zeros
    formatted = f"{length:.3f}".rstrip('0').rstrip('.')

    # If it becomes empty or just a decimal point, it was effectively zero
    if not formatted or formatted == '.':
        return None

    return formatted


def tree_to_newick(tree, format_lengths=True):
    """Convert tree to Newick format with proper branch length formatting"""

    def clade_to_string(clade):
        result = ""

        if clade.clades:
            # Internal node - add children in parentheses
            child_strings = []
            for child in clade.clades:
                child_strings.append(clade_to_string(child))
            result = "(" + ",".join(child_strings) + ")"

        # Add node name if it exists
        if clade.name:
            result += clade.name

        # Add branch length if it exists and is not effectively zero
        if format_lengths and clade.branch_length is not None:
            formatted_length = format_branch_length(clade.branch_length)
            if formatted_length:
                result += ":" + formatted_length

        return result

    return clade_to_string(tree.root) + ";"


class ToNewick:
    def __init__(self):
        self.bi_count = 0
        self.newsets = None

    def buildtree(self, relationfile, outtree):
        """
        Given a file of relationships, return a tree file
        Report error if some node has 2 moms
        """
        relation = {}
        nodes = {}
        rootlabel = None

        try:
            with open(relationfile, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')

                    if parts[0] in relation and parts[1] != relation[parts[0]]:
                        raise ValueError(
                            f"A child {parts[0]} should not have more than 1 mom: {parts[1]}! {relation[parts[0]]}")

                    relation[parts[0]] = parts[1]  # from child to mom

                    for label in parts:
                        if label not in nodes:
                            nodes[label] = Clade(name=label)
        except IOError:
            raise IOError(f"Cannot open file: {relationfile}")

        # Build tree structure
        for child in relation:
            mom = relation[child]
            childnode = nodes[child]
            momnode = nodes[mom]
            momnode.clades.append(childnode)

        # Find root
        n = 0
        for mom in relation.values():
            if mom not in relation:
                n += 1
                if n > 1 and rootlabel != mom:
                    raise ValueError(f"{rootlabel} and {mom} both have no mom?")
                rootlabel = mom

        rootnode = nodes[rootlabel]
        tree = Tree(root=rootnode)

        # Write tree with custom formatting
        newick_string = tree_to_newick(tree)
        try:
            with open(outtree, 'w') as f:
                f.write(newick_string)
        except IOError:
            raise IOError(f"Cannot output to {outtree}")

    def speciestree(self, infile, outfile):
        """Convert tree with species annotations to species tree"""
        try:
            with open(infile, 'r') as f:
                treestring = f.read().strip()
        except IOError:
            raise IOError(f"Cannot open {infile}")

        # Replace species annotations
        treestring = re.sub(r'\[\&\&N[HX]+:S=([\w_-]+)\]', r'\1', treestring)

        try:
            with open(outfile, 'w') as f:
                f.write(f"({treestring});")
        except IOError:
            raise IOError(f"Cannot output to {outfile}")

    def tobifurcate(self, infile, outfile):
        """Convert multifurcating tree to bifurcating tree"""
        self.bi_count = 0

        try:
            with open(infile, 'r') as f:
                line = f.readline().strip()  # Only read the first line (tree string)
        except IOError:
            raise IOError(f"Cannot open {infile}")

        # Clean tree string for BioPython usage
        # Original Perl: $line =~ s/\)(:[0-9\.]+)\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]/\)$2$1/g;
        line = re.sub(r'\)(:[0-9\.]+)\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]', r')\2\1', line)
        # Original Perl: $line =~ s/\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]/$1:0.15/g;
        line = re.sub(r'\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]', r'\1:0.15', line)

        # Create temporary file for BioPython
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.nwk') as tmp:
            tmp.write(line)
            tmp_name = tmp.name

        try:
            # Read tree
            tree = Phylo.read(tmp_name, 'newick')

            # Update tree to bifurcate
            self._updatetree(tree.root)

            # Write result with custom formatting
            newick_string = tree_to_newick(tree)
            try:
                with open(outfile, 'w') as f:
                    f.write(newick_string)
            except IOError:
                raise IOError(f"Cannot output to {outfile}")
        finally:
            os.unlink(tmp_name)

    def removeleaf(self, intree, removeid, outtree):
        """Remove a leaf from the input tree and return resulting tree"""
        try:
            tree = Phylo.read(intree, 'newick')
        except IOError:
            raise IOError(f"Cannot open {intree}")

        # Find node to remove
        node_to_remove = None
        for clade in tree.find_clades():
            if clade.name == removeid:
                node_to_remove = clade
                break

        if not node_to_remove:
            print(f"Cannot find node for ID {removeid}, recheck!!")
            return

        # Find parent
        parent = None
        for clade in tree.find_clades():
            if node_to_remove in clade.clades:
                parent = clade
                break

        if not parent:
            print(f"Cannot find parent for node {removeid}")
            return

        children = list(parent.clades)
        size = len(children)

        if size > 2:
            parent.clades.remove(node_to_remove)
        else:
            # Find grandparent
            grandparent = None
            for clade in tree.find_clades():
                if parent in clade.clades:
                    grandparent = clade
                    break

            if grandparent:
                for child in children:
                    if child.name != removeid:
                        grandparent.clades.append(child)
                grandparent.clades.remove(parent)
            else:  # parent is root
                for child in children:
                    if child.name != removeid:
                        tree.root = child
                        break

        # Write with custom formatting
        newick_string = tree_to_newick(tree)
        try:
            with open(outtree, 'w') as f:
                f.write(newick_string)
        except IOError:
            raise IOError(f"Cannot output to {outtree}")

    def _updatetree(self, node):
        """Recursively update tree to ensure bifurcation"""
        self.newsets = []

        if node.is_terminal():
            return

        children = list(node.clades)
        size = len(children)

        if size < 3:
            for child in children:
                self._updatetree(child)
        else:
            node.clades.clear()
            self._group(children)

            for child in self.newsets:
                node.clades.append(child)
                self._updatetree(child)

    def _group(self, nodes):
        """Group nodes into binary pairs"""
        size = len(nodes)
        di = size // 2
        newnodes = []

        for i in range(1, di + 1):
            new_node = Clade(name=f"BI{self.bi_count}", branch_length=0.15)
            self.bi_count += 1

            for j in range(2):
                index = (i - 1) * 2 + j
                new_node.clades.append(nodes[index])

            newnodes.append(new_node)

        if size % 2:
            newnodes.append(nodes[size - 1])

        newsize = len(newnodes)
        if newsize == 2:
            self.newsets = newnodes
            return
        else:
            self._group(newnodes)

    def tonewick_only_an(self, infile, outfile):
        """Convert tree to Newick format with only AN identifiers"""
        try:
            with open(infile, 'r') as f:
                treestring = f.readline().strip()  # Only read the first line (tree string)
        except IOError:
            raise IOError(f"Cannot open {infile}")

        # Process tree string
        # Original Perl: $treestring =~ s/\)(:[0-9\.]+)\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]/\)$2$1/g;
        treestring = re.sub(r'\)(:[0-9\.]+)\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]', r')\2\1', treestring)
        # Original Perl: $treestring =~ s/\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]/$1\)/g;
        treestring = re.sub(r'\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]', r'\1)', treestring)

        treestring = "(" + treestring

        try:
            with open(outfile, 'w') as f:
                f.write(treestring)
        except IOError:
            raise IOError(f"Cannot output to {outfile}")

    def tonewick(self, infile, outfile):
        """Convert tree to Newick format with species information"""
        try:
            with open(infile, 'r') as f:
                content = f.read()
        except IOError:
            raise IOError(f"Cannot open {infile}")

        lines = content.strip().split('\n')
        treestring = lines[0] if lines else ""

        an_species = {}

        # Extract species information from tree string
        for match in re.finditer(r'(\[[&:=>\_\/\-A-Za-z0-9]+:ID=(AN[0-9]+)\])', treestring):
            string = match.group(1)
            node = match.group(2)

            species_match = re.search(r':S=(.*):ID=', string)
            if species_match:
                species = species_match.group(1)
                an_species[node] = f"{node}|{species}|"
            elif 'Ev=1>0' in string:
                an_species[node] = f"{node}DUPLICATION"

        # Clean tree string
        # Original Perl: $treestring =~ s/\)(:[0-9\.]+)\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]/\)$2$1/g;
        treestring = re.sub(r'\)(:[0-9\.]+)\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]', r')\2\1', treestring)
        # Original Perl: $treestring =~ s/\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]/$1\)/g;
        treestring = re.sub(r'\[[&:=>\/_\-A-Za-z0-9]+:ID=(AN[0-9]+)\]', r'\1)', treestring)

        # Process additional lines for species mapping
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
            match = re.match(r'^(AN[0-9]+):([A-Za-z0-9\-_]+)\|', line)
            if match:
                an = match.group(1)
                species = match.group(2)
                an_species[an] = f"{species}_{an}"

        # Fill in missing AN identifiers
        for match in re.finditer(r'(AN[0-9]+)', treestring):
            an_id = match.group(1)
            if an_id not in an_species:
                an_species[an_id] = an_id

        # Replace AN identifiers with species names
        for an_id, species_name in an_species.items():
            treestring = re.sub(f'{re.escape(an_id)}:', f'{species_name}:', treestring)
            treestring = re.sub(f'{re.escape(an_id)}\\)', f'{species_name})', treestring)

        try:
            with open(outfile, 'w') as f:
                f.write("(" + treestring)
        except IOError:
            raise IOError(f"Cannot output to {outfile}")


# Convenience functions for backward compatibility
def buildtree(relationfile, outtree):
    tn = ToNewick()
    return tn.buildtree(relationfile, outtree)


def speciestree(infile, outfile):
    tn = ToNewick()
    return tn.speciestree(infile, outfile)


def tobifurcate(infile, outfile):
    tn = ToNewick()
    return tn.tobifurcate(infile, outfile)


def removeleaf(intree, removeid, outtree):
    tn = ToNewick()
    return tn.removeleaf(intree, removeid, outtree)


def tonewick(infile, outfile):
    tn = ToNewick()
    return tn.tonewick(infile, outfile)


def tonewickOnlyAN(infile, outfile):
    tn = ToNewick()
    return tn.tonewick_only_an(infile, outfile)