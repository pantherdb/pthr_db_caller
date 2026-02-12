# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**pthr_db_caller** (v2.0.2) is a Python 3.10 library for querying PostgreSQL databases and handling results tailored to PantherDB-related uses. It provides models for PANTHER protein family data, phylogenetic trees, ortholog mapping (OrthoXML), and GO annotations (PAINT IBA/IBD).

## Build & Install

```bash
pip install -e .
```

Dependencies are defined in `requirements.txt` and auto-installed via `setup.py`. Key deps: `psycopg2` (Postgres), `biopython` (sequences/phylo), `ete3` (trees), `networkx` (graph algorithms), `lxml` (XML parsing), `PyYAML` (config).

## Testing

Uses Python's `unittest` framework. All tests are in `test.py` at the project root.

```bash
# Run all tests
python -m unittest test -v

# Run a single test class
python -m unittest test.TestPantherTreeGraph -v

# Run a single test method
python -m unittest test.TestPantherTreeGraph.test_family_tree -v
```

Test data lives in `resources/test/`. No linting or formatting tools are configured.

## Architecture

### Database Layer
- `pthr_db_caller/db_caller.py` — `DBCallerConfig` loads connection settings from YAML (`config/config.yaml` selects a named DB definition). `DBCaller` wraps psycopg2 for query execution with variable substitution and result formatting. Can be run as a CLI: `python -m pthr_db_caller.db_caller <query_file> [-v vars] [-o outfile] [-d delim] [-n]`.

### Models (`pthr_db_caller/models/`)
- **panther.py** — Core data types: `PthrSequence` (parses long IDs like `"SPECIES|Gene=X|UniProtKB=Y"`), `RefProtPantherMapping` (UniProt-to-PANTHER family mapping with lookup methods), `NodeDatFile` (PTN-to-AN node mapping), `GeneDatFile`, `OrganismDatFile`.
- **paint.py** — Largest module. Parses PAINT-curated family XML into `AnnotatedNode`/`Annotation` objects. `PaintIbaWriter` produces GAF 2.2 format output with GO aspect-dependent default qualifiers/relations. `PaintIbaXmlParser.parse()` is the XML entry point.
- **orthoxml.py** — `PthrOrthoXmlParser.parse()` accepts a file or directory of OrthoXML files. Models: `Gene`, `GeneCollection`, `OrthoXmlGroup`.
- **refprot_file.py** — Parsers for UniProt reference proteome files: `RefProtGeneAccFile`, `RefProtIdmappingFile`, `RefProtFastaFile`.
- **protein_class.py** — `ProteinClassGraph` wraps a NetworkX MultiDiGraph for protein classification hierarchies.
- **metadata.py** — `TaxonomyFile`/`TaxonomyRecord` for taxonomy lookups; `PaintIbaFile`/`PaintIbdFile` for output file metadata.

### Tree/Graph Operations
- `pthr_db_caller/panther_tree_graph.py` — `PantherTreeGraph` wraps NetworkX MultiDiGraph for phylogenetic trees parsed from Newick/NHX format. Key operations: `prune_species()`, `subtree()`, `ancestors()`, `descendants()`, `leaves()`. Factory via `PantherTreeGraph.parse()`.

### CLI Scripts (`bin/`)
Installed as console scripts via `setup.py`:
- `format_xml_iba_to_gaf.py` — Primary tool: converts PAINT XML to GAF/IBD annotation files, supports splitting by species.
- `pthrtree2newick.py` — Converts PANTHER tree format to Newick, supports species pruning.
- `etree2orthoxml.py` — ETE3-based tree to OrthoXML conversion.
- `merge_orthoxml.py` — Merges multiple OrthoXML files, filters singletons.
- `align_taxon_term_table_species.py` — Aligns species in taxon-term constraint tables.

### Common Patterns
- **Parser factory methods**: Most model classes use `ClassName.parse(filename)` as the entry point.
- **Lookup collections**: Models maintain both ordered entry lists and fast-lookup dicts (e.g., `find_uniprot()`, `find_long_id()`).
- **Config-driven DB**: DB connections selected by named definition in `config/config.yaml`; extra YAML properties become query variables.
