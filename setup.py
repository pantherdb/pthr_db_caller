from setuptools import setup, find_packages

setup(
    name="pthr_db_caller",
    version='0.0.19',
    packages=find_packages(),
    author="dustine32",
    author_email="debert@usc.edu",
    description="Python library for querying postgresl DBs and handling results tailored to PantherDB-related uses",
    long_description=open("README.md").read(),
    url="https://github.com/dustine32/pthr_db_caller",
    install_requires=[
        "psycopg2==2.7.4",
        "biopython==1.73",
        "networkx==2.3",
        "matplotlib==3.1.1",
        "PyYAML==3.12",
        "ete3>=3.1.2",
        "lxml>=4.6.3",
        "dataclasses"
    ],
    scripts=[
        "bin/align_taxon_term_table_species.py",
        "bin/etree2orthoxml.py",
        "bin/pthrtree2newick.py",
        "bin/taxon_term_tbl_lkp.py",
        "bin/format_xml_iba_to_gaf.py",
        "bin/merge_orthoxml.py"
    ]
)
