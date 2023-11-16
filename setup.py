from setuptools import setup, find_packages

with open("requirements.txt", "r") as FH:
    REQUIREMENTS = FH.readlines()

setup(
    name="pthr_db_caller",
    version='2.0.1',
    packages=find_packages(),
    author="dustine32",
    author_email="debert@usc.edu",
    description="Python library for querying postgresl DBs and handling results tailored to PantherDB-related uses",
    long_description=open("README.md").read(),
    url="https://github.com/pantherdb/pthr_db_caller",
    install_requires=[r for r in REQUIREMENTS if not r.startswith("#")],
    scripts=[
        "bin/align_taxon_term_table_species.py",
        "bin/etree2orthoxml.py",
        "bin/pthrtree2newick.py",
        "bin/taxon_term_tbl_lkp.py",
        "bin/format_xml_iba_to_gaf.py",
        "bin/merge_orthoxml.py"
    ]
)
