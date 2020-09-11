from setuptools import setup, find_packages

setup(
    name="pthr_db_caller",
    version='0.0.11',
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
        "PyYAML==3.12"
    ]
)
