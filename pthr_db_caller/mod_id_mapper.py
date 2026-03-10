import re
import sys
import csv
from typing import Dict, Optional, Tuple, List


class MODIDMapper:
    def __init__(self, gpi_uniprot_mappings: Dict[str, str] = None,
                 tair_mappings: Dict[str, str] = None,
                 araport_mappings: Dict[str, str] = None):
        """
        Initialize the ID mapper with required mapping dictionaries.

        Args:
            gpi_uniprot_mappings: Maps protein IDs to UniProt IDs
            tair_mappings: Maps TAIR gene IDs to locus IDs
            araport_mappings: Maps Araport protein IDs to gene IDs
        """
        self.gpi_uniprot_mappings = gpi_uniprot_mappings or {}
        self.tair_mappings = tair_mappings or {}
        self.araport_mappings = araport_mappings or {}
        self.id_lookup = {}

    @classmethod
    def from_files(MODIDMapper, gpi_uniprot_files: List[str] = None, tair_mapping_file: str = None, araport_mapping_file: str = None):
        gpi_uniprot_mappings = MODIDMapper.parse_gpi_uniprot_file(gpi_uniprot_files)
        tair_mappings = MODIDMapper.parse_tair_file(tair_mapping_file)
        araport_mappings = MODIDMapper.parse_araport_file(araport_mapping_file)
        return MODIDMapper(gpi_uniprot_mappings, tair_mappings, araport_mappings)

    @staticmethod
    def parse_gpi_uniprot_file(gpi_files: List[str]) -> Dict[str, str]:
        gpi_mappings = {}
        for gpi_file in gpi_files:
            with open(gpi_file, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                gpi_version = None
                for row in reader:
                    if not row:
                        continue
                    if row[0].startswith('!gpi-version:'):
                        gpi_version = row[0].split(":")[1].strip()
                    if len(row) < 8 or row[0].startswith('!'):
                        continue
                    if gpi_version and gpi_version == "2.0":
                        entity_id = row[0]
                        encoded_by = row[6]
                        xrefs = row[9]
                    else:
                        entity_id = row[0] + ":" + row[1]
                        encoded_by = row[7]
                        xrefs = row[8]
                    if encoded_by:
                        # Skip GPI line if Encoded_by is set since this is typically used by non-gene entities like isoforms or transcripts
                        continue
                    for xref in xrefs.split("|"):
                        db, ent_id = xref.split(":", maxsplit=1)
                        if db == "UniProtKB":
                            gpi_mappings[xref] = entity_id
        return gpi_mappings

    @staticmethod
    def parse_tair_file(tair_file: str) -> Dict[str, str]:
        with open(tair_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            tair_mappings = {}
            for row in reader:
                if len(row) < 2:
                    continue
                locus_id = row[0].strip()
                agi_id = row[1].strip()
                tair_mappings[agi_id] = locus_id
        return tair_mappings

    @staticmethod
    def parse_araport_file(araport_file: str) -> Dict[str, str]:
        with open(araport_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            araport_mappings = {}
            for row in reader:
                if len(row) < 2:
                    continue
                uniprot_id = row[0].strip()
                agi_id = row[1].strip()
                araport_mappings[uniprot_id] = agi_id
        return araport_mappings

    def parse_long_id(self, long_id: str) -> Tuple[str, str, str]:
        """Parse the long ID into organism, gene ID, and protein ID components."""
        parts = long_id.split('|')
        if len(parts) != 3:
            raise ValueError(f"Invalid long ID format: {long_id}")

        org, gene_id, protein_id = parts

        # Replace '=' with ':' in both IDs
        gene_id = gene_id.replace('=', ':')
        protein_id = protein_id.replace('=', ':')

        return org, gene_id, protein_id

    def handle_uniprot_mapping(self, protein_id: str) -> Optional[str]:
        """Check if protein ID has a UniProt mapping."""
        return self.gpi_uniprot_mappings.get(protein_id)

    def handle_gene_ensembl_case(self, gene_id: str, protein_id: str) -> Optional[str]:
        """Handle Gene or Ensembl prefixed gene IDs."""
        if re.match(r'^Gene|Ensembl', gene_id):
            return protein_id
        return None

    def handle_flybase(self, gene_id: str) -> str:
        """Handle FlyBase gene IDs by replacing prefix with FB."""
        return gene_id.replace('FlyBase', 'FB')

    def handle_wormbase(self, gene_id: str) -> str:
        """Handle WormBase gene IDs by replacing prefix with WB."""
        return gene_id.replace('WormBase', 'WB')

    def handle_tair(self, gene_id: str, protein_id: str) -> Optional[str]:
        """Handle TAIR gene IDs with complex locus mapping logic."""
        # Skip if it's already in the correct TAIR:locus:digit format
        if re.match(r'^TAIR:locus:\d+', gene_id):
            return None

        # Remove prefix (e.g., 'TAIR:')
        clean_gene_id = re.sub(r'^\w+:', '', gene_id)

        # Special case for 'locus' - use protein ID instead
        if clean_gene_id == 'locus':
            clean_protein_id = re.sub(r'^\w+:', '', protein_id)
            if clean_protein_id in self.araport_mappings:
                clean_gene_id = self.araport_mappings[clean_protein_id]

        # Look up the locus ID
        if clean_gene_id in self.tair_mappings:
            locus = self.tair_mappings[clean_gene_id]
            return f"TAIR:locus:{locus}"
        else:
            print(f"TAIR ID {clean_gene_id} has no mapped locus link ID.", file=sys.stderr)
            return None

    def handle_araport(self, gene_id: str) -> Optional[str]:
        """Handle Araport gene IDs."""
        clean_gene_id = re.sub(r'^\w+:', '', gene_id)

        if clean_gene_id in self.tair_mappings:
            locus = self.tair_mappings[clean_gene_id]
            return f"TAIR:locus:{locus}"
        else:
            print(f"Araport ID {clean_gene_id} has no mapped locus link ID.", file=sys.stderr)
            return None

    def handle_hgnc_ecogene(self, protein_id: str) -> str:
        """Handle HGNC and EcoGene cases - return protein ID."""
        return protein_id

    def get_short_id(self, long_id: str) -> Optional[str]:
        """
        Convert a long ID to its corresponding short ID based on the mapping logic.

        Args:
            long_id: The input long ID string

        Returns:
            The mapped short ID or None if mapping fails
        """
        try:
            org, gene_id, protein_id = self.parse_long_id(long_id)
        except ValueError as e:
            print(f"Error parsing ID: {e}", file=sys.stderr)
            return None

        short_id = None

        # Check UniProt mapping first
        short_id = self.handle_uniprot_mapping(protein_id)
        if short_id:
            self.id_lookup[long_id] = short_id
            return short_id

        # Check Gene/Ensembl case
        short_id = self.handle_gene_ensembl_case(gene_id, protein_id)
        if short_id:
            self.id_lookup[long_id] = short_id
            return short_id

        # Handle specific database cases
        if 'FlyBase' in gene_id:
            short_id = self.handle_flybase(gene_id)
        elif 'WormBase' in gene_id:
            short_id = self.handle_wormbase(gene_id)
        elif gene_id.startswith('TAIR') and not re.match(r'^TAIR:locus:\d+', gene_id):
            short_id = self.handle_tair(gene_id, protein_id)
            if short_id is None:  # Failed mapping, skip this ID
                return None
        elif 'Araport' in gene_id:
            short_id = self.handle_araport(gene_id)
            if short_id is None:  # Failed mapping, skip this ID
                return None
        elif 'HGNC' in gene_id:
            short_id = self.handle_hgnc_ecogene(protein_id)
        elif 'EcoGene' in gene_id:
            short_id = self.handle_hgnc_ecogene(protein_id)
        else:
            # Default case - use gene ID
            short_id = gene_id

        if short_id:
            self.id_lookup[long_id] = short_id

        return short_id


# Example usage and testing
def test_id_mapper():
    """Test the ID mapper with the provided examples."""

    # Mock data for testing (you'll need to provide real mappings)
    gpi_uniprot = {
        'UniProtKB:Q90Z46': 'ZFIN:ZDB-GENE-011026-1',
        'UniProtKB:A0A8M9QID4': 'ZFIN:ZDB-GENE-050410-16',
        'UniProtKB:F1QME4': 'ZFIN:ZDB-GENE-070202-8'
    }

    tair_mappings = {
        '2025595': '2025595'  # Mock mapping
    }

    mapper = MODIDMapper(gpi_uniprot, tair_mappings)

    test_cases = [
        ("HUMAN|HGNC=11179|UniProtKB=P00441", "UniProtKB:P00441"),
        ("DROME|FlyBase=FBgn0003462|UniProtKB=P61851", "FB:FBgn0003462"),
        ("ARATH|TAIR=locus=2025595|UniProtKB=P24704", "TAIR:locus:2025595"),
        ("MOUSE|MGI=MGI=98351|UniProtKB=P08228", "MGI:MGI:98351")
    ]

    print("Testing ID Mapper:")
    print("-" * 50)

    for long_id, expected in test_cases:
        result = mapper.get_short_id(long_id)
        status = "✓" if result == expected else "✗"
        print(f"{status} {long_id}")
        print(f"  Expected: {expected}")
        print(f"  Got:      {result}")
        print()


if __name__ == "__main__":
    test_id_mapper()