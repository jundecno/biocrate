import periodictable
import re

MAP_SIZE = 1073741824
elem2num = lambda x: periodictable.elements.symbol(x[:1].upper() + x[1:].lower()).number
elem2mass = lambda x: periodictable.elements.symbol(x[:1].upper() + x[1:].lower()).mass
main_enum = {"N": 0, "CA": 1, "C": 2, "O": 3}
PDB_MIRROR = "https://files.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb"
CIF_MIRROR = "https://files.wwpdb.org/pub/pdb/data/structures/all/mmCIF/"
PDB_BUNDLE_MIRROR = "ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/"
UNIPROT_MIRROR = "https://www.uniprot.org/uniprot/"
PDBE_MIRROR = "https://www.ebi.ac.uk/pdbe-srv/view/entry/"
RCSB_FASTA_MIRROR = "https://www.rcsb.org/fasta/entry/"
RCSB_PDB_MIRROR = "https://www.rcsb.org/pdb/files/"
SCOP_MIRROR = "https://scop.berkeley.edu/downloads/pdbstyle"

