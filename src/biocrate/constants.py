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

pdb_pattern = re.compile(r"^\d\w{3}$")
pdb_chain_pattern = re.compile(r"^\d\w{4,7}$")
accession_pattern = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")
entry_pattern = re.compile(r"^\w{1,5}_\w{1,5}")  # uniprot entry name
PDBe_entry_pattern = re.compile(r"DR\s+PDB;\s*(\d\w{3});[\w\W]+?;[\w\W]+?;\s*([\w\W]+?)\.")
DBREF_pattern = re.compile(r"DBREF[\d]{0,1}\s{1,2}\d\w{3}\s\w{0,1}\s+(\d+)\s+(\d+)\s+UNP\s+(\w+)")
PFAM_pattern = re.compile(r"[Pp][Ff]\d{5}")
pdb_bundle_pattern = re.compile(r"\d\w{3}\-pdb\-bundle\d+\.pdb")
scop_pattern = re.compile(r"^d\d[0-9a-z]{3}[_0-9a-zA-Z]{2,}$")
