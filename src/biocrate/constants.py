import re

import periodictable

MAP_SIZE = 1073741824
elem2num = lambda x: periodictable.elements.symbol(x[:1].upper() + x[1:].lower()).number
elem2mass = lambda x: periodictable.elements.symbol(x[:1].upper() + x[1:].lower()).mass
MAIN_ATOMS = frozenset({"N": 0, "CA": 1, "C": 2, "O": 3})
H_ATOMS = frozenset({"H", "D"})
PDB_MIRROR = "https://files.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb"
CIF_MIRROR = "https://files.wwpdb.org/pub/pdb/data/structures/all/mmCIF/"
PDB_BUNDLE_MIRROR = "ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/"
UNIPROT_MIRROR = "https://www.uniprot.org/uniprot/"
PDBE_MIRROR = "https://www.ebi.ac.uk/pdbe-srv/view/entry/"
RCSB_FASTA_MIRROR = "https://www.rcsb.org/fasta/entry/"
RCSB_PDB_MIRROR = "https://www.rcsb.org/pdb/files/"
RCSB_LIGAND_MIRROR = "https://data.rcsb.org/rest/v1/core/chemcomp/"
SCOP_MIRROR = "https://scop.berkeley.edu/downloads/pdbstyle"
UNIPROT_ACCESSION_PATTERN = re.compile(
    r"^(?:" r"[OPQ][0-9][A-Z0-9]{3}[0-9]" r"|" r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]" r"|" r"[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z0-9]{3}[0-9]" r")$"
)
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
# fmt: off
SEC_STRUC = ["H", "B", "E", "G", "I", "T", "S", "-"]
SEC_STRUC_DICT = {s: i for i, s in enumerate(SEC_STRUC)}
STAND_RESI_DICT = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

SIDECHAIN_HBD_DICT = {
    "A": 0, "R": 3, "N": 1, "D": 0, "C": 1,
    "Q": 1, "E": 0, "G": 0, "H": 1, "I": 0,
    "L": 0, "K": 1, "M": 0, "F": 0, "P": 0,
    "S": 1, "T": 1, "W": 1, "Y": 1, "V": 0,
}

SIDECHAIN_HBA_DICT = {
    "A": 0, "R": 0, "N": 1, "D": 2, "C": 1,
    "Q": 1, "E": 2, "G": 0, "H": 1, "I": 0,
    "L": 0, "K": 0, "M": 1, "F": 0, "P": 0,
    "S": 1, "T": 1, "W": 0, "Y": 1, "V": 0,
}

HYDROPHOBICITY_DICT = {
    "A":  1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C":  2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I":  4.5,
    "L":  3.8, "K": -3.9, "M":  1.9, "F":  2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V":  4.2,
}

FORMAL_CHARGE_DICT = {
    "A":  0.0, "R":  1.0, "N":  0.0, "D": -1.0, "C":  0.0,
    "Q":  0.0, "E": -1.0, "G":  0.0, "H":  0.0, "I":  0.0,
    "L":  0.0, "K":  1.0, "M":  0.0, "F":  0.0, "P":  0.0,
    "S":  0.0, "T":  0.0, "W":  0.0, "Y":  0.0, "V":  0.0,
}
STAND_ONE_RESI = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
]
STAND_RESI_SET = set(STAND_RESI_DICT.keys())
STAND_RESI_ONE_SET = set(STAND_RESI_DICT.values())
POLAR_DICT = {aa: int(aa in {"R", "N", "D", "C", "Q", "E", "H", "K", "S", "T", "Y"}) for aa in HYDROPHOBICITY_DICT}
AROMATIC_DICT = {aa: int(aa in {"F", "W", "Y", "H"}) for aa in HYDROPHOBICITY_DICT}
# fmt: on
