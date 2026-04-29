import os
from .general import tmp_name
from .convertors import PDBParser, MMCIFParser  # type: ignore
from rdkit import Chem


def tranverse_folder(folder):
    filepath_list = []
    for root, _, files in os.walk(folder):
        for file in files:
            filepath_list.append(os.path.join(root, file))
    return filepath_list


def load_structure(file_path):
    ext = os.path.splitext(file_path)[-1].lower()

    if ext == ".pdb":
        return PDBParser(QUIET=True).get_structure(tmp_name(), file_path)
    elif ext in [".cif", ".mmcif"]:
        return MMCIFParser(QUIET=True).get_structure(tmp_name(), file_path)


def uid2path(uid, is_dir=False, format="cif"):
    if is_dir:
        return os.path.join(uid[:2], uid[2:4], uid[4:6], uid)
    return os.path.join(uid[:2], uid[2:4], uid[4:6], f"{uid}.{format}")


def pdbid2path(pdbid, is_dir=False, format="cif"):
    if is_dir:
        return os.path.join(pdbid[1:3], pdbid)
    return os.path.join(pdbid[1:3], f"{pdbid}.{format}")


def parse_sdf(file_path):
    suppl = Chem.SDMolSupplier(file_path, sanitize=True, removeHs=True)
    molecules = []
    for mol in suppl:
        if mol is not None:
            molecules.append(mol)
    if len(molecules) == 0:
        return Chem.SDMolSupplier(file_path, sanitize=False, removeHs=False)  # if no valid molecules found, return unsanitized
    return molecules
