import os
from .general import tmp_name
from .convertors import PDBParser, MMCIFParser, PDBIO, MMCIFIO  # type: ignore
from .convertors import structure_load, structure_dump
from rdkit import Chem


def tranverse_folder(folder):
    filepath_list = []
    for root, _, files in os.walk(folder):
        for file in files:
            filepath_list.append(os.path.join(root, file))
    return filepath_list


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


def remove_solvent(file_path, save_path, solvents={"HOH", "WAT", "GOL", "DOD", "SOL"}):
    structure = structure_load(file_path)
    solvents = {x.upper() for x in solvents}
    for model in structure:  # type: ignore
        for chain in model:
            solvent_residues = [res for res in chain if res.get_resname().strip() in solvents]
            for res in solvent_residues:
                chain.detach_child(res.id)
    structure_dump(structure, save_path)
