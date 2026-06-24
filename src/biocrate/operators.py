import os

import numpy as np

from .constants import STAND_RESI_DICT
from .general import structure_load, structure_dump


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
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(file_path, sanitize=True, removeHs=True)
    molecules = [mol for mol in suppl if mol is not None]
    if not molecules:
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


def remove_hetatm(file_path, save_path):
    structure = structure_load(file_path)
    for model in structure:  # type: ignore
        for chain in model:
            hetatm_residues = [res for res in chain if res.id[0].strip() != ""]
            for res in hetatm_residues:
                chain.detach_child(res.id)
    structure_dump(structure, save_path)


def MolFromSmiles(smiles):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    params = Chem.RemoveHsParameters()
    params.removeDummyNeighbors = True
    return Chem.RemoveHs(mol, params)


# cano smiles
def cano_smiles(smiles: str, isomeric: bool = True):
    from rdkit import Chem

    smiles = smiles.strip()
    mol = MolFromSmiles(smiles)
    if not mol:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric)


# cano rxn
def cano_rxn(rxn, exchange_pos=False, isomeric=True):
    if not isinstance(rxn, str) or not rxn.strip():
        return None
    data = rxn.strip().split(">")
    if len(data) != 3:
        return None

    reactants = tuple(cano_smiles(each, isomeric) for each in data[0].split("."))
    products = tuple(cano_smiles(each, isomeric) for each in data[-1].split("."))

    if None in reactants or None in products:
        return None

    reactants = sorted(reactants)
    products = sorted(products)

    if exchange_pos:
        new_rxn = f"{'.'.join(products)}>>{'.'.join(reactants)}"
    else:
        new_rxn = f"{'.'.join(reactants)}>>{'.'.join(products)}"
    return new_rxn


def get_substrates(rxn_smiles: str):
    # Implementation for extracting substrates from reaction SMILES
    from rdkit import Chem

    reactants = rxn_smiles.strip().split(">", 1)[0].split(".")
    best_reactant = None
    best_count = -1
    for reactant in reactants:
        mol = Chem.MolFromSmiles(reactant)
        if mol is None:
            continue
        count = mol.GetNumHeavyAtoms()
        if count > best_count:
            best_reactant = reactant
            best_count = count
    return best_reactant


def get_molecule_features(smiles: str):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None, None, None, None, None
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.InchiToInchiKey(inchi)  # type: ignore
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)  # type: ignore
    weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)  # type: ignore
    charge = Chem.GetFormalCharge(mol)
    return inchi, inchikey, formula, weight, charge


def smiles2name(smiles: str) -> str | None:
    import pubchempy as pcp

    compounds = pcp.get_compounds(smiles, namespace="smiles")
    if not compounds:
        return None
    compound = compounds[0]
    return compound.iupac_name


# onehot encoder
def onehot_encode(index: str, alphabet: list[str] | dict[str, int]) -> list:
    onehot = [0] * len(alphabet)
    if index not in alphabet:
        return onehot
    idx = alphabet.index(index) if isinstance(alphabet, list) else alphabet[index]
    onehot[idx] = 1
    return onehot


def letter3to1(res_name):
    return STAND_RESI_DICT.get(res_name, "X")


def norm_angle_deg(angle):
    if angle is None or np.isnan(angle):
        return [0.0, 0.0]

    rad = np.deg2rad(angle)
    return [float(np.sin(rad)), float(np.cos(rad))]


def norm_hbond_energy(e, min_energy=-10.0):
    if e is None or np.isnan(e):
        return 0.0

    e = np.clip(e, min_energy, 0.0)
    return float(-e / abs(min_energy))
