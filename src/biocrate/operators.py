import os
from .general import structure_load, structure_dump
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


def remove_hetatm(file_path, save_path):
    structure = structure_load(file_path)
    for model in structure:  # type: ignore
        for chain in model:
            hetatm_residues = [res for res in chain if res.id[0].strip() != ""]
            for res in hetatm_residues:
                chain.detach_child(res.id)
    structure_dump(structure, save_path)


def MolFromSmiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    params = Chem.RemoveHsParameters()
    params.removeDummyNeighbors = True
    return Chem.RemoveHs(mol, params)


# cano smiles
def cano_smiles(smiles: str, isomeric: bool = True):
    smiles = smiles.strip()
    if Chem is None:
        return smiles or None
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

    reactants = [cano_smiles(each, isomeric) for each in data[0].split(".")]
    products = [cano_smiles(each, isomeric) for each in data[-1].split(".")]

    if None in reactants or None in products:
        return None

    reactants = sorted(reactants)  # type: ignore
    products = sorted(products)  # type: ignore

    if exchange_pos:
        new_rxn = f"{'.'.join(products)}>>{'.'.join(reactants)}"
    else:
        new_rxn = f"{'.'.join(reactants)}>>{'.'.join(products)}"
    return new_rxn


def get_substrates(rxn_smiles: str):
    # Implementation for extracting substrates from reaction SMILES
    reactants = rxn_smiles.strip().split(">")[0].split(".")
    if Chem is None:
        return max(reactants, key=len) if reactants else None
    heavy_atoms = {}
    for reactant in reactants:
        heavy_atoms[reactant] = Chem.MolFromSmiles(reactant).GetNumHeavyAtoms()
    sorted_reactants = sorted(heavy_atoms.items(), key=lambda x: x[1], reverse=True)
    return sorted_reactants[0][0] if sorted_reactants else None


def get_molecule_features(smiles: str):
    if Chem is None:
        return None, None, None, None, None
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
