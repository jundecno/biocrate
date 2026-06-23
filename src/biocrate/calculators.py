from .general import structure_load
from .constants import H_ATOMS
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Polypeptide import is_aa


def extract_pocket_residues(file_path, comp_id="", radius=4.5):
    structure = structure_load(file_path)
    model = structure[0]  # type: ignore
    ligand_heavy_atoms = []
    model_heavy_atoms = []
    comp_id = comp_id.strip()

    for res in model.get_residues():
        is_ligand = res.resname.strip() == comp_id and res.id[0].strip() != ""
        for atom in res:
            if atom.element.upper() in H_ATOMS:
                continue
            if is_ligand:
                ligand_heavy_atoms.append(atom)
            else:
                model_heavy_atoms.append(atom)

    if not ligand_heavy_atoms or not model_heavy_atoms:
        return set()

    ns = NeighborSearch(model_heavy_atoms)
    pocket_res_keys = set()

    for l_atom in ligand_heavy_atoms:
        close_resi = ns.search(l_atom.coord, radius, level="R")

        for c_res in close_resi:
            c_chain = c_res.get_parent()  # type: ignore
            if c_res.resname != comp_id and is_aa(c_res, standard=False):  # type: ignore
                pocket_res_keys.add((c_chain.id, c_res.id))  # type: ignore
    return pocket_res_keys


## generate mols
def generate_mol_conformation(smiles, n_repeat=20):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)

    ps = AllChem.ETKDGv3()  # type: ignore
    ps.randomSeed = 42
    ps.pruneRmsThresh = 0.75

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=n_repeat, params=ps)  # type: ignore
    if len(cids) == 0:
        return None

    if AllChem.MMFFHasAllMoleculeParams(mol):  # type: ignore
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=16)  # type: ignore
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(mol, numThreads=16)  # type: ignore

    energies = [r[1] if r[0] == 0 else float("inf") for r in res]
    best_conf_id = energies.index(min(energies))

    conf = Chem.Conformer(mol.GetConformer(best_conf_id))

    new_mol = Chem.Mol(mol)
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(conf, assignId=True)

    return Chem.RemoveHs(new_mol)
