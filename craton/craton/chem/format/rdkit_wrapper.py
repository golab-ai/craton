from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, rdmolops


RDLogger.DisableLog("rdApp.warning")


@dataclass
class RdkitMolecluleInfo:
    molecule: object = None
    smiles: str = ""
    inchi: str = ""
    molecule_noh: object = None
    canonical_smiles: str = ""
    inchi_key: str = ""
    mol_string: str = ""


class RdKitWrapper:
    def __init__(self, mol, chiral=True):
        self.mol = mol
        mol_script = MolFile("normal")
        mol_script.import_moleobj(mol, has3d="yes")
        self.rd_mole_info = RdkitMolecluleInfo(
            molecule := Chem.MolFromMolBlock(mol_script.script, removeHs=False),
            smiles := Chem.MolToSmiles(molecule),
            inchi := Chem.inchi.MolToInchi(molecule),
            molecule_noh=Chem.MolFromMolBlock(mol_script.script),
            canonical_smiles=AllChem.CanonSmiles(smiles),
            inchi_key=Chem.inchi.InchiToInchiKey(inchi),
            mol_string=mol_script.script,
        )
        if chiral:
            Chem.rdmolops.AssignAtomChiralTagsFromStructure(self.rd_mole_info.molecule)

    @property
    def smiles(self):
        return self.rd_mole_info.smiles

    @property
    def inchi(self):
        return self.rd_mole_info.inchi

    @property
    def canonical_smiles(self):
        return self.rd_mole_info.canonical_smiles

    @property
    def inchi_key(self):
        return self.rd_mole_info.inchi_key

    @property
    def mol_string(self):
        return self.rd_mole_info.mol_string

    @property
    def molecule(self):
        return self.rd_mole_info.molecule

    @property
    def molecule_noh(self):
        return self.rd_mole_info.molecule_noh

    def convert2d(self):
        AllChem.Compute2DCoords(self.molecule)
        AllChem.Compute2DCoords(self.molecule_noh)

    def map_to_no_hydrogen(self):
        heavy_atoms_mapper_with_hydreogn_removed = {}
        index = 0
        for i, atom in enumerate(self.molecule.GetAtoms()):
            if atom.GetAtomicNum() != 1:
                heavy_atoms_mapper_with_hydreogn_removed[i] = index
                index += 1
        return heavy_atoms_mapper_with_hydreogn_removed

    def save(self, name=None, output_dir=".", include_hydrogen=False):
        self.convert2d()
        if name is None:
            name = self.mol.name
        if include_hydrogen:
            img = Draw.MolToImage(self.molecule, legend=name)
        else:
            img = Draw.MolToImage(self.molecule_noh, legend=name)
        img.save(str(Path(output_dir) / (name + ".png")))



