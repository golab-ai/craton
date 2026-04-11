import numpy as np
from ...chem.elements import Element
from ...chem.atom import Atom

HYDROGEN_BOND_ACCEPTOR = ["O","N"]
HYDROGEN_BOND_DORNOR = ["H"]
HALOGEN = ["F","Cl","Br","I"]


from dataclasses import dataclass
from typing import Dict, List, Optional
def normalize_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if not norm == 0 else v

@dataclass
class Site:
    type: str
    subtype:str = None
    atom: Atom = None
    adj_atom: Atom = None
    atoms: List[Atom] = None
    group_id: int = 0
    group: str = "LIG"
    chain: str = "S"
    
    @property
    def center(self):
        if self.atoms is not None:
            return np.mean([atom.coor for atom in self.atoms], axis=0)
        elif self.adj_atom is not None:
            return np.array(self.adj_atom.coor)
        else:
            return np.array(self.atom.coor)

    @property
    def h_center(self):
        if self.adj_atom is not None:
            return np.array(self.atom.coor)
        return None

    @property
    def normal(self):
        ring_vec1 = np.array(self.atoms[2].coor) - np.array(self.atoms[0].coor)
        ring_vec2 = np.array(self.atoms[4].coor) - np.array(self.atoms[0].coor)
        return normalize_vector(np.cross(ring_vec1, ring_vec2))

    def __repr__(self):
        if self.atoms is not None:
            return f"{self.group}-{self.group_id}-{self.type}-{self.center.tolist()}"
        if self.adj_atom is not None:
            atom_name = self.adj_atom.name if hasattr(self.adj_atom,"name") else self.adj_atom.element
            return f"{self.group}-{self.group_id}-{atom_name}-{self.type}"
        atom_name = self.atom.name if hasattr(self.atom,"name") else self.atom.element
        return f"{self.group}-{self.group_id}-{atom_name}-{self.type}"
    
    @property
    def site_script(self):
        if self.atoms is not None:
            atom_arr_str = ":".join([str(an.ID) for an in self.atoms])
            ss = f"{self.group}-{self.group_id}-{atom_arr_str}-{self.type}-{self.subtype}"
        else:
            ss = f"{self.group}-{self.group_id}-{self.atom.ID}-{self.type}-{self.subtype}"
        return ss


class Model:
    def __init__(self, molecule,atoms=None):
        self.molecule = molecule
        if atoms is None:
            atoms = [atom.ID for atom in molecule.Atoms]
        self.atoms = atoms

    def run(self):
        __FUNC = {  "hydrophobic": self.find_hydrophobic_atoms,
                "hbond_acceptor": self.find_hba,
                "hbond_donor": self.find_hbd,
               #"charged_atoms": self.find_charged,
                #"rings": self.find_rings,
                #"metal_bindings": self.find_metal_binding,
                "halogen_donor": self.find_hald,
                "halogen_acceptor": self.find_hala,
                }

        interaction_model = {kk:[] for kk in __FUNC.keys()}
        for an in self.atoms:
            atom = self.molecule.Atoms[an]
            for kk,ff in __FUNC.items():
                vv = ff(atom)
                if vv is not None:
                    interaction_model[kk].append(vv)
                    
        interaction_model["charged"] = []
        charge_group = {}
        if hasattr(self.molecule,"charge_group"):
            charge_group = self.molecule.charge_group
        #    if isinstance(self.molecule.charge_group,list):
        #        charge_group = {ii:cg for ii,cg in enumerate(self.molecule.charge_group)}
        #    else:
        #        charge_group = self.molecule.charge_group

        for kk,vv in charge_group.items(): 
            if len(set(vv[:-1]).intersection(set(self.atoms))) > 0:
                interaction_model["charged"].append(self.find_charge(vv))
        
        interaction_model["rings"] = []
        for kk,vv in getattr(self.molecule,"rings",{}).items():
            if len(set(vv[:-1]).intersection(set(self.atoms))) > 0:
                ring = self.find_rings(vv)
                if ring is not None:
                    interaction_model["rings"].append(ring)
        return interaction_model
            
    def find_hydrophobic_atoms(self,atom):
        if atom.element in ["C"]:
            neighbor = set([self.molecule.Atoms[an].element for an in atom.connectivity])
            if neighbor.issubset({"H","C"}):
                return Site(
                            type="hydrophobic",
                            atom = atom,
                            group=getattr(atom,"residue","LIG"),
                            group_id=getattr(atom,"residue_ID",atom.ID),
                            chain=getattr(atom,"chain_name","S")
                            )
        return None

    def find_hba(self,atom):
        if atom.element in HYDROGEN_BOND_ACCEPTOR:
            if len(atom.connectivity) <= Element.get(atom.element).valents[0]:
                return Site(
                            type="hba",
                            atom = atom,
                            group=getattr(atom,"residue","LIG"),
                            group_id=getattr(atom,"residue_ID",atom.ID),
                            chain=getattr(atom,"chain_name","S")
                            )
        return None

    def find_hbd(self,atom):
        if atom.element in HYDROGEN_BOND_DORNOR:
            #adj_an = atom.connectivity[0]
            adj_atom = self.molecule.Atoms[atom.connectivity[0]]
            if adj_atom.element in HYDROGEN_BOND_ACCEPTOR:
                return Site(
                            type="hbd",
                            atom=atom,
                            adj_atom=adj_atom,
                            group=getattr(adj_atom,"residue","LIG"),
                            group_id=getattr(adj_atom,"residue_ID",adj_atom.ID),
                            chain=getattr(adj_atom,"chain_name","S")
                            )
        return None

    def find_hala(self,atom):
        if atom.element in ["N","O","S"]:
            n_atoms = [self.molecule.Atoms[an] for an in atom.connectivity if self.molecule.Atoms[an].element in ["C","N","P","S"]]
            if len(n_atoms) == 1:
                return Site(
                            type="hala",
                            atom = atom,
                            adj_atom=n_atoms[0],
                            group=getattr(atom,"residue","LIG"),
                            group_id=getattr(atom,"residue_ID",atom.ID),
                            chain=getattr(atom,"chain_name","S")
                            )
        return None

    def find_hald(self,atom):
        if atom.element in HALOGEN and len(atom.connectivity) == 1:
            #an = atom.connectivity[0]
            adj_atom = self.molecule.Atoms[atom.connectivity[0]]
            if adj_atom.element == "C":
                return Site(type="hald",
                            atom=atom,
                            adj_atom=adj_atom,
                            group=getattr(adj_atom,"residue","LIG"),
                            group_id=getattr(adj_atom,"residue_ID",adj_atom.ID),
                            chain=getattr(adj_atom,"chain_name","S"))
        return None

    def find_charge(self,vv):
        return Site(type="charge",
                    subtype = vv[-1],
                    atoms=[self.molecule.Atoms[an]for an in vv[:-1]],
                    group=getattr(self.molecule.Atoms[vv[0]],"residue","LIG"),
                    group_id=getattr(self.molecule.Atoms[vv[0]],"residue_ID",self.molecule.Atoms[vv[0]].ID),
                    chain=getattr(self.molecule.Atoms[vv[0]],"chain_name","S")
                    )
        
    def find_rings(self,vv):
        if 4 < len(vv[:-1]) <= 6:
            fg = sum([self.molecule.Atoms[an].formal_charge for an in vv[:-1]])
            if vv[-1] != "nonar" and fg == 0:
                return Site(type="ring",
                            atoms=[self.molecule.Atoms[an] for an in vv[:-1]],
                            group=getattr(self.molecule.Atoms[vv[0]],"residue","LIG"),
                            group_id=getattr(self.molecule.Atoms[vv[0]],"residue_ID",self.molecule.Atoms[vv[0]].ID),
                            chain=getattr(self.molecule.Atoms[vv[0]],"chain_name","S")
                            )
        return None
                


        #rings: List[Ring] = []
        #ring_candidates = self.ob_mol.OBMol.GetSSSR()
        #for ring in ring_candidates:
        #    if 4 < len(ring._path) <= 6:
        #        ring_atoms = [self.atom_dict_reverse[atom] for atom in ring._path]
        #        if (ring.IsAromatic() or ring_is_planar(ring_atoms)) and ions_not_in_ring(ring_atoms):
        #            ring_vec1 = ring_atoms[2].coor - ring_atoms[0].coor
        #            ring_vec2 = ring_atoms[0].coor - ring_atoms[4].coor
        #            center = centroid([atom.coor for atom in ring_atoms])
        #            normal = normalize_vector(np.cross(ring_vec1, ring_vec2))
        #            rings.append(Ring(atoms=ring_atoms, normal=normal, center=center))
        #return rings

    def find_metal_binding_ligand(self):
        metal_binding: List[MetalBinding] = []
        for atom in self.atoms:
            ob_atom = self.atom_dict[atom.atom_orig_idx]
            n_atoms = pybel.ob.OBAtomAtomIter(ob_atom.OBAtom)
            n_atoms_atomicnum = [atom.GetAtomicNum() for atom in n_atoms]
            if ob_atom.atomicnum == 8:
                if n_atoms_atomicnum.count(1) == 1 and len(n_atoms_atomicnum) == 2:  # Oxygen in alcohol (R-[O]-H)
                    metal_binding.append(MetalBinding(atom=atom, type="O", group="alcohol", location="ligand"))
                if True in [n.IsAromatic() for n in n_atoms] and not ob_atom.OBAtom.IsAromatic():  # phenolate oxygen
                    metal_binding.append(MetalBinding(atom=atom, type="O", group="phenolate", location="ligand"))
            if ob_atom.atomicnum == 6:  # carbon atom
                if n_atoms_atomicnum.count(8) == 2 and n_atoms_atomicnum.count(6) == 1:  # It's a carboxylate group
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        atom_p = self.atom_dict_reverse[neighbor.GetIdx()]
                        metal_binding.append(
                            MetalBinding(atom=atom_p, type="O", group="carboxylate", location="ligand")
                        )
            if ob_atom.atomicnum == 15:  # It's a phosphor atom
                if n_atoms_atomicnum.count(8) >= 3:  # It's a phosphoryl
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        atom_p = self.atom_dict_reverse[neighbor.GetIdx()]
                        metal_binding.append(MetalBinding(atom=atom_p, type="O", group="phosphoryl", location="ligand"))
                if n_atoms_atomicnum.count(8) == 2:
                    for neighbor in [n for n in n_atoms if n.GetAtomicNum() == 8]:
                        atom_p = self.atom_dict_reverse[neighbor.GetIdx()]
                        self.metal_binding.append(
                            MetalBinding(atom=atom_p, type="O", group="phosphor.other", location="ligand")
                        )
            if ob_atom.atomicnum == 7:
                if n_atoms_atomicnum.count(6) == 2:
                    metal_binding.append(
                        MetalBinding(atom=atom, type="N", group="imidazole/pyrrole", location="ligand")
                    )
            if ob_atom.atomicnum == 16:
                if True in [n.IsAromatic() for n in n_atoms] and not ob_atom.OBAtom.IsAromatic():  # Thiolate
                    metal_binding.append(MetalBinding(atom=atom, type="S", group="thiolate", location="ligand"))
                if set(n_atoms_atomicnum) == {26}:
                    metal_binding.append(
                        MetalBinding(atom=atom, type="S", group="iron-sulfur.cluster", location="ligand")
                    )
        return metal_binding

    def find_metal_binding_protein(self):
        metal_binding: List[MetalBinding] = []
        for res in self.binding_site_residues:
            if res.name in ("ASP", "GLU", "SER", "THR", "TYR"):
                for atom in pybel.ob.OBResidueAtomIter(res):
                    if atom.GetType().startswith("O") and res.OBResidue.GetAtomProperty(atom.OBAtom, 8):
                        metal_binding.append(
                            MetalBinding(
                                self.protein.atom_dict_reverse[atom.GetIdx()],
                                "O",
                                res.name,
                                "sidechain",
                            )
                        )
            elif res.name == "HIS":  # Look for nitrogen here
                for atom in pybel.ob.OBResidueAtomIter(res):
                    if atom.GetType().startswith("O") and res.OBResidue.GetAtomProperty(atom.OBAtom, 8):
                        metal_binding.append(
                            MetalBinding(
                                self.protein.atom_dict_reverse[atom.GetIdx()],
                                "N",
                                res.name,
                                "sidechain",
                            )
                        )
            elif res.name == "CYS":
                for atom in pybel.ob.OBResidueAtomIter(res):
                    if atom.GetType().startswith("S") and res.OBResidue.GetAtomProperty(atom.OBAtom, 8):
                        metal_binding.append(
                            MetalBinding(
                                self.protein.atom_dict_reverse[atom.GetIdx()],
                                "S",
                                res.name,
                                "sidechain",
                            )
                        )
            for atom in pybel.ob.OBResidueAtomIter(res):
                if atom.GetType().startswith("O") and res.OBResidue.GetAtomProperty(atom.OBAtom, 2):
                    metal_binding.append(
                        MetalBinding(
                            self.protein.atom_dict_reverse[atom.GetIdx()],
                            "O",
                            res.name,
                            "mainchain",
                        )
                    )
        return metal_binding
