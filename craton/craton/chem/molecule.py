import itertools

from copy import deepcopy
from typing import Iterable, List, Any

import networkx as nx
import numpy as np
#import string

from ..utils import logger
#from .atom import VS, Atom

from ..utils.geometry import find_center
from ..property.property import calc_inertia

#from .topology import Angle, Bond, Constrain, Dihedral, Improper, Pair

from .constants import ATOM_DEFAULT_ATTRIBUTES

from ._molecule import _Molecule

class Molecule(_Molecule):

    def __init__(self, style=""):
        self.style = style
        self.steps = []
        super(Molecule,self).__init__(self.style)

    #######暂时用到#####################
    @property
    def scan_term(self):
        return self.torsions

    @scan_term.setter
    def scan_term(self,a):
        self.torsions = a

    @scan_term.deleter
    def scan_term(self):
        del self.torsions

    @property
    def multi(self):
        return self.multiple

    @multi.setter
    def multi(self,a):
        self.multiple = a

    @property
    def id(self):
        return self.ID

    @id.setter
    def id(self,a):
        self.ID = a

    ####################################
    @property
    def atom_count(self):
        if not hasattr(self,"_atom_count"):
            self._atom_count = len(self.Atoms)
        else:
            if self._atom_count != len(self.Atoms):
                self._atom_count = len(self.Atoms)
        return self._atom_count
    
    @atom_count.setter
    def atom_count(self,n):
        self._atom_count = n

    @property
    def bond_count(self):
        if not hasattr(self,"_bond_count"):
            self._bond_count = len(self.Bonds) if hasattr(self,"Bonds") else 0
        else:
            if self._bond_count != len(self.Bonds):
                self._bond_count = len(self.Bonds)
        return self._bond_count
    
    @bond_count.setter
    def bond_count(self,n):
        self._bond_count = n
    
    @property
    def angle_count(self):
        return len(self.Angles) if hasattr(self,"Angles") else 0
    
    @property
    def dihedral_count(self):
        return len(self.Dihedrals) if hasattr(self,"Dihedrals") else 0

    @property
    def improper_count(self):
        return len(self.Impropers) if hasattr(self,"Impropers") else 0

    @property
    def connectivity(self):
        return [aa.connect for aa in self.Atoms]

    @property
    def elements(self):
        return [aa.elem for aa in self.Atoms]

    @property
    def bond_type(self):
        return [aa.bond_type for aa in self.Atoms]

    @property
    def local(self):
        return [aa.local for aa in self.Atoms] 
        
    @property
    def formal_charge(self):
        return [aa.formal_charge for aa in self.Atoms]

    @property
    def coordinates(self):
        return [aa.coor for aa in self.Atoms]

    @property
    def atom_name(self):
        return [aa.atom_name for aa in self.Atoms]

    @property
    def residue(self):
        return [aa.residue for aa in self.Atoms]
    
    @property
    def residue_ID(self):
        return [aa.residue_ID for aa in self.Atoms]
    
    @property
    def plate(self):
        return [aa.plate for aa in self.Atoms]
    
    @property
    def atom_type_name(self):
        return [aa.atom_type_name for aa in self.Atoms]

    @property
    def ff_charge(self):
        """
        Required by conformation DB insertion
        """
        if not hasattr(self.Atoms[0], "point_charge"):
            return None
        return [aa.ff_charge for aa in self.Atoms]

    @property
    def esp_charge(self):
        """
        Required by conformation DB insertion
        """
        if not hasattr(self.Atoms[0], "esp_charge"):
            return None
        return [aa.esp_charge for aa in self.Atoms]



    @property
    def mulliken_charge(self):
        """
        Required by conformation DB insertion
        """
        if not hasattr(self.Atoms[0], "mulliken_charge"):
            return None

        return [aa.mulliken_charge for aa in self.Atoms]

    @property
    def charges(self) -> {str: [float]}:
        d = {
            "point_charge":self.point_charge,
            "ff_charge": self.point_charge,
            "esp_charge": self.esp_charge,
            "mulliken_charge": self.mulliken_charge,
            }
        tmp = []
        for aa,bb in d.items():
            if bb == None:
                tmp.append(aa)
        for aa in tmp:
            del d[aa]
        return d

    @property
    def point_charge(self):
        try:
            return [aa.ff_charge for aa in self.Atoms]
        except:
            return None

    @property
    def point_charge_sum(self):
        return sum([atom.point_charge for atom in self.Atoms])

    @property
    def atom_local(self, ignoreH=True):
        return [aa.local for aa in self.Atoms if aa.elem not in ["H", "F", "Cl", "Br", "I"]]

    @property
    def element_set(self):
        return set(self.elements)

    @property
    def element_count(self):
        return {elem:sum([1 for atom in self.Atoms if atom.element == elem]) for elem in self.element_set}

    @property
    def zelement(self):
        __main_elem = ["C","O","N","S","P",]
        __halogen_elem = ["F","Cl","Br","I",]
        __rare_elem = ["B","Si","As","Se""Te"]
        __order = ["Z","X","R"]
        label = ""
        for nn,rr in enumerate([__main_elem,__halogen_elem,__rare_elem]):
            label += __order[nn]
            for ee in rr:
                if ee in self.element_set:
                    label += ee
            label += "-"
        return label[:-1]

    @property
    def molecule_mass(self):
        if not hasattr(self,"_molecule_mass"):
            self._molecule_mass = sum([atom.mass for atom in self.Atoms])
        return self._molecule_mass
    
    @molecule_mass.setter
    def molecule_mass(self,f):
        self._molecule_mass = f

    @property
    def mass(self):
        if not  hasattr(self,"_molecule_mass"):
            self._molecule_mass = sum([atom.mass for atom in self.Atoms])
        return self._molecule_mass

    @mass.setter
    def mass(self,a):
        self._molecule_mass = a

    @property
    def formula(self):
        return self.get_formula()

    @property
    def net_charge(self):
        if self.style == "protein":
            return sum([atom.ff_charge for atom in self.Atoms])
        return sum([atom.formal_charge for atom in self.Atoms])


    @property
    def heavy_atoms(self):
        return len([1 for aa in self.Atoms if aa.elem != "H"])

    @property
    def torsion_number(self):
        if not hasattr(self,"_torsion_number"):
            if hasattr(self,"torsions"):
                self._torsion_number = len(self.torsions)
            else:
                return None
        return self._torsion_number

    @torsion_number.setter
    def torsion_number(self,n):
        self._torsion_number = n

    @property
    def inchikey_3d_structure_flag(self):
        if not hasattr(self, "_inchikey_3d_structure_flag"):
            return False
        return self._inchikey_3d_structure_flag

    @inchikey_3d_structure_flag.setter
    def inchikey_3d_structure_flag(self,b: bool):
        self._inchikey_3d_structure_flag = b

    @property
    def renew_inchikey_flag(self):
        if not hasattr(self, "_renew_inchikey_flag"):
            return False
        return self._renew_inchikey_flag

    @renew_inchikey_flag.setter
    def renew_inchikey_flag(self,b: bool):
        self._renew_inchikey_flag = b

    @property
    def inchi_key(self):
        if self._renew_inchikey_flag:
            self.get_inchikey()
        return self._inchi_key
    
    @inchi_key.setter
    def inchi_key(self,s: str):
        self._inchi_key = s

    @property
    def inchi(self):
        if self._renew_inchikey_flag:
            self.get_inchikey()
        return self._inchi
    
    @inchi.setter
    def inchi(self,s: str):
        self._inchi = s    

    @property
    def smiles(self):
        if self._renew_inchikey_flag:
            self.get_inchikey()
        return self._smiles
    
    @smiles.setter
    def smiles(self,s: str):
        self._smiles = s

    @property
    def mole_name(self):
        return self._molecule_name

    @mole_name.setter
    def mole_name(self,s: str):
        self._molecule_name = s

    @property
    def name(self):
        return self._molecule_name
    
    @name.setter
    def name(self,s: str):
        self._molecule_name = s

    @property
    def molecule_name(self):
        return self._molecule_name

    @molecule_name.setter
    def molecule_name(self,s):
        self._molecule_name = s


    @property
    def ring_stru(self):
        return self.ring_blocks

    @ring_stru.setter
    def ring_stru(self,a):
        self.ring_blocks = a

    @property
    def ring_dict(self):
        #if not self._rings_defined:
        #    logger.warning("Current molecule's ring_dict is not defined. Returning a dummy one.")
        return self._rings

    @ring_dict.setter
    def ring_dict(self, x):
        self._rings_defined = True
        self._rings = x

    @property
    def rings(self):
        #if not self._rings_defined:
        #    logger.warning("Current molecule's ring_dict is not defined. Returning a dummy one.")
        return self._rings

    @rings.setter
    def rings(self,a):
        self._rings_defined = True
        self._rings = a

    @property
    def _ring_dict(self):
        if not self._rings_defined:
            logger.warning("Current molecule's ring_dict is not defined. Returning a dummy one.")
        return self._rings   
    
    @_ring_dict.setter
    def _ring_dict(self,x):
        self._rings_defined = True
        self._rings = x


    @property
    def ring_number(self):
        if not hasattr(self,"_ring_number"):
            self._ring_number = len(getattr(self, "ring_dict", {}))
        return self._ring_number
        return len(getattr(self, "ring_dict", {}))

    @ring_number.setter
    def ring_number(self,x):
        self._ring_number = x

    @property
    def ring_size(self):
        if not hasattr(self,"_ring_size"):
            self._ring_size =  [(len(vv) - 1) for kk, vv in self._rings.items()]
        return self._ring_size

    @ring_size.setter
    def ring_size(self,a):
        self._ring_size = a

    @property
    def ring_property(self):
        if not hasattr(self,"_ring_property"):
            self._ring_property = [vv[-1] for kk, vv in self._rings.items()]
        return self._ring_property

    @ring_property.setter
    def ring_property(self,a):
        self._ring_property = a 

    @property
    def fused_ring_number(self):
        tmp = []
        for __, frag in self.elem_frag.items():
            if frag["label"][0] == "C":
                tmp.append(0)
            else:
                if frag["label"][1] == "F":
                    tmp.append(int(frag["label"][2]))
                else:
                    tmp.append(1)
        return tmp

    @property
    def largest_fuse_ring_number(self):
        return max(self.fused_ring_number)

    @property
    def ring_stru_number(self):
        return len(self.ring_stru)

    def get_term_by_indices(self, index: list):
        if len(index) == 2:
            for term in self.Bonds:
                if tuple(term.atoms) == tuple(index):
                    return term
        elif len(index) == 3:
            for term in self.Angles:
                if tuple(term.atoms) == tuple(index):
                    return term
        elif len(index) == 4:
            for term in self.Dihedrals:
                if tuple(term.atoms) == tuple(index):
                    return term
        else:
            logger.error("Invalid term index.")
        return None

    @property
    def constrain_term(self):
        if hasattr(self, "constrain"):
            return self.constrain[0].atoms
        else:
            return None

    @property
    def constrain_value(self):
        if hasattr(self, "constrain"):
            return self.constrain[0].fix_value
        else:
            return None
 
    @property
    def ef_number(self):
        return len(self.elem_frag)

    @property
    def ef_number_ignore_side(self):
        return len([0 for __, frag in self.elem_frag.items() if frag["ef_type"] != "side"])

    @property
    def mol_script(self):
        from .format.mol_parse import MolData
        molobj = MolData()
        return molobj._convert(self)

    def write_mol(self, filename):
        from .format.mol_parse import MolData

        molobj = MolData()
        molobj.import_moleobj(self, has3d="yes")
        with open(filename, "w") as f:
            f.write(molobj.script)

    @property
    def ob_mol(self):
        if self._ob_mol is None:
            return self.to_pybel_mol()
        else:
            return self._ob_mol

    @ob_mol.setter
    def ob_mol(self, mol):
        self._ob_mol = mol

    def to_pybel_mol(self):
        from openbabel import pybel

        ob_mol = pybel.ob.OBMol()
        for atom in self.Atoms:
            ob_atom = ob_mol.NewAtom()
            ob_atom.SetAtomicNum(atom.atom_number)
            ob_atom.SetVector(*atom.coor)
        for bond in self.Bonds:
            ###ob_mol.AddBond(bond.a1 + 1, bond.a2 + 1, int(bond.get_type(self)))
            ob_mol.AddBond(bond.a1+1,bond.a2+1, 1)
        ob_mol.PerceiveBondOrders()
        return pybel.Molecule(ob_mol)

    @property
    def sdf_script(self):
        from .format.mol_parse import MolData

        sdfobj = MolData()
        sdfobj.import_moleobj(self, has3d="yes")
        return sdfobj.script

    @property
    def inertia(self):
        if hasattr(self, "inertia_principal"):
            return self.inertia_principal
        inertia = calc_inertia(self)
        if isinstance(inertia, tuple):
            return inertia[0]
        return inertia

    @inertia.setter
    def inertia(self, values):
        self.inertia_principal = np.asarray(values)

    def get_center(self, ct="com"):
        coor = self.coordinates
        mass = [aa.mass for aa in self.Atoms]
        return find_center(coor, mass, center_type=ct)

    def get_energy(self, update_flag=False):
        if update_flag:
            self.update_topol_value()
            self.energy = FFcalculator.single_mole_energy(self)
        if not hasattr(self, "energy"):
            if not hasattr(self.Bonds[0], "value"):
                self.update_topol_value()
            self.energy = FFcalculator.single_mole_energy(self)
        return self.energy

    @property
    def get_force(self, update_flag=False):
        if update_flag:
            self.update_topol_value()
            self.force = FFcalculator.single_mole_force(self)
        if not hasattr(self, "force"):
            if not hasattr(self.Bonds[0], "value"):
                self.update_topol_value()
            self.force = FFcalculator.single_mole_force(self)
        return self.force

    @property
    def get_hessian(self, update_flag=False):
        if update_flag:
            self.update_topol_value()
            self.hessian = FFcalculator.single_mole_hessian(self)
        if not hasattr(self, "hessian"):
            if not hasattr(self.Bonds[0], "value"):
                self.update_topol_value()
            self.hessian = FFcalculator.single_mole_hessian(self)
        return self.hessian

    @property
    def get_freq(self, update_flag=False):
        if update_flag:
            self.update_topol_value()
            self.frep = FFcalculator.single_mole_freq(self)
        if not hasattr(self, "freq"):
            if not hasattr(self.Bonds[0], "value"):
                self.update_topol_value()
            self.freq = FFcalculator.single_mole_freq(self)
        return self.freq

    @property
    def get_frequency(self, update_flag=False):
        return self.get_freq(update_flag=update_flag)

    @property
    def AtomsNumber(self):
        return len(self.Atoms)
