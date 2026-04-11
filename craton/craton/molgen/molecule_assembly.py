from copy import deepcopy
from ..chem.molecule import Molecule
import numpy as np
from ..chem.elements import get_bonded_type_distance
from ..utils.geometry import *

from craton import molxpert as MX
from ..chemistry.molecule import Molecule
#tetrahedron, triangle, line, change_bond, calculate_dihedral, rotation_dihedral,change_angle,line_to_line_angle

#molecule = deepcopy(molecule_target)
#connects = [atom.connect for atom in molecule.Atoms if atom.elem not in ["H"]]

hybrid = {"C4": tetrahedron,"N4": tetrahedron,"N3": tetrahedron,"O2": tetrahedron,
          "S4": tetrahedron,"S2": tetrahedron,
          "C3":triangle,"N2":triangle,"O1":triangle,
          "C2":line,"N1":line,
            "H1": "end",
         }

plate_type = ["C3","N2"]
line_type = ["C2","N1"]


def create_mole_json(molecule):
    pass


class MolAssembly:
    def __init__(self,config) -> None:
        self.config = config



    def change_atom_number(self,ref_Atoms,cut_atoms,fr):
        frag = Molecule("molecule")
        new_number = {}
        ii = 0
        for atom in ref_Atoms:
            if atom.ID not in cut_atoms:
                new_number[atom.ID] = ii
                ii += 1
        Atoms = []
        for atom in ref_Atoms:
            if atom.ID not in cut_atoms:
                atom.ID = new_number[atom.ID]
                connect = []
                bond_type = []
                for an,at in zip(atom.connect,atom.bond_type):
                    if an not in cut_atoms:
                        connect.append(new_number[an])
                        bond_type.append(at)
                atom.connectivity = connect
                atom.bond_type = bond_type
                Atoms.append(atom)
        frag.Atoms = Atoms
        frag.frag_connect_atoms = [new_number[fr[0][0]],new_number[fr[1][0]]]

        return frag

    def assign_atom_number(self):
        fragments = []
        for mn,frag in enumerate(self.config["molecules"]):
            f_cut_atoms = frag.find_side_componend(frag.f_cut_atoms[1],frag.f_cut_atoms[0]) + [frag.f_cut_atoms[1]]
            r_cut_atoms = frag.find_side_componend(frag.r_cut_atoms[1],frag.r_cut_atoms[0]) + [frag.r_cut_atoms[1]]
            fragments.append(self.change_atom_number(deepcopy(frag.Atoms),f_cut_atoms+r_cut_atoms,[frag.f_cut_atoms,frag.r_cut_atoms]))  
            if mn == 0:
                first_frag = self.change_atom_number(deepcopy(frag.Atoms),r_cut_atoms,[frag.f_cut_atoms,frag.r_cut_atoms])
            if mn == len(self.config["molecules"]) -1:
                last_frag = self.change_atom_number(deepcopy(frag.Atoms),f_cut_atoms,[frag.f_cut_atoms,frag.r_cut_atoms])
        
        self.fragments = fragments + [last_frag]
        self.first_frag = first_frag
        self.last_frag = last_frag
        self.fragment_number = [nn if ii not in [0,len(self.fragments)] else nn - 1 for ii,nn in enumerate(self.config["molecule_number"]) ]
        self.fragment_number.append(1)

    def fusion_unit(self,monomers):
        monomer_mol = self.molecule_dict[monomers[0]]
        if len(monomer_mol) > 0:
            for mono in monomers[1:]:
                pass

        return monomer_mol


    def assembly(self):
        #####[[mole1,mole2],[mole3],[mole1,mole4],[]]
        mols = self.config["assembley"]["monomers"]
        mol_number = self.config["assembley"]["monomer_number"]
        tmp_mol = [m for mm in mols for m in mm]
        tmp_mol = list(set(tmp_mol))
        molecules = MX.molecule_create(tmp_mol)
        molecules = MX.molecule_structure(molecules)

        self.molecule_dict = {molecule.mole_name:molecule for molecule in molecules}

        self.molecule = deepcopy(self.molecule_dict(mols[0][0]))


        for monomers,nn in zip(mols,mol_number):
            self.assembly_unit(monomers,nn)


        molecule = self.first_frag

        pre_connectivity = self.first_frag.frag_connect_atoms
        coordinates_type = hybrid(f"{molecule.Atoms[pre_connectivity[0]].element}{len(molecule.Atoms[pre_connectivity[0]].connect)+1}")
        
        for fragment,nn in zip(self.fragments,self.fragment_number):
            connectivity = fragment.frag_connect_atoms
            for ii in range(nn):
                an_shift = len(molecule.Atoms)
                dis = dis = get_bonded_type_distance(molecule.Atoms[pre_connectivity[0]].element,fragment.Atoms[connectivity[0]].element,"1")
                coor_shift = coordinates_type(molecule.Atoms[pre_connectivity[0]].coordinates,dis,[molecule.Atoms[an].coordinates for an in molecule.Atoms[pre_connectivity[0]].connectivity])
                for atom in fragment.Atoms:
                    molecule.Atoms.append(deepcopy(atom))
                    molecule.Atoms[-1].ID = molecule.Atoms[-1].ID + an_shift
                    molecule.Atoms[-1].connectivity = [an + an_shift for an in molecule.Atoms[-1].connectivity]  
                    molecule.Atoms[-1].coordiantes = [coor+coor_shift[kk] for coor,kk in enumerate(molecule.Atoms[-1].coordinates)]         
                    if atom.ID == connectivity[0]:
                        molecule.Atoms[-1].connectivity = [pre_connectivity[0]] + molecule.Atoms[-1].connectivity
                        molecule.Atoms[-1].bond_type = ["1"] + molecule.Atoms[-1].bond_type
                        molecule.Atoms[pre_connectivity[0]].connectivity.append(atom.ID)
                        molecule.Atoms[pre_connectivity[0]].bond_type.append("1")
                
                pre_connectivity = [an + an_shift for an in connectivity]
                coordinates_type = hybrid(f"{molecule.Atoms[pre_connectivity[0]].element}{len(molecule.Atoms[pre_connectivity[0]].connect)+1}")

        return molecule

            

                    
                    

                
