import re
from collections import defaultdict

import numpy as np


class FlexibleBond:
    """
    Object storing flexible torsions for pdbqt string generation

    Attributes:
        _flexible_bonds (list): [[atom_index1(int), atom_index2(int)]], list of
            pairs of indices of atoms forming flexible bonds
    """

    def __init__(self):
        self._flexible_bonds = []

    @property
    def flexible_atom(self):
        """Return unique indices of atoms forming flexible bonds"""
        if self._flexible_bonds == []:
            return set()
        return set.union(*self._flexible_bonds)

    @property
    def flexible_bonds(self):
        return self._flexible_bonds

    @property
    def number_of_bonds(self):
        """Return number of flexible bonds"""
        return len(self._flexible_bonds)

    def add_felxible_bond(self, atom_index_1, atom_index_2):
        """Add pair of atom indices to _flexible_bonds"""
        self._flexible_bonds.append({atom_index_1, atom_index_2})

    def pop(self, atom_index_1):
        """
        Pop indices of all the atoms connect to atom_index_1 by a flexible bond
        """
        indice = set()
        new_flexible_bonds = []
        for pair in self._flexible_bonds:
            if atom_index_1 not in pair:
                new_flexible_bonds.append(pair)
                continue
            pair.remove(atom_index_1)
            indice.add(pair.pop())
        self._flexible_bonds = new_flexible_bonds
        return indice


class PdbqtData:
    """Object converting CFL.molecule to pdbqt string block"""

    _name = "PDBQT"

    def __init__(self, style="normal"):
        self.style = style

    def convert_info(self,molecule,extra_var=None):
        text = self._convert(molecule,extra_var=extra_var)
        return text, self.atom_index_mapping

    def _convert(self,molecule,extra_var=None):
        self.import_moleobj(molecule)
        if extra_var is not None:
            self.igonre_connect = extra_var["ignore_connect"]
        else:
            self.igonre_connect = False
        return self.write_pdbqt_file()
    

    def import_moleobj(self, m):
        """
        Import CFL.molecule object to pdbqt
        Attributes:
            mole_name (string): name of molecule
            coor (list): [[x(float), y(float), z(float)]], coordinates of atoms
            elem (list): [string], element of atoms
            connect (list): [[atom_index(int)]], connected atom indecs of each
                atom
            bond_type (list): [[string]], bond types of bonds of each atom, the
                order should be the same as connect
            name (list): [string], atom type of atoms in PDB format
            residu (list): [string], name of residue to which atoms belong
            residu_number (list): [int], index of residue to which atoms belong
            chain_name (list): [string], label of chain to which atoms belong
            atom_type_name (list): [string], atom type of atoms in force field
            ff_charge (list): [float], point charge of atoms of atoms in force
                field
            formal_charge (list): [int], formal charge of atoms
            scan_term (list): [[atom_1(int), atom_2(int), atom_3(int),
                atom_4(int)]], list of indices of atoms in flexible dihedrals.
                The middle two atoms are atoms forming flexible bonds

        Raises:
            AttributeError: scan_term not available in CFL.molecule object.
                Currently only scan_term is required even for rigid molecule.
        """
        self.mole_style = m.style
        if hasattr(m, "mole_name"):
            self.mole_name = m.mole_name
        if m.style not in ["pdb","protein","template","dna","rna","DNA","RNA","Protein"]:
            self.smiles = m.smiles
        self.coor = []  # ->   每个原子的坐标
        self.elem = []  # ->
        self.connect = []  # ->
        self.bond_type = []  # ->
        self.name = []  # ->
        self.residu = []
        self.chain_name = []
        self.residu_number = []
        self.ff_charge = []
        self.formal_charge = []
        self.atom_type_name = []
        __label_dict = {
            "coor": [self.coor, [0.000, 0.000, 0.000]],
            "elem": [self.elem, "zz"],
            "connect": [self.connect, []],
            "bond_type": [self.bond_type, []],
            "name": [self.name, "XX"],
            "residu": [self.residu, "UNK"],
            "residu_number": [self.residu_number, 1],
            "chain_name": [self.chain_name, "A"],
            "atom_type_name": [self.atom_type_name, "XX"],
            "ff_charge": [self.ff_charge, "0.000"],
            "formal_charge": [self.formal_charge, 0],
        }

        for aa in m.Atoms:
            for attr in __label_dict.keys():
                if hasattr(aa, attr):
                    if attr == "residu_number":
                        if len(getattr(aa,attr)) > 4:
                            __label_dict[attr][0].append(getattr(aa, attr)[:4])
                        else:
                            __label_dict[attr][0].append(getattr(aa, attr))
                    else:
                        __label_dict[attr][0].append(getattr(aa, attr))
                else:
                    __label_dict[attr][0].append(__label_dict[attr][1])

        if hasattr(m, "scan_term"):
            self.scan_term = m.scan_term
        else:
            raise AttributeError(
                "Attribute scan_term missed. Currently the class is designed "
                "for flexible structures. Even for rigid molecule scan_term is "
                "also required. If non-flexible pdbqt required, please use Pdb "
                "object with ff_charge=True"
            )

    def dump_atom(self, connect_atom):
        """Dump the atom line of pdbqt for each atom"""
        self.dumped.append(connect_atom)
        self.atom_index += 1
        
        self.atom_index_mapping[connect_atom] = self.atom_index

        return ("ATOM  %5d %4s %3s %1s%4s    %8.3f%8.3f%8.3f  " "1.00  0.00    %6.3f %2s \n") % (
            self.atom_index,
            self.name[connect_atom],
            self.residu[connect_atom],
            self.chain_name[connect_atom],
            str(self.residu_number[connect_atom]),
            self.coor[connect_atom][0],
            self.coor[connect_atom][1],
            self.coor[connect_atom][2],
            self.ff_charge[connect_atom],
            self.atom_type_name[connect_atom],
        )

    def get_bond_flexible(self):
        self.flexible_bonds = FlexibleBond()
        self.branch_connection = defaultdict(set)

        for dihedral in self.scan_term:
            self.flexible_bonds.add_felxible_bond(dihedral[1], dihedral[2])
        for atom in self.flexible_bonds.flexible_atom:
            self.branch_connection[atom].add(atom)
            self.recursive_find_branch(atom, atom)

    def recursive_find_branch(self, center_atom, initial_atom):
        """Identify rigid branches by flexible bonds"""
        for connect_atom in self.connect[center_atom]:
            if {connect_atom, center_atom} in self.flexible_bonds.flexible_bonds:
                continue
            if connect_atom not in self.branch_connection[initial_atom]:
                self.branch_connection[initial_atom].add(connect_atom)
                self.recursive_find_branch(connect_atom, initial_atom)

    def recursive_write_flexible_branch(self, center_atom, content):
        """Iterate though branches and dump atoms by branches"""
        for anchor_atom in self.branch_connection[center_atom]:
            if anchor_atom not in self.flexible_bonds.flexible_atom:
                continue
            for branch_atom in self.flexible_bonds.pop(anchor_atom):
                if branch_atom in self.dumped:
                    continue
                atom_index_tmp = self.atom_index + 1
                content += "BRANCH%4d%4d\n" % (self.atom_index_mapping[anchor_atom], atom_index_tmp)
                content += self.dump_atom(branch_atom)
                for connect_atom in self.branch_connection[branch_atom]:
                    if connect_atom == branch_atom or connect_atom in self.dumped:
                        continue
                    content += self.dump_atom(connect_atom)
                content = self.recursive_write_flexible_branch(branch_atom, content)
                content += "ENDBRANCH%4d%4d\n" % (
                    self.atom_index_mapping[anchor_atom],
                    atom_index_tmp,
                )
        return content

    def generate_pdbqt_string(self, mole_name=None):
        """
        Generate pdbqt string content from a CFL.molecule object. The order of
        atom would be different from the original order of atoms in
        CFL.molecule

        Args:
            mole_name (string): name of molecule. Would try to get
                molecule.mole_name if not specified

        Attributes:
            atom_index (int): current dumping index of atom in pdbqt
            atom_index_mapping (dict): {int: int}, mapping of original index of
                atom to new index in pdbqt
            flexible_bonds (FlexibleBond): object containing flexible torsions
                for pdbqt string generaton
            branch_connection (dict): {center_atom_index(int):
                [branch_atom_index(int)]}, contain indices of branch atoms of
                each rigid branch connected to other branches via
                center_atom_index
        """
        self.get_bond_flexible()
        self.atom_index_mapping = {}
        self.atom_index = 0
        number_of_bonds = self.flexible_bonds.number_of_bonds
        self.dumped = []

        if mole_name is None:
            mole_name = self.mole_name
        #content = f"TITLE    {mole_name}\n"
        #content += "AUTHOR    cpy\n"
        content = ""
        if self.mole_style not in ["pdb","protein","template","dna","rna","DNA","RNA","Protein"]:
            content += f"REMARK SMILES {self.smiles}\n"

        if number_of_bonds == 0:
            #content += "ROOT\n"
            for atom_tmp in range(len(self.elem)):
                content += self.dump_atom(atom_tmp)
            #content += "ENDROOT\n"
            #content += "TORSDOF 0\n"
        else:
            content += "REMARK%3d active torsions:\n" % number_of_bonds
            content += "REMARK  status: ('A' for Active; 'I' for Inactive)\n"
            for bond_index, atoms in enumerate(self.flexible_bonds.flexible_bonds, start=1):
                atom_1, atom_2 = atoms
                content += "REMARK%5d  A    between atoms: %s_%d  and  %s_%d\n" % (
                    bond_index,
                    self.name[atom_1],
                    atom_1,
                    self.name[atom_2],
                    atom_2,
                )
            content += "ROOT\n"
            for center_atom, connect_atom_list in self.branch_connection.items():
                if 0 not in connect_atom_list:
                    continue
                for connect_atom in connect_atom_list:
                    content += self.dump_atom(connect_atom)
                break
            content += "ENDROOT\n"
            content = self.recursive_write_flexible_branch(center_atom, content)
            content += "TORSDOF%2d\n" % number_of_bonds

        if not self.igonre_connect:
            for center_atom, connect_atom_list in enumerate(self.connect):
                content += "CONECT%5d" % self.atom_index_mapping[center_atom]
                for connect_atom in connect_atom_list:
                    content += "%5d" % self.atom_index_mapping[connect_atom]
                content += "\n"
            content += "END\n\n"
        return content

    def write_pdbqt_file(self, mole_name=None):
        """Write pdbqt string to file"""
        #if mole_name is None:
        #    mole_name = outputf.strip(".pdbqt")
        #with open(outputf, "w") as outf:
        #    outf.write(self.generate_pdbqt_string(mole_name))
        return self.generate_pdbqt_string(mole_name)

    def split_connect_info_of_pdbqt(self, mole_name=None):
        """Split connect block from pdbqt for autodock vina"""
        pdbqt_block = self.generate_pdbqt_string(mole_name)
        pdbqt_heads = ["REMARK", "ROOT", "ATOM", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"]
        new_pdbqt = []
        connect_info = []
        for line in pdbqt_block.splitlines():
            if np.any([re.match(i, line) for i in pdbqt_heads]):
                new_pdbqt.append(line)
            elif re.match("CONNECT", line):
                connect_info.append(line)
        return "\n".join(new_pdbqt), "\n".join(connect_info)

    def write_part_of_pdbqt_file(self, outputf, mole_name=None):
        """Write pdbqt string to file"""
        if mole_name is None:
            mole_name = outputf.strip(".pdbqt")

        new_pdbqt, connect_info = self.split_connect_info_of_pdbqt(mole_name)
        with open(outputf, "w") as outf:
            outf.write(new_pdbqt)
        return connect_info
