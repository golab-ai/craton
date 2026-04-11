import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from shlex import split
from typing import Dict, List



@dataclass
class Atom:
    x: float = None
    y: float = None
    z: float = None
    atomic_number: int = None
    formal_charge: float = 0


@dataclass
class Bond:
    a1: int = None
    a2: int = None
    order: int = None


@dataclass
class MaeMol:
    property: Dict[str, str] = None
    atoms: List[Atom] = field(default_factory=list)
    bonds: List[Bond] = field(default_factory=list)


class MaeFile:
    def __init__(self, input):
        self.input = input
        self.mae_mols = []
        self.fmct_block = []

    def open_mae(self, filename):
        func_dict = {
            ".mae": open(filename, "r"),
            ".maegz": gzip.open(filename, "rt"),
        }
        return func_dict[Path(filename).suffix]

    def read_file(self):
        with self.open_mae(self.input) as f:
            lines = f.read().splitlines()
        idx_fmct = []
        for i, line in enumerate(lines):
            if line.strip().startswith("f_m_ct"):
                idx_fmct.append(i)
        idx_fmct.append(len(lines))
        for i in range(len(idx_fmct) - 1):
            idx_start = idx_fmct[i]
            idx_end = idx_fmct[i + 1] - 1
            self.fmct_block.append(lines[idx_start:idx_end])
        self.read_fmct()

    def read_fmct(self):
        for block in self.fmct_block:
            mae = MaeMol()
            for i, line in enumerate(block):
                if line.strip().startswith("m_"):
                    break
            self.read_property(mae, block[1:i])
            m_idx = []
            m_block = block[i:]
            for i, line in enumerate(m_block):
                if line.strip().startswith("m_"):
                    m_idx.append(i)
            # m_idx.append(len(m_block) - 2)
            m_idx.append(len(m_block))
            for i in range(len(m_idx) - 1):
                if m_block[m_idx[i]].strip().startswith("m_atom"):
                    self.read_m_atoms(mae, m_block[m_idx[i] : m_idx[i + 1]])
                if m_block[m_idx[i]].strip().startswith("m_bond"):
                    self.read_m_bonds(mae, m_block[m_idx[i] : m_idx[i + 1]])
            self.mae_mols.append(mae)

    def read_property(self, mae, lines):
        keys = []
        values = []
        key_fininsh_flag = False

        for line in lines:
            if line.strip().startswith(":::"):
                key_fininsh_flag = True
                continue
            if line.strip().startswith("#") or line.strip() == "":
                continue
            if not key_fininsh_flag:
                keys.append(line.strip())
            else:
                values.append(line.strip())
        mae.property = {key: value for key, value in zip(keys, values)}

    def read_m_atoms(self, mae, lines):
        key_finish_flag = False
        keys = []
        values = []
        for line in lines[1:-1]:
            if line.strip().startswith("#") or line.strip() == "":
                continue
            if line.strip().startswith(":::"):
                key_finish_flag = True
                continue
            if not key_finish_flag:
                keys.append(line.strip())
            else:
                values.append(line.strip())
        formal_charge_idx = None
        for i, key in enumerate(keys):
            if "m_x_coord" in key:
                x_idx = i
            elif "m_y_coord" in key:
                y_idx = i
            elif "m_z_coord" in key:
                z_idx = i
            elif "atomic_number" in key:
                atomic_number_idx = i
            elif "formal_charge" in key:
                formal_charge_idx = i
        for value in values:
            atom = Atom()
            token = split(value)[1:]
            atom.x = float(token[x_idx])
            atom.y = float(token[y_idx])
            atom.z = float(token[z_idx])
            atom.atomic_number = int(token[atomic_number_idx])
            if formal_charge_idx is not None:
                atom.formal_charge = int(token[formal_charge_idx])
            mae.atoms.append(atom)

    def read_m_bonds(self, mae, lines):
        key_finish_flag = False
        keys = []
        values = []
        for line in lines[1:-1]:
            if line.strip().startswith("#") or line.strip() == "":
                continue
            if line.strip().startswith(":::"):
                key_finish_flag = True
                continue
            if not key_finish_flag:
                keys.append(line.strip())
            else:
                values.append(line.strip())

        for i, key in enumerate(keys):
            if key.strip() == "i_m_from":
                from_idx = i
            elif key.strip() == "i_m_to":
                to_idx = i
            elif key.strip() == "i_m_order":
                order_idx = i

        for value in values:
            bond = Bond()
            token = split(value)[1:]
            bond.a1 = int(token[from_idx])
            bond.a2 = int(token[to_idx])
            bond.order = int(token[order_idx])
            mae.bonds.append(bond)

    def export_molobj(self):
        inf_mols = []
        num_to_elem = partial(get_elem_property, "number", "elem")
        for mol in self.mae_mols:
            bond_dict = defaultdict(list)
            inf_mol_dict = defaultdict(list)
            for bond in mol.bonds:
                bond_dict[bond.a1 - 1].append((bond.a2 - 1, bond.order))
                bond_dict[bond.a2 - 1].append((bond.a1 - 1, bond.order))
            for i, atom in enumerate(mol.atoms):
                bond_to_list = []
                type_list = []
                for bond_to, type in bond_dict[i]:
                    bond_to_list.append(bond_to)
                    type_list.append(str(type))
                inf_mol_dict["connect"].append(bond_to_list)
                inf_mol_dict["bond_type"].append(type_list)
                inf_mol_dict["elem"].append(num_to_elem(atom.atomic_number))
                inf_mol_dict["coor"].append([atom.x, atom.y, atom.z])
                if atom.formal_charge != 0:
                    inf_mol_dict["formal_charge"].append([i, atom.formal_charge])
            m = Molecule("normal")
            m = get_mole_info(m,self.__dict__)
            m = get_atoms_info(m,self.__dict__)
            #m.get_mole_info(inf_mol_dict)
            #m.create_atoms(len(mol.atoms))
            #m.get_atoms_info(inf_mol_dict)
            m.create_topols()
            m.mole_name = mol.property.get("s_m_title", "UNK")
            inf_mols.append(m)
        return inf_mols
