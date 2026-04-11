import re
from string import digits

from rdkit import Chem

from .mol_parse import MolData


class Mol2Data:
    """Represents a .mol2 file"""

    def __init__(self, style=""):
        pass

    def _parse(self,input_script,extra_var=None):
        
        return self._parse_mol2_to_mol(input_script)
        return {key:value for key,value in molobj.__dict__.items()}

    def _parse_mol2_to_mol(self,input_script):
        rdkm = Chem.MolFromMol2Block(input_script, removeHs=False)
        mol_script = Chem.MolToMolBlock(rdkm)
        molobj = MolData("normal")
        return molobj._parse(mol_script)

    def _convert(self,molecule,extra_var=None):
        return self._convert_molecule(molecule)

    def _convert_molecule(self, molecule):
        if not hasattr(molecule.Atoms[0],"atom_name"):
            molecule.get_null_atom_type()

        sss = "@<TRIPOS>MOLECULE\n"
        sss += "LIG\n"
        sss += "%d %d\n" % (len(molecule.Atoms), len(molecule.Bonds))
        sss += "SMALL\n"
        sss += "USER_CHARGES\n"
        sss += "@<TRIPOS>ATOM\n"
        for aa in molecule.Atoms:
            
            atom_type = f"{aa.elem}.{getattr(aa,'atom_type_name')}"
            if hasattr(aa, "ff_charge"):
                charge = aa.ff_charge
            else:
                charge = 0.0000
            sss += f"{0:>7} {1:<8}{2:>11.3f}{3:>11.3f}{4:>11.3f} {5:<5}{6:>10} {7:<9}{8:>10.6f}\n".format(
                    aa.ID + 1, aa.atom_name, aa.coor[0], aa.coor[1], aa.coor[2], atom_type, "1", "LIG", charge)
            #sss += f"{aa.ID + 1:>7} {aa.atom_name:<8}{aa.coor[0]:>11.3f}{aa.coor[1]:>11.3f}{aa.coor[2]:>11.3f} {atom_type:<5}{'1':>10} {'LIG':<9}{charge:>10.6f}\n"
        sss += "@<TRIPOS>BOND\n"
        n = 0
        using_bond_type = "bond_type_aromatic" if hasattr(molecule.Atoms[0],"bond_type_aromatic") else "bond_type"
        for bb in molecule.Bonds:
            n += 1
            
            sss += "        %d  %d  %d  %s\n" % (
                n,
                bb.a1 + 1,
                bb.a2 + 1,
                getattr(molecule.Atoms[bb.a1],using_bond_type)[molecule.Atoms[bb.a1].connect.index(bb.a2)],
            )
        sss += "@<TRIPOS>SUBSTRUCTURE\n"
        sss += "       1 LIG             1 RESIDUE    0 **** ROOT      0\n"
        return sss


    ####下面的测试后删除
    def read_file(self, input_script):
        _tt = [[], []]
        for i in range(0, len(input_script)):
            if re.match("^@", input_script[i]):
                _tt[0].append(input_script[i].strip())
                _tt[1].append(i)
        _tt[1].append(i + 1)

        for i in range(0, len(_tt[0])):
            if _tt[0][i] == "@<TRIPOS>ATOM":
                self.read_coor(input_script, _tt[1][i], _tt[1][i + 1])
            elif _tt[0][i] == "@<TRIPOS>BOND":
                self.read_connect(input_script, _tt[1][i], _tt[1][i + 1])
        self.script = "".join(input_script)

    def read_coor(self, input_script, a, b):
        self.elem = []
        self.coor = []
        self.ff_charge = []
        self.atom_type_name = []
        self.name = []
        for i in range(a + 1, b):
            ss = input_script[i].strip().split()
            # self.elem.append(ss[5].split('.')[0])
            remove_digits = str.maketrans("", "", digits)
            self.elem.append(ss[1].translate(remove_digits))
            if ss[5].find(".") != -1:
                self.atom_type_name.append(ss[5].split(".")[1])
                # self.name.append(ss[5].split('.')[0])
            else:
                self.atom_type_name.append(ss[5].strip())
            self.name.append(ss[1].strip())
            self.ff_charge.append(float(ss[8]))
            self.coor.append([float(ss[2].split("#")[0]), float(ss[3].split("#")[0]), float(ss[4].split("#")[0])])
        self.mole_n = len(self.elem)

    def read_connect(self, input_script, a, b):
        arr = []
        for i in range(a + 1, b):
            ss = input_script[i].strip().split()
            arr.append([int(ss[1]) - 1, int(ss[2]) - 1, ss[3]])
        self.connect = [[] for i in range(self.mole_n)]
        self.bond_type_aromatic = [[] for i in range(self.mole_n)]
        for rr in arr:
            if rr[1] not in self.connect[rr[0]]:
                self.connect[rr[0]].append(rr[1])
                self.bond_type_aromatic[rr[0]].append(rr[2])
            if rr[0] not in self.connect[rr[1]]:
                self.connect[rr[1]].append(rr[0])
                self.bond_type_aromatic[rr[1]].append(rr[2])

    def import_moleobj_old(self, m, has3d="no", vsflag=False):
        sss = "@<TRIPOS>MOLECULE\n"
        sss += "Molecule Name\n"
        if vsflag:
            sss += "%d %d\n" % (len(m.Atoms), len(m.Bonds))
        else:
            if hasattr(m, "Vss"):
                sss += "%d %d\n" % (len(m.Atoms) - len(m.Vss), len(m.Bonds))
            else:
                sss += "%d %d\n" % (len(m.Atoms), len(m.Bonds))
        sss += "SMALL\n"
        sss += "NO_CHARGES\n\n"
        sss += "@<TRIPOS>ATOM\n"
        if vsflag:
            for aa in m.Atoms:
                if has3d == "yes":
                    sss += "%d %s%d %15.4f %15.4f %15.4f" % (
                        aa.No + 1,
                        aa.elem,
                        aa.No,
                        aa.coor[0],
                        aa.coor[1],
                        aa.coor[2],
                    )
                else:
                    sss += "%d %s%d %15.4f %15.4f %15.4f" % (aa.No + 1, aa.elem, aa.No, 0.0000, 0.0000, 0.0000)
                if aa.atom_type_name is not None:
                    sss += " %10s.%s\n" % (aa.elem, aa.atom_type_name)
                else:
                    sss += " %10s\n" % (aa.elem)
        else:
            for aa in m.Atoms:
                if aa.elem != "EP":
                    if has3d == "yes":
                        sss += "%d %s%d %15.4f %15.4f %15.4f" % (
                            aa.No + 1,
                            aa.elem,
                            aa.No,
                            aa.coor[0],
                            aa.coor[1],
                            aa.coor[2],
                        )
                    else:
                        sss += "%d %s%d %15.4f %15.4f %15.4f" % (aa.No + 1, aa.elem, aa.No, 0.0000, 0.0000, 0.0000)
                    if aa.atom_type_name is not None:
                        sss += " %10s.%s\n" % (aa.elem, aa.atom_type_name)
                    else:
                        sss += " %10s\n" % (aa.elem)
        sss += "@<TRIPOS>BOND\n"
        n = 0
        for bb in m.Bonds:
            n += 1
            sss += "%d  %d  %d  %s\n" % (
                n,
                bb.a1 + 1,
                bb.a2 + 1,
                m.Atoms[bb.a1].bond_type_aromatic[m.Atoms[bb.a1].connect.index(bb.a2)],
            )
        self.script = sss


