"""
测试：名称带有%等特殊符号
    V3000 和 V2000
    带电荷的分子
    原子数目，键数目不对的文件

"""
import re
from urllib.parse import quote
from typing import List, Dict
from ...utils import  logger


class MolData:
    """.mol format"""

    def __init__(self, style="normal"):
        self.style = style

    def _parse(self, input_script : str,extra_var=None) -> Dict:
        script = input_script.splitlines()
        return self._parse_script(script)

    def _parse_v3000_script(self, input_script : str) -> Dict:
        atom_count = int(input_script[5].split()[3])
        bond_count = int(input_script[5].split()[4])
        datas = {
            "molecule_name"  :quote(input_script[0].strip()),
            "atom_count"     :atom_count,
            "bond_count"     :bond_count,
            "elements"       : [],
            "coordinates"    : [],
            "formal_charge"  : [0 for i in range(atom_count)],
            "connectivity"   : [[] for i in range(atom_count)],
            "bond_type"      : [[] for i in range(atom_count)],
            "associated_data": {},
            }
        
        #datas["formal_charge"] = [0 for i in range(datas["atom_count"])]
        for ii, line in enumerate(input_script[7:datas["atom_count"] + 7]):
            rr = line.strip().split()
            datas["elements"].append(rr[3])
            datas["coordinates"].append([float(rr[4]), float(rr[5]), float(rr[6])])
            for rrr in rr[6:]:
                if rrr[:3] == "CHG":
                    #datas["formal_charge"].append(int(rrr.split("=")[1].strip()))
                    datas["formal_charge"][ii] = int(rrr.split("=")[1].strip())
                    break

        for ii,rr in enumerate(input_script[datas["atom_count"] + 9 : datas["atom_count"] + 9 + datas["bond_count"]]):
            rrr = rr.split()
            atom1 = int(rrr[4]) - 1
            atom2 = int(rrr[5]) - 1
            bond_type = rrr[3]
            if bond_type == "4":
                bond_type = "ar"
            if atom2 not in datas["connectivity"][atom1]:
                datas["connectivity"][atom1].append(atom2)
                datas["bond_type"][atom1].append(bond_type)
            if atom1 not in datas["connectivity"][atom2]:
                datas["connectivity"][atom2].append(atom1)
                datas["bond_type"][atom2].append(bond_type)
        
        # regex for the pattern "> <>"
        regexp = re.compile(">\s+<(.*)>")
        for i in range(len(input_script)):
            result = regexp.match(input_script[i])
            if result:
                datas["associated_data"][result.group(1).lower()] = input_script[i + 1].strip()
        return datas

    def _parse_v2000_script(self, input_script : str) -> Dict:
        #atom_count = int(input_script[3].split()[0])
        atom_count = int(input_script[3][:3])
        bond_count = int(input_script[3][3:6])
        datas = {
            "molecule_name"  :quote(input_script[0].strip()),
            "atom_count"     :atom_count,
            "bond_count"     :bond_count,
            "elements"       : [],
            "coordinates"    : [],
            "formal_charge"  : [0 for i in range(atom_count)],
            "connectivity"   : [[] for i in range(atom_count)],
            "bond_type"      : [[] for i in range(atom_count)],
            "associated_data": {},
            }

        for ii, line in enumerate(input_script[4:datas["atom_count"] + 4]):
            rr = line.strip().split()
            datas["elements"].append(rr[3])
            datas["coordinates"].append([float(rr[0]), float(rr[1]), float(rr[2])])
        
        for ii,rr in enumerate(input_script[datas["atom_count"] + 4 : datas["atom_count"] + 4 + datas["bond_count"]]):
            atom1 = int(rr[:3]) - 1
            atom2 = int(rr[3:6]) - 1
            bond_type = rr[6:9].strip()
            if bond_type == "4":
                bond_type = "ar"
            if atom2 not in datas["connectivity"][atom1]:
                datas["connectivity"][atom1].append(atom2)
                datas["bond_type"][atom1].append(bond_type)
            if atom1 not in datas["connectivity"][atom2]:
                datas["connectivity"][atom2].append(atom1)
                datas["bond_type"][atom2].append(bond_type)
        
        if input_script[datas["atom_count"] + 4 + datas["bond_count"]].find("M  CHG") != -1:
            ss = input_script[datas["atom_count"] + 4 + datas["bond_count"]].strip().split()
            for ii in range(3,len(ss),2):
                datas["formal_charge"][int(ss[ii]) - 1] = int(ss[ii + 1])

        regexp = re.compile(">\s+<(.*)>")
        for i in range(len(input_script)):
            result = regexp.match(input_script[i])
            if result:
                datas["associated_data"][result.group(1).lower()] = input_script[i + 1].strip()
        return datas

    def _parse_script(self, script: str) -> Dict:
        if script is None or len(script) < 4:
            logger.warning("file script error: mol file is empty")
            return None
        script_type = script[3].strip().split()[-1]
        # # `script` is a list of lines from `splitlines()`. SDF may include empty blocks.
        # if len([line for line in script if line.strip()]) < 4:
        #     logger.warning("file script error: the script of this file is null")
        #     return None
        # count_line = script[3].strip().split()
        # if len(count_line) == 0:
        #     logger.warning("file script error: invalid counts line")
        #     return None
        # script_type = count_line[-1]
        if script_type == "V2000":
            return self._parse_v2000_script(script)
        elif script_type == "V3000":
            return self._parse_v3000_script(script)
        else:
            logger.warning(f"file version error: the version must V2000 or V3000, {script_type} is illegal")
            return None

    def _convert(self, molecule, extra_var=None) -> str:
        structure_3d = extra_var["structure_3d"] if extra_var is not None and "structure_3d" in extra_var \
            else True
        return self._convert_molecule(molecule,structure_3d)

    def _convert_molecule(self, molecule, structure_3d) -> str:
        """
        create the script of mol file format
        structure_3d = False, mean the coordinate are zero, usually for generate the smiles of a moleucle
        """
        
        molecule_name = getattr(molecule, "molecule_name")
        script = "%s\n%10s%12s\n\n" % (molecule_name, "DY", "3D") if structure_3d \
            else "%s\n%10s%12s\n\n" % (molecule_name, "DY", "2D")
        script += "%3d%3d  0  0  0  0  0  0  0  0999 V2000\n" % (molecule.atom_count, molecule.bond_count if hasattr(molecule,"Bonds") else 0)

        charges = []
        for atom in molecule.Atoms:
            element = "D" if hasattr(atom, "atom_type_name") and atom.atom_type_name == "_D" else atom.element
            x,y,z = atom.coordinates if structure_3d else [0.0000,0.0000,0.0000]
            script += "%10.4f%10.4f%10.4f%3s  0  0  0  0  0  0  0  0  0  0  0  0\n" % (x,y,z,element)
            if getattr(atom,"formal_charge",0) != 0:
                charges.append([atom.ID, atom.formal_charge])

        if hasattr(molecule,"Bonds") and molecule.bond_count > 0:
            for bond in molecule.Bonds:
                atom1 = bond.a1
                atom2 = bond.a2
                script += "%3d%3d%3s  0\n" % (atom1 + 1, atom2 + 1, 
                                            bond.bond_type if hasattr(bond,"bond_type") else bond.get_type(molecule)) 
        if len(charges) > 0:
            script += "M  CHG%3d" % len(charges)
            for cc in charges:
                script += "%4d%4d" % (cc[0] + 1, cc[1])
            script += "\n"
        script += "M  END\n"

        for k, v in getattr(molecule, "associated_data", {}).items():
            script += f"> <{k}>\n{v}\n\n"
        return script


class SdfData:
    """.mol or .sdf format"""

    def __init__(self, style="normal"):
        self.style = style

    def _parse(self,input_script: str,extra_var=None) -> Dict:
        return self._parse_script(input_script)
        
    def _parse_script(self, input_script : str) -> Dict:
        blocks = input_script.split("$$$$\n")
        datas_arr = []

        molf = MolData()
        for block in blocks:
            datas = molf._parse(block)
            if datas is not None:
                datas_arr.append(datas)
        return datas_arr

    def _convert(self,molecules,extra_var=None) -> str:
        if not isinstance(molecules,list):
            molecules = [molecules]
        structure_3d = extra_var["structure_3d"] if extra_var is not None and "structure_3d" in extra_var \
            else True
        return self._convert_molecules(molecules,structure_3d)

    def _convert_molecules(self,molecules,structure_3d) -> str:
        molf = MolData()
        scripts = []
        for molecule in molecules:
            script = molf._convert(molecule, extra_var={"structure_3d":structure_3d})
            if script is not None:
                scripts.append(script)

        return "$$$$\n".join(scripts) + "$$$$\n"
