import re

class PdbData:
    _name = "PDB"

    def __init__(self, style="pdb"):
        self.style = style

    def _parse(self, input_script,extra_var=None):
        script = input_script.splitlines()
        return self.read_coord(script)

    def read_file(self, input_script):
        self.read_coord(input_script)
        if self.style != "pdb":
            self.read_connect(input_script)

    def read_coord(self, input_script):
        datas = {
            "style": self.style,
            "molecule_name"  : "protein",
            "atom_count"     : 0,
            "elements"       : [],
            "coordinates"    : [],
            "formal_charge"  : [],
            "atom_name"      : [],
            "residue_ID"     : [],
            "residue"        : [],
            "chain_name"     : [],
            }

        # self.residu_subnumber=[]
        for line in input_script:
            if re.match("^TITLE", line):
                datas["molecule_name"] = line.split()[1].strip()

            # if re.match("^HETATM", line) or re.match("^ATOM", line):
            if re.match("^ATOM", line) and line[16:17] != "B":
                atom_coord = []
                atom_coord.append(float(line[30:38].strip()))
                atom_coord.append(float(line[38:46].strip()))
                if self.style == "charmm":
                    atom_coord.append(float(line[46:53].strip()))
                else:
                    atom_coord.append(float(line[46:54].strip()))
                datas["coordinates"].append(atom_coord)
                datas["elements"].append(line[76:78].strip())
                datas["residue"].append(line[17:20].strip())
                datas["atom_name"].append(line[12:16].strip())
                datas["residue_ID"].append(line[22:27].strip())
                # self.residu_subnumber.append(line[26:27].strip())
                datas["chain_name"].append(line[21:22].strip())
                if len(line) > 79 and line[79] == "+":
                    datas["formal_charge"].append(int(line[78].strip()))
                elif len(line) > 79 and line[79] == "-":
                    datas["formal_charge"].append(int(line[78].strip()) * -1)
                else:
                    datas["formal_charge"].append(0)
        datas["atom_count"] = len(datas["elements"])
        return datas

    def read_connect(self, input_script):
        self.connect_dict = {}
        for line in input_script:
            if re.match("^CONECT", line):
                atom_center = line[6:11].strip()
                if atom_center not in self.connect_dict.keys():
                    self.connect_dict[atom_center] = []
                for i in range(0, 4):
                    a = line[11 + i * 5 : 11 + (i + 1) * 5].strip()
                    if len(a) != 0:
                        self.connect_dict[atom_center].append(int(a))

    def _convert(self, molecule, extra_var=None) -> str:
        structure_3d = extra_var["structure_3d"] if extra_var is not None and "structure_3d" in extra_var \
            else True
        return self._convert_molecule(molecule)

    def _convert_molecule(self,molecule):
        self.import_moleobj(molecule)
        return self.generate_pdb_string(molecule.mole_name)

    def import_moleobj(self, m):
        """
        导入一个分子到pdb的对象中
        输入：m    ->   分子对象
        """
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
        self.atom_type = []
        __label_dict = {
            "coor": [self.coor, [0.000, 0.000, 0.000]],
            "elem": [self.elem, "zz"],
            "connect": [self.connect, []],
            "bond_type": [self.bond_type, []],
            "name": [self.name, "XX"],
            "residu": [self.residu, "UNK"],
            "residu_number": [self.residu_number, 1],
            "chain_name": [self.chain_name, "A"],
            "atom_type_name": [self.atom_type, "XX"],
            "ff_charge": [self.ff_charge, "0.000"],
            "formal_charge": [self.formal_charge, 0],
        }

        for aa in m.Atoms:
            for attr in __label_dict.keys():
                if hasattr(aa, attr):
                    __label_dict[attr][0].append(getattr(aa, attr))
                else:
                    __label_dict[attr][0].append(__label_dict[attr][1])
        # self.Coor_dict = moleobj.Coor_dict
        # self.Elem_dict = moleobj.Elem_dict
        # self.Connect_dict = moleobj.Connect_dict

    def generate_pdb_string(self, name, connect_flag=False, ff_charge_flag=False):
        content = f"COMPND    Molecule {name}\n"
        content += "AUTHOR    cpy\n"
        n = len(self.elem)
        for i in range(n):
            content += "ATOM  %5d %4s %3s %1s%4s    %8.3f%8.3f%8.3f  1.00  0.00    " % (
                i + 1,
                self.name[i],
                self.residu[i],
                self.chain_name[i],
                self.residu_number[i],
                self.coor[i][0],
                self.coor[i][1],
                self.coor[i][2],
            )
            if ff_charge_flag:
                content += "%6.3f %2s \n" % (self.ff_charge[i], self.atom_type[i])
            else:
                content += "      %2s" % self.elem[i]
                if self.formal_charge[i] != 0:
                    if self.formal_charge[i] > 0:
                        content += "%d+\n" % abs(self.formal_charge[i])
                    else:
                        content += "%d-\n" % abs(self.formal_charge[i])
                else:
                    content += "  \n"
        if connect_flag:
            for nn, vv in enumerate(self.connect):
                content += "CONECT%5d" % (nn + 1)
                for jj in vv:
                    content += "%5d" % (jj + 1)
                content += "\n"
        content += "END\n\n"
        return content

    def write_pdb_file(
        self,
        outputf,
        connect_flag=False,
        ff_charge_flag=False,
    ):
        """
        生成pdb格式的文件
        输入：outputf         ->    pdb文件的名称
              connect_flag    -> 是否生成连接性部分
              ff_charge_flag  -> 是否添加原子电荷。AutoDock软件pdbqt文件的专用
        """
        # with open(outputf, 'w') as outf:
        #    outf.write(
        #        self.generate_pdb_string(
        #            connect_flag=connect_flag, ff_charge_flag=ff_charge_flag
        #        )
        #    )
        outf = open(outputf, "w")
        outf.write("COMPND    Molecule for %s\n" % outputf.strip(".pdb"))
        outf.write("AUTHOR    cpy\n")
        n = len(self.elem)
        for i in range(n):
            outf.write(
                "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00    "
                % (
                    i + 1,
                    self.name[i],
                    self.residu[i],
                    self.chain_name[i],
                    self.residu_number[i],
                    self.coor[i][0],
                    self.coor[i][1],
                    self.coor[i][2],
                )
            )
            if ff_charge_flag:
                outf.write("%6.3f %2s \n" % (self.ff_charge[i], self.atom_type[i]))
            else:
                outf.write("      %2s" % self.elem[i])
                if self.formal_charge[i] != 0:
                    if self.formal_charge[i] > 0:
                        outf.write("%d+\n" % abs(self.formal_charge[i]))
                    else:
                        outf.write("%d-\n" % abs(self.formal_charge[i]))
                else:
                    outf.write("  \n")
        if connect_flag:
            for nn, vv in enumerate(self.connect):
                outf.write("CONECT%5d" % (nn + 1))
                for jj in vv:
                    outf.write("%5d" % (jj + 1))
                outf.write("\n")
        outf.write("END\n\n")
