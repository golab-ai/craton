import math
import time

import numpy as np



class CifData:
    """.cif file"""

    def __init__(self, style):
        self.s = style

        __Head_default = [  # noqa
            "_symmetry_space_group_name",
            "_sysmmetry_Int_Tables_number",
            "_symmetry_cell_setting",
            "_cell_length_a",
            "_cell_length_b",
            "_cell_length_c",
            "_cell_angle_alpha",
            "_cell_angle_beta",
            "_cell_angle_gama",
        ]
        self.Head_arr = []
        self.Atom_Site = []
        self.Geom_Bond = []

    def cif_extra_loop(self, n, script):
        Arr = [[]]
        exit_flag = 0
        for i in range(n + 1, len(script)):
            if script[i].find("loop_") != -1:
                break
            elif script[i][0] == "_":
                if exit_flag == 1:
                    break
                string = script[i].strip().split()
                for j in range(0, len(string)):
                    Arr[0].append(string[j])
            else:
                Arr.append([s for s in script[i].strip().split() if len(script[i].strip().split()) > 0])
                exit_flag = 1
        return Arr

    def read_info(self, input_script):
        for term in self.__Head_default:
            if term not in input_script:
                self.Head_arr.append("None")
            else:
                self.Head_arr.append([line for line in input_script if line.find(term) != -1][0].strip().split()[1])
        section_arr = [i for i in range(0, len(input_script)) if "loop_" in input_script[1]]
        n = [i for i in section_arr if "_atom_site" in input_script[i + 1]]
        self.Atom_Site = self.cif_extra_loop(n[0], input_script)
        if "connect" in self.s:
            m = [i for i in section_arr if "_geom_bond" in input_script[i + 1]]
            self.Geom_Bond = self.cif_extra_loop(m[0], input_script)

    def write_cif(self):
        out_script = "data_from_cpy\n"
        out_script += "_audit_creation_date %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
        out_script += "_audit_creation_method %s\n" % '"cpy"'
        for i in range(0, len(self.__Head_default)):
            if self.Head_arr[i] != "None":
                out_script += "%s %.8f\n" % (self.__Head_default[i], self.Head_arr[i])
        out_script += "\nloop_\n_space_group_symop_operation_xyz\nx,y,z\n"
        out_script += "loop_\n"
        for term in self.Atom_Site[0]:
            out_script += "%s\n" % term
        for record in self.Atom_Site[1:]:
            out_script += "   ".join(record)
            out_script += "\n"
        if len(self.Geom_Bond) > 0:
            for term in self.Geom_Bond[0]:
                out_script += "%s\n" % term
            for record in self.Geom_Bond[1:]:
                out_script += "   ".join(record)
                out_script += "\n"

    def transfer_matrix(self):
        a = float(self.Head_arr[3])
        b = float(self.Head_arr[4])
        c = float(self.Head_arr[5])
        alpha = math.radians(float(self.Head_arr[6]))
        beta = math.radians(float(self.Head_arr[7]))
        gamma = math.radians(float(self.Head_arr[8]))
        n2 = (math.cos(alpha) - math.cos(gamma) * math.cos(beta)) / math.sin(gamma)
        n3 = (math.sin(beta) * math.sin(beta) - n2 * n2) ** 0.5
        AA = np.array([[a, b * math.cos(gamma), c * math.cos(beta)], [0, b * math.sin(gamma), c * n2], [0, 0, c * n3]])
        return AA

    def get_connectivity(self):
        connect_dict = {}
        arr = [-1] + [term[0] for term in self.Atom_Site]
        for record in self.Geom_Bond:
            if str(arr.index(record[0])) not in connect_dict.keys():
                connect_dict[str(arr.index(self.recrd[0]))] = [arr.index(record[1])]
            else:
                if arr.index(record[1]) not in connect_dict[str(arr.index(record[0]))]:
                    connect_dict[str(arr.index(self.recrd[0]))].append(arr.index(record[1]))
            if str(arr.index(record[1])) not in connect_dict.keys():
                connect_dict[str(arr.index(self.recrd[1]))] = [arr.index(record[0])]
            else:
                if arr.index(record[0]) not in connect_dict[str(arr.index(record[1]))]:
                    connect_dict[str(arr.index(self.recrd[1]))].append(arr.index(record[0]))
        return connect_dict

    def make_atoms(self):
        atoms = {}
        if "connect" in self.s:
            Connect_dict = self.get_connectivity()
        AA = self.transfer_matrix()
        for i in range(1, len(self.Atom_Site)):
            atoms[str(i)] = Atom("atomic")
            atoms[str(i)].number = i
            atoms[str(i)].elem = self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_type_symbol")]
            atoms[str(i)].mass = get_elem_property("elem", "mass", atoms[str(i)].elem)
            atoms[str(i)].atom_number = get_elem_property("elem", "number", atoms[str(i)].elem)
            if "connect" in self.s:
                atoms[str(i)].connectivity = Connect_dict[str(i)]
            atoms[str(i)].name = self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_label")]
            if "fract" in self.s:
                coor_f = [
                    self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_fract_x")],
                    self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_fract_y")],
                    self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_fract_z")],
                ]
                atoms[str(i)].x, atoms[str(i)].y, atoms[str(i)].z = np.matmul(AA, coor_f)
            else:
                atoms[str(i)].x = self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_Cartn_x")]
                atoms[str(i)].y = self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_Cartn_y")]
                atoms[str(i)].z = self.Atom_Site[i][self.Atom_Site[0].index("_atom_site_Cartn_z")]
        self.atoms = atoms

    def creat_Cif(self, system):
        arr = system.transfer_matrix_lattic()
        self.Head_arr = ["'P1'", "1", "None", arr[0], arr[1], arr[2], arr[3], arr[4], arr[5]]
        Atom_Site = [
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
            "_atom_site_occupancy",
        ]
        coord_f = system.transfer_cartn_fract()
        n = -1
        for record in system.System:
            for i in range(0, record[1]):
                for j in range(1, len(range(0, system.Molecules[record[0]]).Atoms) + 1):
                    n += 1
                    name = system.Molecules[record[0]].Atoms[str(j)].name + "_" + str(i)
                    elem = system.Molecules[record[0]].Atoms[str(j)].elem
                    Atom_Site.append([name, elem, coord_f[n][0], coord_f[n][1], coord_f[n][2], "1.0000"])
