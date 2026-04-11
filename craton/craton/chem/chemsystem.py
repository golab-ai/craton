import copy
import math
#from ..force_field import MolForceField as MFF #checkout_force_field

import numpy as np


class System:
    __Head_default = ["lx", "ly", "lz", "xy", "xz", "yz"]
    __defaultattrs = ["lattics", "molecules", "molecule_number", "coordinates", "ff", "md_para", "name"]
    __defterms = [
        "Atoms",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "Pair12",
        "Pair13",
        "Pair14",
        "Pair1n",
        "InterPair",
    ]
    __rest2attrs = [
        "rest2_scale_factor",
        "rest2_ligand_hot_atoms",
        "rest2_ligand_hot_torsions",
        "rest2_protein_hot_atoms",
        "rest2_protein_hot_torsions",
    ]

    def __init__(self, style=None):
        self.Coord = None
        self.style = style
        self.lattics = []
        self.molecules = []
        self.molecule_number = []
        self.coordinates = []
        self.ff = {}

    def set_info(self, dd):
        for aa, bb in dd.items():
            if aa in self.__defaultattrs + self.__rest2attrs:
                setattr(self, aa, bb)

    def assign_at(self, at_system):
        tmp = []
        for m in self.Molecules:
            m.get_at(at_system)
            for a in m.Atoms:
                tmp.append(a.atom_type_name)
        self.atom_types = {}
        for a in set(tmp):
            self.atom_types[a] = force_field.at_equ(a)
            self.atom_types[a].get_equ(at_system.equ_table)
        self.improper_types = set(tmp).intersection(at_system.improper)

    def constructe_topol(self):
        for m in self.Molecules:
            if m.hasattr(self.Bonds) is False:
                m.creat_topol()
            m.create_improper(self.improper)
            m.addATtoTerm(self.atom_types)

    # TODO Caofl
    def create_cross_term(self):
        pass

    # TODO Caofl
    def assign_hydrogen_bond(self):
        pass

    # TODO Caofl
    def add_virtual_point(self):
        pass

    # TODO Caofl
    def create_cg_bead(self):
        pass

    # TODO Caofl
    def grasp_para(self, at_system, ff_repo):
        ff = force_field.Force_Field(at_system, ff_repo)  # noqa
        # ff.grasp_ff(self.Molecules, terms) # noqa

    def create_topols(self):
        self.Atoms = []
        self.Bonds = []
        self.Angles = []
        self.Dihedrals = []
        self.Impropers = []
        self.tmp_Pair12 = []
        self.tmp_Pair13 = []
        self.tmp_Pair14 = []
        self.tmp_Pair1n = []
        __tt = {
            "Bonds": self.Bonds,
            "Angles": self.Angles,
            "Dihedrals": self.Dihedrals,
            "Impropers": self.Impropers,
            "Pair12": self.tmp_Pair12,
            "Pair13": self.tmp_Pair13,
            "Pair14": self.tmp_Pair14,
            "Pair1n": self.tmp_Pair1n,
        }
        mole_id = 0
        atom_id_start = 0
        for i in range(len(self.molecules)):
            mole_n = len(self.molecules[i].Atoms)
            for j in range(self.molecule_number[i]):
                this_atom_id_start = atom_id_start + mole_n * (j - 1)
                for k in range(mole_n):
                    aa = copy.deepcopy(self.molecules[i].Atoms[k])
                    aa.mole_id = mole_id
                    aa.mole_type = self.molecules[i].mole_name
                    aa.mole_xuhao = i
                    aa.coor = self.coor[this_atom_id_start + k]
                    for i in len(aa.connect):
                        aa.connect[i] = aa.connect[i] + this_atom_id_start
                    self.Atoms.append(aa)
                this_terms = [term for term in self.molecules[i].__dict__.keys() if term[0].isupper()]
                del this_terms[this_terms.index("Atoms")]
                for term in this_terms:
                    term_n = len(getattr(self.molecules[i], term))
                    for k in range(term_n):
                        tt = copy.deepcopy(getattr(self.molecules[i], term)[k])
                        tt.mole_id = mole_id
                        tt.mole_type = self.molecules[i].mole_name
                        tt.mole_xuhao = i
                        for ai in range(4):
                            a_term = "a%s" % ai
                            if hasattr(tt, a_term):
                                new_id = getattr(tt, a_term) + this_atom_id_start
                                setattr(tt, a_term, new_id)
                        __tt[term].append(tt)
                mole_id += 1
            atom_id_start = this_atom_id_start

    def create_atoms(self):
        self.Atoms = []
        mole_id = 0
        atom_id_start = 0
        for i in range(len(self.molecules)):
            mole_n = len(self.molecules[i].Atoms)
            for j in range(self.molecule_number[i]):
                this_atom_id_start = atom_id_start + mole_n * (j - 1)
                for k in range(mole_n):
                    aa = copy.deepcopy(self.molecules[i].Atoms[k])
                    aa.mole_id = mole_id
                    aa.mole_type = self.molecules[i].mole_name
                    aa.mole_xuhao = i
                    aa.No = this_atom_id_start + k
                    aa.No_offset = k
                    aa.coor = self.coor[this_atom_id_start + k]
                    for ii in range(len(aa.connect)):
                        aa.connect[ii] = aa.connect[ii] + this_atom_id_start
                    if not hasattr(aa, "para"):
                        aa.para = self.ff["atomtypes"][aa.atom_type_name]["para"]
                    aa.style = "LJ12_6"
                    self.Atoms.append(aa)
                mole_id += 1
            atom_id_start = this_atom_id_start

    def determine_pair(self, value, cutoff, pbc):
        if cutoff == "infinite":
            if pbc == "None":
                return True
            else:
                pass
        else:
            if pbc == "None":
                if value < cutoff:
                    return True
                else:
                    return False
            else:
                if value < cutoff:
                    return True
                else:
                    return False

    def create_intra_pairs(self, cutoff="infinite", pbc="all"):
        self.Pair12 = []
        self.Pair13 = []
        self.Pair14 = []
        self.Pair1n = []
        __tt = {
            "Pair12": [self.Pair12, self.tmp_Pair12],
            "Pair13": [self.Pair13, self.tmp_Pair13],
            "Pair14": [self.Pair14, self.tmp_Pair14],
            "Pair1n": [self.Pair1n, self.tmp_Pair1n],
        }
        for term in __tt.keys():
            tmp_arr = __tt[term][1]
            for rr in tmp_arr:
                rr.calc_value(self.Atoms[rr.a1].coor, self.Atoms[rr.a2].coor)
                flag = self.determine_pair(rr.value, cutoff, pbc)
                if flag:
                    __tt[term][0].append(rr)

    def create_inter_pairs(self, atom_arr1, atom_arr2, cutoff="infinite", pbc="all"):
        InterPairs = []
        for i in atom_arr1:
            for j in atom_arr2:
                this_pair = Pair("inter_pair", i, j)
                this_pair.calc_value([self.Atoms[i].coor, self.Atoms[j].coor])
                flag = self.determine_pair(this_pair.value, cutoff, pbc)
                if flag:
                    this_pair.charge_para = [self.Atoms[i].ff_charge, self.Atoms[j].ff_charge]
                    this_pair.style = self.Atoms[i].style
                    sigma = self.Atoms[i].para[0] * 0.5 + self.Atoms[j].para[0] * 0.5
                    espi = (self.Atoms[i].para[1] * self.Atoms[j].para[1]) ** 0.5
                    this_pair.para = [sigma, espi]
                    this_pair.mole_id_1 = self.Atoms[i].mole_id
                    this_pair.mole_type_1 = self.Atoms[i].mole_type
                    this_pair.mole_xuhao_1 = self.Atoms[i].mole_xuhao

                    this_pair.mole_id_2 = self.Atoms[j].mole_id
                    this_pair.mole_type_2 = self.Atoms[j].mole_type
                    this_pair.mole_xuhao_2 = self.Atoms[j].mole_xuhao
                    InterPairs.append(this_pair)
        return InterPairs

    def group_group_inter_energy(self, arr1, arr2, cutoff="infinite"):
        calculator = ff_calc.calculator("normal")
        # atom_arr1 = [self.Atoms[i] for i in arr1]
        # atom_arr2 = [self.Atoms[i] for i in arr2]
        InterPairs = self.create_inter_pairs(arr1, arr2, cutoff=cutoff)
        energy = calculator.inter_energy(InterPairs)
        return energy

    def mole_start(self, mn):
        atom_n = -1
        mole_n = -1
        flag = False
        for i in range(len(self.molecule_number)):
            for j in range(self.molecule_number[i]):
                mole_n += 1
                if mole_n == mn:
                    flag = True
                    break
                else:
                    atom_n += len(self.molecules[i].Atoms)
            if flag:
                break
        return atom_n, i

    @property
    def total_mole_number(self):
        return np.array(self.molecule_number).sum()

    @property
    def inter_energy(self):
        inter_energy = [0.0, 0.0]
        mn = -1
        for ii in range(len(self.mole)):
            for jj in range(self.mole_number[ii]):
                mn += 1
                start_n, __ = self.mole_start(mn)
                arr1 = [i for i in range(start_n, start_n + self.mole[ii].mole_n)]
                arr2 = [i for i in range(start_n + self.mole[ii].mole_n, len(self.Atoms))]
                energy = self.group_group_inter_energy(arr1, arr2)
            inter_energy[0] += energy[0]
            inter_energy[1] += energy[1]
        return inter_energy

    @property
    def intra_energy(self):
        calculator = ff_calc.calculator("normal")
        mn = -1
        intra_energy = []
        for ii in range(len(self.molecules)):
            for jj in range(self.molecule_number[ii]):
                mn += 1
                start_n, __ = self.mole_start(mn)
                for aa in range(start_n, start_n + self.molecules[ii].mole_n):
                    self.molecules[ii].Atoms[aa - start_n].coor = self.coor[aa]
                self.molecules[ii].update_topol_value()
                intra_energy.append(calculator.single_mole_energy(self.molecules[ii]))
        return intra_energy

    @property
    def intra_structure(self):
        pass

    def divide_by_topol(self, atoms, pp):
        __topol_level = [
            "mole_type",
            "mole_id",
            "chain_name",
            "residu",
            "residu_number",
            "ring",
            "name",
            "atom_type_name",
            "formal_charge",
        ]
        pp_level = []
        for p in __topol_level:
            if p in pp:
                pp_level.append(p)
        groups = {}
        for an in atoms:
            key_name = ""
            for p in pp_level:
                key_name += str(getattr(self.Atoms[an], p))
                key_name += "$"
            key_name = key_name[:-1]
            if key_name not in groups:
                groups[key_name] = []
            groups[key_name].append(an)
        return groups

    def get_scopes(self, scopes):
        """
        得到指定范围的原子编号
        """
        if scopes[0] in ["atom", "a"]:
            if ">" in scopes:
                return [i for i in range(int(scopes[1].strip(">")), len(self.Atoms))]
            elif "<" in scopes:
                return [i for i in range(0, int(scopes[1].strip("<")))]
            elif "_" in scopes:
                return [i for i in range(int(scopes[1].split("_")[0]), int(scopes[1].split("_")[1]))]
            else:
                return scopes[1:]
        elif scopes[0] in ["molecule_number", "mn"]:
            if ">" in scopes:
                start_n, mole_i = self.mole_start(int(scopes[1].strip(">")))
                return [i for i in range(start_n, len(self.Atoms))]
            elif "<" in scopes:
                start_n, mole_i = self.mole_start(int(scopes[1].strip("<")))
                start_n += len(self.molecules[mole_i].Atoms)
                return [i for i in range(0, start_n)]
            elif "_" in scopes:
                start_n, __ = self.mole_start(int(scopes[1].split("_")[0]))
                end_n, mole_i = self.mole_start(int(scopes[1].split("_")[1]))
                end_n += len(self.molecules[mole_i].Atoms)
                return [i for i in range(start_n, end_n)]
            else:
                arr = []
                for ii in scopes[1:]:
                    start_n, mole_i = self.mole_start(ii)
                    end_n = start_n + len(self.molecules[mole_i].Atoms)
                    arr.extend([i for i in range(start_n, end_n)])
                return arr
        elif scopes[0] in ["molecule_type", "mt"]:
            arr = []
            mole_start_n = 0
            for ii in range(len(self.molecules)):
                if self.molecules[ii].mole_name in scopes[1:]:
                    start_n, __ = self.mole_start(mole_start_n)
                    mole_start_n += self.molecule_number[ii]
                    end_n, mole_i = self.mole_start(mole_start_n - 1)
                    end_n += len(self.molecules[mole_i].Atoms)
                    arr.extend([i for i in range(start_n + 1, end_n + 1)])
                else:
                    mole_start_n += self.molecule_number[ii]
            return arr

    def get_groups(self, scopes, pp):
        """
        把指定的体系分割与满足要求的groups
        输入：scopes  ->  list 进行分割的范围。
                         scopes[0]是作用范围的标方式：atom or a: 按原子编号划定
                                                   molecule_number or mn:按分子编号进行划定
                                                   molecule_type or mt:按分子类型进行划定
                         scopes[1]，可以包含>,<,_，分别表示大于等于某一编号，小于等于某一编号，或在一个范围内
                         如果len(scopes) >= 2表示枚举
             pp      ->  list or string， group分割的条件
        """
        atoms = self.get_scopes(scopes)
        if isinstance(pp, str):
            pp = [pp]
        return self.divide_by_topol(atoms, pp)

    def get_energy(self):
        this_terms = [term for term in self.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]
        energy = 0.0
        for term in this_terms:
            for rr in getattr(self, term):
                energy += rr.energy
        self.energy = energy

    def energy_analysis(self, term=None, mole_type=None, mole_id=None):
        this_terms = [term for term in self.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]

        if term is None:
            if mole_type is None:
                if mole_id is None:
                    return self.get_energy
                else:
                    pass

    def transfer_lattic_matrix(self, arr):
        a = float(arr[0])
        b = float(arr[1])
        c = float(arr[2])
        alpha = math.radians(float(arr[3]))
        beta = math.radians(float(arr[4]))
        gamma = math.radians(float(arr[5]))
        n2 = (math.cos(alpha) - math.cos(gamma) * math.cos(beta)) / math.sin(gamma)
        n3 = (math.sin(beta) * math.sin(beta) - n2 * n2) ** 0.5
        Head = np.array(
            [[a, b * math.cos(gamma), c * math.cos(beta)], [0, b * math.sin(gamma), c * n2], [0, 0, c * n3]]
        )
        return Head

    def transfer_matrix_lattic(self, Head):
        a = Head[0]
        b = (Head[3] ** 2 + Head[1] ** 2) ** 0.5
        c = (Head[4] ** 2 + Head[5] ** 2 + Head[2] ** 2) ** 0.5
        alpha = math.degrees(math.acos((Head[3] * Head[4] + Head[1] * Head[5]) / b / c))
        beta = math.degrees(math.acos(Head[4] / c))
        gamma = math.degrees(math.acos(Head[3] / b))
        AA_t = [a, b, c, alpha, beta, gamma]
        return AA_t

    def transfer_cartn_fract(self, Head):
        BB = np.matrix([[Head[0], Head[3], Head[4]], [0.0, Head[1], Head[5]], [0.0, 0.0, Head[2]]]).I
        coord_f = []
        for term in self.Coord:
            coord_f.append(np.matmul(BB, term))
        return coord_f

    def transfer_fract_cartn(self, arr, Head):
        AA = np.array([[Head[0], Head[3], Head[4]], [0.0, Head[1], Head[5]], [0.0, 0.0, Head[2]]])
        for term in arr:
            self.Coord.append(np.matmul(AA, term))
        return self.Coord
