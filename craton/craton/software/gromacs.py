import math
import os
from copy import deepcopy

import numpy as np

from ..utils import logger
from ..chem.atom import Atom
from ..chem.topology import Angle, Bond, Dihedral, Improper, Pair
from ..chem.molecule import Molecule
from ..chem.chemsystem import System
from ..chem.ensemble import Ensemble
from ..chem import constants


from .read_file_utils import extra_info, extra_section, search_line

from ..utils.common.utils import add_two_dict
_XTC_IMPORTED = False
#try:
    #from compuchem.chemistry.software.gmx_traj_parser.libxdr import XTCFile

#    _XTC_IMPORTED = True
#except ImportError:
#    _XTC_IMPORTED = False
#    logger.warning("Cython extension for XTCFile not found")


current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)


def format_str(element):
    """
    Always add a leading space in case a string is too long and kiss the butt of previous string
    """
    if isinstance(element, str):
        element = element.strip()
        if len(element) <= 2 and not element.startswith("_"):
            return " %2s" % element
        return " %12s" % element
    elif isinstance(element, (int, np.integer)):
        return " %6i" % element
    elif isinstance(element, (float, np.floating)):
        return " %14.6f" % element
    else:
        raise Exception("Invalid type for GROMACS output: %s %s" % (element, type(element)))


def coul_vdw_value(params, rr=3.5):
    sigma = params[3] / 2.0 + params[5] / 2.0
    espi = (params[4] * params[6]) ** 0.5
    coul = params[2] * 138.935458 * params[0] * params[1] / rr
    vdw = params[7] * 4 * espi * ((sigma / rr) ** 12 - (sigma / rr) ** 6)
    return coul + vdw


class TableFunction:
    def __init__(self, style, params, interval=0.002, steps=1500, fname="table.xvg"):
        self.rr = [interval * n for n in range(steps + 1)]
        self.style = style
        params.append(self.rr)
        self.params = params
        self.fname = fname
        arr = self.__Func_ID[self.style](self, *self.params)
        with open(self.fname, "w") as outf:
            for rr in arr:
                outf.write("%.10e %.10e %.10e\n" % (rr[0], rr[1], rr[2]))

    def LJ12_6(self, sigma1, espi1, sigma2, espi2, factor, rr):
        sigma = sigma1 / 2.0 + sigma2 / 2.0
        espi = (espi1 * espi2) ** 0.5
        arr = []
        for r in rr[:10]:
            arr.append([r, 0, 0])
        for r in rr[10:]:
            energy = factor * 4 * espi * ((sigma / r) ** 12 - (sigma / r) ** 6)
            force = factor * 4 * espi * (12 * sigma**12 / r**13 - 6 * sigma**6 / r**7)
            arr.append([r, energy, force])
        return arr

    def Coul(self, q1, q2, factor, rr):
        arr = []
        for r in rr[:10]:
            arr.append([r, 0, 0])
        for r in rr[10:]:
            energy = factor * 138.935458 * q1 * q2 / r
            force = factor * 138.935458 * q1 * q2 / r**2
            arr.append([r, energy, force])
        return arr

    def Coul_LJ12_6(self, q1, q2, qfactor, sigma1, espi1, sigma2, espi2, vdwfactor, rr, conver_gro_units=True):
        if conver_gro_units:
            sigma1 = sigma1 / 10.0
            sigma2 = sigma2 / 10.0
            espi1 = espi1 * 4.184
            espi2 = espi2 * 4.184
        qarr = self.Coul(q1, q2, qfactor, rr)
        vdwarr = self.LJ12_6(sigma1, espi1, sigma2, espi2, vdwfactor, rr)
        arr = []
        for i in range(len(qarr)):
            energy = qarr[i][1] + vdwarr[i][1]
            force = qarr[i][2] + vdwarr[i][2]
            arr.append([qarr[i][0], energy, force])
        return arr

    def Soft_Bond(self, k, r0, alpha, rr):
        pass

    __Func_ID = {
        "coul": Coul,  # coul interaction
        "vdw": LJ12_6,  # vdw interaction
        "q_vdw": Coul_LJ12_6,  # full nonbonded
    }

    def __call__(self):
        func = self.__Func_ID[self.simu_type]
        func(self, *self.params)


class GroInputFile:
    # special_para: object
    _name = "Gromacs"
    __Parameters_keys = [
        "atomtypes",
        "bondtpes",
        "pairtypes",
        "angletypes",
        "dihedraltypes",
        "constrainttypes",
        "nonbond_params",
    ]
    __Mole_Define_keys = [
        "atoms",
        "bonds",
        "pairs",
        "pairs_nb",
        "angles",
        "dihedrals",
        "impropers",
        "exclusions",
        "constraints",
        "virtual_sites2",
        "virtual_sites3",
        "virtual_sitesn",
        "settles",
    ]
    __Mole_Define_keys_Speical = [
        "virtual_sites2",
        "virtual_sites3",
        "virtual_sitesn",
        "position_restraints",
        "distane_restraints",
        "dihedral_restraints",
        "orientation_restraints",
        "angle_restraints",
        "angle_restraints_z",
    ]
    __System_keys = ["system", "molecules"]
    __Key_styles = {
        "atoms": ["No", "atom_type_name", "residue_number", "residue", "element", "charge_group", "charge", "mass"],
        "bonds": ["a1", "a2", "pstyle"],
        "pairs": ["a1", "a2", "pstyle"],
        "pairs_nb": ["a1", "a2", "pstyle"],
        "angles": ["a1", "a2", "a3", "pstyle"],
        "dihedrals": ["a1", "a2", "a3", "a4", "pstyle"],
        "impropers": [
            "a1",
            "a2",
            "a3",
            "a4",
            "pstyle",
        ],
        "exclusions": ["ai", "a2", "......"],
        "constraints": ["ai", "aj", "ak", "al", "fun_type"],
        "settles": ["settles for water"],
        "virtual_sites2": ["ai", "aj", "ak", "fun_type"],
        "virtual_sites3": ["ai", "aj", "ak", "al", "fun_type"],
        "virtual_sites4": ["ai", "aj", "ak", "al", "am", "fun_type"],
        "virtual_sitesn": ["ai", "aj", "ak", "al", "am", "fun_type"],
        "position_restraints": ["ai", "......"],
        "distance_restraints": ["ai", "aj", "fun_type"],
        "dihedral_restraints": ["ai", "aj", "ak", "al", "fun_type"],
        "orientation_restraints": ["ai", "aj", "fun_type"],
        "angle_restraints": ["ai", "aj", "ak", "al", "fun_type"],
        "angle_restraints_z": ["ai", "aj", "fun_type"],
    }
    __vstype_transfer = {
        "style1": [1, [1]],
        "style2": [2, [0.1]],
        "style3": [1, [1, 1]],
        "style4": [2, [1, 0.1, 1]],
        "style5": [3, [1, 0.1]],
        "style6": [4, [1, 1, 10]],
        "style7": [2, [1, 1, 0.1]],
        "style8": [4, [1, 1, 10]],
    }

    def __init__(self, style=""):
        pass

    def read_mdp(self, inscript, parent_dir):
        __param_dict = {vv: kk for kk, vv in constants.md_para["gromacs"].items()}
        __term_dict = {
            "ref_t": "float",
            "ref_p": "float",
            "rvdw": "float",
            "rcoulomb": "float",
            "coul_lambdas": "list_float",
            "vdw_lambdas": "list_float",
            "init_lambda_state": "int",
            "dt": "float",
            "nsteps": "int",
        }
        params = {}
        for line in inscript:
            ss = line.strip().split("=")
            if ss[0] in __term_dict:
                if __term_dict[ss[0]] == "float":
                    params[__param_dict[ss[0]]] = float(ss[1])
                elif __term_dict[ss[0]] == "int":
                    params[__param_dict[ss[0]]] = int(ss[1])
                elif __term_dict[ss[0]] == "list_float":
                    params[__param_dict[ss[0]]] = [float(rr) for rr in ss[1].split()]
        self.params = params
        return params

    def read_gro(self, inscript, parent_dir):
        self.Total_Atom = int(inscript[1].strip())
        self.lattics = [float(s) * 10.0 for s in inscript[-1].strip().split()]
        self.coor = []
        self.volicity = []
        volicity_flag = True
        for line in inscript[2:-1]:
            self.coor.append(
                [
                    float(line[20:28]) * 10.0,
                    float(line[28:36]) * 10.0,
                    float(line[36:44]) * 10.0,
                ]
            )
            try:
                # 单位需要处理
                self.volicity.append(
                    [
                        float(line[44:52]),
                        float(line[52:60]),
                        float(line[60:68]),
                    ]
                )
            except:  # noqa
                volicity_flag = False
                self.volicity.append(
                    [
                        "nan",
                        "nan",
                        "nan",
                    ]
                )
        if volicity_flag is False:
            delattr(self, "volicity")
            return (self.lattics, self.coor)
        else:
            return (self.lattics, self.coor, self.volicity)

    def call_read_file(self, inscript, parent_dir):
        include_f = []
        for line in inscript:
            if line.startswith("#include"):
                include_f.append(line.strip().split()[1].strip('"'))
        for f in include_f:
            this_script = []
            strlist = f.split("/")
            if len(strlist) == 2 and (
                strlist[0].startswith("gromos") or strlist[0].startswith("charmm") or strlist[0].startswith("amber")
            ):
                gmx_data = os.environ.get("GMXDATA")
                if not gmx_data:
                    raise RuntimeError(
                        "Cannot find GMXDATA enviroment variable, maybe GROMACS is not installed properly!"
                    )  # noqa
                with open(os.path.join(gmx_data, "top", f), "r") as inf:
                    for line in inf:
                        if line.startswith("#include"):
                            sub_this_script = []
                            sub_gmx_top = line.strip().split()[1].strip('"')
                            with open(os.path.join(gmx_data, "top", strlist[0], sub_gmx_top), "r") as sub_f:
                                for line in sub_f:
                                    sub_this_script.append(line)
                            self.file_scripts.append(sub_this_script)
                        else:
                            this_script.append(line)
                self.file_scripts.append(this_script)
            else:
                inf = open(f"{parent_dir}/{f}")
                for line in inf:
                    this_script.append(line)
                self.file_scripts.append(this_script)
                self.call_read_file(this_script, parent_dir)

    def read_top(self, inscript, parent_dir):
        self.file_scripts = []
        self.file_scripts.append(inscript)
        self.call_read_file(inscript, parent_dir)
        self.Moles_dict = {}
        self.System_dict = {}
        self.Force_dict = {}
        break_arr = self.__Parameters_keys + self.__System_keys + ["moleculetype"]
        for term in self.file_scripts:
            molestartline = search_line("moleculetype", term, multi="yes", annotationSymbol=";")
            for mole in molestartline:
                outscript = extra_section(mole[1] + 1, term, break_arr=break_arr, annotationSymbol=";")
                section = extra_section(0, outscript, break_arr=self.__Mole_Define_keys[1:], annotationSymbol=";")
                mole_name = section[0].split()[0]
                self.Moles_dict[mole_name] = extra_info(
                    self.__Mole_Define_keys, outscript, self.__Mole_Define_keys, annotationSymbol=";"
                )
                self.Moles_dict[mole_name]["nrexcl"] = int(section[0].split()[1])
            this_force_dict = extra_info(self.__Parameters_keys, term, break_arr, annotationSymbol=";")
            this_system_dict = extra_info(self.__System_keys, term, break_arr, annotationSymbol=";")
            self.Force_dict = add_two_dict(this_force_dict, self.Force_dict)
            self.System_dict = add_two_dict(this_system_dict, self.System_dict)

    def read_files(self, files, parent_dir):
        run_read = {
            "gro": self.read_gro,
            "mdp": self.read_mdp,
            "top": self.read_top,
        }
        for ff in files:
            f_type = str(ff)[-3:]
            if f_type in ["top", "gro", "mdp"]:
                inscript = open(ff).readlines()
                run_read[f_type](inscript, parent_dir)

    def make_system(self):
        sysobj = System("normal")
        if hasattr(self, "Force_dict"):
            sysobj.ff = self.make_force_field()

        self.this_ff = sysobj.ff
        if hasattr(self, "System_dict"):
            for rr in self.System_dict["molecules"]:
                for mm in self.Moles_dict:
                    if mm == rr[0]:
                        sysobj.mole.append(self.make_mole(mm, self.Moles_dict[mm]))
                        break
                sysobj.mole_number.append(int(rr[1]))
        if hasattr(self, "lattics"):
            sysobj.lattics = self.lattics
        if hasattr(self, "coor"):
            sysobj.coor = self.coor
        if hasattr(self, "volicity"):
            sysobj.volicity = self.volicity
        if hasattr(self, "params"):
            sysobj.params = self.params

        return sysobj

    def make_force_field(self):
        this_ff = {}
        for term, vv in self.Force_dict.items():
            this_ff[term] = {}
            if term == "atomtypes":
                for rr in self.Force_dict[term]:
                    this_ff[term][rr[0]] = {}
                    this_ff[term][rr[0]]["name"] = rr[0]
                    if len(rr) == 6:
                        this_ff[term][rr[0]]["mass"] = float(rr[1])
                        this_ff[term][rr[0]]["parameter"] = [round(float(rr[4]) * 10.0, 4), round(float(rr[5]) / 4.184, 4)]
                    elif len(rr) == 7:
                        this_ff[term][rr[0]]["mass"] = float(rr[2])
                        this_ff[term][rr[0]]["parameter"] = [round(float(rr[5]) * 10.0, 4), round(float(rr[6]) / 4.184, 4)]
        return this_ff

    def make_mole(self, mm, dicts):
        imole = Molecule("mole")
        imole.name = mm
        imole.mole_name = mm
        imole.Atoms = []
        for atom in dicts["atoms"]:
            iatom = Atom("aa")
            iatom.atom_type_name = atom[1]
            iatom.elem = atom[4]
            iatom.ff_charge = float(atom[6])
            try:
                iatom.mass = float(atom[7])
            except:  # noqa
                pass
            iatom.residu = atom[3]
            iatom.residu_number = int(atom[2])
            iatom.charge_group = int(atom[5])
            iatom.No = int(atom[0]) - 1
            if len(atom) > 8:
                iatom.atom_type_name_m2 = atom[8]
                iatom.mass_m2 = float(atom[10])
                iatom.ff_charge_m2 = float(atom[9])
            if hasattr(self, "this_ff"):
                iatom.parameter = self.this_ff["atomtypes"][iatom.atom_type_name]["parameter"]
            imole.Atoms.append(iatom)
        imole.mole_n = len(imole.Atoms)
        if "bonds" in dicts:
            if dicts["bonds"] != "NONE":
                imole.Bonds = []
                for bb in dicts["bonds"]:
                    if bb[2] == "9":
                        continue
                    a1 = int(bb[0]) - 1
                    a2 = int(bb[1]) - 1
                    ibond = Bond("bond", a1, a2)
                    for kk, vv in constants.gromacs_ftype_trans["bond"].items():
                        if vv[0] == int(bb[2]):
                            bond_type = vv
                            ibond.pstyle = kk
                            break
                    if len(bb) > 3:
                        ibond.parameter = [float(bb[3]) / bond_type[1][1], float(bb[4]) / bond_type[2][1]]
                    imole.Bonds.append(ibond)
        if "angles" in dicts:
            if dicts["angles"] != "NONE":
                imole.Angles = []
                for bb in dicts["angles"]:
                    a1 = int(bb[0]) - 1
                    a2 = int(bb[1]) - 1
                    a3 = int(bb[2]) - 1
                    iangle = Angle("angle", a1, a2, a3)
                    for kk, vv in constants.gromacs_ftype_trans["angle"].items():
                        if vv[0] == int(bb[3]):
                            angle_type = vv
                            iangle.pstyle = kk
                            break
                    if len(bb) > 4:
                        iangle.parameter = [float(bb[4]) / angle_type[1][1], float(bb[5]) / angle_type[2][1]]
                    imole.Angles.append(iangle)
        if "dihedrals" in dicts:
            if dicts["dihedrals"] != "NONE":
                imole.Dihedrals = []
                imole.impropers = []
                name1 = ""
                for dihe in dicts["dihedrals"]:
                    a1 = int(dihe[0]) - 1
                    a2 = int(dihe[1]) - 1
                    a3 = int(dihe[2]) - 1
                    a4 = int(dihe[3]) - 1
                    name = f"{a1}-{a2}-{a3}-{a4}"
                    if dihe[4] != "4":
                        if name != name1:
                            try:
                                imole.Dihedrals.append(idihe)
                            except:  # noqa
                                pass
                            idihe = Dihedral("Dihedral", a1, a2, a3, a4)
                            for kk, vv in constants.gromacs_ftype_trans["dihedral"].items():
                                if vv[0] == int(dihe[4]):
                                    dihe_type = vv
                                    idihe.pstyle = kk
                                    break
                            if len(dihe) > 5:
                                idihe.parameter = [0.0, 0.0, 0.0, 0.0]
                                ii = int(dihe[7])
                                idihe.parameter[ii - 1] = float(dihe[6]) / dihe_type[ii][1]
                            name1 = name
                        else:
                            if len(dihe) > 5:
                                ii = int(dihe[7])
                                idihe.parameter[ii - 1] = float(dihe[6]) / dihe_type[ii][1]
                    else:
                        iimproper = Improper("Improper", a3, a1, a2, a4)
                        for kk, vv in constants.gromacs_ftype_trans["improper"].items():
                            if vv[0] == int(dihe[4]):
                                improper_type = vv
                                iimproper.pstyle = kk
                                break
                        if len(dihe) > 5:
                            iimproper.parameter = [float(dihe[6]) / improper_type[0][1]]
        if "pairs" in dicts:
            if dicts["pairs"] != "NONE":
                imole.Pair14 = []
                for pair in dicts["pairs"]:
                    a1 = int(pair[0]) - 1
                    a2 = int(pair[1]) - 1
                    ipair = Pair("pair", a1, a2)
                    for kk, vv in constants.gromacs_ftype_trans["pair"].items():
                        if vv[0] == int(pair[2]):
                            pair_type = vv
                            ipair.pstyle = kk
                            break
                    if len(pair) > 3:
                        ipair.parameter = [float(pair[3]) / pair_type[1][1], float(pair[4]) / pair_type[2][1]]

                    imole.Pair14.append(ipair)
        imole.connectivity_from_bonds()
        if imole.mole_name not in ["water", "Cl-", "Na+", "AION"] and imole.mole_n <= 1000:
            exclusions = []
            if "exclusions" in dicts:
                if dicts["exclusions"] != "NONE":
                    for pair in dicts["exclusions"]:
                        a1 = int(pair[0]) - 1
                        a2 = int(pair[1]) - 1
                        exclusions.append(f"{a1}-{a2}")
            tmp_mol = deepcopy(imole)
            tmp_mol.create_intra_nonbond()
            Pair1n = []
            for pair in tmp_mol.Pair1n:
                a1 = int(pair.a1) - 1
                a2 = int(pair.a2) - 1
                name1 = f"{a1}-{a2}"
                name2 = f"{a2}-{a1}"
                if name1 not in exclusions and name2 not in exclusions:
                    pair.pstyle = "LJ12_6"
                    Pair1n.append(pair)
            if len(Pair1n) > 0:
                imole.Pair1n = Pair1n
        imole.assign_charge_to_pair()
        imole.assign_vdw_to_pair()

        return imole

    def write_input_file(self, this_path="./", write_gro=True, write_top=True):
        if write_gro:
            self.write_gro(this_path)
        if write_top:
            self.write_mole_itp(this_path)
            self.write_ff_itp(this_path)
            self.write_top(this_path)

    def write_gro(self, this_path=".", opened_file=None):
        outf = opened_file if opened_file else open(f"{this_path}/conf.gro", "w")
        outf.write("cpy generated gromacs coordination file\n")
        outf.write("%d\n" % len(self.Atoms_info))
        i = 0
        for atom in self.Atoms_info:
            i += 1
            an = atom[3] % 100000
            __dit = "".join([s for s in atom[0] if s.isdigit()])
            mn = int(__dit) % 100000
            outf.write("%5s%5s%5s%5d%8.3f%8.3f%8.3f\n" % (mn, atom[1], atom[2], an, atom[4], atom[5], atom[6]))
            #outf.write("%5s%5s%5s%5d%8.3f%8.3f%8.3f\n" % (atom[0], atom[1], atom[2], an, atom[4], atom[5], atom[6]))
        for s in self.lattics:
            outf.write(f"{s / 10:10.4f}")
        outf.write("\n")
        if not opened_file:
            outf.close()

    def write_top(self, this_path, exclusions=False):
        exclusion_str = "_exclusions" if exclusions else ""
        fname = f"{this_path}/topol{exclusion_str}.top"

        with open(fname, "w") as outf:
            if hasattr(self, "protein_force_field"):
                outf.write('#include "%s.ff/forcefield.itp"\n' % self.protein_force_field)
            outf.write('#include "ff.itp"\n')
            outf.write(";cpy generated gromacs top file\n")
            tmp = []
            for mm in self.molecules:
                if mm[0] not in tmp:
                    if mm[0] == "water":
                        outf.write('#include "%s.itp"\n' % mm[0])
                    else:
                        if exclusions is False:
                            outf.write('#include "%s.itp"\n' % mm[0])
                        else:
                            outf.write('#include "%s_exclusions.itp"\n' % mm[0])
                    tmp.append(mm[0])
            outf.write("\n[ system ]\n\n\n[ molecules ]\n")
            for vv in self.molecules:
                outf.write("%s %s\n" % (vv[0], vv[1]))
            if hasattr(self, "intermolecular_interaction"):
                outf.write("\n[ intermolecular_interactions ]\n")
                for term in ["bonds", "angles", "dihedrals"]:
                    if term in self.intermolecular_interaction:
                        outf.write("[ %s ]\n" % term)
                        for tt in self.intermolecular_interaction[term]:
                            for vv in tt:
                                outf.write("%15s" % vv)
                            outf.write("\n")
                        outf.write("\n")
                    outf.write("\n")
            outf.write("\n\n")

    def write_small_molecule_system(self, top_file_name="topol.top", coor_file_name="conf.gro"):
        """For vacuum md_simulation and energy calculation"""
        with open(top_file_name, "w") as f:
            if hasattr(self, "protein_force_field"):
                f.write('#include "amber99sb.ff/forcefield.itp"\n\n')
            self.write_ff_itp(opened_file=f)
            self.write_mole_itp(opened_file=f)
            f.write("\n[ system ]\n\n\n[ molecules ]\n")
            for vv in self.molecules:
                f.write("%s %s\n\n" % (vv[0], vv[1]))
        with open(coor_file_name, "w") as f:
            self.write_gro(opened_file=f)

    def write_mole_itp(self, this_path=".", exclusions=False, opened_file=None):
        for mole_name, info in self.mole_dict.items():
            #if mole_name not in ["water", "Cl-", "Na+", "AION"]:
            if mole_name not in ["QQQQQQQQ"]:
                exclusion_str = "_exclusions" if exclusions else ""
                fname = f"{this_path}/{mole_name}{exclusion_str}.itp"
                outf = opened_file if opened_file else open(fname, "w")
                outf.write(";cpy generated gromacs molecule top file for %s\n" % mole_name)
                if "gmx_itp_script" in info:
                    outf.write(info["gmx_itp_script"])
                    outf.close()
                    continue
                
                outf.write("[ moleculetype ]\n")
                outf.write("; Name        nrexcl\n")
                outf.write(" %s   3\n\n" % mole_name)
                for item in self.__Mole_Define_keys:
                    if item in info.keys() and item not in ["pairs", "exclusions"]:
                        if item == "impropers":
                            outf.write("[ dihedrals ]\n; ")
                        else:
                            outf.write("[ %s ]\n; " % item)
                        outf.write(" ".join(self.__Key_styles[item]))
                        outf.write("\n")
                        for record in info[item]:
                            for rr in record:
                                outf.write(format_str(rr))
                            outf.write("\n")
                        outf.write("\n")
                if exclusions is False:
                    if "pairs" in info.keys():
                        outf.write("[ pairs ]\n;")
                        outf.write(" ".join(self.__Key_styles["pairs"]))
                        outf.write("\n")
                        for record in info["pairs"]:
                            for rr in record:
                                outf.write(format_str(rr))
                            outf.write("\n")
                        outf.write("\n")
                    if "exclusions" in info.keys():
                        outf.write("[ exclusions ]\n;")
                        outf.write(" ".join(self.__Key_styles["exclusions"]))
                        outf.write("\n")
                        for record in info["exclusions"]:
                            for rr in record:
                                outf.write(format_str(rr))
                            outf.write("\n")
                        outf.write("\n")
                else:
                    if "exclusions" in info.keys():
                        outf.write("[ exclusions ]\n;")
                        outf.write(" ".join(self.__Key_styles["exclusions"]))
                        outf.write("\n")
                        for record in info["exclusions"]:
                            for rr in record:
                                outf.write(format_str(rr))
                            outf.write("\n")
                        outf.write("\n")
                if not opened_file:
                    outf.close()

    def write_ff_itp(self, this_path=".", opened_file=None):
        __defaults = ["nbfunc", "comb-rule", "gen-pairs", "fudgeLJ", "fudgeQQ"]
        outf = opened_file if opened_file else open(f"{this_path}/ff.itp", "w")
        if not hasattr(self, "protein_force_field"):
            outf.write("[ defaults ]\n;")
            outf.write(" ".join(__defaults))
            outf.write("\n")
            for aa in __defaults:
                outf.write(format_str(self.ff_dict["defaults"][aa]))
            outf.write("\n\n")

        for item in self.__Parameters_keys:
            if item in self.ff_dict.keys():
                outf.write("[ %s ]\n" % item)
                for record in self.ff_dict[item]:
                    for rr in record:
                        outf.write(format_str(rr))
                    outf.write("\n")
                outf.write("\n")
        outf.write("\n")
        if not opened_file:
            outf.close()

    @staticmethod
    def write_mdp(mdpara,this_path,free_energy=False):
        def _write_section(i,text,data):
            for kk,vv in data.items():
                if vv is None:
                    continue
                if isinstance(vv,list):
                    for ii,vvv in enumerate(vv):
                        if vvv is True:
                            vv[ii] = "yes"
                        elif vvv is False:
                            vv[ii] = "no"
                else:
                    if vv is True:
                        vv = "yes"
                    elif vv is False:
                        vv = "no"
                if isinstance(vv,list):
                    text += f"{constants.md_para['gromacs'][kk]:{format_string}}= {vv[i]}\n"
                else:
                    text += f"{constants.md_para['gromacs'][kk]:{format_string}}= {vv}\n"
            return text

        velocity_flag = True
        format_string = ">15s"
        jobs = mdpara["jobs"]
        for i in range(len(jobs)):
            with open(f"{this_path}/_{jobs[i]}.mdp", "w") as f:
                text = ""
                for kk in ["integrator","nsteps","timestep"]:
                    text += f"{constants.md_para['gromacs'][kk]:{format_string}}= {mdpara[kk][i]}\n"
                sections = ["energy","output"]
                if "mini" in jobs[i]:
                    sections.append("minimize")
                if "nvt" in jobs[i] or "npt" in jobs[i]:
                    sections.append("temperature")
                    if velocity_flag:
                        sections.append("velocity")
                        velocity_flag = False
                if "npt" in jobs[i]:
                    sections.append("pressure")
                if mdpara["constraint_used"][i]:
                    sections.append("constraint")
                if free_energy:
                    sections.append("free_energy")
                for section in sections:
                    text = _write_section(i,text,mdpara[section])
                f.write(text)

    @staticmethod
    def old_write_mdp(mdpara, this_path, free_energy=False):

        velocity_flag = True
        format_string = ">15s"

        #jobs = mdpara.get("md", "jobs").split()
        jobs = mdpara["jobs"]
        #integrator = mdpara.get("md", "integrator").split()
        integrator = mdpara["integrator"]
        #nsteps = mdpara.get("md", "nsteps").split()
        nsteps = mdpara["nsteps"]
        #dt = mdpara.get("md", "timestep").split()
        dt = mdpara["timestep"]
        coulomb_type = mdpara.get("md", "coulomb_type").split()
        Pcoupl = mdpara.get("md", "pressure_coupl").split()
        constraint = mdpara.get("md", "constraint_used").split()

        for i in range(len(jobs)):
            with open(f"{this_path}/_{jobs[i]}.mdp", "w") as f:
                f.write(f"{'integrator':{format_string}}= {mdpara['integrator'][i]}\n")
                f.write(f"{'nsteps':{format_string}}= {mdpara['nsteps'][i]}\n")
                f.write(f"{'dt':{format_string}}= {mdpara['timestep'][i]}\n")
                #f.write(f"{'coulomb_type':{format_string}}= {coulomb_type[i]}\n")

                # energy cutoff
                for key, value in mdpara["energy"].items():
                    f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")

                # output frequency
                for key, value in mdpara["output"].items():
                    f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")

                # minimize parameters
                if "mini" in jobs[i]:
                    f.write("define = -DFLEXIBLE\n")
                    for key, value in mdpara["minimize"].items():
                        f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")

                # temperature control
                if "nvt" in jobs[i] or "npt" in jobs[i]:
                    for key, value in mdpara["temperature"].items():
                        f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")
                    # only generate velocity in first eq_nvt md_simulation
                    if velocity_flag:
                        for key, value in mdpara["velocity"].items():
                            f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")
                        velocity_flag = False

                # pressure control
                if "npt" in jobs[i]:
                    f.write(f"{'Pcoupl':{format_string}}= {Pcoupl[i]}\n")
                    for key, value in mdpara["pressure"].items():
                        f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")

                # constraint
                if constraint[i] == "yes":
                    for key, value in mdpara["constraint"].items():
                        f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")

                # free_energy
                if free_energy:
                    for key, value in mdpara["free_energy"].items():
                        f.write(f"{constants.md_para['gromacs'][key]:{format_string}}= {value}\n")

    def import_systemobj(self, sm, include_ff=True, mdp=False, exclusions=False):
        self.name = sm.name
        self.lattics = sm.lattics
        # if hasattr(sm,"special_para"):
        #    self.special_para = sm.special_para
        if hasattr(sm, "intermolecular_interaction"):
            _gromacs_ii = {
                "bonds": ["6", 0.1, 100000.0],
                "angles": ["1", 1.0, 1500.0],
                "dihedrals": ["2", 1.0, 200.0],
            }
            self.intermolecular_interaction = {}
            for term, vv in sm.intermolecular_interaction.items():
                self.intermolecular_interaction[term] = []
                for atoms in vv:
                    self.intermolecular_interaction[term].append([an + 1 for an in atoms[:-1]])
                    self.intermolecular_interaction[term][-1].append(_gromacs_ii[term][0])
                    self.intermolecular_interaction[term][-1].append(round(atoms[-1] * _gromacs_ii[term][1], 2))
                    self.intermolecular_interaction[term][-1].append(_gromacs_ii[term][2])
        self.import_systemobj_mole(sm, include_ff=include_ff, exclusions=exclusions)
        self.import_systemobj_ff(sm, include_ff=include_ff)
        self.import_systemobj_Atoms(sm)
        self.import_systemobj_top(sm)
        if mdp:
            self.import_systemobj_mdp(sm.mdpara)
        if hasattr(sm, "protein_force_field"):
            self.protein_force_field = sm.protein_force_field

    def old_import_systemobj_mole(self, sm, include_ff=True, exclusions=False):
        self.mole_dict = {}
        __term_dict = {
            "Bonds": ["bonds", 2, "bond"],
            "Angles": ["angles", 3, "angle"],
            "Dihedrals": ["dihedrals", 4, "dihedral"],
            "Impropers": ["impropers", 4, "improper"],
            "Pair14": ["pairs", 2, "pair"],
        }
        pre_include_ff = include_ff
        for m in sm.molecules:
            if m.style == "protein":
                include_ff = False
            else:
                include_ff = pre_include_ff
            self.mole_dict[m.name] = {}
            self.mole_dict[m.name]["atoms"] = []
            for atom in m.Atoms:
                if hasattr(atom, "name"):
                    self.mole_dict[m.name]["atoms"].append(
                        [
                            atom.No + 1,
                            atom.atom_type_name,
                            atom.residu_number,
                            atom.residu,
                            atom.name,
                            atom.charge_group,
                            atom.ff_charge,
                            atom.mass,
                        ]
                    )
                else:
                    self.mole_dict[m.name]["atoms"].append(
                        [
                            atom.No + 1,
                            atom.atom_type_name,
                            atom.residu_number,
                            atom.residu,
                            atom.elem,
                            atom.charge_group,
                            atom.ff_charge,
                            atom.mass,
                        ]
                    )
                if hasattr(atom, "atom_type_name_m2"):
                    self.mole_dict[m.name]["atoms"][-1].extend(
                        [atom.atom_type_name_m2, atom.ff_charge_m2, atom.mass_m2],
                    )

            if hasattr(m, "Vss"):
                for vs in m.Vss:
                    nn = len(vs.patoms)
                    stylename = "virtual_sites" + str(nn)
                    if stylename not in self.mole_dict[m.name].keys():
                        self.mole_dict[m.name][stylename] = []
                    tmp_arr = [vs.a1]
                    for rra in vs.patoms:
                        tmp_arr.append(rra)
                    tmp_arr.append(self.__vstype_transfer[vs.style])
                    for rra in vs.psets:
                        tmp_arr.append(rra)
                    self.mole_dict[m.name][stylename].append(tmp_arr)
            if include_ff:
                for kk in ["Bonds", "Angles"]:
                    if hasattr(m, kk):
                        vv = __term_dict[kk]
                        self.mole_dict[m.name][vv[0]] = self.import_mole_term_ff(getattr(m, kk), vv[2], vv[1])
                if hasattr(m, "Dihedrals"):
                    if len(m.Dihedrals) > 0 and m.Dihedrals[0].pstyle in ["fourier", "amber", "opls"]:
                        self.mole_dict[m.name]["dihedrals"] = self.import_mole_dihedral_fourier(m.Dihedrals)
                    else:
                        self.mole_dict[m.name]["dihedrals"] = self.import_mole_term_ff(
                            getattr(m, "Dihedrals"), "dihedral", 4
                        )
                if hasattr(m, "Impropers"):
                    if len(m.Impropers) > 0 and m.Impropers[0].pstyle in ["fourier", "amber", "opls"]:
                        self.mole_dict[m.name]["impropers"] = self.import_mole_improper_fourier(m.Impropers)
                    else:
                        self.mole_dict[m.name]["impropers"] = self.import_mole_term_ff(
                            getattr(m, "impropers"), "improper", 4
                        )
                if hasattr(m, "Pair14"):
                    self.mole_dict[m.name]["pairs"] = self.import_mole_term(getattr(m, "Pair14"), "pair", 2)
                if hasattr(m, "AlteredPairs"):
                    self.mole_dict[m.name]["pairs"].extend(self.import_altered_pairs(getattr(m, "AlteredPairs")))
            else:
                for kk, vv in __term_dict.items():
                    if hasattr(m, kk):
                        self.mole_dict[m.name][vv[0]] = self.import_mole_term(getattr(m, kk), vv[2], vv[1])
            if exclusions:
                self.mole_dict[m.name]["exclusions"] = self.import_systemobj_exclusions(m)
            elif hasattr(m, "rfe_exclusions"):
                self.mole_dict[m.name]["exclusions"] = []
                for pp in m.rfe_exclusions:
                    self.mole_dict[m.name]["exclusions"].append([pp[0] + 1, pp[1] + 1])

    def import_systemobj_mole(self, sm, include_ff=True, exclusions=False):
        self.mole_dict = {}
        __term_dict = {
            "Bonds": ["bonds", 2, "bond"],
            "Angles": ["angles", 3, "angle"],
            "Dihedrals": ["dihedrals", 4, "dihedral"],
            "Impropers": ["impropers", 4, "improper"],
            "Pair14": ["pairs", 2, "pair"],
        }
        for m in sm.molecules:
            self.mole_dict[m.name] = {}
            self.mole_dict[m.name]["atoms"] = []
            if hasattr(m,"gmx_itp_script"):
                self.mole_dict[m.name]["gmx_itp_script"] = m.gmx_itp_script

            for atom in m.Atoms:
                if hasattr(atom, "name"):
                    self.mole_dict[m.name]["atoms"].append(
                        [
                            atom.No + 1,
                            #### change to atom.atom_type_used_name by CFL on 2024-05-24
                            ###atom.atom_type_name,
                            atom.atom_type_used_name,
                            atom.residu_number,
                            atom.residu,
                            atom.name,
                            atom.charge_group,
                            atom.ff_charge,
                            atom.mass,
                        ]
                    )
                else:
                    self.mole_dict[m.name]["atoms"].append(
                        [
                            atom.No + 1,
                            ### change to atom.atom_type_used_name by CFL on 2024-05-24
                            ###atom.atom_type_name,
                            atom.atom_type_used_name,
                            atom.residu_number,
                            atom.residu,
                            atom.elem,
                            atom.charge_group,
                            atom.ff_charge,
                            atom.mass,
                        ]
                    )
                if hasattr(atom, "atom_type_name_m2"):
                    self.mole_dict[m.name]["atoms"][-1].extend(
                        [atom.atom_type_name_m2, atom.ff_charge_m2, atom.mass_m2],
                    )

            if hasattr(m, "Vss"):
                for vs in m.Vss:
                    nn = len(vs.patoms)
                    stylename = "virtual_sites" + str(nn)
                    if stylename not in self.mole_dict[m.name].keys():
                        self.mole_dict[m.name][stylename] = []
                    tmp_arr = [vs.a1 + 1]
                    for rra in vs.patoms:
                        tmp_arr.append(rra + 1)
                    tmp_arr.append(self.__vstype_transfer[vs.style][0])

                    for index, rra in enumerate(vs.psets):
                        tmp_arr.append(rra * self.__vstype_transfer[vs.style][1][index])
                    self.mole_dict[m.name][stylename].append(tmp_arr)

            for kk in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
                if hasattr(m, kk):
                    vv = __term_dict[kk]
                    self.mole_dict[m.name][vv[0]] = self.import_mole_term_ff(getattr(m, kk), vv[2], vv[1])
            if hasattr(m, "Pair14"):
                self.mole_dict[m.name]["pairs"] = self.import_mole_term(getattr(m, "Pair14"), "pair", 2)
            if hasattr(m, "AlteredPairs"):
                self.mole_dict[m.name]["pairs"].extend(self.import_altered_pairs(getattr(m, "AlteredPairs")))
            if exclusions:
                self.mole_dict[m.name]["exclusions"] = self.import_systemobj_exclusions(m)
            elif hasattr(m, "rfe_exclusions"):
                self.mole_dict[m.name]["exclusions"] = []
                for pp in m.rfe_exclusions:
                    self.mole_dict[m.name]["exclusions"].append([pp[0] + 1, pp[1] + 1])

    ###暂时不用的函数
    def close_pair_parameters(self,atom1,atom2,fterm,dummy_tag):
        for aa in [atom1,atom2]:
            if not hasattr(aa,"ff_charge_m2"):
                setattr(aa,"ff_charge_m2",aa.ff_charge)
            if not hasattr(aa,"parameter_m2"):
                setattr(aa,"parameter_m2",aa.parameter)
        if fterm == "pair14_1n" and dummy_tag == "PP":
            tag = "-14+1n"
        elif fterm == "pair1n_14" and dummy_tag == "PP":
            tag = "-1n+14"
        elif (fterm == "pair23_1n" and dummy_tag in ("PP", "DP")) or (fterm == "pair14_1n" and dummy_tag == "DP"):
            tag = "+1n"
        elif (fterm == "pair23_14" and dummy_tag in ("PP", "DP")) or (fterm == "pair1n_14" and dummy_tag == "DP"):
            tag = "+14"
        elif (fterm == "pair1n_23" and dummy_tag in ("PP", "PD")) or (fterm == "pair1n_14" and dummy_tag == "PD"):
            tag = "-1n"
        elif (fterm == "pair14_23" and dummy_tag in ("PP", "PD")) or (fterm == "pair14_1n" and dummy_tag == "PD"):
            tag = "-14"
        else:
            tag = None
    
        if tag == "-14+1n":
            parameters1 = [atom1.ff_charge,atom2.ff_charge,0.8333,atom1.parameter[0],atom1.parameter[1],atom2.parameter[0],atom2.parameter[1], 0.5] # noqa
            parameters2 = [atom1.ff_charge_m2,atom2.ff_charge_m2,1.0, atom1.parameter_m2[0],atom1.parameter_m2[1],atom2.parameter_m2[0],atom2.parameter_m2[1], 1.0,] # noqa
        elif tag == "-1n+14":
            parameters1 = [atom1.ff_charge,atom2.ff_charge,1.0,atom1.parameter[0],atom1.parameter[1],atom2.parameter[0],atom2.parameter[1], 1.0,] # noqa
            parameters2 = [atom1.ff_charge_m2,atom2.ff_charge_m2,0.8333,atom1.parameter_m2[0],atom1.parameter_m2[1],atom2.parameter_m2[0],atom2.parameter_m2[1], 0.5,] # noqa
        elif tag == "+14":
            parameters1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            parameters2 = [atom1.ff_charge_m2,atom2.ff_charge_m2,0.8333, atom1.parameter_m2[0],atom1.parameter_m2[1],atom2.parameter_m2[0],atom2.parameter_m2[1], 0.5,] # noqa
        elif tag == "-14":
            parameters1 = [atom1.ff_charge,atom2.ff_charge,0.8333,atom1.parameter[0],atom1.parameter[1],atom2.parameter[0],atom2.parameter[1], 0.5,] # noqa
            parameters2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        elif tag == "-1n":
            parameters1 = [atom1.ff_charge,atom2.ff_charge,1.0,atom1.parameter[0],atom1.parameter[1],atom2.parameter[0],atom2.parameter[1], 1.0,] # noqa
            parameters2 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        elif tag == "+1n":
            parameters1 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
            parameters2 = [atom1.ff_charge_m2,atom2.ff_charge_m2,1.0,atom1.parameter_m2[0],atom1.parameter_m2[1],atom2.parameter_m2[0],atom2.parameter_m2[1], 1.0,] # noqa
        else:
            return False
        return parameters1, round(coul_vdw_value(parameters1),2), parameters2, round(coul_vdw_value(parameters2),2)

    def close_convert_pair_to_bond(self, pair_dict, Atoms):
        arr = []
        bn = -1
        parameter_dicts = {}
        for dummy_tag, group in pair_dict.items():
            for fterm, terms in group.items():
                for term in terms:
                    if not self.pair_parameters(Atoms[term[0]], Atoms[term[1]], fterm, dummy_tag):
                        continue
                    param1, label1, param2, label2 = self.pair_parameters(Atoms[term[0]], Atoms[term[1]], fterm, dummy_tag)
                    if label1 not in parameter_dicts and label1 != 0.0:
                        bn += 1
                        parameter_dicts[label1] = [bn, param1]
                    if label2 not in parameter_dicts and label2 != 0.0:
                        bn += 1
                        parameter_dicts[label2] = [bn, param2]
                    if label1 != 0.0:
                        arr.append([term[0] + 1,term[1] + 1,9,parameter_dicts[label1][0],1.0,parameter_dicts[label1][0],0.0]) # noqa
                    if label2 != 0.0:
                        arr.append([term[0] + 1,term[1] + 1,9,parameter_dicts[label2][0],0.0,parameter_dicts[label2][0],1.0]) # noqa
        for lable, term in parameter_dicts.items():
            fname = self.special_para["table_path"] + "/t_b%d.xvg" % term[0]
            Tablen = TableFunction("q_vdw",term[1],fname=fname) # noqa
        return arr

    ######
    def import_systemobj_ff(self, sm, include_ff=True):
        __term_dict = {
            "bondterm": ["bondtypes", 2, "bond"],
            "angleterm": ["angletypes", 3, "angle"],
            "dihedralterm": ["dihedraltypes", 4, "dihedral"],
            "improperterm": ["dihedraltypes", 4, "improper"],
        }
        default_para = {"nbfunc": 1, "comb-rule": 2, "gen-pairs": "yes", "fudgeLJ": 0.500, "fudgeQQ": 0.8333}
        self.ff_dict = {}
        self.ff_dict["defaults"] = {}
        for aa, bb in default_para.items():
            self.ff_dict["defaults"][aa] = bb
        self.ff_dict["atomtypes"] = []
        for at,para in sm.ff["atomtype"].items():
            self.ff_dict["atomtypes"].append([at,para["mass"],0.0,"A",para["parameter"][0] * 0.1, para["parameter"][1] * 4.184])
        if "pairwise" in sm.ff:
            self.ff_dict["nonbond_params"] = []
            for at, para in sm.ff["pairwise"].items():
                ss = at.split("$")
                self.ff_dict["nonbond_params"].append([ss[0], ss[1],1, para["parameter"][0] * 0.1, para["parameter"][1] * 4.184])
        if not include_ff:
            for kk in ["bondterm", "angleterm", "improperterm"]:
                vv = __term_dict[kk]
                self.ff_dict[bb[0]] = self.import_mole_term_ff(sm.ff[kk], vv[2], vv[1], ff_flag=True)
            for m in sm.molecules:
                if len(m.Dihedrals) > 0 and m.Dihedrals[0].pstyle in ["fourier", "amber", "opls"]:
                    self.ff_dict[m.name]["dihedrals"] = self.import_mole_dihedral_fourier(
                        sm.ff["dihedralterm"], ff_flag=True
                    )
                else:
                    self.ff_dict[m.name]["dihedrals"] = self.import_mole_term_ff(sm.ff[kk], "dihedral", 4, ff_flag=True)

    def import_systemobj_ff_old(self, sm, include_ff=True):
        __term_dict = {
            "bondterm": ["bondtypes", 2, "bond"],
            "angleterm": ["angletypes", 3, "angle"],
            "dihedralterm": ["dihedraltypes", 4, "dihedral"],
            "improperterm": ["dihedraltypes", 4, "improper"],
        }
        default_para = {"nbfunc": 1, "comb-rule": 2, "gen-pairs": "yes", "fudgeLJ": 0.500, "fudgeQQ": 0.8333}
        self.ff_dict = {}
        self.ff_dict["defaults"] = {}
        for aa, bb in default_para.items():
            self.ff_dict["defaults"][aa] = bb
        self.ff_dict["atomtypes"] = []
        tmp = []
        water_exist = False
        for m in sm.molecules:
            for aa in m.Atoms:
                if hasattr(aa, "parameter") and aa.atom_type_name not in tmp:
                    tmp.append(aa.atom_type_name)
                    if aa.atom_type_name == "o_2w":
                        water_exist = True
                        o_2w_para = aa.parameter
                    self.ff_dict["atomtypes"].append(
                        [aa.atom_type_name, aa.mass, 0.0, "A", aa.parameter[0] * 0.1, aa.parameter[1] * 4.184]
                    )
                if hasattr(aa, "parameter_m2") and aa.atom_type_name_m2 not in tmp:
                    tmp.append(aa.atom_type_name_m2)
                    if aa.atom_type_name_m2 == "o_2w":
                        water_exist = True
                        o_2w_para = aa.parameter
                    self.ff_dict["atomtypes"].append(
                        [aa.atom_type_name_m2, aa.mass_m2, 0.0, "A", aa.parameter_m2[0] * 0.1, aa.parameter_m2[1] * 4.184]
                    )
        if water_exist:
            self.ff_dict["nonbond_params"] = []
            tmp = []
            for m in sm.molecules:
                for aa in m.Atoms:
                    if (
                        hasattr(aa, "parameter")
                        and aa.atom_type_name not in tmp
                        and aa.atom_type_name != "h_1w"
                        and len(aa.parameter) > 2
                        and aa.parameter[2] != 1.0
                    ):
                        tmp.append(aa.atom_type_name)
                        sigma_1 = o_2w_para[0]
                        epsilon_1 = o_2w_para[1]
                        sigma_2 = aa.parameter[0]
                        epsilon_2 = aa.parameter[1]
                        sigma_mix = 0.5 * (sigma_1 + sigma_2)
                        epsilon_mix = math.pow(epsilon_1 * epsilon_2, 0.5)
                        sigma = sigma_mix / math.pow(aa.parameter[2], 1.0 / 6.0)
                        epsilon = epsilon_mix * math.pow(aa.parameter[2], 2.0)
                        self.ff_dict["nonbond_params"].append(
                            ["o_2w", aa.atom_type_name, 1, sigma * 0.1, epsilon * 4.184]
                        )
                    if (
                        hasattr(aa, "parameter_m2")
                        and aa.atom_type_name_m2 not in tmp
                        and aa.atom_type_name_m2 != "h_1w"
                        and len(aa.parameter_m2) > 2
                        and aa.parameter_m2[2] != 1.0
                    ):
                        tmp.append(aa.atom_type_name_m2)
                        sigma_1 = o_2w_para[0]
                        epsilon_1 = o_2w_para[1]
                        sigma_2 = aa.parameter_m2[0]
                        epsilon_2 = aa.parameter_m2[1]
                        sigma_mix = 0.5 * (sigma_1 + sigma_2)
                        epsilon_mix = math.pow(epsilon_1 * epsilon_2, 0.5)
                        sigma = sigma_mix / math.pow(aa.parameter_m2[2], 1.0 / 6.0)
                        epsilon = epsilon_mix * math.pow(aa.parameter_m2[2], 2.0)
                        self.ff_dict["nonbond_params"].append(
                            ["o_2w", aa.atom_type_name_m2, 1, sigma * 0.1, epsilon * 4.184]
                        )
            if self.ff_dict["nonbond_params"] == []:
                del self.ff_dict["nonbond_params"]
        if not include_ff:
            for kk in ["bondterm", "angleterm", "improperterm"]:
                vv = __term_dict[kk]
                self.ff_dict[bb[0]] = self.import_mole_term_ff(sm.ff[kk], vv[2], vv[1], ff_flag=True)
            for m in sm.molecules:
                if len(m.Dihedrals) > 0 and m.Dihedrals[0].pstyle in ["fourier", "amber", "opls"]:
                    self.ff_dict[m.name]["dihedrals"] = self.import_mole_dihedral_fourier(
                        sm.ff["dihedralterm"], ff_flag=True
                    )
                else:
                    self.ff_dict[m.name]["dihedrals"] = self.import_mole_term_ff(sm.ff[kk], "dihedral", 4, ff_flag=True)

    def import_rest2_system_para(self, sm, window_idx):
        if not hasattr(sm, "rest2_scale_factor"):
            return
        scale_factor = sm.rest2_scale_factor[window_idx]
        for m in sm.molecules:
            if m.style == "protein":
                for dihedral in m.Dihedrals:
                    if (
                        torsion_key := (dihedral.a1, dihedral.a2, dihedral.a3, dihedral.a4)
                    ) in sm.rest2_protein_hot_torsions:
                        dihedral.parameter = [0.0, 0.0, 0.0, 0.0]
                        for item in sm.rest2_protein_hot_torsions[torsion_key]:
                            dihedral.parameter[item.pn] = scale_factor * item.kd / 4.184
                for idx, atom in sm.rest2_protein_hot_atoms.items():
                    m.Atoms[idx].ff_charge *= math.sqrt(scale_factor)
                    # m.Atoms[idx].atom_type_name += "!!"
            if m.style == "ligand":
                for dihedral in m.Dihedrals:
                    if (dihedral.a1, dihedral.a2, dihedral.a3, dihedral.a4) in sm.rest2_ligand_hot_torsions:
                        dihedral.parameter = [item * scale_factor for item in dihedral.parameter]
                for atom in m.Atoms:
                    if atom.No in sm.rest2_ligand_hot_atoms:
                        atom.ff_charge *= math.sqrt(scale_factor)

        # SKIP this now
        # atom_types_data = {}
        # for atom_type in sm.rest2_protein_hot_atoms.values():
        #     atom_types_data[atom_type.name] = atom_type
        # for atom_type in atom_types_data.values():
        #     self.ff_dict['atomtypes'].append([atom_type.name + "!", atom_type.mass, atom_type.charge, atom_type.ptype,
        #                                       atom_type.sigma,  atom_type.epsilon*scale_factor])
        self.import_systemobj_mole(sm)

    def import_mole_term(
        self,
        terms,
        fterm,
        tn,
    ):
        arrs = []
        for term in terms:
            arr = []
            if fterm != "improper":
                for i in range(1, tn + 1):
                    arr.append(getattr(term, "a%d" % i) + 1)
            else:
                arr.append(term.a2 + 1)
                arr.append(term.a3 + 1)
                arr.append(term.a1 + 1)
                arr.append(term.a4 + 1)
            arr.append(constants.gromacs_ftype_trans[fterm][term.pstyle][0])
            arrs.append(arr)
        return arrs

    def import_altered_pairs(self, terms):
        arrs = []
        for term in terms:
            arr = []
            arr.append(term.a1 + 1)
            arr.append(term.a2 + 1)
            arr.append(2)
            arr.extend(term.ff_charge)
            arr.extend([term.parameter[0] * 0.1, term.parameter[1] * 4.184])
            arrs.append(arr)
        return arrs

    def import_single_mole_term_ff(self, arr, fterm, term, first_para, second_para):
        if not hasattr(term, first_para):
            return
        arr.append(constants.gromacs_ftype_trans[fterm][term.pstyle][0])
        for pp in constants.gromacs_ftype_trans[fterm][term.pstyle][1:]:
            arr.append(getattr(term, first_para)[pp[0]] * pp[1])
        # perturbed parameters if exist
        if second_para is not None:
            for pp in constants.gromacs_ftype_trans[fterm][term.pstyle][1:]:
                try:
                    arr.append(getattr(term, second_para)[pp[0]] * pp[1])
                except AttributeError:
                    pass
        return arr

    def import_single_dihedral_fourier(self, arr, fterm, term, first_para, second_para):
        if not hasattr(term, first_para):
            return
        arrs = []
        pp = getattr(term, first_para)
        pp_m2 = getattr(term, second_para, pp)
        for i in range(0,len(pp),2):
            tmp_arr = deepcopy(arr)
            #aa = 180.0 * (i % 2)
            if second_para is not None:
                if pp[i] != 0.0 or pp_m2[i] != 0.0:
                    tmp_arr.extend([9, pp[i+1], pp[i] * 4.184, (i // 2) + 1])
                    if pp_m2[i] != pp[i] or pp_m2[i+1] != pp[i+1]:
                        tmp_arr.extend([pp_m2[i+1], pp_m2[i] * 4.184, (i //2) + 1])
            else:
                if pp[i] != 0.0:
                    tmp_arr.extend([9, pp[i+1], pp[i] * 4.184, (i // 2) + 1])
            if len(tmp_arr) > 4:
                arrs.append(tmp_arr)
        if len(arrs) == 0:
            arrs.append(arr + [9, 180, 0.0, 2])
        return arrs

    def import_single_improper_fourier(self, arr, fterm, term, first_para, second_para):
        if not hasattr(term, first_para):
            return
        arrs = []
        for i in range(len(term.parameter)):
            tmp_arr = deepcopy(arr)
            tmp_arr.extend([4, 180.0, getattr(term, first_para)[i] * 4.184, 2])

            if second_para is not None:
                try:
                    tmp_arr.extend([180.0, getattr(term, second_para)[i] * 4.184, 2])
                except AttributeError:
                    pass
            arrs.append(tmp_arr)
        return arrs

    __get_para_func = {
        "normal_bonded": import_single_mole_term_ff,
        "dihedral_fourier": import_single_dihedral_fourier,
        "improper_fourier": import_single_improper_fourier,
    }

    def import_mole_term_ff(self, terms, fterm, tn, ff_flag=False):
        arrs = []
        for term in terms:
            arr = []
            if fterm != "improper":
                for i in range(1, tn + 1):
                    if ff_flag:
                        arr.append(getattr(term, "a%d_atom_type" % i))
                    else:
                        arr.append(getattr(term, "a%d" % i) + 1)
            else:
                if ff_flag:
                    arr.append(term.a2_atom_type)
                    arr.append(term.a3_atom_type)
                    arr.append(term.a1_atom_type)
                    arr.append(term.a4_atom_type)
                else:
                    arr.append(term.a2 + 1)
                    arr.append(term.a3 + 1)
                    arr.append(term.a1 + 1)
                    arr.append(term.a4 + 1)
            if fterm == "dihedral" and term.pstyle in ["fourier", "amber", "opls"]:
                __func_type = "dihedral_fourier"
            elif fterm == "improper" and term.pstyle in ["fourier", "amber", "opls"]:
                __func_type = "improper_fourier"
            else:
                __func_type = "normal_bonded"

            __flag = True
            if hasattr(term, "core_hopping_flag"):
                if fterm == "bond":
                    if term.core_hopping_flag:
                        # bond_style = self.special_para["core_hopping_bond"][0]
                        # if bond_style == "morse":
                        beta_scale_factor = term.parameter[2] if len(term.parameter) == 3 else 1.0
                        if beta_scale_factor == 0.0:
                            arr.append(constants.gromacs_ftype_trans[fterm]["harmonic"][0])
                            para_scale = [pp[1] for pp in constants.gromacs_ftype_trans[fterm]["harmonic"][1:]]
                            arr.append(term.parameter[0] * para_scale[0])
                            arr.append(term.parameter[1] * para_scale[1])
                        else:
                            arr.append(constants.gromacs_ftype_trans[fterm]["morse"][0])
                            morse_alpha = float(constants.gromacs_ftype_trans[fterm]["morse"][1]) * beta_scale_factor
                            para_scale = [pp[1] for pp in constants.gromacs_ftype_trans[fterm]["morse"][2:]]
                            arr.append(term.parameter[0] * para_scale[0])
                            arr.append(term.parameter[1] / morse_alpha**2 * para_scale[1])
                            arr.append(morse_alpha * para_scale[2])
                        arr_term = arr
                        __flag = False
                    else:
                        arr_term = None
                        __flag = False
                else:
                    if term.core_hopping_flag:
                        args = [self, arr, fterm, term, "parameter", None]
                        # arr = self.import_single_mole_term_ff(arr,fterm,term,"para","para_m2")
                    else:
                        arr_term = None
                        __flag = False
            elif hasattr(term, "dummy_bonded_m1") or hasattr(term, "dummy_bonded_m2") or hasattr(term, "use_parameter_m2"):
                for bonded_attr in ["dummy_bonded_m1", "dummy_bonded_m2", "use_parameter_m2"]:
                    if hasattr(term, bonded_attr):
                        break
                if getattr(term, bonded_attr) == "m2":
                    args = [self, arr, fterm, term, "parameter_m2", None]
                    # arr = self.import_single_mole_term_ff(arr,fterm,term,"para_m2",None)
                elif getattr(term, bonded_attr) == "m1-m2":
                    args = [self, arr, fterm, term, "parameter", "parameter_m2"]
                    # arr = self.import_single_mole_term_ff(arr,fterm,term,"para","para_m2")
                elif getattr(term, bonded_attr) == "m1":
                    args = [self, arr, fterm, term, "parameter", None]
                    # arr = self.import_single_mole_term_ff(arr,fterm,term,"para",None)
            else:
                args = [self, arr, fterm, term, "parameter", "parameter_m2"]
                # arr = self.import_single_mole_term_ff(arr,fterm,term,"para","para_m2")
            if __flag:
                arr_term = self.__get_para_func[__func_type](*args)

            if isinstance(arr_term, list):
                if np.array(arr_term).ndim == 1:
                    arrs.append(arr_term)
                else:
                    arrs.extend(arr_term)
            else:
                arr.append(constants.gromacs_ftype_trans[fterm][term.pstyle][0])
                arrs.append(arr)
        return arrs

    def import_mole_dihedral_fourier(self, terms, include_ff=True, ff_flag=False):
        arrs = []
        if include_ff:
            for term in terms:
                pp = term.parameter
                pp_m2 = getattr(term, "parameter_m2", pp)
                flag = True
                for i in range(0,len(pp),2):
                    if pp[i] != 0.0 or pp_m2[i] != 0.0:
                        flag = False
                        #aa = 180.0 * (i % 2)
                        arrs.append([term.a1 + 1, term.a2 + 1, term.a3 + 1, term.a4 + 1, 9, pp[i+1], pp[i] * 4.184, (i // 2) + 1])
                        if pp_m2[i] != pp[i] or pp_m2[i+1] != pp[i+1]:
                            arrs[-1].extend([pp_m2[i+1], pp_m2[i] * 4.184, (i // 2) + 1])
                if flag:
                    arrs.append([term.a1 + 1, term.a2 + 1, term.a3 + 1, term.a4 + 1, 9, 0.0, 0.0000, 1])
        else:
            for term in terms:
                arrs.append([term.a1 + 1, term.a2 + 1, term.a3 + 1, term.a4 + 1, 9])
        return arrs

    def import_mole_improper_fourier(self, terms, include_ff=True, ff_flag=False):
        arrs = []
        for term in terms:
            pp = term.parameter
            pp_m2 = getattr(term, "parameter_m2", pp)
            for i in range(len(pp)):
                if include_ff:
                    arrs.append([term.a2 + 1, term.a3 + 1, term.a1 + 1, term.a4 + 1, 4, 180.0, pp[i] * 4.184, 2])
                    if pp_m2[i] != pp[i]:
                        arrs[-1].extend([180.0, pp_m2[i] * 4.184, 2])
                else:
                    arrs.append([term.a2 + 1, term.a3 + 1, term.a1 + 1, term.a4 + 1, 4])
        return arrs

    def import_systemobj_exclusions(self, m):
        arrs = []
        if hasattr(m, "Pair14"):
            for tt in m.Pair14:
                arrs.append([tt.a1 + 1, tt.a2 + 1])
        if hasattr(m, "Pair1n"):
            for tt in m.Pair1n:
                arrs.append([tt.a1 + 1, tt.a2 + 1])
        return arrs

    def import_systemobj_Atoms(self, sm):
        self.Atoms_info = []
        n = 0
        nn = 0
        for i in range(len(sm.molecules)):
            m = sm.molecules[i]
            for _ in range(sm.molecule_number[i]):
                for k in range(len(m.Atoms)):
                    n += 1
                    residue_ID = m.Atoms[k].residue_ID
                    if isinstance(residue_ID,int) or residue_ID.isdigit():
                        resn = str(nn + int(residue_ID))
                    else:
                        digit = [d for d in residue_ID if d.isdigit()]
                        string = [d for d in residue_ID if d.isalpha()]
                        resn = nn + int("".join(digit))
                        resn = str(resn) + "".join(string)
                    arr = [
                        resn,
                        m.Atoms[k].residu,
                        m.Atoms[k].atom_name if m.style == "protein" else m.Atoms[k].elem,
                        n,
                        sm.coordinates[n - 1][0] / 10.0,
                        sm.coordinates[n - 1][1] / 10.0,
                        sm.coordinates[n - 1][2] / 10.0,
                    ]
                    self.Atoms_info.append(arr)
                nn += m.residu_n

    def import_systemobj_top(self, sm):
        self.molecules = []
        for i in range(len(sm.molecules)):
            self.molecules.append([sm.molecules[i].name, sm.molecule_number[i]])

    def import_systemobj_mdp(self, mdpara):
        self.md_para = {}
        # for kword, para in sm.md_para.items():
        for kword, para in mdpara.items():
            self.md_para[constants.md_para["gromacs"][kword]] = para

    def get_connectivity(self, bond_arr):
        connect_dict = {}
        for record in bond_arr:
            if record[0] not in connect_dict.keys():
                connect_dict[record[0]] = [int(record[1])]
            else:
                if int(record[1]) not in connect_dict[record[0]]:
                    connect_dict[record[0]].append(int(record[1]))
            if record[1] not in connect_dict.keys():
                connect_dict[record[1]] = [int(record[0])]
            else:
                if int(record[0]) not in connect_dict[record[1]]:
                    connect_dict[record[1]].append(int(record[0]))
        return connect_dict


class GroOutputFile:
    def __init__(self, style):
        self.style = style
        self.thermo = {}
        self.trj = {"coor": [], "volicity": [], "force": [], "lattics": []}
        self.property = {}
        self.mole_n = 0

    def read_files(self, file_dicts):
        if "thermo" in file_dicts:
            for tt, ff in file_dicts["thermo"].items():
                inscript = open(ff).readlines()
                self.read_thermo(inscript, tt)
        if "trj" in file_dicts:
            f_type = file_dicts["trj"][-3:]
            if f_type == "xtc":
                self.read_trj_xtc(file_dicts["trj"])
            else:
                inscript = open(file_dicts["trj"]).readlines()
                self.read_trj(inscript, f_type)
        if "property" in file_dicts:
            for tt, ff in file_dicts["property"].items():
                inscript = open(ff).readlines()
                self.read_property(inscript, tt)
        if "mole_n" in file_dicts:
            inscript = open(file_dicts["mole_n"]).readlines()
            for i in range(len(inscript)):
                if inscript[i].find("molecules") != -1:
                    self.mole_n = int(inscript[i + 1].split()[1].strip())

    def read_energy(self, inscript):
        __energy = [
            "bond",
            "angle",
            "proper",
            "improper",
            "lj-14",
            "coulomb-14",
            "lj",
            "disper.",
            "coulomb",
            "coul.",
            "potential",
            "kinetic",
            "total",
            "dvcoul/dl",
            "dvvdw/dl",
        ]
        __volume = ["Volume"]
        __label = {
            "bond": "Bond",
            "angle": "Angle",
            "proper": "Dihedral",
            "improper": "Improper",
            "lj-14": "LJ-14",
            "coulomb-14": "Coul-14",
            "lj": "LJ",
            "disper.": "LJ-long",
            "coulomb": "Coul",
            "coul.": "Coul-long",
            "potential": "Potential",
            "kinetic": "Kinetic",
            "total": "Energy",
        }
        pps = []
        for line in inscript:
            if line[0] != "#":
                ss = line.strip().split()
                if line[0] == "@":
                    if ss[1][0] == "s":
                        pre_term = ss[3].strip('"').lower()
                        if pre_term in __label:
                            term = __label[pre_term]
                        else:
                            term = pre_term
                        self.thermo[term] = []
                        pps.append(pre_term)
                else:
                    for i in range(1, len(ss)):
                        if pps[i - 1] in __label:
                            term = __label[pps[i - 1]]
                        else:
                            term = pps[i - 1]
                        if pps[i - 1] in __energy:
                            self.thermo[term].append(float(ss[i]) / 4.184)
                        elif pps[i - 1] in __volume:
                            self.thermo[term].append(float(ss[i]) * 1000.0)
                        else:
                            self.thermo[term].append(float(ss[i]))

    def read_hov(self, inscript):
        self.thermo["conhesive"] = []
        for line in inscript:
            if line[0] != "#" and line[0] != "@":
                ss = line.strip().split()
                v = 0
                for s in ss[1:]:
                    v += float(s.strip())
                self.thermo["conhesive"].append(v)

    def read_dhdl(self, inscript):
        self.thermo["dhdl"] = []
        for line in inscript:
            if line[0] != "#":
                ss = line.strip().split()
                if line[0] == "@":
                    if ss[1][0] == "s":
                        self.thermo["dhdl"].append([])
                else:
                    for i in range(1, len(ss)):
                        self.thermo["dhdl"][i - 1].append(float(ss[i]))

    def read_thermo(self, inscript, f_type):
        __run_read = {
            "energy": self.read_energy,
            "hov": self.read_hov,
            "dhdl": self.read_dhdl,
        }
        __run_read[f_type](inscript)

    def read_epsilon(self, inscript):
        self.property["dielectronic"] = []
        for line in inscript:
            if line[0] != "#" and line[0] != "@":
                ss = line.strip().split()
                self.property["dielectronic"].append(float(ss[1]))

    def read_property(self, inscript, f_type):
        __run_read = {
            "dielectronic": self.read_epsilon,
        }
        __run_read[f_type](inscript)

    def read_trj(self, inscript, f_type):
        __run_read = {
            "gro": self.read_trj_gro,
        }
        __run_read[f_type](inscript)

    def read_trj_gro(self, inscript):
        groinput = GroInputFile("read_gro")
        Total_Atom = int(inscript[1].strip())
        n = Total_Atom + 3
        for ii in range(0, len(inscript), n):
            this_inscript = inscript[ii : ii + n]
            record = groinput.read_gro(this_inscript, "./")
            self.trj["lattics"].append(record[0])
            self.trj["coor"].append(record[1])
            if len(record) > 2:
                self.trj["volicity"].append(record[2])

    def gro_to_xtc(self, filename):
        if not _XTC_IMPORTED:
            raise ImportError("XTCFile cython extension not compiled")

        time = 1.0
        with XTCFile(filename, "w") as xtc:
            for i in range(len(self.trj["coor"])):
                xtc.write(np.array(self.trj["coor"][i]) / 10, np.diag(self.trj["lattics"][i]), i + 1, time * (i + 1))

    def read_trj_xtc(self, filename):
        if not _XTC_IMPORTED:
            raise ImportError("XTCFile cython extension not compiled")

        with XTCFile(filename) as xtc:
            for frame in xtc:
                self.trj["coor"].append(frame.x * 10)
                self.trj["lattics"].append([frame.box[0, 0] * 10, frame.box[1, 1] * 10, frame.box[2, 2] * 10])

    def make_ensemble(self):
        ense = Ensemble("normal")
        ense.thermo = self.thermo
        ense.trj = self.trj
        ense.property = self.property
        return ense
