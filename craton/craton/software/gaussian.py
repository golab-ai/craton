import os
from io import StringIO
from typing import List, Union
from copy import deepcopy

from ..chemkit.conformation.conformation import ConformType, get_scan_curve, assign_scan_conf_type
###from ..chemistry.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
#from compuchem.chemistry.constants import const
#from compuchem.chemistry.molecule import Molecule
#from compuchem.chemistry.molecule_mechanics.result_analyze import ResultAnaly
#from compuchem.chemistry.topology import Constrain
#from compuchem.chemistry.format.assign import get_atoms_info, get_mole_info

class GauInputFile:
    _name = "Gaussain"

    def __init__(self, style=""):
        self.s = style

    def _convert(self,molecule,extra_var=None):
        text = "# hf/3-21g\n\n"
        #text += "CPY Job: %s %s\n\n" %(molecule.inchi_key,molecule.smiles)
        text += "CPY Job: %s %s\n\n" %(molecule.mole_name,molecule.smiles)
        text += "%d %d\n"%(molecule.net_charge,molecule.multi if hasattr(molecule,"multi") else 1)
        for atom in molecule.Atoms:
            text += "%s %8.4f %8.4f %8.4f\n"%(atom.element,atom.coor[0],atom.coor[1],atom.coor[2])
        text += "\n\n\n"
        return text

    def _parse(self,input_script,extra_var=None):
        script = input_script.splitlines()
        self.read_file(script)
        datas = {
            "molecule_name"  :"",
            "atom_count"     :self.mole_n,
            "elements"       : self.elem,
            "coordinates"    : self.coor,
            "formal_charge"  : [0 for i in range(self.mole_n)],
            "connectivity"   : self.connect,
            "bond_type"      : self.bond_type,
            "net_charge"     : self.net_charge,
            "multi"          : self.multi,
            }
        
        return datas

    def determine_style(self):
        if len(self.script[self.anchor_point[1] + 2]) == 1:
            self.style = "intracoor"
        else:
            self.style = "cartesian"

    def read_file(self, script):
        self.read_script(script)
        
    def read_script(self, script):
        self.script = script
        self.anchor_point = []
        for i in range(len(script)):
            if script[i].strip() == "":
                self.anchor_point.append(i)
        self.determine_style()
        if self.style == "cartesian":
            self.read_info_cartesian()
        elif self.style == "intracoor":
            self.read_info_intracoor()

    def read_info_cartesian(self):
        self.elem = []
        self.coor = []
        self.connect = []
        self.bond_type = []
        # self.constrain=[]
        for i in range(0, self.anchor_point[0]):
            if self.script[i].find("connectivity") != -1:
                connectivity = "yes"
                break
        self.net_charge = int(self.script[self.anchor_point[1] + 1].strip().split()[0])
        self.multi = int(self.script[self.anchor_point[1] + 1].strip().split()[1])
        for i in range(self.anchor_point[1] + 2, self.anchor_point[2]):
            string = self.script[i].strip().split()
            self.elem.append(string[0])
            self.coor.append([float(string[1]), float(string[2]), float(string[3])])
        self.mole_n = len(self.elem)
        if connectivity == "yes":
            for i in range(len(self.elem)):
                self.connect.append([])
                self.bond_type.append([])
            for i in range(self.anchor_point[2] + 1, self.anchor_point[3]):
                string = self.script[i].strip().split()
                if len(string) != 1:
                    a = int(string[0]) - 1
                    for j in range(1, len(string), 2):
                        b = int(string[j]) - 1
                        if float(string[j + 1]) == 1.5:
                            bn = "ar"
                        else:
                            bn = str(int(float(string[j + 1])))
                        self.connect[a].append(b)
                        self.bond_type[a].append(bn)
                        if a not in self.connect[b]:
                            self.connect[b].append(a)
                            self.bond_type[b].append(bn)

    def read_info_intracoor(self):
        pass

    def read_info(self):
        if self.style == "cartesian":
            self.read_info_cartesian()
        elif self.style == "intracoor":
            self.read_info_intracoor()

    def import_moleobj(self, m, zmatrix=None):
        self.elem = []
        self.coor = []
        self.connect = []
        self.bond_type = []
        self.charge = m.net_charge
        self.multi = 1
        self.smiles = m.smiles
        self.inchi_key = m.inchi_key
        self.mole_name = m.mole_name
        #if hasattr(m, "conform_id"):
        self.conform_id = getattr(m,"confID","None")
        #else:
        #    self.conform_id = "None"
        #if hasattr(m, "charge"):
        #    self.charge = m.charge
        if hasattr(m, "multi"):
            self.multi = m.multi
        for aa in m.Atoms:
            self.elem.append(aa.elem)
            self.coor.append(aa.coor)
        if hasattr(m, "scan_term"):
            self.scan_term = m.torsions
        if hasattr(m, "elem_set"):
            self.elem_set = m.elem_set
        if hasattr(m, "normal_elem_arr"):
            self.normal_elem_arr = m.normal_elem_arr
        if hasattr(m, "special_elem_arr"):
            self.special_elem_arr = m.special_elem_arr

        self.zmatrix = zmatrix

    def write_gjf(self, fpath, step="", qmpara=None) -> list:
        fname = os.path.basename(fpath)
        gjf_content, addcontext = StringIO(), StringIO()

        basisset_append = qmpara["basisset_append"]
        radii_append = qmpara["radii_append"]
        fpaths = []

        job_type = qmpara["job_type"]
        if job_type in ["sp", "freq", "charge", "freqcharge"]:
            qmpara["method_basisset_level"] = qmpara["sp_method_basisset_level"]
        else:
            qmpara["method_basisset_level"] = qmpara["opt_method_basisset_level"]

        if "charge_model" in qmpara:
            qmodel = qmpara["charge_model"]
        else:
            qmodel = "chelpg"


        gjf_content.write(f"%nproc={qmpara['cpucores']}\n" f"%mem=4gb\n")
        if job_type.find("charge") != -1:
            gjf_content.write(f"%chk={fname}.chk\n")

        if job_type == "scan":
            if qmpara["scan_term"][0] != "all-flexible":
                for aa in qmpara["scan_term"]:
                    addcontext.write(f" {aa}")
                    fpath += "-%d" %aa
                addcontext.write(" S")
                if "scan_setting" not in qmpara:
                    addcontext.write(" 12 30.0\n")
                else:
                    for aa in qmpara["scan_setting"]:
                        addcontext.write(f" {aa}")
                    addcontext.write("\n")
            else:
                if "keyword" in qmpara:
                    gjf_content.write("# %s opt=addred %s \n\n" %(qmpara["method_basisset_level"]," ".join(qmpara["keyword"])))
                else:
                    gjf_content.write("# %s opt=addred \n\n" % qmpara["method_basisset_level"])
                scan_setting = qmpara.get("scan_setting", None) or " 12 30.0"
                #gjf_content.write(f"CPY Job: {self.inchi_key} {self.smiles} {self.conform_id} {step}\n\n")
                gjf_content.write(f"CPY Job: {self.mole_name} {self.smiles} {self.conform_id} {step}\n\n")
                gjf_content.write(f"{self.charge} {self.multi}\n")
                if self.zmatrix is None:
                    for i in range(len(self.elem)):
                        gjf_content.write(
                            "%5s%12.4f%12.4f%12.4f\n"
                            % (
                                self.elem[i],
                                self.coor[i][0],
                                self.coor[i][1],
                                self.coor[i][2],
                            )
                        )
                else:
                    gjf_content.write(self.zmatrix.to_str())

                gjf_content.write("\n")
                
                fpaths += self.write_multi_scan_gjf(
                    fpath, gjf_content, scan_setting, qmpara["basisset_append"], start=1, step=step,
                )
                return fpaths

        # opt and charge have to be separated into two link1

        if job_type == "sp":
            gjf_content.write(f"# {qmpara['method_basisset_level']} force ")

        elif job_type == "scan":
            gjf_content.write("# %s opt=addred " % qmpara["method_basisset_level"])

        elif job_type == "opt":
            gjf_content.write(f"# {qmpara['method_basisset_level']} opt ")

        elif job_type == "charge":
            if len(self.special_elem_arr) > 0:
                gjf_content.write("# %s pop=(%s,readradii) " % (qmpara["method_basisset_level"], qmodel))
            else:
                gjf_content.write("# %s pop=%s " % (qmpara["method_basisset_level"], qmodel))

        elif job_type == "freqcharge":
            if len(self.special_elem_arr) > 0:
                gjf_content.write("# %s freq pop=(%s,readradii) " % (qmpara["method_basisset_level"], qmodel))
            else:
                gjf_content.write("# %s freq pop=%s " % (qmpara["method_basisset_level"], qmodel))

        elif job_type == "optfreq":
            gjf_content.write("# %s opt freq " % qmpara["method_basisset_level"])

        elif job_type == "optcharge":
            if len(self.special_elem_arr) > 0:
                gjf_content.write("# %s opt " % qmpara["method_basisset_level"])
            else:
                gjf_content.write("# %s opt pop=%s " % (qmpara["method_basisset_level"], qmodel))

        elif job_type == "optfreqcharge":
            if len(self.special_elem_arr) > 0:
                gjf_content.write("# %s opt freq " % qmpara["method_basisset_level"])
            else:
                gjf_content.write("# %s opt freq pop=%s " % (qmpara["method_basisset_level"], qmodel))

        elif job_type == "fixopt":
            gjf_content.write("# %s opt=addred " % qmpara["method_basisset_level"])
            if "fix_term" not in qmpara:
                raise Exception("Fixed term not found")

            for aa in qmpara["fix_term"]:
                fpath += "_%s" % round(aa, 0)
                addcontext.write(f" {aa}")
            addcontext.write(" F \n")

        else:
            raise Exception("Invalid job type: %s" % job_type)

        if "keyword" in qmpara:
            gjf_content.write("%s \n\n" %" ".join(qmpara["keyword"]))
        else:
            gjf_content.write(" \n\n")

        #gjf_content.write(f"CPY Job: {self.inchi_key} {self.smiles} {self.conform_id} {step}\n\n")
        gjf_content.write(f"CPY Job: {self.mole_name} {self.smiles} {self.conform_id} {step}\n\n")
        gjf_content.write(f"{self.charge} {self.multi}\n")
        if self.zmatrix is None:
            for i in range(len(self.elem)):
                gjf_content.write(
                    "%5s%12.4f%12.4f%12.4f\n" % (self.elem[i], self.coor[i][0], self.coor[i][1], self.coor[i][2])
                )
        else:
            gjf_content.write(self.zmatrix.to_str())

        gjf_content.write("\n")
        if addcontext.tell() > 0:
            gjf_content.write(f"{addcontext.getvalue()}\n")

        gjf_content.write(f"{basisset_append}")

        if job_type.find("charge") != -1 and len(self.special_elem_arr) > 0:
            if job_type.find("opt") != -1:
                gjf_content.write(
                    f"--Link1--\n"
                    f"%chk={fname}.chk\n"
                    f"# {qmpara['method_basisset_level']} pop=({qmodel},readradii) geom=AllCheck guess=Read\n\n"
                )
            gjf_content.write(f"{radii_append}")

        gjf_content.write("\n")

        fpath = f"{fpath}.gjf"
        with open(fpath, "w") as outf:
            outf.write(gjf_content.getvalue())
        fpaths.append(fpath)

        return fpaths


    def write_multi_scan_gjf(
        self, fpath_base, gjf_content, scan_setting, basisset_append, start=1, step=""
    ) -> List[str]:
        fpaths = []
        for scan in self.scan_term:
            addcontext = StringIO()
            fpath = f"{fpath_base}_Scan"
            for s in scan:
                addcontext.write(f" {s + start} ")
                fpath = f"{fpath}_{s + start}"
            addcontext.write("S ")
            for aa in scan_setting:
                addcontext.write(f"{aa} ")
            addcontext.write("\n")
            content = StringIO()
            content.write(gjf_content.getvalue())
            content.write(addcontext.getvalue())
            content.write("\n")
            content.write(basisset_append)
            content.write("\n")
            fpath = f"{fpath}.gjf"
            fname = os.path.basename(fpath)
            with open(fpath, "w") as outf:
                outf.write(content.getvalue())
            fpaths.append(fpath)

        return fpaths


class GauOutputLink1:
    """
    Parse Gaussian 09 log file.
    """
    

    def __init__(self, lines, read_optimizing_flag=False):
        self.read_optimizing_flag = read_optimizing_flag
        self.script = lines[:]
        self.command: str = ""
        self.title: str = ""
        self.delimiter = ""
        self.result_list: [str] = []
        self.style: str = ""

        self.qm_method: str = ""
        self.qm_basis_set: str = ""
        self.energy: float = 0.0
        self.dipole: [float] = []
        self.mole_n: int = 0
        self.net_charge: int = 0
        self.multi: int = 1
        self.elem: [str] = []
        self.coor: [[float]] = []
        self.charge: {str: [float]} = {}
        self.extra_data: {str: str} = {}
        self.internal_coor: {str: []} = {}  # redundant internal coordinate

        # for vibrational analysis
        self.freq: [float] = []
        self.force: [[float]] = []
        self.hessian: [float] = []
        self.thermochem = {}  # not used for now

        # for constrained opt and torsion scan
        self.constrained_internal_dof = []  # ['D12', ...]

        # for 1-D torsion scan
        self.scan_energy: [float] = []
        self.scan_coor: [[[float]]] = []
        self.scan_internal_coor: [{str: []}] = []

        self.search_method = "unknow"
        self.confID = None
        self.inchi_key = None
        self.mole_name = None
        self.smiles = None
        self.force_flag = False
        self.coor_search_tag = "Standard orientation:"
        for line in lines:
            if line.find("Input orientation:") != -1:
                self.force_flag = True
                self.coor_search_tag = "Input orientation:"
                break

        self._content_parse()

        self.error_log = None
        try:
            func = self.__Func_ID[self.style]
        except:
            self.error_log = "Unknown Gaussian job style: %s" % self.style
            #raise Exception("Unknown Gaussian job style: %s" % self.style)

        func(self)


        self.conform_type = self._CONFORM_TYPE[self.style]

    def _content_parse(self):
        r"""
        First get the command title and line number of final result.
        Every line may or may not starts with a space.

        Command startswith # and is wrapped by two dash lines
        -------------------------------------------
        # b3lyp def2svp opt freq pop=chelpg
        -------------------------------------------

        Title is also wrapped by two dash lines
        -------------------------------------------
        Title may contain any character
        -------------------------------------------

        Final results starts with 1\1\ or 1|1| and ends with @ and occupies multiple lines.
        Below the final result it is a empty line.
        For Windows, the terms in result are separated by '|';
        For Linux, the terms in result are separated by '\'.
        Title is also included in final result, which may contains separator also.
        """

        lines = self.script[:]
        n_line = len(lines)
        i_last = n_line - 1
        while i_last > 0:
            if lines[i_last].strip() != "":
                break
            i_last -= 1
        if not lines[i_last].startswith("Normal termination"):
            self.error_log = "Gaussian job not terminate normally"
            return
            ####raise Exception("Gaussian job not terminate normally")

        ### command and title
        i_cmd = -1
        for i, line in enumerate(lines):
            if line.startswith("*" * 12) and lines[i + 1].startswith("Gaussian"):
                i_cmd = i + 1

        self.command = ""
        for i in range(i_cmd + 1, len(lines)):
            if lines[i].startswith("-" * 12):
                i_cmd = i + 1
                break

        while True:
            if lines[i_cmd].startswith("-" * 12):
                break
            self.command += lines[i_cmd]
            i_cmd += 1
        self.title = ""
        for i in range(i_cmd + 1, len(lines)):
            if lines[i].startswith(("-" * 12)):
                i_cmd = i + 1
                break

        while True:
            if lines[i_cmd].startswith("-" * 12):
                break
            self.title += lines[i_cmd]
            i_cmd += 1
        ###

        title_ss = self.title.split(":")
        if title_ss[0] == "CPY Job":
            sub_title_ss = title_ss[1].split()
            self.search_method = sub_title_ss[-1].strip()
            #self.inchi_key = sub_title_ss[0].strip()
            self.mole_name = sub_title_ss[0].strip()
            self.smiles = sub_title_ss[1].strip()
            if len(sub_title_ss) > 2:
                if sub_title_ss[2] != "None":
                    self.confID = sub_title_ss[2]


        for i in range(i_last, 1, -1):
            if (lines[i].strip() == "" or lines[i].strip() == "The archive entry for this job was punched.") and lines[i - 1].endswith("@"):
                i_result_end = i - 1
                break
        else:
            # assume it's a scan calculation if no result block
            self.style = "Scan"

            i_atom = -1
            for ii, line in enumerate(lines):
                if line.startswith("Charge ="):
                    words = line.split()
                    self.net_charge, self.multi = int(words[2]), int(words[5])
                    i_atom = ii + 1
                    break

            for ii, line in enumerate(lines):
                if line.startswith("Standard basis:"):
                    self.qm_basis_set = line.split()[2].strip()
                    break
            for ii, line in enumerate(lines):
                if line.startswith("SCF Done:"):
                    self.qm_method = line.split()[2].split("(")[1].strip(")")
                    break

            for ii in range(i_atom, len(lines)):
                if lines[ii] == "":
                    break
                words = lines[ii].split()
                self.elem.append(words[0])

            self.mole_n = len(self.elem)
            self.coor = [[0, 0, 0] for elem in self.elem]

            return

        if "|" in lines[i_result_end] or "|" in lines[i_result_end - 1]:
            self.delimiter = "|"
        else:
            self.delimiter = "\\"

        for i in range(i_result_end, 1, -1):
            if lines[i].strip().startswith(f"1{self.delimiter}1{self.delimiter}"):
                i_result_start = i
                break
        else:
            raise Exception("Result block not found")

        str_result = "".join(lines[i_result_start : i_result_end + 1])
        if self.delimiter in self.title:
            str_result = str_result.replace(f"{self.delimiter}{self.title}{self.delimiter}", f"{self.delimiter * 2}")

        self.result_list = str_result.split(self.delimiter)
        self.style, self.qm_method, self.qm_basis_set = self.result_list[3:6]
        self.net_charge, self.multi = map(int, self.result_list[15].split(","))

        self.mole_n = 0
        while True:
            line = self.result_list[16 + self.mole_n]
            if "," not in line:
                break
            words = line.split(",")
            self.elem.append(words[0])
            self.coor.append(list(map(float, words[-3:])))
            self.mole_n += 1

        for result in self.result_list[16 + self.mole_n :]:
            if result.find("=") > -1:
                i_cmd, v = result.split("=")
                self.extra_data[i_cmd] = v

    def _read_redundant_internal_coor(self):
        """
        Get redundant internal coordinate after optimization
        Useful for determining bond forming and breaking
        """
        i = 0
        while True:
            if i >= len(self.script) - 1:
                break
            if self.script[i].find("Optimized Parameters") != -1:
                i += 5
                internal_coor = {}
                while True:
                    if self.script[i].startswith("-" * 12):
                        break
                    ss = self.script[i].split()
                    vv = [int(ii) for ii in ss[2][2:-1].split(",")]
                    vv.append(float(ss[3]))
                    key = ss[1]
                    internal_coor[key] = vv
                    i += 1
                if self.style == "Scan":
                    self.scan_internal_coor.append(internal_coor)
                else:
                    self.internal_coor = internal_coor
            i += 1

    def _read_charge(self):
        i = 0
        charge = {}
        while True:
            if i >= len(self.script) - 1:
                break
            elif self.script[i].startswith("Mulliken charges:"):
                arr = []
                for j in range(i + 2, i + 2 + self.mole_n):
                    arr.append(float(self.script[j].split()[2]))
                    i = j
                charge["mulliken"] = arr
            elif self.script[i].startswith("APT charges:"):
                arr = []
                for j in range(i + 2, i + 2 + self.mole_n):
                    arr.append(float(self.script[j].split()[2]))
                    i = j
                charge["apt"] = arr
            elif self.script[i].startswith("ESP charges:"):
                arr = []
                for j in range(i + 2, i + 2 + self.mole_n):
                    arr.append(float(self.script[j].split()[2]))
                    i = j
                charge["esp"] = arr
            i += 1
        self.charge = charge

    def _read_frequencies(self):
        for i in range(0, len(self.script)):
            if self.script[i].find("Frequencies --") != -1:
                ss = self.script[i].split()
                for s in ss[2:]:
                    freq = float(s)
                    self.freq.append(freq)

    def _read_forces(self):
        from ..chem.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
        forces = []
        for i in range(len(self.script) - 1, 0, -1):
            if self.script[i].find("Forces (Hartrees/Bohr)") != -1:
                i_force = i + 3
                break
        else:
            raise Exception("Force section not found")

        for i in range(i_force, len(self.script)):
            if self.script[i].startswith("-" * 12):
                break
            words = self.script[i].split()
            forces.append([float(i) * HF_TO_KCAL_PER_MOL / BOHR_TO_ANGSTROM for i in words[2:]])
        # change by CFL
        # force 不用二维数组，而是用一维数组
        self.force = [f for force in forces for f in force]

    def _read_thermochemistry(self):
        from ..chem.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
        thermochem = {}
        for i in range(0, len(self.script)):
            if self.script[i].find("Zero-point correction=") != -1:
                thermochem["ZPE"] = float(self.script[i].split("=")[1].split("(")[0]) * HF_TO_KCAL_PER_MOL
            if self.script[i].find("Sum of electronic and zero-point Energies=") != -1:
                thermochem["total_e"] = float(self.script[i].split("=")[1]) * HF_TO_KCAL_PER_MOL
            if self.script[i].find("Sum of electronic and thermal Energies=") != -1:
                thermochem["thermal_e"] = float(self.script[i].split("=")[1]) * HF_TO_KCAL_PER_MOL
            if self.script[i].find("Sum of electronic and thermal Enthalpies=") != -1:
                thermochem["enthalpy"] = float(self.script[i].split("=")[1]) * HF_TO_KCAL_PER_MOL
            if self.script[i].find("Sum of electronic and thermal Free Energies=") != -1:
                thermochem["free_energy"] = float(self.script[i].split("=")[1]) * HF_TO_KCAL_PER_MOL
            if self.script[i].find("E (Thermal)             CV                S") != -1:
                ss = self.script[i + 2].split()
                thermochem["capacity"] = float(ss[2])
                thermochem["entropy"] = float(ss[3])
        self.thermochem = thermochem

    def old_read_scan_coor_energy(self):
        from ..chem.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
        i = 0
        self.scan_coor = []
        self.scan_energy = []
        while 1:
            if i >= len(self.script) - 1:
                break
            if self.script[i].find("Standard orientation:") != -1:
                coor = []
                energy = 0
                i += 5
                for j in range(i, i + self.mole_n):
                    coor.append([float(a) for a in self.script[j].split()[3:]])
                    i = j
            if self.script[i].startswith("SCF Done:"):
                energy = float(self.script[i].split()[4]) * HF_TO_KCAL_PER_MOL

            if self.script[i].find("Optimized Parameters") != -1:
                self.scan_coor.append(coor)
                self.scan_energy.append(energy)
            i += 1

    def _read_scan_coor_energy(self):
        i = 0
        self.scan_coor = []
        self.scan_energy = []
        self.scan_force = []
        if self.read_optimizing_flag:
            self.optimizing_coor = []
            self.optimizing_energy = []
            self.optimizing_force = []
        optimizing_range = [
            0,
        ]
        for i in range(len(self.script)):
            if self.script[i].find("Optimized Parameters") != -1:
                optimizing_range.append(i)

        for nn in range(1, len(optimizing_range)):
            aa, bb, cc = self._read_optimizing_coor_energy(self.script[optimizing_range[nn - 1] : optimizing_range[nn]])
            self.scan_coor.append(aa[-1])
            self.scan_energy.append(bb[-1])
            self.scan_force.append(cc[-1])
            if self.read_optimizing_flag:
                self.optimizing_coor.append(aa[:-1])
                self.optimizing_energy.append(bb[:-1])
                self.optimizing_force.append(cc[:-1])

    def _read_optimizing_coor_energy(self, this_script):
        from ..chem.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
        i = 0
        optimizing_coor = []
        optimizing_energy = []
        optimizing_force = []
        force = []
        coor = []
        while 1:
            if i >= len(this_script) - 1:
                break
            ###if this_script[i].find("Standard orientation:") != -1:
            if this_script[i].find(self.coor_search_tag) != -1:
                coor = []
                i += 5
                for j in range(i, i + self.mole_n):
                    coor.append([float(a) for a in this_script[j].split()[3:]])
                    i = j
            if this_script[i].startswith("SCF Done:"):
                energy = float(this_script[i].split()[4]) * HF_TO_KCAL_PER_MOL
            if this_script[i].find("Forces (Hartrees/Bohr)") != -1:
                i_force = i + 3
                force = []
                for i in range(i_force, len(this_script)):
                    if this_script[i].startswith("-" * 12):
                        optimizing_coor.append(coor)
                        optimizing_energy.append(energy)
                        optimizing_force.append([f for force in force for f in force])
                        break
                    words = this_script[i].split()
                    force.append([float(i) * HF_TO_KCAL_PER_MOL / BOHR_TO_ANGSTROM for i in words[2:]])
            i += 1
        return optimizing_coor, optimizing_energy, optimizing_force

    def _get_info_SP(self):
        from ..chem.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
        self._read_charge()
        #self.search_method = self.title.split()[-1].strip()
        self.energy = float(self.extra_data["HF"]) * HF_TO_KCAL_PER_MOL
        self.dipole = list(map(float, self.extra_data["Dipole"].split(",")))
        if self.command.find("opt") != -1 and self.style != "Scan" and self.read_optimizing_flag:
            this_script = []
            for line in self.script:
                if line.find("@") != -1:
                    break
                else:
                    this_script.append(line)
            self.optimizing_coor, self.optimizing_energy, self.optimizing_force = self._read_optimizing_coor_energy(
                this_script
            )
            self.optimizing_coor = self.optimizing_coor[:-1]
            self.optimizing_energy = self.optimizing_energy[:-1]
            self.optimizing_force = self.optimizing_force[:-1]

    def _get_element(self):
        for ii, line in enumerate(self.script):
            if line.find("Multiplicity") != -1:
                break
        self.elements = []
        for ii in range(ii+1, len(self.script)):
            if self.script[ii].strip() == "":
                break
            else:
                self.elements.append(self.script[ii].split()[0].strip())

    def _get_info_Force(self):
        self._get_info_SP()
        self._read_forces()

    def _get_info_Fopt(self):
        self._get_info_SP()
        self._read_forces()
        self._read_redundant_internal_coor()
        
        Grad_n = 0
        for line in self.script:
            if Grad_n == 2:
                break
            if line.find("GradGradGradGradGradGradGrad") != -1:
                Grad_n += 1
                continue
            if line.find("Frozen") != -1:
                dof = line.split()[1]
                self.constrained_internal_dof.append(dof)

    def _get_info_Freq(self):
        from ..chem.constants import HF_TO_KCAL_PER_MOL, BOHR_TO_ANGSTROM
        self._get_info_SP()

        self._read_frequencies()
        self._read_thermochemistry()
        self._read_forces()
        self.hessian = [
            float(i) * HF_TO_KCAL_PER_MOL / BOHR_TO_ANGSTROM**2 for i in self.result_list[-6].split(",")
        ]

    def _get_info_Scan(self):
        self._read_redundant_internal_coor()
        self._read_scan_coor_energy()
        #self.search_method = self.title.split()[-1].strip()
        i_start = 0
        for i, line in enumerate(self.script):
            if line.find("Initial Parameters") != -1:
                i_start = i
                break
        for i in range(i_start + 5, len(self.script)):
            if self.script[i].startswith("-" * 12):
                break
            if self.script[i].find("Scan") != -1:
                dof = self.script[i].split()[1]
                self.constrained_internal_dof.append(dof)

    def get_data(self):

        def get_constrain(data_dict, ss, ns=0):
            data_dict["constrain"] = []
            data_dict["scan_term"] = []
            for tt in ss:
                tmp = []
                for i in range(len(tt) - 1):
                    tmp.append(tt[i] - ns)
                data_dict["constrain"].append([tmp, tt[-1]])
                data_dict["scan_term"].append(tmp)
            return data_dict

        if self.error_log is not None:
            print(self.error_log)
            return []

        data_dict = {}
        data_dict["force_flag"] = self.force_flag
        data_dict["search_method"] = self.search_method
        data_dict["mole_n"] = self.mole_n
        data_dict["net_charge"] = self.net_charge
        data_dict["multi"] = self.multi
        data_dict["qm_method"] = self.qm_method
        data_dict["qm_basis_set"] = self.qm_basis_set
        data_dict["energy"] = self.energy
        data_dict["coordinates"] = self.coor
        data_dict["elements"] = self.elem
        data_dict["conform_type"] = self.conform_type
        #if self.inchi_key is not None:
        if self.mole_name is not None:
            data_dict["inchi_key"] = self.mole_name
            data_dict["mole_name"] = self.mole_name
        if self.confID is not None:
            data_dict["confID"] = self.confID

        if self.smiles is not None:
            data_dict["smiles"] = self.smiles

        if self.style == "Force":
            data_dict["force"] = self.force
        
        if self.style == "FOpt":
            data_dict["force"] = self.force

        if self.style == "Freq":
            data_dict["force"] = self.force
            data_dict["hessian"] = self.hessian
            data_dict["freq"] = self.freq
            data_dict["thermochem"] = self.thermochem

        if self.style != "Scan":
            data_dict["charge"] = self.charge
            data_dict["dipole"] = self.dipole
            if self.style == "FOpt":
                data_dict = get_constrain(data_dict, [self.internal_coor[dof] for dof in self.constrained_internal_dof], ns=1)
            if hasattr(self, "optimizing_energy"):
                data_dit_arr = [deepcopy(data_dict)]
                for attr in ["charge", "dipole", "hessian", "freq", "thermochem"]:
                    if attr in data_dict:
                        del data_dict[attr]
                data_dict["search_method"] = data_dict["search_method"] + "optimizing"
                data_dict["conform_type"] = self._CONFORM_TYPE["Optimizing"]
                data_dict["mole_n"] = self.mole_n
                for i in range(len(self.optimizing_coor)):
                    data_dict["coordinates"] = self.optimizing_coor[i]
                    data_dict["energy"] = self.optimizing_energy[i]
                    data_dict["force"] = self.optimizing_force[i]
                    if self.style == "FOpt":
                        data_dict = get_constrain(data_dict, [self.internal_coor[dof] for dof in self.constrained_internal_dof], ns=1)
                    data_dit_arr.append(deepcopy(data_dict))
                return data_dit_arr
            return [data_dict]

        else:
            data_dict_arr = []
            for i in range(len(self.scan_energy)):
                data_dict["energy"] = self.scan_energy[i]
                data_dict["coordinates"] = self.scan_coor[i]
                data_dict["force"] = self.scan_force[i]
                
                data_dict = get_constrain(data_dict, [self.scan_internal_coor[i][dof] for dof in self.constrained_internal_dof], ns=1)
                data_dict_arr.append(deepcopy(data_dict))

            scan_curve = get_scan_curve(data_dict_arr)
            assign_scan_conf_type(scan_curve)

            if hasattr(self, "optimizing_energy"):
                data_dict["search_method"] = data_dict["search_method"] + "optimizing"
                data_dict["conform_type"] = self._CONFORM_TYPE["Optimizing"]
                for i in range(len(self.optimizing_energy)):
                    for j in range(len(self.optimizing_energy[i])):
                        data_dict["coordinates"] = self.optimizing_coor[i][j]
                        data_dict["energy"] = self.optimizing_energy[i][j]
                        data_dict["force"] = self.optimizing_force[i][j]
                        get_constrain(
                            data_dict, [self.scan_internal_coor[i][dof] for dof in self.constrained_internal_dof], ns=1
                        )
                        data_dict_arr.append(deepcopy(data_dict))
            return data_dict_arr

    __Func_ID = {
        "SP": _get_info_SP,
        "Force": _get_info_Force,
        "FOpt": _get_info_Fopt,
        "Freq": _get_info_Freq,
        "Scan": _get_info_Scan,
    }

    _CONFORM_TYPE = {
        "SP": ConformType.SINGLE_POINT,
        "Force": ConformType.SINGLE_POINT,
        "FOpt": ConformType.LOCAL_MINIMUM,
        "Freq": ConformType.LOCAL_MINIMUM,
        "Scan": ConformType.SCAN,
        "Optimizing": ConformType.OPTIMIZING,
    }


class GauOutputFile:
    """
    Parse Gaussian 09 log file
    """
    def __init__(self, string=None):

        if string is not None:
            self.read_file(string)

    def _parse(self,input_script,extra_var=None):
        self.read_optimizing = extra_var["read_optimizing"] if extra_var is not None and "read_optimizing" in extra_var \
                               else False
        self.element_flag = extra_var["element_flag"] if extra_var is not None and "element_flag" in extra_var \
                               else False


        self.link1s = []
        lines = [line[1:] for line in input_script.splitlines()]
 
        if not lines[-1].startswith("Normal termination"):
            #raise Exception("Gaussian job not terminate normally")
            return []

        idx_link1 = [i for i, line in enumerate(lines) if line == "Initial command:"]
        idx_link1.append(len(lines))
        for i, j in zip(idx_link1[:-1], idx_link1[1:]):
            self.link1s.append(GauOutputLink1(lines[i:j],read_optimizing_flag=self.read_optimizing))

        if len(self.link1s) < 1:
            ####raise Exception("Invalid Gaussian log file")
            return []

        # TODO Only consider charge for now
        if len(self.link1s) > 1:
            self.link1s[0].charge.update(self.link1s[1].charge)

        ####self.read_file(input_script)
        return self.get_data()

    def read_file(self, string):
        self.link1s = []
        lines = [line[1:] for line in string.splitlines()]
        idx_link1 = [i for i, line in enumerate(lines) if line == "Initial command:"]
        idx_link1.append(len(lines))
        for i, j in zip(idx_link1[:-1], idx_link1[1:]):
            self.link1s.append(GauOutputLink1(lines[i:j],read_optimizing_flag=self.read_optimizing))

        if len(self.link1s) < 1:
            raise Exception("Invalid Gaussian log file")

        # TODO Only consider charge for now
        if len(self.link1s) > 1:
            self.link1s[0].charge.update(self.link1s[1].charge)

    def get_data(self):
        return self.link1s[0].get_data()
