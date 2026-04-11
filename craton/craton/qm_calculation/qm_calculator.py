import os
from typing import List

#from ..software.gaussian import GauInputFile
from ..chem.elements import Element


def parse_qm_setting(file):
    qmpara = {}
    with open(file) as f:
        for line in f:
            if line[0] != "#" and line.strip() != "":
                if "#" in line:
                    ss = line[: line.index("#")].strip().split(":")
                else:
                    ss = line.strip().split(":")
                qmpara[ss[0].strip()] = ss[1].strip().split()
    return qmpara


class QmCalculator:
    def __init__(self, path="./"):
        self.path = path
    def special_elem_radii(self, special_elem_arr, qmpara: dict) -> dict:
        """
        处理Br I等元素的原子半径的问题。
        """
        qmpara["radii_append"] = ""
        for ee in special_elem_arr:
            qmpara["radii_append"] += "%s %s\n" % (ee, Element.get(ee).vdw_radius)
        qmpara["radii_append"] += "\n"
        return qmpara

    def create_qm_input_files(self, m, qmpara=None, step="", index=None, zmatrix=None,fpath_pre=""):
        """
        Returns List of .gjf file paths
        """
        from ..software.gaussian import GauInputFile
        m.elem_set = set([aa.elem for aa in m.Atoms])
        m.normal_elem_arr = list(m.elem_set & {"H", "C", "O", "N", "S", "P", "F", "Cl", "Si", "B"})
        m.special_elem_arr = list(m.elem_set - set(m.normal_elem_arr))
        
        if qmpara is None:
            qmpara={"qmjobs":["sp"]}
        # qm_input_file_list = []
        engine = qmpara["engine"] if "engine" in qmpara else "g09"
        
        if engine == "g09":
            qmengineobj = GauInputFile()
            qmengineobj.import_moleobj(m, zmatrix=zmatrix)
        fpaths = []
        for job in qmpara["qmjobs"]:
            qmpara.update({"job_type": job})
            job_type = ""
            if job in ["sp", "freq", "charge", "freqcharge"]:
                job_type = "sp"
            elif job in ["opt", "optfreq", "optcharge", "optfreqcharge", "scan", "fixopt"]:
                job_type = "opt"
            if job_type != "":
                if f"{job_type}_method_basisset" in qmpara:
                    for level in qmpara[f"{job_type}_method_basisset"]:
                        qmpara.update(
                            {
                                f"{job_type}_method_basisset_level": level,
                                "basisset_append": "",
                                "radii_append": "",
                            }
                        )

                        if len(m.special_elem_arr) != 0:
                            qmpara = self.special_elem_radii(m.special_elem_arr, qmpara)
                        level1 = level.replace("/", "_").replace("(", "").replace(")", "").replace("*", "").replace(",","")
                        #filename = "_".join([m.inchi_key, job, level1])
                        filename = "_".join([m.mole_name, job, level1])
                        filename += fpath_pre
                        if index is not None:
                            filename += f"_{index}"
                        full_path = os.path.join(self.path, filename)
                        fpaths += qmengineobj.write_gjf(full_path, step=step, qmpara=qmpara)
                else:
                    qmpara[f"{job_type}_method_basisset_level"] = "hf"
                    qmpara["basisset_append"] = ""
                    qmpara["radii_append"] = ""
                    #filename = "_".join([m.inchi_key, job, "hf"])
                    filename = "_".join([m.mole_name, job, "hf"])
                    filename += fpath_pre
                    if index is not None:
                        filename += f"_{index}"
                    full_path = os.path.join(self.path, filename)
                    fpaths += qmengineobj.write_gjf(full_path, step=step, qmpara=qmpara)
        #return fpaths
