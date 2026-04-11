import gzip
import json
import multiprocessing as mp
from pathlib import Path
from copy import deepcopy
from typing import *

import munch
import numpy as np

from ...utils import logger
#from ...mm_calculator.mm.mm_calculator import MMCalc
#from ...mm_calculator.optimizer import optimize
####from ._fitting import st 
from .result_analyze import ResultAnaly
from ...chemkit.conformation.conformation import ConformType


def initi_mole_info_count(molecules) -> Dict:
    """
    统计用来拟合的分子的信息，包括分子数目，构象的数目，势能面的数目，local minimum的数目等
    输入：
        moles: List[Molecule] 拟合分子信息
    输出：
        datas: Dict, 包含上面信息的字典
    """
    datas = {
        "total": {
            "pes_number": [
                0,
                0,
            ],
            "conformations": 0,
            "miniminum_number": 0,
            "molecules_number": 0,
        },
    }
    inchi_key_dict = {}
    for mol in molecules:
        if mol.inchi_key not in inchi_key_dict:
            datas["total"]["molecules_number"] += 1
            datas[mol.inchi_key] = {"conformations": 0, "pes_number": [0, 0], "miniminum_number": 0}
            inchi_key_dict[mol.inchi_key] = []
        datas[mol.inchi_key]["conformations"] += 1
        datas["total"]["conformations"] += 1
        if hasattr(mol, "constrain"):
            #scan_term_name = "-".join([str(n) for n in mol.scan_term[0]])
            scan_term_name = mol.constrain[0].name
            if scan_term_name not in inchi_key_dict[mol.inchi_key]:
                datas[mol.inchi_key]["pes_number"][0] += 1
                datas[mol.inchi_key]["pes_number"][1] += 1
                datas["total"]["pes_number"][0] += 1
                datas["total"]["pes_number"][1] += 1
                inchi_key_dict[mol.inchi_key].append(scan_term_name)
            else:
                datas[mol.inchi_key]["pes_number"][1] += 1
                datas["total"]["pes_number"][1] += 1
        else:
            datas[mol.inchi_key]["miniminum_number"] += 1
            datas["total"]["miniminum_number"] += 1
    return datas

def initi_ff_info_count(this_ff: Dict) -> Dict[str, Any]:
    """
    统计用来拟合的力场的信息，包括每种类型有多少参数等
    输入：
        this_ff: Dict, 一些统计上的数据，一般由validation方法生成
    输出：
        datas: Dict, 包含上面信息的字典
    """
    __tt = {"atomtype", "bondterm", "angleterm", "dihedralterm", "improperterm"}
    datas = dict()
    fix_datas = dict()
    for key, value in this_ff.items():
        if key in __tt:
            for term, para in value.items():
                if "isfitting" in para:
                    if key not in datas:
                        datas[key] = {}
                    datas[key][term] = [round(v, 4) for v in para["para"]]
                else:
                    if key not in fix_datas:
                        fix_datas[key] = {}
                    fix_datas[key][term] = [round(v, 4) for v in para["para"]]
    datas["fix_datas"] = fix_datas
    return datas


def _validation(
    moles,
    target_moles,
    terms: List[str] = [
        "rmse",
        "pes",
        "charge",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Pair1n",
        "energy",
        "force",
        "hessian",
        "freq",
    ],
    is_fitting_param: List[str] = [],
    param_validation_flag: bool = False,
    save_path: str = "./",
    write_gjf_energy_threshold: float = -1,
) -> Tuple[Dict, Dict, Optional[Dict]]:
    """
    对计算或拟合结果进行验证。验证的对象包含键长、键角、二面角，能量、力，频率、势能面等方面
    The attributes (energy, force, hessian, freq) of moles will be modified

    输入：
        moles: List[Molecule], 拟合的分子集，进行力场参数的拟合和结构的优化
        target_molecules: List[Molecule], 目标分子集，记录了QM构象和性质数据
        terms: List[str], 要验证的性质，如pes, Bonds, Angles, Dihedrals, energy, freq等
        is_fitting_param：List[str], 记录拟合的参数
        generate_figure_flag: True or False, 是否生成检查的图片
        para_validation_flag: True or False，是否进行参数的验证
        save_path: path, 文件生成的路径
    输出：
        aa: Dict, 包含有验证数据的字典
        bb: Dict, 验证数据的统计信息，如平均值，RMSE,R2等
        param_validation_data: Dict, 参数验证信息字典
    """
    from ...mm_calculator.mm.mm_calculator import MMCalc
    # recalculate properties
    calc = MMCalc("intra")

    if "pes" in terms or "energy" in terms:
        energy = calc.intra_mole_energy(moles)
        for i in range(len(energy)):
            moles[i].energy = energy[i]["total"]

    if "force" in terms:
        moles_force = [m for m in moles if hasattr(m, "force")]
        force = calc.mole_force(moles_force)
        for i in range(len(moles_force)):
            moles_force[i].force = force[i]

    if "hessian" in terms or "freq" in terms:
        moles_freq = [m for m in moles if hasattr(m, "freq")]
        hessian = calc.mole_hessian(moles_freq)
        for i in range(len(moles_freq)):
            moles_freq[i].hessian = hessian[i]
        freq = calc.mole_freq(moles_freq)
        for i in range(len(moles_freq)):
            moles_freq[i].freq = freq[i]

    result = ResultAnaly("fitting", save_path=save_path)
    aa, bb, bad_moles = result.get_ff_fitting_result(
        moles,
        target_moles,
        terms=terms,
        isfittingparas=is_fitting_param,
        energy_diff_threshold=write_gjf_energy_threshold,
    )
    if bad_moles != []:
        moles_qm, moles_mm = zip(*bad_moles)
        from ...chem import FormatMolecule as FM
        Path(f"{save_path}/bad_moles").mkdir(exist_ok=True)
        FM._convert(moles_qm,otype="gjf",opath= f"{save_path}/bad_moles",)
        FM._convert(moles_mm,otype="gjf",opath=f"{save_path}/bad_moles")
        #MMCalc(moles_qm).save_gjf(save_path + "/bad_moles", zmatrix=True, merge=True, suffix="_QM")
        #MMCalc(moles_mm).save_gjf(save_path + "/bad_moles", zmatrix=True, merge=True, suffix="_MM")

    param_validation_data = None
    if param_validation_flag:
        param_validation_data = result.get_para_fitting_result(moles, target_moles, is_fitting_params=is_fitting_param)
    #if generate_figure_flag:
    #    result.show_figure_fitting(aa, bb)
    #    if param_validation_flag:
    #        result.show_figure_violin(param_validation_data)

    return (aa, bb, param_validation_data)

def find_pseudo_fitting(this_ff: Dict, pes_data: Dict):
    """
    检查没有scan数据，而被拟合的dihedral参数，将其标记为pseudo_fitting。
    `this_ff` will be modified

    输入：
        this_ff: Dict, 力场参数
        pes_data: Dict, 势能面的数据。由validation方法生成
    """
    scanned_dihe = []
    for inchi_key, bb in pes_data.items():
        for scan, bbb in bb.items():
            for name in bbb[2][1]:
                scanned_dihe.append(name)
                scanned_dihe.append("$".join(reversed(name.split("$"))))

    fitting_dihe = [
        name for name, para in this_ff["dihedralterm"].items() if para["tag"] == "isfitting" and para["isfitting"]
    ]
    pseudo_dihe = set(fitting_dihe).difference(set(scanned_dihe))

    for name in pseudo_dihe:
        # ignore endocyclic torsion
        if "@" in name.split("$")[1]:
            continue
        this_ff["dihedralterm"][name]["tag"] = "pseudofitting"


def intra_fitting_validation(
    target_moles,
    this_ff: Dict,
    init_moles: Dict,
    init_param: Dict,
    output_dir: str,
    optimizer: str = "numpy",
    #figure_flag: bool = True,
    hessian_flag: bool = True,
    fitting_info={},
) -> Dict[str, Any]:
    """
    分子内参数validation，包括energy,force, hessian, frequency, PES, charge, bond, angle, dihedrals, pairs, et al.
    输入：
        target_molecules: 目标分子，通常是QM计算出来的构象
        this_ff: Dict,字典形式记录的力场参数,标记了需要拟合的参数
        init_moles: Dict, 参与拟合 的分子的信息
        init_param: Dict, 力场参数的信息


        moles_ts: List[Molecule], 目标分子集，记录了QM构象和性质数据
        this_ff: Dict,字典形式记录的力场参数,标记了需要拟合的参数
        init_moles: Dict, 参与拟合 的分子的信息
        init_param: Dict, 力场参数的信息
        fitting_info: Dict, 拟合过程的一些信息
        out_put_dir: path, 结果文件保存的路径
        optimizer: str, the optimizer used when optimizing
        figure_flag: bool, 是否生成图片
    输出：
        results_data: Dict, same as results.json.gz
    """

    # separate conformations into groups w/wo optimization
    moles_ts_opt = [
        mol for mol in target_moles if mol.conform_type in ConformType.TORSION_SCAN_TYPES + [ConformType.LOCAL_MINIMUM]
    ]
    moles_ts_stretch = [mol for mol in target_moles if mol.conform_type == ConformType.STRETCH]
    moles_opt = [deepcopy(mol) for mol in moles_ts_opt]
    moles_stretch = [deepcopy(mol) for mol in moles_ts_stretch]

    #logger.debug(
    #    f"Validation: Number of total/opt/stretch conformations: "
    #    f"{len(target_moles)}/{len(moles_opt)}/{len(moles_stretch)}"
    #)

    is_fitting_params = []
    for term, v0 in this_ff.items():
        if term in [
            "atomtype",
            "bondterm",
            "angleterm",
            "dihedralterm",
            "improperterm",
        ]:
            for para, v1 in v0.items():
                if v1["tag"] == "isfitting":
                    is_fitting_params.append(para)

    # compare force on same structures from QM
    data_f, stat_f, _ = _validation(
        moles_opt + moles_stretch,
        moles_ts_opt + moles_ts_stretch,
        is_fitting_param=is_fitting_params,
        #generate_figure_flag=figure_flag,
        param_validation_flag=False,
        save_path=output_dir,
        terms=["force"],
    )
    from ...mm_calculator.optimizer import optimize
    # compare other properties on separately optimized structures
    logger.info(f"Validation: Optimizing {len(moles_opt)} conformations ...")
    moles_opt = optimize(moles_opt, optimizer=optimizer)
    #mm = MMCalc([mol for mol in moles_ts_opt + moles_opt if mol.conform_type == ConformType.LOCAL_MINIMUM])
    
    #mm.save_gjf(dir=output_dir, zmatrix=True, merge=True)

    terms = ["energy", "pes", "charge", "Bonds", "Angles", "Dihedrals", "rmse", "Pair1n"]
    if hessian_flag:
        terms += ["hessian", "freq"]
    data, statics, param_val_data = _validation(
        moles_opt,
        moles_ts_opt,
        is_fitting_param=is_fitting_params,
        #generate_figure_flag=figure_flag,
        param_validation_flag=True,
        save_path=output_dir,
        terms=terms,
        write_gjf_energy_threshold=5.0,
    )
    data["force"] = data_f["force"]
    statics["force"] = stat_f["force"]
    find_pseudo_fitting(this_ff, data.get("pes", {}))

    results_datas = {
        "initi_info": {
            "init_moles": init_moles,
            "init_para": init_param,
        },
        "process": fitting_info,
        "results": {
            "properties": data,
            "property_statics":statics,
            "parameters": param_val_data,
        },
    }

    # output as gzip file
    with gzip.open(f"{output_dir}/fitting_results.json.gz", "wt") as f:
        f.write(json.dumps(results_datas, indent=2))

    return munch.munchify(results_datas)

def validation(
        molecules,
        this_ff,
        output_dir,
        optimizer="openmm",
        hessian_flag=False,
        fitting_info=None,
        init_this_ff=None,
):
    init_moles = initi_mole_info_count(molecules)
    if init_this_ff is None:
        init_this_ff = this_ff
    init_para = initi_ff_info_count(init_this_ff)
    logger.info("Get validation results ......")
    return intra_fitting_validation(
        molecules,
        this_ff,
        init_moles,
        init_para,
        output_dir,
        optimizer=optimizer,
        hessian_flag=hessian_flag,
        fitting_info=fitting_info,
    )