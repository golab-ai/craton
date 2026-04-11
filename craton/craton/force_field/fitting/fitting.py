from copy import deepcopy
from typing import *

import munch
import numpy as np

from ...utils import logger
#from ..calculator.optimizer import optimize
#from ...mm_calculator import Calculator as Calc


from ...chemkit.conformation.conformation import ConformType
from ...chemkit import Structure as Stru
from ..force_field import ForceField as FF
from ..forcefield_manager import ForceFieldManager as FFM
#from ..force_field import _special_dihedrals_infos, assign_force_field_parameter

from ._fitting import st, FFFConstructor

def _pseudo_fitting_dihedralterm(this_ff, molecules):
    """
    检查没有scan数据，而被拟合的dihedral参数，将其标记为pseudo_fitting。
    `this_ff` will be modified

    输入：
        this_ff: Dict, 力场参数
        pes_data: Dict, 势能面的数据。由validation方法生成
    """
    FFM.get_special_infos(molecules,this_ff)
    _tmp = [para for para,vv in this_ff["dihedralterm"].items() if vv["ptag"] == "isfitting" and para not in this_ff["torsion_para"]]
    _tmp = [para for para in _tmp if "@" not in para.split("$")[1] and "@" not in para.split("$")[2]]
    for para in _tmp:
        this_ff["dihedralterm"][para]["ptag"] = "pseudofitting"

def bonded_fitting_func(
    moles,
    target_moles,
    this_ff: Dict,
    target_properties: List[str] = ["energy", "force", "hessian"],
    ):
    for i in range(len(moles)):
        moles[i].energy = target_moles[i].energy
        if hasattr(target_moles[i], "force"):
            moles[i].force = target_moles[i].force
        if hasattr(target_moles[i], "hessian"):
            moles[i].hessian = target_moles[i].hessian
    weight_factor = [WEIGHT[pp] for pp in target_properties]
    fitt1 = para_fitting_func(
        "intra",
        targets=target_properties,
        weight_factor=weight_factor,
    )
    fitt1.import_claculator()

    w0, w1, w2, w3 = fitt1.fitting_qmdata(moles, this_ff)
    return w0, w1, w2, w3

def bonded_do_fitting(
    molecules,
    this_ff: Dict,
    w0: List[Callable],
    w1: Dict,
    w2: List[float],
    w3: List[str],
    bup: np.ndarray,
    bdown: np.ndarray,
    ):
    sst = st("intra", w0, w1, w2, w3, bup, bdown)
    fit_ff = sst.calc()
    label = {
        "Angles": "angleterm",
        "Bonds": "bondterm",
        "Dihedrals": "dihedralterm",
        "Impropers": "improperterm",
    }
    for aa, bb in w1.items():
        for aaa, bbb in bb.items():
            for aaaa, bbbb in bbb.items():
                this_ff[label[aa]][aaa]["parameter"][aaaa] = fit_ff.x[bbbb]
    for m in molecules:
        m.assign_ff_para(this_ff)

def update_parameter(molecules,this_ff,fit_ff,fitting_parameter):
    pass
    #return this_ff

def bonded_fitting(
    target_molecules,
    this_ff,
    fitting_terms=["bondterm", "angleterm", "dihedralterm", "improperterm"],
    target_prop = ["energy", "force", "hessian"],
    weight_factor = None,
    optimizer = "numpy",
    torsion_constraint_step = None,
    ):
    from ...mm_calculator import Calculator as Calc

    WEIGHT: Dict[str, float] = {
    "energy": 1.0,
    "force": 0.0001,
    "hessian": 0.0001,
    "penalty_torsion": 0.1,
    }

    if weight_factor is None:
        weight_factor = [WEIGHT[term] for term in target_prop]
    if not torsion_constraint_step:
        torsion_constraint_step = [1500.0]

    molecules = [deepcopy(molecule) for molecule in target_molecules]
    opt_molecules = [molecule for molecule in molecules if molecule.conform_type in ConformType.TORSION_SCAN_TYPES + [ConformType.LOCAL_MINIMUM]]
    relax_molecules = [deepcopy(molecule) for molecule in opt_molecules]

    results = {"target_function": [-1]}
    for i, k in enumerate(torsion_constraint_step):
        logger.info(f"Bonded fitting: Step {i + 1}/{len(torsion_constraint_step)} ...")
        #for molecule in relax_molecules:
        #    molecule.assign_ff_para(this_ff)
        relax_molecules = FFM.assign_force_field_parameter(relax_molecules,this_ff)
        relax_molecules = Calc._optimize(relax_molecules, all_torsion_constraint=k,optimizer=optimizer)

        for molecule, molecule_relaxed in zip(opt_molecules, relax_molecules):
            for atom, atom_relaxed in zip(molecule.Atoms, molecule_relaxed.Atoms):
                atom.coor = atom_relaxed.coor
        opt_molecules = Stru._update_topol_values(opt_molecules)
        #for molecule in opt_molecules:
        #    molecule.update_topol_value()

        logger.info("Bonded fitting: fitting function generation ......")
        
        FFFC = FFFConstructor(molecules,this_ff,fitting_terms=fitting_terms,targets=target_prop,weight_factor=weight_factor)
        funcs,parameter_init,fitting_parameter,boundary = FFFC.run()
        logger.info("Bonded fitting: fitting funciton generation Done")
        if len(parameter_init) > 0:
            logger.info("Bonded fitting run ......")

            ##### fitting parameter
            sst = st(funcs,parameter_init,fitting_parameter,boundary)
            fit_ff = sst.run()
            
            #####update parameter
            for item_name, data in fitting_parameter.items():
                ss = item_name.split("-")
                for kk,vv in data.items():
                    this_ff[ss[0]][ss[1]]["parameter"][kk] = fit_ff.x[vv]
            molecules = FFM.assign_force_field_parameter(molecules,this_ff)
            #for molecule in molecules:
            #    molecule.assign_ff_para(this_ff)
            
            results["target_function"].append(-1)
            logger.info("Bonded fitting run Done")

    return results

def binc_fitting(moles_ts, this_ff: Dict, target: str = "esp"):
    """
    Fit binc parameters to ESP/am1-bcc charge
    `this_ff` will be modified
    Set ff_charge of atoms in `moles_ts` based on formal charge and binc

    Parameters
    ----------
    moles_ts : list of Molecule
    this_ff : dict
    """
    flexible_binc_terms = [
        k
        for k, v in this_ff["binc"].items()
        if k.split("$")[0] != k.split("$")[1] and v["ptag"] == "isfitting" and 0 in v["isfitting"]
    ]
    n_binc = len(flexible_binc_terms)

    mol_list = []
    charge_list = []
    n_row = 0
    charge_name = f"{target}_charge"
    for mol in moles_ts:
        if charge_name not in mol.charges:
            continue

        mol_list.append(mol)
        #n_row += mol.mole_n + mol.ring_number
        n_row += mol.atom_count + mol.ring_number
        charge_list += getattr(mol,charge_name) + [0] * mol.ring_number

    if not mol_list:
        logger.error(f"Target charge data not found for any molecule: {target}")
        raise Exception("BINC fitting failed")

    def get_binc_direction(binc, this_ff: Dict) -> Tuple[Any, int]:
        name = "$".join(binc)
        if name in this_ff["binc"]:
            return name, 1

        name = "$".join(reversed(binc))
        if name in this_ff["binc"]:
            return name, -1

        raise Exception(f"BINC term not found in this_ff: {name}")

    A = np.zeros([n_row, n_binc])
    b = np.array(charge_list)[:, np.newaxis]
    res = np.zeros([n_row, 1])
    i_row = 0
    for mol in mol_list:
        for atom in mol.Atoms:
            res[i_row][0] += atom.ff_charge_base
            for neigh in atom.connect:
                binc_atoms = [atom.binc_atom_type, mol.Atoms[neigh].binc_atom_type]
                name, direction = get_binc_direction(binc_atoms, this_ff)
                try:
                    idx = flexible_binc_terms.index(name)
                except ValueError:
                    res[i_row][0] += this_ff["binc"][name]["parameter"][0] * direction
                else:
                    A[i_row][idx] += direction
            i_row += 1
        for _, ring in mol.ring_dict.items():
            for i in range(len(ring) - 1):
                a1 = mol.Atoms[ring[i]]
                if i != len(ring) - 2:
                    a2 = mol.Atoms[ring[i + 1]]
                else:
                    a2 = mol.Atoms[ring[0]]
                binc_atoms = [a1.binc_atom_type, a2.binc_atom_type]
                name, direction = get_binc_direction(binc_atoms, this_ff)
                try:
                    idx = flexible_binc_terms.index(name)
                except ValueError:
                    res[i_row][0] += this_ff["binc"][name]["parameter"][0] * direction
                else:
                    A[i_row][idx] += direction
            i_row += 1

    result, rsq, _, _ = np.linalg.lstsq(A, b - res, rcond=None)

    X = np.hstack([A, res, b, np.matmul(A, result) + res])

    for term, val in zip(flexible_binc_terms, result):
        this_ff["binc"][term]["parameter"][0] = float(val[0])
        this_ff["binc"][term]["ptag"] = "Fit"

    for mol in moles_ts:
        FF.assign_charge_para(mol,this_ff)
        #mol.assign_charge_para(this_ff)

def fitting(
        this_ff,
        molecules,
        fitting_terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
        target_prop = ["energy", "force", "hessian", "penalty_torsion"],
        torsion_constraint_step=None,
        optimizer="openmm",
    ):
    molecules = Stru._update_topol_values(molecules)
    molecules = FFM.assign_force_field_parameter(molecules,this_ff)
    #for molecule in molecules:
    #    molecule.update_topol_value()
    #    molecule.assign_para(this_ff)
    
    results=[]
    try:
        for key in set(fitting_terms) - {"binc"}:
            for term in this_ff[key].values():
                if term["ptag"] == "isfitting" and term["isfitting"]:
                    raise StopIteration
    except StopIteration:
        results = bonded_fitting(
            molecules,
            this_ff,
            fitting_terms=fitting_terms,
            target_prop=target_prop,
            optimizer=optimizer,
            torsion_constraint_step=torsion_constraint_step,
            )
    else:
        logger.info("No intra parameters to be fit ...")

    _pseudo_fitting_dihedralterm(this_ff,molecules)    

    return results

