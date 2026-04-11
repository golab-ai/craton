import copy
import sys
from pathlib import Path
from ..utils import logger
from .mm.mm_calculator import MMCalc, NumpyMMCalc, OpenMMCalc
#from ..chem import FormatMolecule as FM
from ..utils.commons import parallel_run
from ..chemkit import Structure as Stru



def save_gjf_file(molecules,dir):
    from ..chem import FormatMolecule as FM
    name_filename = dict()
    for molecule in molecules:
        finename = molecule.name
        for cons in getattr(molecule, "constrain", []):
            for idx in cons.atoms:
                filename += "_%i" % (idx + 1)
            filename += "_%.1f" % cons.fix_value
        if filename in name_filename:
            name_filename[filename] += 1
            filename += "%s" %name_filename[filename]
        else:
            name_filename[filename] = 0
        
        Path(dir).mkdir(exist_ok=True)
        FM._convert(molecules,otype="gjf",ofilename=filename,opath=dir)


def optimize(molecules, optimizer="openmm", all_torsion_constraint: float = 0.0, write_mol=None, **kwargs):
    """
    对分子进行力场的结构优化，可以按分子进行并行计算。调用do_optimize完成
    输入：
        moles: List[Molecule], 优化的分子集
        optimizer: str, can be openmm, numpy, multip or None
        all_torsion_constraint: float, force constant for cosine constrain on all torsion, in kcal/mol/rad^2
    输出：
        moles: List[Molecule], 优化后的分子集
    """
    if all_torsion_constraint and optimizer not in ["openmm", "numpy"]:
        logger.error("Non-zero constrain_all_torsion requires openmm or numpy optimizer")
        sys.exit(1)

    if write_mol is not None:
        save_gjf_file(molecules,write_mol + "/gjf_before")

    if optimizer == "numpy":  # use numpy to vectorize the energy and force calculation
        mm = NumpyMMCalc(molecules)
        mm.optimize(all_torsion_constraint=all_torsion_constraint, **kwargs)

    elif optimizer == "openmm":  # use openmm md_simulation engine
        mm = OpenMMCalc(molecules)
        mm.optimize(all_torsion_constraint=all_torsion_constraint, **kwargs)

    elif optimizer == "multip": 
         # multi-process: get better performance
        mm = MMCalc()
        molecules = parallel_run("optimize",molecules,objs=mm,single_args_flag=False,keep_order=True,return_result=True)
        

        #molecules = [None] * len(new_molecules)
        #for ii,tt in enumerate(new_molecules):
        #    molecules[ii] = tt
        #    molecules.
        #parallel_run("update_topol_value",None,objs=molecules,single_args_flag=True,keep_order=False,return_result=False)
        molecules = Stru._update_topol_values(molecules)
        #for m in molecules:
        #    m.update_topol_value()
    else:
        mm = MMCalc()  # single process - for debug
        for i in range(len(molecules)):
            logger.debug(f"Optimizing molecule {i + 1} / {len(molecules)} ...")
            molecules[i] = mm.optimize(molecules[i], i)[0]
            molecules[i].update_topol_value()

    if write_mol is not None:
        save_gjf_file(molecules,write_mol + "/gjf_after")

    return molecules


def empirical_optimize(moles, optimizer="openmm", two_step=False):
    """
    TODO Slow code. Avoid calling it
    """

    mols_opt = [copy.deepcopy(m) for m in moles]
    MM.molecule_structure(mols_opt, create_intra_nonbond=True, update_mole_info=False)
    #MoleculeManager.prepare_molecule(mols_opt, create_intra_nonbond=True, update_mole_info=False)
    ff = ForceFieldManager.prepare_ff(mols_opt)
    for m in mols_opt:
        m.assign_ff(ff)
        m.update_topol_value()

    if two_step:
        optimize(mols_opt, optimizer=optimizer, all_torsion_constraint=1500.0, tol_force=1.0, max_iter=100)
    optimize(mols_opt, optimizer=optimizer, tol_force=1.0, max_iter=100)

    for mole, mol_opt in zip(moles, mols_opt):
        for i, atom in enumerate(mole.Atoms):
            atom.coor = mol_opt.Atoms[i].coor[:]
