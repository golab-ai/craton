from .conformation import assign_scan_conf_type
from .conformation import get_conformation_RMSD
from .conformation import get_scan_curve
from .conformation import create_lm_by_combine_scan_rlm
from .conformation import local_minimum_pes
from .conformation import conformation_id_hash
from .conformation import ignore_alkane_torsion
from .conformation import find_stablest_molecule
from .conformation import remove_similar_conformer
from .conformation import bond_angle_extend_conformer
from .conformation import get_bond_angle_scan_term
from .conformation import ConformType

from ...utils.commons import parallel_run

class MolConformer:
    def __init__(self) -> None:
        pass

    @staticmethod
    def _scan_curve(molecules):
        return get_scan_curve(molecules)

    @staticmethod
    def _scan_curve_data(scan_curve):
        datas = {f"{inchikey}_{term}":[[m.constrain[0].fix_value for m in dd],[m.energy for m in dd]] 
                    for inchikey,pes in scan_curve.items() for term,dd in pes.items()}
        for name,data in datas.items():
            e = min(data[1])
            data[1] = [round(rr-e,4) for rr in data[1]]
        return datas

    @staticmethod
    def _scan_conf_type(scan_curve):
        assign_scan_conf_type(scan_curve)

    @staticmethod
    def _pes_local_minimum(scan_curve):
        return local_minimum_pes(scan_curve)

    @staticmethod
    def _lm_by_combine_scan_curve(molecules,rlm_dicts,n=64, create_constrain=False,parallel=True):
        if not isinstance(molecules,list):
            molecules = [molecules]
        if not isinstance(rlm_dicts,list):
            rlm_dicts=[rlm_dicts]
        if parallel:
            args = [[molecules[ii],rlm_dicts[ii]] for ii in range(len(molecules))]
            kwds = {"n":n,"create_constrain":create_constrain}
            return parallel_run(create_lm_by_combine_scan_rlm,args,kwds=kwds)
        else:
            return [create_lm_by_combine_scan_rlm(molecules[ii],rlm_dicts[ii],n=n,create_constrain=create_constrain) for ii in range(len(molecules)) ]

    @staticmethod
    def _conformer_RMSD(molecule1,molecule2):
        return get_conformation_RMSD(molecule1,molecule2)

    @staticmethod
    def _conformation_id_hash_(molecule):
        return conformation_id_hash(molecule)
    
    @staticmethod
    def _ignore_alkane_torsion_(molecules):
        if not isinstance(molecules,list):
            ignore_alkane_torsion(molecules)
        else:
            for molecule in molecules:
                ignore_alkane_torsion(molecule)
        return molecules

    @staticmethod
    def _find_stablest_molecule(molecules):
        return find_stablest_molecule(molecules)

    @staticmethod
    def _remove_similar_conformer(molecules,target_molecule=None):
        return remove_similar_conformer(molecules,target_molecule=target_molecule)
    
    @staticmethod
    def _bond_angle_extend_conformer(molecule,ignore_alkane=True):
        return bond_angle_extend_conformer(molecule,ignore_alkane=ignore_alkane)
    
    @staticmethod
    def _get_bond_angle_scan_term(molecule,inter_val=[0.1,0.2,5.0,5.0],ignore_ring=True,exists_type=None):
        return get_bond_angle_scan_term(molecule,inter_val=inter_val,ignore_ring=ignore_ring,exists_type=exists_type)