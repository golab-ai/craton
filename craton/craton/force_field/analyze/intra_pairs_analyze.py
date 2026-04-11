#!/usr/bin/env python
"""

"""
from copy import deepcopy

from ...utils import logger
#from ...chemkit import MolConformer as MConf  

from ...utils.numerical_algorithm import rmse_calculate, linear_fitting

DEL_ENERGY = 5.0


class IntraAnalyze:
    def __init__(self, style="", save_path="./"):
        self.style = style
        self.save_path = save_path

    _intra_term ={
            "energy":"energy",
            "force":"energy_dev",
            "hessian":"energy_dev",
            "freq": "energy_dev",
            "frequency":"energy_dev",
            "pes": "pes",
            "torsion_scan":"pes",
            "rmsd":"rmsd",
            "Bonds":"bonded",
            "Angles":"bonded",
            "Dihedrals":"bonded",
            "Impropers":"bonded",
            "Pair14":"bonded",
            "Pair1n":"bonded",
            "Pair13":"bonded",
            "Pair12":"bonded",
            "esp_charge":"charge",
            "mulliken_charge":"charge",
            "point_charge":"charge",
            "Bond":"bonded",
            "Angle":"bonded",
            "Dihedral":"bonded",
            "Improper":"bonded",
            "ff_charge":"charge",
        }
    _default_ = ["energy","force","hessian","freq","rmsd","pes",
                 "Bonds","Angles","Dihedrals","Impropers","pair1n","esp_charge"]
    @staticmethod
    def get_same_mole(moles):
        same_mole = {}
        for i in range(len(moles)):
            if moles[i].mole_name not in same_mole.keys():
                same_mole[moles[i].mole_name] = []
            same_mole[moles[i].mole_name].append(i)
        return same_mole

    @staticmethod
    def _get_statics(X,Y):
        rmse, mae = rmse_calculate(X, Y)
        a, b, r2 = linear_fitting(X, Y)
        
        return [rmse,mae,a,b,r2]

    @staticmethod
    def _energy(same_molecules,term,molecules):
        values = {"total":[]}
        for name,idxs in same_molecules.items():
            tmp = [getattr(molecules[idx],term) for idx in idxs if hasattr(molecules[idx],term)]
            min_e = min(tmp)
            values[name] = [round(rr-min_e,4) for rr in tmp]
            values["total"].extend(values[name])
        return values

    @staticmethod
    def _energy_pairs(same_molcules, term, molecules, ts_molecules):
        values = {"total":[[],[]],"check_conformer":[]}
        statics = {}
        for name, idxs in same_molcules.items():
            _tmp0 = [[getattr(ts_molecules[idx],term),idx] for idx in idxs if hasattr(ts_molecules[idx],term)]
            _tmp1 = [[getattr(molecules[idx],term),idx] for idx in idxs if hasattr(molecules[idx],term)]

            tmp0 = [rr[0] for rr in _tmp0]
            idx0 = [rr[1] for rr in _tmp0]

            tmp1 = [rr[0] for rr in _tmp1]
            idx1 = [rr[1] for rr in _tmp1]

            if len(tmp0) > 0:
                min0 = min([rr for rr in tmp0])
                idx_min = tmp0.index(min0)
                min1 = tmp1[idx_min]
                values[name] = [[round(rr-min0,4) for rr in tmp0],[round(rr-min1,4) for rr in tmp1]]
                for ii,rr in enumerate(values[name][0]):
                    if abs(rr-values[name][1][ii]) >= DEL_ENERGY:
                        values["check_conformer"].append(idx0[ii])
                print("#########################################")
                print(values[name])
                print("#########################################")
                statics[name] = IntraAnalyze._get_statics(*values[name])
                values["total"][0].extend(values[name][0])
                values["total"][1].extend(values[name][1])
        statics["total"] = IntraAnalyze._get_statics(*values["total"])

        return values, statics

    @staticmethod
    def _energy_dev(same_molecule,term,molecules):
        values = {"total":[]}
        for name,idxs in same_molecule.items():
            tmp = [round(v,4) for idx in idxs if hasattr(molecules[idx],term) for v in getattr(molecules[idx],term)]
            if len(tmp) > 0:
                values[name] = tmp

                values["total"].extend(tmp)
        return values

    @staticmethod
    def _energy_dev_pairs(same_molecule, term, molecules, ts_molecules):
        values = {"total":[[],[]]}
        statics = {}
        for name,idxs in same_molecule.items():
            tmp0 = [getattr(ts_molecules[idx],term,None) for idx in idxs]
            tmp1 = [getattr(molecules[idx],term,None) for idx in idxs]
            vv0 = []
            vv1 = []
            for ii,rr in enumerate(tmp0):
                if rr is not None and tmp1[ii] is not None:
                    vv0.extend([round(v,4) for v in rr])
                    vv1.extend([round(v,4) for v in tmp1[ii]])
            if len(vv0) > 0:
                values[name] = [vv0,vv1]
                statics[name] = IntraAnalyze._get_statics(vv0,vv1)
                values["total"][0].extend(vv0)
                values["total"][1].extend(vv1)
        statics["total"] = IntraAnalyze._get_statics(*values["total"])
        return values,statics

    @staticmethod
    def _pes_pairs(same_molecule, term, molecules, ts_molecules):
        from ...chemkit import MolConformer as MConf  
        values = {}
        statics = {}
        ts_scan_curve = MConf._scan_curve(ts_molecules)
        scan_curve = MConf._scan_curve(molecules)

        ts_scan_data = MConf._scan_curve_data(ts_scan_curve)
        scan_data = MConf._scan_curve_data(scan_curve)
        
        for name,vv in ts_scan_curve.items():
            values[name] = {}
            statics[name] = {}
            for torsion, vvv in vv.items():
                torsion_name = f"{name}_{torsion}"
                tmp0 = ts_scan_data[torsion_name]
                tmp1 = scan_data[torsion_name]

                min0 = min(tmp0[1])
                min0_index = tmp0[1].index(min0)
                min1 = tmp1[1][min0_index]
                tmp1[1] = [rr-min1 for rr in tmp1[1]]
                
                
                rmsd = [tmp0[0],[round(MConf._conformer_RMSD(scan_curve[name][torsion][ii],mm),4) for ii,mm in enumerate(vvv)]]
                values[name][torsion] = [tmp0,tmp1,rmsd]
                statics[name][torsion] = IntraAnalyze._get_statics(tmp0[1],tmp1[1])
        return values, statics

    @staticmethod
    def _rmsd_pairs(same_molecule, term, molecules, ts_molecules):
        from ...chemkit import MolConformer as MConf  
        values = {"total":[]}
        statics = {}
        for name,idxs in same_molecule.items():
            values[name] = [round(MConf._conformer_RMSD(molecules[idx],ts_molecules[idx]),4) for idx in idxs]
            values["total"].extend(values[name])
        return values,statics

    @staticmethod
    def _bonded_pairs(same_molecule, term, molecules, ts_molecules):
        values = {"total":[[],[]]}
        statics = {}
        _convert = {"angles":"","dihedrals":"",}
        _label = {
                    "Angles":"value",
                    "Dihedrals":"value",
                    "Impropers":"value",
                    "Bonds":"value",
                    "Pair12":"value",
                    "Pair13":"value",
                    "Pair14":"value",
                    "Pair1n":"value",
                  }
        attr = _label[term]
        for name,idxs in same_molecule.items():
            this_idxs = [idx for idx in idxs if hasattr(molecules[idx],term)]
            if len(this_idxs) > 0:
                if term in ["dihedrals","dihedral"]:
                    tmp0 = [round(getattr(topol,attr),4) for idx in this_idxs for topol in getattr(ts_molecules[idx],term) if not topol.is_linear]
                    tmp1 = [round(getattr(topol,attr),4) for idx in this_idxs for topol in getattr(molecules[idx],term) if not topol.is_linear]
                    for ii,v in enumerate(tmp0):
                        if tmp1[ii] - v > 180:
                            tmp1[ii] -= 360
                        elif tmp1[ii] - v < -180:
                            tmp1[ii] += 360
                else:
                    tmp0 = [round(getattr(topol,attr),4) for idx in this_idxs for topol in getattr(ts_molecules[idx],term)]
                    tmp1 = [round(getattr(topol,attr),4) for idx in this_idxs for topol in getattr(molecules[idx],term)]
                values[name] = [tmp0,tmp1]
                statics[name] = IntraAnalyze._get_statics(tmp0,tmp1)
                values["total"][0].extend(tmp0)
                values["total"][1].extend(tmp1)
        statics["total"] = IntraAnalyze._get_statics(*values["total"])
        return values,statics

    @staticmethod
    def _charge_pairs(same_molecule, term, molecules, ts_molecules):
        values = {"total":[[],[]]}
        statics = {}
        for name,idxs in same_molecule.items():
            this_idxs = [idx for idx in idxs if hasattr(ts_molecules[idx].Atoms[0],term)]
            if len(this_idxs) > 0:
                tmp0 = [round(v,4) for idx in this_idxs for v in getattr(ts_molecules[idx],term)]
                tmp1 = [round(v,4) for idx in this_idxs for v in getattr(molecules[idx],"ff_charge")]
                values[name] = [tmp0,tmp1]
                statics[name] = IntraAnalyze._get_statics(tmp0,tmp1)
                values["total"][0].extend(tmp0)
                values["total"][1].extend(tmp1)
        statics["total"] = IntraAnalyze._get_statics(*values["total"])
        return values, statics

    @staticmethod
    def intra_molecule_pairs(molecules, ts_molecules, terms=None):
        _func = {
            "energy": IntraAnalyze._energy_pairs,
            "energy_dev": IntraAnalyze._energy_dev_pairs,
            "charge": IntraAnalyze._charge_pairs,
            "bonded": IntraAnalyze._bonded_pairs,
            "pes": IntraAnalyze._pes_pairs,
            "rmsd": IntraAnalyze._rmsd_pairs,
        }
        
        if terms is None:
            terms = list(IntraAnalyze._default_) 
        for ii,term in enumerate(terms):
            if term in ["bond","bonds","angle","angles","dihedral","dihedrals","improper","impropers",
                        "pair12","pair12s","pair13","pair13s","pair14","pair14s","pair1n","pair1ns"]:
                term = term.capitalize()
            if term in ["Bond","Angle","Dihedral","Improper"]:
                term += "s"
            if term in ["Pair12s","Pair13s","Pair14s","Pair1ns"]:
                term = term[:-1]
            terms[ii] = term


        same_mole = IntraAnalyze.get_same_mole(molecules)  # {inchi_key: [id in moles, ...]}

        values = {}
        statics = {}
        for term in terms:
            term_type = IntraAnalyze._intra_term[term]
            tmp_values,tmp_statics = _func[term_type](same_mole,term,molecules,ts_molecules)

            values[term],statics[term] = tmp_values,tmp_statics
            
        return values, statics


