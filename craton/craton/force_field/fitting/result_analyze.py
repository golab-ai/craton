#!/usr/bin/env python
"""

"""
from copy import deepcopy

import numpy as np
from rdkit.Chem.AllChem import AlignMol
from scipy.optimize import leastsq

from ...utils import logger
from ...chemkit.conformation import get_conformation_RMSD

from ...utils.numerical_algorithm import rmse_calculate, linear_fitting

R_N = 4
R_N_ANGLE = 1


class ResultAnaly:
    def __init__(self, style, save_path="./"):
        self.style = style
        self.save_path = save_path

    @staticmethod
    def get_same_mole(moles):
        same_mole = {}
        for i in range(len(moles)):
            if moles[i].mole_name not in same_mole.keys():
                same_mole[moles[i].mole_name] = []
            same_mole[moles[i].mole_name].append(i)
        return same_mole

    @staticmethod
    def get_scan_curve(moles):
        scan_curve = {}
        for i in range(len(moles)):
            if hasattr(moles[i], "constrain"):
                if moles[i].mole_name not in scan_curve.keys():
                    scan_curve[moles[i].mole_name] = {}
                scan_term_name = moles[i].constrain[0].name
                if scan_term_name not in scan_curve[moles[i].mole_name].keys():
                    scan_curve[moles[i].mole_name][scan_term_name] = []
                scan_curve[moles[i].mole_name][scan_term_name].append(i)
        return scan_curve

    def get_ff_fitting_t1(self, mms, term, moles, moles_ts):
        datas = {}
        tmp = [[], []]
        for inchi_key, indexes in mms.items():
            datas[inchi_key] = [[], []]
            for idx in indexes:
                a = getattr(moles[idx], term)
                b = getattr(moles_ts[idx], term)
                datas[inchi_key][0].append(b)
                datas[inchi_key][1].append(a)
            min_v0 = min(datas[inchi_key][0])
            idx_min = datas[inchi_key][0].index(min_v0)
            min_v1 = datas[inchi_key][1][idx_min]
            for i in range(len(datas[inchi_key][0])):
                datas[inchi_key][0][i] = round(datas[inchi_key][0][i] - min_v0, R_N)
                datas[inchi_key][1][i] = round(datas[inchi_key][1][i] - min_v1, R_N)
                tmp[0].append(datas[inchi_key][0][i])
                tmp[1].append(datas[inchi_key][1][i])
        datas["total"] = tmp
        return datas

    def get_ff_fitting_t2(self, mms, term, moles, moles_ts, chargetype="esp"):
        datas = {}
        tmp = [[], []]
        _term = f"{chargetype}_{term}"
        for aa, bb in mms.items():
            datas[aa] = [[], []]
            for nn in bb:
                if hasattr(moles_ts[nn].Atoms[0],_term):
                    for ii,atom in enumerate(moles[nn].Atoms):
                        a = atom.ff_charge
                        b = getattr(moles_ts[nn].Atoms[ii],_term)
                        datas[aa][0].append(round(b, R_N))
                        datas[aa][1].append(round(a, R_N))
                        tmp[0].append(round(b, R_N))
                        tmp[1].append(round(a, R_N))

                #if term in moles[nn].Atoms[0].__dict__.keys():
                #    for i in range(len(moles[nn].Atoms)):
                #        a = moles[nn].Atoms[i].ff_charge
                #        b = moles_ts[nn].Atoms[i].charge[chargetype]
                #        datas[aa][0].append(round(b, R_N))
                #        datas[aa][1].append(round(a, R_N))
                #        tmp[0].append(round(b, R_N))
                #        tmp[1].append(round(a, R_N))
        datas["total"] = tmp
        return datas

    def get_ff_fitting_t3(self, mms, term, moles, moles_ts):
        datas = {}
        tmp = [[], []]
        for aa, bb in mms.items():
            datas[aa] = [[], []]
            for nn in bb:
                if hasattr(moles[nn], term):
                    for i in range(len(getattr(moles[nn], term))):
                        if term not in ["Angles", "Dihedrals", "Impropers"]:
                            a = getattr(moles[nn], term)[i].value
                            b = getattr(moles_ts[nn], term)[i].value
                            a = round(a, R_N)
                            b = round(b, R_N)
                        else:
                            a = getattr(moles[nn], term)[i].value_a
                            b = getattr(moles_ts[nn], term)[i].value_a
                            if term == "Dihedrals":
                                # ignore linear dihedral
                                if getattr(moles[nn], term)[i].is_linear:
                                    continue
                                # adjust dihedral for 180 and -180
                                if a - b > 180:
                                    a -= 360
                                elif a - b < -180:
                                    a += 360
                            a = round(a, R_N_ANGLE)
                            b = round(b, R_N_ANGLE)
                        datas[aa][0].append(b)
                        datas[aa][1].append(a)
                        tmp[0].append(b)
                        tmp[1].append(a)
        datas["total"] = tmp
        dd = []
        for aa, bb in datas.items():
            if len(bb[0]) == 0:
                dd.append(aa)
        for aa in dd:
            del datas[aa]
        return datas

    def get_ff_fitting_t4(self, mms, term, moles, moles_ts):
        datas = {}
        tmp = [[], []]
        for aa, bb in mms.items():
            datas[aa] = [[], []]
            for nn in bb:
                if hasattr(moles[nn], term):
                    for i in range(len(getattr(moles[nn], term))):
                        a = getattr(moles[nn], term)[i]
                        b = getattr(moles_ts[nn], term)[i]
                        datas[aa][0].append(round(b, R_N))
                        datas[aa][1].append(round(a, R_N))
                        tmp[0].append(round(b, R_N))
                        tmp[1].append(round(a, R_N))
        datas["total"] = tmp
        dd = []
        for aa, bb in datas.items():
            if len(bb[0]) == 0:
                dd.append(aa)
        for aa in dd:
            del datas[aa]
        return datas

    def get_ff_fitting_t5(self, mms, term, moles, moles_ts, isfittingparas, rmse_datas=None):
        datas = {}
        for inchi_key, bb in mms.items():
            datas[inchi_key] = {}
            for name_scan, indexes in bb.items():
                datas[inchi_key][name_scan] = [[[], []], [[], []], [[], [], [], []], [[], []]]
                scan = datas[inchi_key][name_scan]
                """
                [[[value],
                  [qm_energy]
                 ],
                 [[value],
                  [mm_energy]
                 ],
                 [[dihedral_list, e.g. [0,2,3,4], [1,2,3,4]],
                  [unique_dihedral_para, e.g. "c_3$c_3$o_2$h_1", "h_1$c_3$o_2$h_1"],
                  [times_unique_dihedral_para],
                  [dihedral_para_is_fitting?]
                 ],
                 [[value],
                  [rmsd]
                 ]
                ]
                """

                m0 = moles[indexes[0]]
                dihedral_para_dict = {}
                for dd in m0.Dihedrals:
                    dname1 = f"{str(dd.a1)}-{str(dd.a2)}-{str(dd.a3)}-{str(dd.a4)}"
                    dname2 = f"{str(dd.a4)}-{str(dd.a3)}-{str(dd.a2)}-{str(dd.a1)}"
                    pname1 = (
                        f"{dd.a1_atom_type_used}${dd.a2_atom_type_used}${dd.a3_atom_type_used}${dd.a4_atom_type_used}"
                    )
                    pname2 = (
                        f"{dd.a4_atom_type_used}${dd.a3_atom_type_used}${dd.a2_atom_type_used}${dd.a1_atom_type_used}"
                    )
                    dihedral_para_dict[dname1] = [pname1, pname2]
                    dihedral_para_dict[dname2] = [pname1, pname2]

                arrs = []
                this_scan_term = m0.constrain[0].atoms
                a0_connect = [this_scan_term[0]]
                a1_connect = [this_scan_term[3]]
                for ii in m0.Atoms[this_scan_term[1]].connect:
                    if ii != this_scan_term[0] and ii != this_scan_term[3] and ii != this_scan_term[2]:
                        a0_connect.append(ii)
                for ii in m0.Atoms[this_scan_term[2]].connect:
                    if ii != this_scan_term[0] and ii != this_scan_term[3] and ii != this_scan_term[1]:
                        a1_connect.append(ii)
                for ii in a0_connect:
                    for jj in a1_connect:
                        arrs.append([ii, this_scan_term[1], this_scan_term[2], jj])
                for rr in arrs:
                    scan[2][0].append(rr)
                    dihedral_name = f"{str(rr[0])}-{str(rr[1])}-{str(rr[2])}-{str(rr[3])}"
                    names = dihedral_para_dict[dihedral_name]
                    dihe_name = list(set(names).intersection(set(scan[2][1])))
                    if len(dihe_name) == 1:
                        scan[2][2][scan[2][1].index(dihe_name[0])] += 1
                    elif len(dihe_name) == 0:
                        scan[2][1].append(names[0])
                        scan[2][2].append(1)
                        if len(set(names).intersection(set(isfittingparas))) == 1:
                            scan[2][3].append(1)
                        else:
                            scan[2][3].append(0)

                for idx in indexes:
                    a0 = getattr(moles[idx], "energy")
                    a1 = moles[idx].constrain[0].fix_value
                    b0 = getattr(moles_ts[idx], "energy")
                    b1 = moles_ts[idx].constrain[0].fix_value
                    scan[0][0].append(round(b1, R_N))
                    scan[0][1].append(b0)
                    scan[1][0].append(round(a1, R_N))
                    scan[1][1].append(a0)

                    scan[3][0].append(round(a1, R_N))
                    if rmse_datas is None:
                        rmse = get_conformation_RMSD(deepcopy(moles[idx]), deepcopy(moles_ts[idx]))
                    else:
                        rmse = rmse_datas[idx]
                    scan[3][1].append(round(rmse, R_N))
                min_v0 = min(scan[0][1])
                idx_min = scan[0][1].index(min_v0)
                min_v1 = scan[1][1][idx_min]
                for i in range(len(scan[0][1])):
                    scan[0][1][i] = round(scan[0][1][i] - min_v0, R_N)
                    scan[1][1][i] = round(scan[1][1][i] - min_v1, R_N)
        return datas

    def get_ff_fitting_t6(self, mms, term, moles, moles_ts):
        datas = {}
        tmp = [[], []]
        for aa, bb in mms.items():
            datas[aa] = [[], []]
            for nn in bb:
                rmse = get_conformation_RMSD(moles[nn], moles_ts[nn])
                datas[aa][0].append(0)
                datas[aa][1].append(round(rmse, R_N))
                tmp[0].append(0)
                tmp[1].append(round(rmse, R_N))
        datas["total"] = tmp
        return datas

    def get_ff_fitting_result(
        self, moles, moles_ts, terms="all", isfittingparas=[], chargetype="esp", energy_diff_threshold=0
    ):
        term1 = ["energy"]
        term2 = ["charge"]
        term3 = ["Bonds", "Angles", "Dihedrals", "Impropers", "Pair14", "Pair1n"]
        term4 = ["force","hessian","freq",]
        term5 = ["pes"]
        term6 = ["rmse"]
        if terms == "all":
            terms = term1 + term2 + term3 + term4 + term6 + term5
        same_mole = self.get_same_mole(moles)  # {inchi_key: [id in moles, ...]}
        scan_curve = self.get_scan_curve(moles)  # {inchi_key: {'0-1-2-6': [id in moles, ...]}}
        datas = {}
        result = {}
        for tt in terms:
            if tt in term1:
                datas[tt] = self.get_ff_fitting_t1(same_mole, tt, moles, moles_ts)
            if tt in term2:
                datas[tt] = self.get_ff_fitting_t2(same_mole, tt, moles, moles_ts, chargetype=chargetype)
            if tt in term3:
                datas[tt] = self.get_ff_fitting_t3(same_mole, tt, moles, moles_ts)
            if tt in term4:
                datas[tt] = self.get_ff_fitting_t4(same_mole, tt, moles, moles_ts)
            if tt in term5:
                if "rmse" in datas.keys():
                    datas[tt] = self.get_ff_fitting_t5(
                        scan_curve, tt, moles, moles_ts, isfittingparas, rmse_datas=datas["rmse"]["total"][1]
                    )
                else:
                    datas[tt] = self.get_ff_fitting_t5(scan_curve, tt, moles, moles_ts, isfittingparas)
            if tt in term6:
                datas[tt] = self.get_ff_fitting_t6(same_mole, tt, moles, moles_ts)
        for tt in terms:
            if tt not in term5:
                result[tt] = {}
                if len(datas[tt]) > 0:
                    for aa, bb in datas[tt].items():
                        a, b = rmse_calculate(bb[0], bb[1])
                        c, d, e = linear_fitting(bb[0], bb[1])

                        result[tt][aa] = [a, b, c, d, e]
                else:
                    del datas[tt]
                    del result[tt]
            else:
                result[tt] = {}
                for aa, bb in datas[tt].items():
                    if len(bb) > 0:
                        result[tt][aa] = {}
                        for aaa, bbb in bb.items():
                            a, b = rmse_calculate(bbb[0][1], bbb[1][1])
                            c, d, e = linear_fitting(bbb[0][1], bbb[1][1])
                            result[tt][aa][aaa] = [a, b, c, d, e]
                    else:
                        del datas[tt][aa]
                if len(datas[tt]) == 0:
                    del datas[tt]
                    del result[tt]

        bad_curves = {}
        if energy_diff_threshold >= 0:
            for inchi_key, d in scan_curve.items():
                for curve, indexes in d.items():
                    min_energy = min([moles_ts[i].energy for i in indexes])
                    idx_min = next(i for i in indexes if moles_ts[i].energy == min_energy)
                    for idx in indexes:
                        energy = moles[idx].energy - moles[idx_min].energy
                        energy_ts = moles_ts[idx].energy - moles_ts[idx_min].energy
                        if abs(energy - energy_ts) > energy_diff_threshold:
                            if inchi_key not in bad_curves:
                                bad_curves[inchi_key] = []
                            bad_curves[inchi_key].append(curve)
                            break

        bad_moles = []
        for inchi_key, curve_list in bad_curves.items():
            for curve in curve_list:
                for idx in scan_curve[inchi_key][curve]:
                    bad_moles.append((moles_ts[idx], moles[idx]))

        return (datas, result, bad_moles)

    def get_para_fitting_result(self, moles, moles_ts, is_fitting_params=[]):
        datas = {
            "Bonds": {},
            "Angles": {},
            "Dihedrals": {},
            "Impropers": {},
        }
        for i in range(len(moles)):
            for tt in datas.keys():
                if hasattr(moles[i], tt):
                    for j in range(len(getattr(moles[i], tt))):
                        if tt == "Bonds":
                            a = getattr(moles[i], tt)[j].value
                            b = getattr(moles_ts[i], tt)[j].value
                            a = round(a, R_N)
                            b = round(b, R_N)
                        else:
                            a = getattr(moles[i], tt)[j].value_a
                            b = getattr(moles_ts[i], tt)[j].value_a
                            a = round(a, R_N_ANGLE)
                            b = round(b, R_N_ANGLE)
                        if getattr(getattr(moles[i], tt)[j], "atom_type_used_name") not in datas[tt].keys():
                            if getattr(getattr(moles[i], tt)[j], "name") in is_fitting_params:
                                datas[tt][getattr(getattr(moles[i], tt)[j], "atom_type_used_name")] = [[b], [a], 1]
                            else:
                                datas[tt][getattr(getattr(moles[i], tt)[j], "atom_type_used_name")] = [[b], [a], 0]
                        else:
                            datas[tt][getattr(getattr(moles[i], tt)[j], "atom_type_used_name")][0].append(b)
                            datas[tt][getattr(getattr(moles[i], tt)[j], "atom_type_used_name")][1].append(a)
        return datas


