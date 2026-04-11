#!/usr/bin/env python
"""

"""

import numpy as np
from copy import deepcopy

from .ff_calculator.nonbond import charge_calculator, vdw_calculator, nonbond_calculator
from .ff_calculator.bonded import bond_calculator, angle_calculator, dihedral_calculator, improper_calculator, constrain_calculator


from ...utils import logger
from ...utils.commons import parallel_run
from ...force_field.force_field import _topol_to_ff_term, _Pair_terms

Calc = {
    "Bonds": bond_calculator,
    "Angles": angle_calculator,
    "Dihedrals": dihedral_calculator,
    "Impropers": improper_calculator,
    "Pair1n": nonbond_calculator,
    "Pair12": nonbond_calculator,
    "Pair13": nonbond_calculator,
    "Pair14": nonbond_calculator,
    "coul": charge_calculator,
    "vdw": vdw_calculator,
    "constrain": constrain_calculator,
}

class Calculator:

    def __init__(self, style="normal"):
        self.style = style

    def run_molecule_energy(self, molecule, ff_parameters=None):
        """
        calculation the energy of a molecule:
            include: bond, angle, dihedral, improper, intra pair14, pair1n, pair13, pair12
        """

        energy = {"total":{"operator":[],"value":0.0}}
        this_terms = [term for term in molecule.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]
        for term in this_terms:
            if term in _Pair_terms:
                if term in ["Pair14","Pair1n"]:
                    energy[f"{term}-charge"] = {"operator":[],"value":0.0}
                    energy[f"{term}-vdw"] = {"operator":[],"value":0.0}
                    for item in getattr(molecule, term):
                        wow = Calc[term]("coul_vdw",["coul", item.pstyle],item.value,[item.charge_parameter, item.parameter],
                                                combination_rule=item.combination_rule,scale_factor=item.scale_factor,)
                        wow()
                        energy[f"{term}-charge"]["value"] += wow.charge_value
                        energy[f"{term}-vdw"]["value"] += wow.vdw_value
                        energy["total"]["value"] += wow.charge_value + wow.vdw_value
            else:
                energy[term] = {"operator":[],"value":0.0}
                for item in getattr(molecule, term):
                    if ff_parameters is not None:
                        item_name = f"{_topol_to_ff_term[term]}-{item.atom_type_used_name}"
                        if  item_name in ff_parameters:
                            energy[term]["operator"].append([term,item.pstyle,item.value,item.parameter,ff_parameters[item_name],1.0])
                            energy["total"]["operator"].append([term,item.pstyle,item.value,item.parameter,ff_parameters[item_name],1.0])
                            continue

                    wow = Calc[term](item.pstyle, item.value, item.parameter)
                    wow()
                    energy[term]["value"] += wow.value
                    energy["total"]["value"] += wow.value
        if ff_parameters is None:
            _energy = deepcopy(energy)
            energy = {}
            for aa,bb in _energy.items():
                energy[aa] = bb["value"]
        return energy

    def run_molecule_force(self, molecule, method="numerical",ff_parameters=None):
        """
        calculate the force of molecule.
            each atom in Molecule.Atoms has three force value (x, y, z)
        """
        _total_force = []
        if method == "numerical":
            for i in range(len(molecule.Atoms)):
                for j in range(0, 3):
                    vv = molecule.Atoms[i].coor[j]
                    molecule.Atoms[i].coor[j] = vv - 0.005
                    molecule.update_topol_value()
                    e1 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                    molecule.Atoms[i].coor[j] = vv + 0.005
                    molecule.update_topol_value()
                    e2 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                    _total_force.append([e1,e2])
                    #total_force.append((e1 - e2) / 0.01)
                    molecule.Atoms[i].coor[j] = vv
            molecule.update_topol_value()
        if ff_parameters is None:
            total_force = [(rr[0] - rr[1]) / 0.01 for rr in _total_force]
        else:
            total_force = []
            for rr in _total_force:
                vv = (rr[0]["value"] - rr[1]["value"]) / 0.01
                _tmp = [[rrr[0],rrr[1],rrr[2],rrr[3],rrr[4],100.0] for rrr in rr[0]["operator"]] + [[rrr[0],rrr[1],rrr[2],rrr[3],rrr[4],-100.0] for rrr in rr[1]["operator"]]
                total_force.append([vv,_tmp])

        return total_force

    def run_molecule_hessian(self, molecule, method="numerical",ff_parameters=None):
        """
        calculate the hessian of molecule:
            each atom in Molecule.Atoms has three*three value (xx,xy,xz,yx,yy,yz,zx,zy,zz)
        """
        _total_hessian = []
        if method == "numerical":
            e0 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
            for i in range(len(molecule.Atoms)):
                for j in range(0, 3):
                    for k in range(0, i + 1):
                        for l in range(0, 3):  # noqa
                            if i == k and j == l:
                                vv = molecule.Atoms[i].coor[j]
                                molecule.Atoms[i].coor[j] = vv - 0.0001
                                molecule.update_topol_value()
                                e1 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                                molecule.Atoms[i].coor[j] = vv + 0.0001
                                molecule.update_topol_value()
                                e2 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                                #_total_hessian.append((e1 + e2 - 2 * e0) / 0.0001 / 0.0001)
                                _total_hessian.append([e1,e2,e0,e0])
                                molecule.Atoms[i].coor[j] = vv
                                break
                            else:
                                vv1 = molecule.Atoms[i].coor[j]
                                vv2 = molecule.Atoms[k].coor[l]
                                molecule.Atoms[i].coor[j] = vv1 + 0.0001
                                molecule.update_topol_value()
                                e2 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                                molecule.Atoms[k].coor[l] = vv2 + 0.0001
                                molecule.update_topol_value()
                                e1 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                                molecule.Atoms[i].coor[j] = vv1
                                molecule.update_topol_value()
                                e3 = self.run_molecule_energy(molecule,ff_parameters=ff_parameters)["total"]
                                #_total_hessian.append((e0 + e1 - e2 - e3) / 0.0001 / 0.0001)
                                _total_hessian.append([e0,e1,e2,e3])
                                molecule.Atoms[k].coor[l] = vv2
            molecule.update_topol_value()
        if ff_parameters is None:
            total_hessian = [(rr[0] + rr[1] - rr[2] - rr[3]) / 0.0001 / 0.0001 for rr in _total_hessian]
        else:
            total_hessian = []
            for rr in _total_hessian:
                vv = (rr[0]["value"] + rr[1]["value"] - rr[2]["value"] - rr[3]["value"]) / 0.0001 / 0.0001 
                _tmp = [[rrr[0],rrr[1],rrr[2],rrr[3],rrr[4],1.0e8] for rrr in rr[0]["operator"]] 
                _tmp += [[rrr[0],rrr[1],rrr[2],rrr[3],rrr[4],1.0e8] for rrr in rr[1]["operator"]]
                _tmp += [[rrr[0],rrr[1],rrr[2],rrr[3],rrr[4],-1.0e8] for rrr in rr[2]["operator"]]
                _tmp += [[rrr[0],rrr[1],rrr[2],rrr[3],rrr[4],-1.0e8] for rrr in rr[3]["operator"]]
                total_hessian.append([vv,_tmp])
        return total_hessian

    def run_molecule_freq(self, molecule):
        """
        calculate frequency from hessian
        """
        n = int(len(molecule.Atoms) * 3)
        arr = []
        for i in range(0, n):
            e = int((i * i + i) / 2)
            arr.append([])
            for j in range(0, i + 1):
                arr[i].append(molecule.hessian[e + j])
            for j in range(i + 1, n):
                arr[i].append(0)
        for i in range(0, n):
            for j in range(i + 1, n):
                arr[i][j] = arr[j][i]
        mass_arr = []
        for ai in molecule.Atoms:
            mass_arr.append(ai.mass)
        for i in range(0, n):
            ii = int(i / 3)
            for j in range(0, n):
                jj = int(j / 3)
                wei = (mass_arr[ii] * mass_arr[jj]) ** 0.5
                arr[i][j] = arr[i][j] / wei

        mx = np.matrix(arr)
        tt, qq = np.linalg.eig(mx)
        for i in range(len(tt)):
            if tt[i] < 0:
                tt[i] = -1 * (abs(tt[i])) ** 0.5 * 108.77  # *5140.15
            else:
                tt[i] = (tt[i]) ** 0.5 * 108.77  # *5140.15
        tt = sorted(tt)
        tt = tt[6:]
        return tt

    def _run_molecule_energy(self,molecule,ff_parameters=None,idx=None):
        return self.run_molecule_energy(molecule,ff_parameters=ff_parameters), idx

    def _run_molecule_force(self,molecule,ff_parameters=None,method="numerical",idx=None):
        return self.run_molecule_force(molecule,ff_parameters=ff_parameters,method=method), idx

    def _run_molecule_hessian(self,molecule,ff_parameters=None,method="numerical",idx=None):
        return self.run_molecule_hessian(molecule,ff_parameters=ff_parameters,method=method), idx

    def _run_molecule_freq(self,molecule,idx=None):
        return self.run_molecule_freq(molecule), idx

    def molecule_energy(self, molecules,ff_parameters=None,parallel=True):
        """
        一系列分子，分子内能量的计算
        输入：
            moles: List[Molecule]
        输出：
            intra_energy： List[float],能量
        """
        if parallel:
            energys = parallel_run("_run_molecule_energy",molecules,kwds={"ff_parameters":ff_parameters},objs=self,single_args_flag=False)
        else:
            energys = []
            for molecule in molecules:
                energys.append(self.run_molecule_energy(molecule),ff_parameters=ff_parameters)
        return energys

    def molecule_force(self, molecules, ff_parameters=None, method="numerical",parallel=True):
        """
        forces
        """
        if parallel:
            forces = parallel_run("_run_molecule_force",molecules,kwds={"method":method,"ff_parameters":ff_parameters},objs=self,single_args_flag=False)
        else:
            forces=[]
            for molecule in molecules:
                forces.append(self.run_molecule_force(molecule,method=method,ff_parameters=ff_parameters))

        return forces

    def molecule_hessian(self, molecules, ff_parameters=None,method="numerical",parallel=True):
        """
        hessians
        """
        if parallel:
            hessians = parallel_run("_run_molecule_hessian",molecules,kwds={"method":method,"ff_parameters":ff_parameters},objs=self,single_args_flag=False)
        else:
            hessians=[]
            for molecule in molecules:
                hessians.append(self.run_molecule_hessian(molecule,method=method))

        return hessians

    def molecule_freq(self, molecules,parallel=True):
        """
        frequency
        """
        if parallel:
            freqs = parallel_run("_run_molecule_freq",molecules,objs=self,single_args_flag=False)
        else:
            freqs = []
            for molecule in molecules:
                freqs.append(self.run_molecule_freq(molecule))
        return freqs

    def single_mole_pes(self, m, x):
        # TODO CFL 还有没完善
        this_terms = [term for term in m.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]
        energy = 0
        total_e = 0.0
        for term in this_terms:
            if term in ["Pair1n", "Pair12", "Pair13", "Pair14"]:
                energy[term] = [0.0, 0.0]
                for b in getattr(m, term):
                    # this_coord = [[], []]
                    # value = b.calc_value()
                    oo = self.__CalcTerm[term](
                        "coul_vdw", ["coul", b.pstyle], "energy", b.value, [b.charge_parameter, b.parameter]
                    )
                    oo()
                    energy[term][0] += oo.charge_value
                    energy[term][1] += oo.vdw_value
                    total_e += oo.charge_value
                    total_e += oo.vdw_value
            else:
                energy[term] = 0
                for b in getattr(m, term):
                    oo = self.__CalcTerm[term](b.pstyle, "energy", b.value, b.parameter)
                    oo()
                    energy[term] += oo.value
                    total_e += oo.value
        energy["total"] = total_e

    def system_energy(self, m):
        """
        系统中某一个分子的能量计算。包含有分子间的贡献。
        输入：
            m: Molecule，
        输出：

        """

        this_terms = [term for term in m.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]
        for term in this_terms:
            if term in ["Pair1n", "Pair12", "Pair13", "Pair14", "InterPair"]:
                for b in getattr(m, term):
                    oo = self.__CalcTerm[term](
                        "coul_vdw",
                        ["coul", b.pstyle],
                        "energy",
                        b.value,
                        [b.charge_parameter, b.parameter],
                        outtype=["value", "value"],
                        combination_rule=b.combination_rule,
                        scale_factor=b.scale_factor,
                    )
                    oo()
                    b.charge_energy = oo.charge_value
                    b.charge_energy = oo.vdw_value
                    b.energy = oo.charge_value + oo.vdw_value
            else:
                for b in getattr(m, term):
                    oo = self.__CalcTerm[term](b.pstyle, "energy", b.value, b.parameter)
                    oo()
                    b.energy = oo.value

