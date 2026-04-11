import numpy as np
from copy import deepcopy

from scipy.optimize import least_squares, leastsq
from scipy.optimize.optimize import OptimizeResult

from ...mm_calculator.mm.ff_calculator.nonbond import charge_calculator, vdw_calculator, nonbond_calculator
from ...mm_calculator.mm.ff_calculator.bonded import bond_calculator, angle_calculator, dihedral_calculator, improper_calculator, constrain_calculator
#from ...mm_calculator.mm.calculator import Calculator


class torsion_penalty:
    def __init__(self, funcid, x, para, style="energy"):
        self.para = para
        self.x = x
        self.style = "normal"

    def calculation(self):
        return self.x * (self.para[0]*self.para[0] + self.para[2]*self.para[2] + self.para[4]*self.para[4] + self.para[6]*self.para[6])

    __Func_ID = {
        "normal": calculation
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)

This_Calc = {
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
    "torsion_penalty":torsion_penalty,
}

class FFFConstructor:
    def __init__(  self,
                   molecules,
                   this_ff,
                   fitting_terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
                   targets=["energy", "force", "hessian","penalty_torsion"],
                   weight_factor=[1.0, 0.0001, 0.0001,0.1],
                   ):
        from ...mm_calculator.mm.calculator import Calculator
        self.molecules = molecules
        self.this_ff = this_ff
        self.fitting_terms = fitting_terms
        self.targets = targets
        self.Calc = Calculator("fitting")
        self.weight_factor = weight_factor


    def parameter_constructor(self):
        parameter_init = []
        fitting_parameter = {}

        bup = []
        bdown = []
        nn = 0
        for term in self.fitting_terms:
            if term in self.this_ff:
                for name,item in self.this_ff[term].items():
                    if item["ptag"] == "isfitting":
                        fitting_parameter[f"{term}-{name}"] = {}
                        for ii in item["isfitting"]:
                            fitting_parameter[f"{term}-{name}"][ii] = nn
                            parameter_init.append(item["parameter"][ii])
                            if term == "dihedralterm":
                                bup.append(max(6.0, abs(parameter_init[-1])))
                                bdown.append(min(-6.0, -abs(parameter_init[-1])))
                            elif term in ["bondterm","angleterm"]:
                                bup.append(parameter_init[-1]*1.05)
                                bdown.append(parameter_init[-1]*0.95)
                            else:
                                bup.append(max(10,abs(parameter_init[-1])))
                                bdown.append(max(0.0,abs(parameter_init[-1])))
                            nn += 1
        self.parameter_init = parameter_init
        self.fitting_parameter = fitting_parameter
        self.bup = bup
        self.bdown = bdown
        #return parameter_init,fitting_parameter,bup,bdown

    def molecule_constructor(self):
        _gm_molecule = []
        mm = {}
        for molecule in self.molecules:
            if molecule.mole_name not in mm.keys():
                mm[molecule.mole_name] = []
            mm[molecule.mole_name].append(molecule)
        for aa,bb in mm.items():
            bb = sorted(bb,key=lambda m:m.energy)
            min_v = bb[0].energy
            for molecule in bb:
                molecule.energy = molecule.energy - min_v
            _gm_molecule.append(deepcopy(bb[0]))
        
        energys = self.Calc.molecule_energy(_gm_molecule,ff_parameters=self.fitting_parameter)
        self.gm_energy = {}
        for ii,molecule in enumerate(_gm_molecule):
            if len(energys[ii]["total"]["operator"]) != 0:
                self.gm_energy[molecule.mole_name] = [energys[ii]["total"]["value"],energys[ii]["total"]["operator"]]

    def get_energy_funcs(self,molecules,weight):
        energys = self.Calc.molecule_energy(molecules,ff_parameters=self.fitting_parameter)
        funcs = []
        for ii,molecule in enumerate(molecules):
            func = [energys[ii]["total"]["value"],energys[ii]["total"]["operator"]]
            mfunc = self.gm_energy[molecule.mole_name]
            func[0] = (func[0] - mfunc[0] - molecule.energy) * weight
            func[1] = [[rr[0],rr[1],rr[2],rr[3],rr[4],rr[5]*weight] for rr in func[1]]
            func[1] += [[rr[0],rr[1],rr[2],rr[3],rr[4],-1*rr[5]*weight] for rr in mfunc[1]]
            funcs.append(func)
        return funcs

    def get_force_hessian_funcs(self,molecules,ttype,weight):
        if ttype == "force":
            total_funcs = self.Calc.molecule_force(molecules,ff_parameters=self.fitting_parameter)
        if ttype == "hessian":
            total_funcs = self.Calc.molecule_hessian(molecules,ff_parameters=self.fitting_parameter)
        for ii,molecule in enumerate(molecules):
            funcs = total_funcs[ii]
            target_funcs = getattr(molecule,ttype)
            for jj,vv in enumerate(funcs):
                if len(vv) != 0:
                    vv[0] = (vv[0] - target_funcs[jj]) * weight
                    for rr in vv[1]:
                        rr[5] = rr[5] * weight
            
        #return total_funcs
        return [r for rr in total_funcs for r in rr ]

    def get_penalty_torsion_func(self,weight):
        pp = []
        for item in self.this_ff["dihedralterm"]:
            if "isfitting" in item:
                ss = item["name"].splie("$")
                xx = 1.0
                if ss[0] == "h_1" or ss[3] == "h_1":
                    xx = 5.0
                elif ss[0].startswith("h_1") or ss[3].startswith("h_1"):
                    xx = 2.0
                pp.append(["torsion_penalty","normal",xx,item["parameter"],self.fitting_parameter[f"dihedralterm-{item['name']}"],weight])
        return [[0.0, pp]]

    def funcs_constructor(self):
        self.molecule_constructor()
        used_molecules = {"energy":[],"force":[],"hessian":[]}
        for molecule in self.molecules:
            if molecule.mole_name in self.gm_energy:
                used_molecules["energy"].append(molecule)
                if hasattr(molecule,"force"):
                    used_molecules["force"].append(molecule)
                if hasattr(molecule,"hessian"):
                    used_molecules["hessian"].append(molecule)
        funcs = []
        for ii, target in enumerate(self.targets):
            if target in ["energy"]:
                funcs.extend(self.get_energy_funcs(used_molecules[target],self.weight_factor[ii]))
            elif target in ["force","hessian"]:
                funcs.extend(self.get_force_hessian_funcs(used_molecules[target],target,self.weight_factor[ii]))
            elif target == "penalty_torsion":
                funcs.extend(self.get_penalty_torsion_func(self.weight_factor[ii]))
        self.funcs = funcs

    def run(self):
        self.parameter_constructor()
        self.funcs_constructor()
        return self.funcs, self.parameter_init,self.fitting_parameter,(self.bdown,self.bup)

class st:
    """
    least_squares方法拟合参数
    """

    def __init__(
        self,
        funcs,
        parameter_init,
        fitting_parameter,
        boundary,
        ):
        """
        molecules: fitting molecules
        parameter_init: inital parameter
        fitting_parameter: the fitting parameters
        gm_energy: the local minimum structure for each molecule
        bup, bdown: the up or down limit of parameters

        """

        self.funcs = funcs
        self.parameter_init = np.asarray(parameter_init)
        self.fitting_parameter = fitting_parameter
        self.boundary = boundary
        

    def residue_intra(self,para):
        x = para[:]
        residue = np.zeros(len(self.funcs), dtype=np.float64)
        for ii,func in enumerate(self.funcs):
            residue[ii] += func[0]
            for rr in func[1]:
                #[term,item.style,item.value,item.parameter,ff_parameters[item_name],1.0]
                this_para = [vv if kk not in rr[4] else x[rr[4][kk]] for kk,vv in enumerate(rr[3])]
                wow = This_Calc[rr[0]](rr[1], rr[2], this_para)
                wow()
                v = wow.value * rr[5]
                residue[ii] += v
        return residue

    def run(self) -> OptimizeResult:
        """
        运行least_squares
        """
        #jac=self.jac_func_intra,
        root = least_squares(
                self.residue_intra,
                self.parameter_init,
                bounds=self.boundary,
                ftol=1e-5,
                max_nfev=100,
                verbose=2,
            )
        return root