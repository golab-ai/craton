import math
#import scipy.optimize
import scipy.optimize as opt
try:
    from simtk import openmm as mm
    from simtk import unit

    _OPENMM_IMPORTED = True
except ImportError:
    _OPENMM_IMPORTED = False

from ...utils import logger
from .intra_kernel_np import *
from .calculator import Calculator, Calc


class MMCalc(Calculator):
    def __init__(self, style="normal"):
        super().__init__(style)

    def energy_pes(self,x):
        #__CalcTerm = {
        #    "Bonds": ff_calc.bond_calculator,
        #    "Angles": ff_calc.angle_calculator,
        #    "Dihedrals": ff_calc.dihedral_calculator,
        #    "Impropers": ff_calc.improper_calculator,
        #    "Pair1n": ff_calc.nonbond_calculator,
        #    "Pair12": ff_calc.nonbond_calculator,
        #    "Pair13": ff_calc.nonbond_calculator,
        ##    "Pair14": ff_calc.nonbond_calculator,
        #    "coul": ff_calc.charge_calculator,
        #    "vdw": ff_calc.vdw_calculator,
        #    "constrain": ff_calc.constrain_calculator,
        #}
        m = MOLE
        this_terms = [term for term in m.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]
        total_e = 0
        for term in this_terms:
            if term in ["Pair1n", "Pair12", "Pair13", "Pair14"]:
                for b in getattr(m, term):
                    this_coord = [
                        [x[b.a1 * 3], x[b.a1 * 3 + 1], x[b.a1 * 3 + 2]],
                        [x[b.a2 * 3], x[b.a2 * 3 + 1], x[b.a2 * 3 + 2]],
                    ]
                    b.calc_value(this_coord)
                    oo = Calc[term](
                        "coul_vdw",
                        ["coul", b.style],
                        b.value,
                        [b.charge_parameter, b.parameter],
                        style="energy"
                    )
                    oo()
                    total_e += oo.charge_value
                    total_e += oo.vdw_value
            else:
                for b in getattr(m, term):
                    this_coor = []
                    if hasattr(b, "a1"):
                        this_coor.append([x[b.a1 * 3], x[b.a1 * 3 + 1], x[b.a1 * 3 + 2]])
                    if hasattr(b, "a2"):
                        this_coor.append([x[b.a2 * 3], x[b.a2 * 3 + 1], x[b.a2 * 3 + 2]])
                    if hasattr(b, "a3"):
                        this_coor.append([x[b.a3 * 3], x[b.a3 * 3 + 1], x[b.a3 * 3 + 2]])
                    if hasattr(b, "a4"):
                        this_coor.append([x[b.a4 * 3], x[b.a4 * 3 + 1], x[b.a4 * 3 + 2]])
                    b.calc_value(this_coor)
                    oo = Calc[term](b.pstyle, b.value, b.parameter,style="energy")
                    oo()
                    total_e += oo.value

        if hasattr(m, "constrain"):
            for b in m.constrain:
                this_coor = []
                if hasattr(b, "a1"):
                    this_coor.append([x[b.a1 * 3], x[b.a1 * 3 + 1], x[b.a1 * 3 + 2]])
                if hasattr(b, "a2"):
                    this_coor.append([x[b.a2 * 3], x[b.a2 * 3 + 1], x[b.a2 * 3 + 2]])
                if hasattr(b, "a3"):
                    this_coor.append([x[b.a3 * 3], x[b.a3 * 3 + 1], x[b.a3 * 3 + 2]])
                if hasattr(b, "a4"):
                    this_coor.append([x[b.a4 * 3], x[b.a4 * 3 + 1], x[b.a4 * 3 + 2]])
                b.calc_value(this_coor)
                oo = Calc["constrain"](b.style, b.value, b.fix_value)
                oo()
                total_e += oo.value
        return total_e

    def minimize(self,x, m, method="fmin_bfgs"):
        global MOLE
        MOLE = m
        fmin_func = opt.__dict__[method]
        if method in ["fmin", "fmin_powell"]:
            result = fmin_func(self.energy_pes, x, maxiter=200)
        elif method in ["fmin_cg", "fmin_bfgs", "fmin_l_bfgs_b", "fmin_tnc"]:
            result = fmin_func(self.energy_pes, x, maxiter=200, disp=False)
        elif method in ["fmin_cobyla"]:
            result = fmin_func(self.energy_pes, x, [], maxiter=200)
        else:
            logger.error("fmin function not found")
            return
        return result

    def optimize(self,molecule,index=None):
        """
        优化结构实际运行的方法
        输入：
            mole: Molecule, 优化的分子
            index: int, 顺序编号，用于并行后，不打乱顺序
        输出：
            moles: Molecule, 优化后的分子
            index: int, 顺序编号，用于并行
        """

        x = [a for aa in molecule.Atoms for a in aa.coor]
        xx = self.minimize(x, molecule)
        for ii,atom in enumerate(molecule.Atoms):
            atom.coor = [xx[ii * 3], xx[ii * 3 + 1], xx[ii * 3 + 2]]
        # mole.update_topol_value()
        # logger.info(f"{mole.Atoms[-1].coor}, {mole.Atoms[0].coor}, {index}"
        return molecule, index

class MMCalc_:
    def __init__(self, moles):
        self.moles = moles
        self.n_atom_list = [len(m.Atoms) for m in moles]
        self.n_atom_offset = [sum(self.n_atom_list[:i]) for i in range(len(moles))]

        positions = []
        for mol in moles:
            positions += mol.coordinates
            mol.update_topol_value()

        self._positions = np.array(positions, dtype=np.float64)
        self._energy = 0
        self._forces = np.zeros(self._positions.shape, dtype=np.float64)

    def set_positions(self, positions):
        self._positions = positions[:]
        for i, mol in enumerate(self.moles):
            offset = self.n_atom_offset[i]
            for j, atom in enumerate(mol.Atoms):
                atom.coor = positions[offset + j].tolist()
            mol.update_topol_value()

        self._update_internal_positions()

    def _update_internal_positions(self):
        pass

    def calc_energy(self):
        """
        Extremely slow. Only for test purpose.
        Use subclasses OpenMMCalc and NumpyMMCalc instead in real scenarios.
        """
        energy = 0
        calc = Calculator()
        for mol in self.moles:
            energy += calc.single_mole_energy(mol)["total"]
        return energy

    def calc_forces(self):
        raise Exception("Not implemented")

    def calc_hessian(self, output_for_fitting=False):
        """
        TODO Haven't been tested

        """
        dof = self._positions.shape[0] * 3
        self._hessian = np.zeros([dof, dof], dtype=np.float64)

        delta = 0.001
        forces0 = self.calc_forces().copy().flatten()
        n_atom = self._positions.shape[0]
        for i in range(n_atom):
            for k in range(3):
                self._positions[i][k] += delta
                self._update_internal_positions()
                forces1 = self.calc_forces().flatten()
                self._hessian[i * 3 + k] = (forces1 - forces0) / delta * -1
                self._positions[i][k] -= delta

        if not output_for_fitting:
            return self._hessian

        hessian = []
        for n_atom, offset in zip(self.n_atom_list, self.n_atom_offset):
            for i1 in range(n_atom):
                for k1 in range(3):
                    for i2 in range(i1 + 1):
                        for k2 in range(3):
                            hessian.append(self._hessian[3 * (offset + i1) + k1][3 * (offset + i2) + k2])
                            if i1 == i2 and k1 == k2:
                                break
        return hessian

    def optimize(self, constrain_all_torsion=False, **kwargs):
        raise Exception("Not implemented")


class NumpyMMCalc(MMCalc_):
    def __init__(self, moles):
        super().__init__(moles)
        self._kernels = self._setup_kernels(moles)

    def _setup_kernels(self, moles):
        kernels = []

        ###bond interaction
        indexes = []
        parameters = []
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for bond in getattr(mol, "Bonds", []):
                indexes.append([bond.a1 + offset, bond.a2 + offset])
                parameters.append([bond.parameter[0], bond.parameter[1]])
        if indexes != []:
            kernel = HarmonicBondKernel(indexes, parameters)
            kernels.append(kernel)

        ###angle interaction
        indexes.clear()
        parameters.clear()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for angle in getattr(mol, "Angles", []):
                indexes.append([angle.a1 + offset, angle.a2 + offset, angle.a3 + offset])
                parameters.append([angle.parameter[0], angle.parameter[1]])
        if indexes != []:
            kernel = HarmonicAngleKernel(indexes, parameters)
            kernels.append(kernel)

        ###dihedral interaction
        indexes.clear()
        parameters.clear()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for dihedral in getattr(mol, "Dihedrals", []):
                indexes.append([dihedral.a1 + offset, dihedral.a2 + offset, dihedral.a3 + offset, dihedral.a4 + offset])
                parameters.append(dihedral.parameter)
        if indexes != []:
            kernel = OplsTorsionKernel(indexes, parameters)
            kernels.append(kernel)

        ###improper interaction
        indexes.clear()
        parameters.clear()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for improper in getattr(mol, "Impropers", []):
                indexes.append([improper.a1 + offset, improper.a2 + offset, improper.a3 + offset, improper.a4 + offset])
                parameters.append([0, improper.parameter[0], 0, 0])
        if indexes != []:
            kernel = OplsTorsionKernel(indexes, parameters)
            kernels.append(kernel)

        ###nonbonded interaction
        indexes.clear()
        parameters.clear()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for pair in getattr(mol, "Pair14", []) + getattr(mol, "Pair1n", []):
                indexes.append([pair.a1 + offset, pair.a2 + offset])
                parameters.append([pair.parameter[1], pair.parameter[0], pair.charge_parameter[0] * pair.charge_parameter[1]])
        if indexes != []:
            kernel = NonbondedKernel(indexes, parameters)
            kernels.append(kernel)

        ###constrain interaction
        indexes.clear()
        parameters.clear()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for cons in getattr(mol, "constrain", []):
                indexes.append([n + offset for n in cons.atoms])
                parameters.append([cons.fix_value, 1500])
        if indexes != []:
            kernel = ConstrainedTorsionKernel(indexes, parameters)
            kernels.append(kernel)

        return kernels

    def calc_energy(self):
        self._energy = 0
        for kernel in self._kernels:
            v, e, f = kernel.evaluate(self._positions)
            self._energy += e.sum()
        return self._energy

    def calc_forces(self):
        self._forces.fill(0)
        for kernel in self._kernels:
            v, e, f = kernel.evaluate(self._positions, force=True)
            self._forces += f
        return self._forces

    def optimize(self, all_torsion_constraint=0.0, tol_force=0.1, max_iter=200, method="CG", verbose=True):
        if verbose:
            logger.debug("Potential energy before optimization: %.2f kcal/mol" % self.calc_energy())

        if all_torsion_constraint:
            indexes = []
            parameters = []
            for i, mol in enumerate(self.moles):
                offset = self.n_atom_offset[i]
                for dihedral in getattr(mol, "Dihedrals", []):
                    if dihedral.is_linear:
                        continue
                    indexes.append(
                        [dihedral.a1 + offset, dihedral.a2 + offset, dihedral.a3 + offset, dihedral.a4 + offset]
                    )
                    parameters.append([dihedral.value_a, all_torsion_constraint])
            if indexes != []:
                kernel = ConstrainedTorsionKernel(indexes, parameters)
                self._kernels.append(kernel)

        self._force_tolerance = tol_force
        self._max_iter = max_iter
        self._n_iter = 0
        self._n_evaluation = 0
        try:
            #result = scipy.optimize.minimize(
            result = opt.minimize(
                self._func_opt,
                self._positions.flatten(),
                method=method,
                jac=True,
                callback=self._func_callback,
            )
        except StopIteration:
            pass
        else:
            self._positions = result.x.reshape(-1, 3)

        if all_torsion_constraint:
            self._kernels.pop()

        self.set_positions(self._positions)

        if verbose:
            logger.debug("Potential energy after optimization: %.2f kcal/mol" % self.calc_energy())

    def _func_opt(self, positions):
        self._n_evaluation += 1
        self._positions = positions.reshape(-1, 3)
        energy = self.calc_energy()
        forces = self.calc_forces()
        return energy, -forces.flatten()

    def _func_callback(self, xk):
        self._n_iter += 1
        rmsf = np.sqrt((self._forces * self._forces).mean())
        logger.debug(
            "Number of iter/eval: %i/%i  Energy: %.2f kcal/mol  RMS force: %.2f kcal/mol/A"
            % (self._n_iter, self._n_evaluation, self._energy, rmsf)
        )
        if rmsf < self._force_tolerance or self._n_iter >= self._max_iter:
            self._positions = xk.reshape(-1, 3)
            raise StopIteration()


class OpenMMCalc(MMCalc_):
    def __init__(self, moles):
        if not _OPENMM_IMPORTED:
            raise Exception("OpenMM not found. Install it with `conda install -c omnia openmm`")

        super().__init__(moles)
        self._system, self._context = self._setup_openmm(moles)

    @property
    def omm_system(self):
        return self._system

    def _setup_openmm(self, moles):
        system = mm.System()
        for _, mol in enumerate(moles):
            for atom in mol.Atoms:
                system.addParticle(atom.mass)

        ###bond interaction
        bforce = mm.HarmonicBondForce()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for bond in getattr(mol, "Bonds", []):
                bforce.addBond(
                    bond.a1 + offset,
                    bond.a2 + offset,
                    bond.parameter[0] / 10,  # A -> nm
                    bond.parameter[1] * 2 * 4.184 * 100,  # kcal/mol/A^2 -> kJ/mol/nm^2
                )
        system.addForce(bforce)

        ###angle interaction
        aforce = mm.HarmonicAngleForce()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for angle in getattr(mol, "Angles", []):
                aforce.addAngle(
                    angle.a1 + offset,
                    angle.a2 + offset,
                    angle.a3 + offset,
                    math.radians(angle.parameter[0]),  # radian
                    angle.parameter[1] * 2 * 4.184,  # kcal/mol -> kJ/mol
                )
        system.addForce(aforce)

        ###dihedral interaction
        dforce = mm.PeriodicTorsionForce()
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for dihedral in getattr(mol, "Dihedrals", []):
                for n in range(4):
                    k = dihedral.parameter[n*2]
                    if k == 0:
                        continue
                    #dforce.addTorsion(
                    #    dihedral.a1 + offset,
                    #    dihedral.a2 + offset,
                    #    dihedral.a3 + offset,
                    #    dihedral.a4 + offset,
                    #    n + 1,
                    #    (n % 2) * math.pi,
                    #    k * 4.184,  # kcal/mol -> kJ/mol
                    #)
                    dforce.addTorsion(
                        dihedral.a1 + offset,
                        dihedral.a2 + offset,
                        dihedral.a3 + offset,
                        dihedral.a4 + offset,
                        n + 1,
                        math.radians(dihedral.parameter[n*2+1]),
                        k * 4.184 # kcal/mol -> kJ/mol
                    )
        system.addForce(dforce)

        ###improper interaction
        iforce = mm.CustomTorsionForce("k*(1-cos(2*theta))")
        iforce.addPerTorsionParameter("k")
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for improper in getattr(mol, "Impropers", []):
                iforce.addTorsion(
                    improper.a1 + offset,
                    improper.a2 + offset,
                    improper.a3 + offset,
                    improper.a4 + offset,
                    [improper.parameter[0] * 4.184],  # kcal/mol -> kJ/mol
                )
        system.addForce(iforce)

        #### nonbond interaction
        nforce = mm.CustomBondForce("4*eps*(c^2-c) + 138.935458*qq/r;" "c=(sig/r)^6;")
        nforce.addPerBondParameter("eps")
        nforce.addPerBondParameter("sig")
        nforce.addPerBondParameter("qq")
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for pair in getattr(mol, "Pair14", []) + getattr(mol, "Pair1n", []):
                nforce.addBond(
                    pair.a1 + offset,
                    pair.a2 + offset,
                    [
                        pair.parameter[1] * 4.184,  # kcal/mol -> kJ/mol
                        pair.parameter[0] / 10,  # A -> nm
                        pair.charge_parameter[0] * pair.charge_parameter[1],  #
                    ],
                )
        system.addForce(nforce)

        ###constrain
        cforce = mm.CustomTorsionForce("k*(1-cos(theta-theta0))")
        cforce.addPerTorsionParameter("theta0")
        cforce.addPerTorsionParameter("k")
        system.addForce(cforce)
        for i, mol in enumerate(moles):
            offset = self.n_atom_offset[i]
            for cons in getattr(mol, "constrain", []):
                cforce.addTorsion(
                    *[n + offset for n in cons.atoms],
                    [math.radians(cons.fix_value), 1500 * 4.184],  # kcal/mol -> kJ/mol
                )

        integrator = mm.VerletIntegrator(0.001)
        try:
            platform = mm.Platform.getPlatformByName("CPU")
        except:  # paracloud has no CPU platform
            platform = mm.Platform.getPlatformByName("Reference")
        context = mm.Context(system, integrator, platform)
        context.setPositions(self._positions / 10)  # A -> nm

        return system, context

    def _update_internal_positions(self):
        self._context.setPositions(self._positions / 10)  # A -> nm

    def calc_energy(self):
        state = self._context.getState(getEnergy=True)
        self._energy = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        return self._energy

    def calc_forces(self):
        state = self._context.getState(getForces=True)
        self._forces = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole / unit.angstrom)
        return self._forces

    def optimize(self, all_torsion_constraint: float = 0, tol_force=0.1, max_iter=200, verbose=False):
        if verbose:
            logger.debug("Potential energy before optimization: %.2f kcal/mol" % self.calc_energy())

        idx_cforce = -1
        if all_torsion_constraint:
            cforce = mm.CustomTorsionForce("k*(1-cos(theta-theta0))")
            cforce.addPerTorsionParameter("theta0")
            cforce.addPerTorsionParameter("k")
            idx_cforce = self._system.addForce(cforce)
            for i, mol in enumerate(self.moles):
                offset = self.n_atom_offset[i]
                for dihedral in getattr(mol, "Dihedrals", []):
                    if dihedral.is_linear:
                        continue
                    cforce.addTorsion(
                        dihedral.a1 + offset,
                        dihedral.a2 + offset,
                        dihedral.a3 + offset,
                        dihedral.a4 + offset,
                        [dihedral.value, all_torsion_constraint * 4.184],  # kcal/mol -> kJ/mol
                    )
            self._context.reinitialize(preserveState=True)

        mm.LocalEnergyMinimizer.minimize(
            self._context, tolerance=tol_force * 41.84, maxIterations=max_iter  # kcal/mol/A -> kJ/mol/nm
        )

        if all_torsion_constraint:
            self._system.removeForce(idx_cforce)
            self._context.reinitialize(preserveState=True)

        state = self._context.getState(getEnergy=True, getPositions=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        self.set_positions(positions)

        if verbose:
            logger.debug("Potential energy after optimization: %.2f kcal/mol" % self.calc_energy())

