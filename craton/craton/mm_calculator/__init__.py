from copy import deepcopy
import numpy as np
from .am1bcc_calculator import get_am1bcc
from .optimizer import optimize
from .mm.calculator import Calculator as mmcalc
from .molecule_geometry import calculate_moment_of_inertia, calculate_multipole_moments, find_center
from .molecule_geometry import MolecularSurfaceCalculator as MSC
from ..utils.commons import parallel_run

class Calculator:
    def __init__(self):
        pass

    @staticmethod
    def am1bcc_charge(molecule):
        get_am1bcc(molecule)

    @staticmethod
    def _optimize(molecules, optimizer="openmm", all_torsion_constraint = 0.0, write_mol=None,):
        
        optimize(molecules,optimizer=optimizer,all_torsion_constraint=all_torsion_constraint,write_mol=write_mol)
        return molecules

    @staticmethod
    def _energy(molecules,prop="energy",parallel=True):
        cal = mmcalc()

        energys = cal.molecule_energy(molecules,parallel=parallel)
        for ii,molecule in enumerate(molecules):
            molecule.energy = energys[ii]["total"]
            molecule._mm_energy = energys[ii]

        if prop == "freq":
            forces = cal.molecule_force(molecules,parallel=parallel)
            for ii,molecule in enumerate(molecules):
                molecule.force = forces[ii]

            hessians = cal.molecule_hessian(molecules,parallel=parallel)
            for ii,molecule in enumerate(molecules):
                molecule.hessian = hessians[ii]

            freqs = cal.molecule_freq(molecules,parallel=parallel)
            for ii,molecule in enumerate(molecules):
                molecule.frequency = freqs[ii]

        if prop == "hessian":
            forces = cal.molecule_force(molecules,parallel=parallel)
            for ii,molecule in enumerate(molecules):
                molecule.force = forces[ii]

            hessians = cal.molecule_hessian(molecules,parallel=parallel)
            for ii,molecule in enumerate(molecules):
                molecule.hessian = hessians[ii]

        if prop == "force":
            forces = cal.molecule_force(molecules,parallel=parallel)
            for ii,molecule in enumerate(molecules):
                molecule.force = forces[ii]
        return molecules

    def _mix_energy(ts_molecules,terms=None,parallel=True,):
        calc = mmcalc()
        if terms is None:
            terms = ["energy"]

        molecules = deepcopy(ts_molecules)
        if "pes" in terms or "energy" in terms:
            energy = calc.molecule_energy(molecules)
            for i in range(len(energy)):
                molecules[i].energy = energy[i]["total"]

        if "force" in terms:
            molecules_force = [m for m in molecules if hasattr(m, "force")]
            if len(molecules_force) > 0:
                force = calc.molecule_force(molecules_force)
                for i in range(len(molecules_force)):
                    molecules_force[i].force = force[i]

        if "hessian" in terms or "freq" in terms:
            molecules_freq = [m for m in molecules if hasattr(m, "freq")]
            if len(molecules_freq) > 0:
                hessian = calc.molecule_hessian(molecules_freq)
                for i in range(len(molecules_freq)):
                    molecules_freq[i].hessian = hessian[i]
                freq = calc.molecule_freq(molecules_freq)
                for i in range(len(molecules_freq)):
                    molecules_freq[i].freq = freq[i]

        return molecules

    @staticmethod
    def _torsion_scan(molecules, scan_interval=30,parallel=True):
        from ..chemkit import MolEdit as ME
        from ..chemkit import MolConformer as MConf
        from copy import deepcopy
        if not isinstance(molecules,list):
            molecules = [molecules]

        total_molecule = []
        for molecule in molecules:
            for torsion in molecule.torsions:
                for angle in range(-180, 180, scan_interval):
                    total_molecule.append(deepcopy(molecule))
                    ME._structure_change(total_molecule[-1],torsion,angle)
                    total_molecule[-1].create_constrain([torsion + [angle]])
        total_molecule = Calculator._optimize(total_molecule)
        total_molecule = Calculator._energy(total_molecule,parallel=parallel)
        
        return MConf._scan_curve(total_molecule)
        
class MolGeo:
    def __init__(self,molecules,parallel=True):
        self.molecules = molecules
        self.parallel = parallel


    def run_volume_surface(self,molecule,idx=None):
        molecule = MSC(molecule).run()
        if idx is None:
            return molecule
        else:
            return molecule,idx

    def volume_surface(self):
        if self.parallel:
            molecules = parallel_run(self.run_volume_surface,self.molecules)
        else:
            molecules = []
            for molecule in self.molecules:
                molecules.append(self.run_volume_surface(molecule))
        return molecules

    def run_moment_of_inertia(self,molecule,idx=None):
        molecule = calculate_moment_of_inertia(molecule)
        if idx is None:
            return molecule
        else:
            return molecule,idx
            
    def moment_of_inertia(self):
        if self.parallel:
            molecules = parallel_run(self.run_moment_of_inertia,self.molecules)
        else:
            molecules = []
            for molecule in self.molecules:
                molecules.append(self.run_moment_of_inertia(molecule))
        return molecules
    
    def run_multipole_moments(self,molecule,idx=None):
        charges = [atom.ff_charge for atom in molecule.Atoms]
        molecule = calculate_multipole_moments(molecule, charges)
        if idx is None:
            return molecule
        else:
            return molecule,idx
        
    def multipole_moments(self):
        if self.parallel:
            molecules = parallel_run(self.run_multipole_moments,self.molecules)
        else:
            molecules = []
            for molecule in self.molecules:
                molecules.append(self.run_multipole_moments(molecule))
        return molecules

    def run_molecule_center(self,molecule,idx=None):
        molecule = find_center(molecule)
        if idx is None:
            return molecule
        else:
            return molecule,idx
        
    def molecule_center(self):
        if self.parallel:
            molecules = parallel_run(self.run_molecule_center,self.molecules)
        else:
            molecules = []
            for molecule in self.molecules:
                molecules.append(self.run_molecule_center(molecule))
        return molecules

        
