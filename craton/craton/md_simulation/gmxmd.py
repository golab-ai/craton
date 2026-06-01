import shutil
import os
import numpy as np
from copy import deepcopy
from itertools import chain
from pathlib import Path

from ..chem.molecule import Molecule
from ..software.gromacs import GroInputFile
from ..software.utils import set_gromacs_atom_info
import getpass


#WATER_DIR = share_template.water
#ION_DIR = share_template.ion
MINIMUM_FRAME_NUMBER = 50
MINIMUM_DHDL_NUMBER = 50

class GmxEngine:
    def __init__(self,system):
        self.system = system
        for molecule in self.system.molecules:
            if molecule.style != "protein":
                set_gromacs_atom_info(molecule,molecule.gmx_residue_name)

        if self.system.style in ["rbfe","rhfe","abfe","ahfe","hfe","mutation","rlogs","rlogp","alogp","mem-rbfe","cov-rbfe","pep-rbfe"]:
            eng  = GmxFEPEachWindow(self.system)
            eng.write_para()
        else:
            self._write_gromacs_prepare_file()
    
    def _load_experiment_data(self, sm):
        if not hasattr(sm, "inchi_key"):
            return
        exp_data = read_md_experiment_file()
        if mol_exp_data := exp_data.get(sm.mole[0].inchi_key.split("-")[0]):
            self.gromacs_md_parameters["temperature"]["temperature"] = str(mol_exp_data.temperature)
            self.gromacs_md_parameters["velocity"]["velocity_generate_temp"] = str(mol_exp_data.temperature)
            self.gromacs_md_parameters["pressure"]["pressure"] = str(mol_exp_data.pressure)

    def _write_gromacs_prepare_file(self):
        dir_output = self.system.output_dir
        md_parameters = self.system.md_setting
        GroInputFile.write_mdp(md_parameters, dir_output)
        grofile = GroInputFile("normal")
        grofile.import_systemobj(self.system)
        grofile.write_input_file(dir_output, write_gro=True)
        if "hov" in md_parameters["property"]:
            grofile = GroInputFile("normal")
            grofile.import_systemobj(self.system, exclusions=True)
            grofile.write_mole_itp(self.path, exclusions=True)
            grofile.write_top(self.path, exclusions=True)
            shutil.copy(f"{MD_SETTING_DIR}/hov_output.in", dir_output)
        if "density" in md_parameters["property"]:
            shutil.copy(f"{MD_SETTING_DIR}/energy_output.in", dir_output)
        
class GmxSlurm:
    def __init__(self,sm):
        self._write_slurm_file(sm)
        
        #self.usr = getpass.getuser()

    def _get_slure_header_string(slef,sm):
        hpc = sm.env_setting["hpc_resource"]
        try:
            hpc_loads = sm.env_setting[hpc]
        except:
            hpc_loads = []
        tasks = [int(dd) for dd in os.listdir(sm.output_dir) if os.path.isdir(f"{sm.output_dir}/{dd}") and dd != "job_info"]
        tasks = sorted(tasks)
        if len(tasks) == 0:
            ntasks = 1
        else:
            ntasks = len(tasks)

        ncpu = sm.env_setting["ncpu"]
        ngpu = sm.env_setting["ngpu"]
        cpu_per_task = ncpu // ntasks
        if ntasks == 1:
            if cpu_per_task > 64:
                cpu_per_task = 64

        text = "#!/bin/bash\n"
        text += f"#SBATCH --job-name={sm.name}\n"
        text += f"#SBATCH --output=_job.out\n"
        text += f"#SBATCH --error=_job.err\n"
        text += f"#SBATCH --partition={sm.env_setting['partition']}\n"
        text += f"#SBATCH --nodes={sm.env_setting['nodes']}\n"
        text += f"#SBATCH --ntasks={ntasks}\n"
        text += f"#SBATCH --cpus-per-task={cpu_per_task}\n"
        if ngpu != 0:
            text += f"#SBATCH --gres=gpu:{ngpu}\n"
        text += "\n\n\n"

        # shiyunfei modified
        # for ss in hpc_loads:
        #    text += f"{ss}\n"
        if hpc == "CFFF":
            # text += "source /cpfs01/xingyun/hpc/software/CentOS/7.9/gromacs+fep+mpi/gromacs.sh\n"        
            text += "source /cpfs01/projects-HDD/cfff-4405968bce88_HDD/public/craton.sh\n"
        text += "env\n"
        text += "\n\n"

        return text, ntasks, cpu_per_task,tasks

    def _write_slurm_file(self,sm):
        
        hpc = sm.env_setting["hpc_resource"]
        ngpu = sm.env_setting["ngpu"]
        text,ntasks,ntomp,tasks = self._get_slure_header_string(sm)
        ntasks_string = " ".join([str(n) for n in tasks])
        if ntasks > 1:
            has_conf = os.path.isfile(f"{sm.output_dir}/conf.gro") 
            has_top = os.path.isfile(f"{sm.output_dir}/topol.top") 
            if has_conf:
                jobs = ["../conf"]
            else:
                jobs = ["/conf"]
            jobs.extend(sm.md_setting["jobs"])
            
            if has_top:
                top_file = "../topol.top"
            else:
                top_file = "topol.top"
        else:
            jobs = ["conf"]
            top_file = "topol.top"
            jobs.extend(sm.md_setting["jobs"])
        fep_hrex = False
        if "free_energy_auixed" in sm.md_setting:
            if "fep_hrex" in sm.md_setting["free_energy_auixed"]:
                fep_hrex = sm.md_setting["free_energy_auixed"]["fep_hrex"]

        for ii,job in enumerate(jobs[1:]):
            if ntasks > 1:
                text += f"for i in {ntasks_string}; do cd $i; "
            text += f"gmx_mpi -quiet grompp -f _{job}.mdp -p {top_file} -c {jobs[ii]}.gro -maxwarn 3 -o {job}"
            if ntasks > 1:
                text += f"; cd ..; done"
            text += "\n"
            if hpc == "CFFF":
                if ntasks == 1:
                     text += f"srun --mpi=pmi2 -n {ntomp} -c 1 gmx_mpi mdrun -pin on -quiet -deffnm {job} -ntomp 1"  
                     if job != "mini":
                         text += " -dlb yes"
                else:
                     text += f"srun --mpi=pmi2 -n {ntasks} -c {ntomp} gmx_mpi mdrun -pin on -quiet -deffnm {job} -ntomp {ntomp}"
            else:  
                text += f"mpirun -np {ntasks} gmx_mpi mdrun -quiet -deffnm {job} -ntomp {ntomp}"
                if ngpu == 0:
                    text += f" -nb cpu -pme cpu -pmefft cpu -bonded cpu -update cpu"
                elif ngpu == 1:
                    text += f" -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu"
                else:
                    text += f" -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu"
            if ntasks > 1:
                text += f" -multidir {ntasks_string}"
                if ii >= 2 and fep_hrex:
                    text += f" -replex 500 -fephrex -dlb no -notunepme "
            text += "\n"
        text += "\n\n"
        with open(f"{sm.output_dir}/job.sh","w") as f:
            f.write(text)


        if hpc == "CFFF":
            self.write_cfff_batch(sm.output_dir,ntasks,ntomp)
        # elif hpc == "local":
        #     self.write_local_batch(sm.output_dir,ntasks,ntomp)

    def write_cfff_batch(self,output_dir,ntasks,ntomp):
        usr = getpass.getuser()
        if not os.path.isabs(output_dir):
            output_dir = os.path.abspath(output_dir)
        try:
            re_path = output_dir.split(f"{usr}/")[1]
        except:
            re_path = output_dir.split("public/")[1]

        batch_text = f"sbatch -N 1 -n 1 -c 128 -D {re_path} --exclusive --mem=300GB {output_dir}/job.sh \n"
        # batch_text = f"sbatch -N 1 -n {ntasks} -c {ntomp} -D {re_path} --exclusive --mem=300GB {output_dir}/job.sh \n"
        if os.path.exists("./batch_0.txt"):
            outf = open("./batch_0.txt",'a')
        else:
            outf = open("./batch_0.txt",'w')
        outf.write(batch_text)
        outf.close()


    def write_local_batch(self,output_dir,ntasks,ntomp):
        if not os.path.isabs(output_dir):
            output_dir = os.path.abspath(output_dir)

        batch_text = f"{output_dir}/job.sh\n"
        if os.path.exists("./batch_0.txt"):
            outf = open("./batch_0.txt",'a')
        else:
            outf = open("./batch_0.txt",'w')
        outf.write(batch_text)
        outf.close()

    def _old_write_multitask_slurm_file(self, sm, df_lambda, mixed_lambda):
        try:
            df_lambda = sm.md_para["free_energy_auixed"]["lambdas"]
        except:
            df_lambda = []
        try:
            is_relative = sm.md_para["free_energy_auixed"]["is_relative"]
        except:
            is_relative = False
        n_window = len(df_lambda)
        if not is_relative:
            self._run_write_multitask_slurm_file(np.arange(n_window), sm.name, "job.sh", fep_hrex=False)
        elif not mixed_lambda:
            split_idx = df_lambda.perturbB_vdw.argmin()
            left_window_number = np.arange(n_window)[:split_idx]
            right_window_number = np.arange(n_window)[split_idx:]
            self._run_write_multitask_slurm_file(left_window_number, sm.name + "1", "job1.sh", fep_hrex=False)
            self._run_write_multitask_slurm_file(right_window_number, sm.name + "2", "job2.sh", fep_hrex=False)
        else:
            self._run_write_multitask_slurm_file(np.arange(n_window), sm.name, "job.sh", fep_hrex=True)


class FEPMolPara:
    @staticmethod
    def vdw_para(atom, stage1, stage2, onlyA, onlyB):
        if stage1:  # perturbB_vdw:
            if onlyB:
                atom.atom_type_name = "_D"
                atom.atom_type_name_m2 = atom.atom_type_name_m2
            elif onlyA:
                atom.atom_type_name = atom.atom_type_name + ":_D"
                atom.atom_type_name_m2 = atom.atom_type_name
            else:
                atom.atom_type_name = atom.atom_type_name
                atom.atom_type_name_m2 = atom.atom_type_name_m2
        elif stage2:
            if onlyB:
                atom.parameter, atom.parameter_m2 = atom.parameter_m2, atom.parameter
                atom.mass, atom.mass_m2 = atom.mass_m2, atom.mass
                atom.atom_type_name = atom.atom_type_name_m2
                atom.atom_type_name_m2 = atom.atom_type_name_m2
            elif onlyA:
                atom.atom_type_name = atom.atom_type_name + ":_D"
                atom.atom_type_name_m2 = "_D"
            else:
                atom.parameter, atom.parameter_m2 = atom.parameter_m2, atom.parameter
                atom.mass, atom.mass_m2 = atom.mass_m2, atom.mass
                atom.atom_type_name = atom.atom_type_name_m2
                atom.atom_type_name_m2 = atom.atom_type_name_m2

    @staticmethod
    def coul_para(atom, stage1, stage2, onlyA, onlyB, coul, in_place_para=True):
        if not in_place_para:
            if stage2:
                if onlyB:  # B has grown
                    atom.ff_charge = atom.ff_charge_m2
                elif onlyA:  # A need to disappear
                    atom.ff_charge, atom.ff_charge_m2 = 0, 0
                else:
                    atom.ff_charge = atom.ff_charge_m2
        else:
            if stage1:
                if onlyB:
                    atom.ff_charge = 0
                    atom.ff_charge_m2 = 0
                elif onlyA:
                    atom.ff_charge = atom.ff_charge * (1 - coul) + atom.ff_charge_m2 * coul
                    atom.ff_charge_m2 = atom.ff_charge
                else:
                    # atom.ff_charge = atom.ff_charge
                    atom.ff_charge_m2 = atom.ff_charge
            elif stage2:
                if onlyB:
                    atom.ff_charge = atom.ff_charge * (1 - coul) + atom.ff_charge_m2 * coul
                    atom.ff_charge_m2 = atom.ff_charge
                elif onlyA:
                    atom.ff_charge = 0
                    atom.ff_charge_m2 = 0
                else:
                    atom.ff_charge = atom.ff_charge * (1 - coul) + atom.ff_charge_m2 * coul
                    atom.ff_charge_m2 = atom.ff_charge

    @staticmethod
    def bond_para(molecule, bond, in_place_para=True):
        if not in_place_para:
            return
        for attr in Molecule.attrs_topol:
            for term in getattr(molecule, attr, []):
                if hasattr(term, "parameter_m2"):
                    term.parameter = [(term.parameter[j] * (1 - bond) + term.parameter_m2[j] * bond) for j in range(len(term.parameter))]
                    del term.parameter_m2

    @staticmethod
    def delete_unnecessary_para(atom):
        if atom.atom_type_name == atom.atom_type_name_m2 and atom.ff_charge == atom.ff_charge_m2:
            del atom.ff_charge_m2
            del atom.atom_type_name_m2
            del atom.parameter_m2
            del atom.mass_m2

class GmxFEPEachWindow:
    def __init__(self, system):
        self.system = system
        self.path = system.output_dir
        self.md_parameters = system.md_setting
        self.is_relative = self.md_parameters["free_energy_auixed"]["is_relative"]
        if "lambdas" in self.md_parameters["free_energy_auixed"]:
            self.lambda_df = self.md_parameters["free_energy_auixed"]["lambdas"]
        self.mixed_lambda = self.md_parameters["free_energy_auixed"]["mixed_lambda"]
        self.absolute_intra_flag = self.md_parameters["free_energy_auixed"]["absolute_intra_flag"]
        if hasattr(self,"lambda_df"):
            self.n_window = len(self.lambda_df)
        else:
            self.n_window = len(self.md_parameters["free_energy"]["vdw_lambdas"])
        self.coul_para_inplace = self.bond_para_inplace = self.mixed_lambda
        
        
    def write_para(self):
        if self.is_relative:
            self._write_modified_mole()
            self._write_mdp_file()
        else:
            if self.absolute_intra_flag:
                self.md_parameters["free_energy"]["couple-moltype"] = self.system.molecules[0].name
            else:  # when charge change
                if "couple-lambda0" in self.md_parameters["free_energy"]:
                    self.md_parameters["free_energy"].pop("couple-lambda0")
                if "couple-lambda1" in self.md_parameters["free_energy"]:
                    self.md_parameters["free_energy"].pop("couple-lambda1")
                if "couple-intramol" in self.md_parameters["free_energy"]:
                    self.md_parameters["free_energy"].pop("couple-intramol")
            self._write_gromacs_gro_file(self.system)
            self._write_gromacs_top_file(self.system, self.path)  
            self._write_mdp_file()
            
            #shutil.copy(f"{ION_DIR}/Cl-.itp", self.path)
            #shutil.copy(f"{ION_DIR}/Na+.itp", self.path)
            
    def _write_modified_mole(self):
        sm_copy = deepcopy(self.system)
        self._write_gromacs_gro_file(sm_copy)
        has_replaced = False
        for i_window, row in self.lambda_df.iterrows():
            sm_copy.molecules[0] = deepcopy(self.system.molecules[0])
            # ====== check the ALW exists
            idx = -1
            target_atoms = []
            for i, mol in enumerate(sm_copy.molecules):
                if mol.name == "ALW":
                    idx = i
                    break
            if idx != -1:
                sm_copy.molecules[idx] = deepcopy(self.system.molecules[-2])
                target_atoms = sm_copy.molecules[idx].Atoms
            # finish check
            for atom in chain(sm_copy.molecules[0].Atoms, target_atoms):
                if not hasattr(atom, "atom_type_name_m2"):
                    continue
                only_a = atom.atom_type_name_m2 == "_D"  # need to eliminated
                only_b = atom.atom_type_name == "_D"  # need to form
                FEPMolPara.vdw_para(atom, row.perturbB_vdw, row.perturbA_vdw, only_a, only_b)
                FEPMolPara.coul_para(
                    atom, row.perturbB_vdw, row.perturbA_vdw, only_a, only_b, row.coul, self.coul_para_inplace
                )
                FEPMolPara.delete_unnecessary_para(atom)
            FEPMolPara.bond_para(sm_copy.molecules[0], row.bonded, self.bond_para_inplace)
            dir_output = Path(self.path) / str(i_window)
            dir_output.mkdir(parents=True, exist_ok=True)
            ######FEPMolPara.neural_system(sm_copy, dir_output, has_replaced)
            has_replaced = True
            self._write_gromacs_top_file(sm_copy, dir_output)
        
    def _write_gromacs_gro_file(self, sm):
        grofile = GroInputFile("normal")
        grofile.import_systemobj(sm)
        grofile.write_input_file(self.path, write_top=False)

    def _write_gromacs_top_file(self, sm, dir_output):
        grofile = GroInputFile("normal")
        grofile.import_systemobj(sm)
        if self.is_relative:
            grofile.import_rest2_system_para(sm, int(dir_output.name))
        grofile.write_input_file(dir_output, write_gro=False)
        #shutil.copy(f"{WATER_DIR}/{self.md_parameters.defaults().get('water_model', 'tip3p')}" + ".itp", dir_output)

    def _write_mdp_file(self):
        if self.is_relative:
            if self.coul_para_inplace:
                self.lambda_df = self.lambda_df.drop(columns="coul")
            if self.bond_para_inplace:
                self.lambda_df = self.lambda_df.drop(columns="bonded")

        if hasattr(self,"lambda_df"):
            self.md_parameters["free_energy"]["vdw_lambdas"] = " ".join(map(str, self.lambda_df.vdw.to_list()))
            if "coul" in self.lambda_df:
                self.md_parameters["free_energy"]["coul_lambdas"] = " ".join(map(str, self.lambda_df.coul.to_list()))
            if "bonded" in self.lambda_df:
                self.md_parameters["free_energy"]["bonded_lambdas"] = " ".join(map(str, self.lambda_df.bonded.to_list()))
        else:
            for attr in ["vdw_lambdas","coul_lambdas","bonded_lambdas"]:
                if attr in self.md_parameters["free_energy"]:
                    self.md_parameters["free_energy"][attr] = " ".join([str(a) for a in self.md_parameters["free_energy"][attr]])


        def write_mdp_for_window(i_window):
            window_dir = Path(self.path, f"{i_window}")
            window_dir.mkdir(exist_ok=True)
            self.md_parameters["free_energy"]["init_lambda_state"] = str(i_window)
            GroInputFile.write_mdp(self.md_parameters, window_dir, free_energy=True)

        for i_window in range(1, self.n_window - 1):
            write_mdp_for_window(i_window)

        # extra setting for initial and end window
        nstxout_compressed = self.md_parameters["output"]["nstxout-compressed"]
        prod_steps = self.md_parameters["nsteps"][-1]
        nstdhdl = self.md_parameters["free_energy"]["nstdhdl"]
        if (prod_steps // nstxout_compressed) < MINIMUM_FRAME_NUMBER:
            self.md_parameters["output"]["nstxout-compressed"] = str(prod_steps // MINIMUM_FRAME_NUMBER)
        else:
            self.md_parameters["output"]["nstxout-compressed"] = str(nstxout_compressed // 5)
        if nstdhdl != 0 and (prod_steps // nstdhdl) < MINIMUM_DHDL_NUMBER:
            self.md_parameters["free_energy"]["nstdhdl"] = str(prod_steps // MINIMUM_DHDL_NUMBER)

        for i_window in [0, self.n_window - 1]:
            write_mdp_for_window(i_window)

        ## restore default value
        self.md_parameters["output"]["nstxout-compressed"] = str(nstxout_compressed)
