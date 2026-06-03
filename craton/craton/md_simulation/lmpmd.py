import shutil
import os
import numpy as np
from copy import deepcopy
from itertools import chain
from pathlib import Path

from ..software.lammps import LmpInputFile

#WATER_DIR = share_template.water
#ION_DIR = share_template.ion
MINIMUM_FRAME_NUMBER = 50
MINIMUM_DHDL_NUMBER = 50

class LmpEngine:
    def __init__(self,system):
        self.system = system
        lmp = LmpInputFile(self.system)
        lmp.write_data()
        
        lmp.write_in()
   
class LmpSlurm:
    def __init__(self,sm):
        self._write_slurm_file(sm)

    def _get_slure_header_string(slef,sm):
        hpc = sm.env_setting["hpc_resource"]
        try:
            hpc_loads = sm.env_setting[hpc]
        except:
            hpc_loads = []
        
        #tasks = [int(dd) for dd in os.listdir(sm.output_dir) if os.path.isdir(f"{sm.output_dir}/{dd}")]
        #tasks = sorted(tasks)
        #if len(tasks) == 0:
        #    ntasks = 1
        #else:
        #    ntasks = len(tasks)
        
        # ntasks = 1
        cpu_per_task = 1
        ncpu = sm.env_setting["ncpu"]
        ngpu = sm.env_setting["ngpu"]
        # cpu_per_task = ncpu // ntasks
        ntasks = ncpu

        # text = "#!/bin/bash\n"
        # text += f"#SBATCH --job-name={sm.name}\n"
        # text += f"#SBATCH --output=_job.out\n"
        # text += f"#SBATCH --error=_job.err\n"
        # text += f"#SBATCH --partition={sm.env_setting['partition']}\n"
        # text += f"#SBATCH --nodes={sm.env_setting['nodes']}\n"
        # text += f"#SBATCH --ntasks={ntasks}\n"
        # text += f"#SBATCH --cpus-per-task={cpu_per_task}\n"
        # text += f"#SBATCH --gres=gpu:{ngpu}\n"
        # text += "\n\n\n"

        # for ss in hpc_loads:
        #     text += f"{ss}\n"
        # text += "\n\n"

        text = """#!/bin/bash

module load CentOS/7.9/gcc/14.0.0
module load CentOS/7.9/LAMMPS/20210630
source /cpfs01/xingyun/hpc/VASP/Basekit/setvars.sh --force

"""

        return text, ntasks, cpu_per_task ###,tasks

    def _write_slurm_file(self,sm):
        text,ntasks,ntomp = self._get_slure_header_string(sm)
        hpc = sm.env_setting["hpc_resource"]
        if hpc == "CFFF":
            text += f"srun --mpi=pmi2 -n {ntasks} -c 1 lmp_mpi -in lmp.in"
        else:
            text += "mpirun -np {ntasks} lmp_mpi < lmp.in > lmp.out"
        text += "\n\n"
        with open(f"{sm.output_dir}/job.sh","w") as outf:
            outf.write(text)


