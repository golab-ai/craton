import os
import subprocess
from pathlib import Path

CRATON_DIR = Path(__file__).parent
#ROOT_DIR = CRATON_DIR.parent
ROOT_DIR = CRATON_DIR

#protein_file = f"{ROOT_DIR}/data/fep/r_group/tyk2_protein.pdb"
ligand_file = f"{ROOT_DIR}/data/liquid_test.csv"
sdf_file = f"{ROOT_DIR}/data/ligand1prepped.sdf"

Path("./runjob").mkdir(exist_ok=True)
Path("./runjob/ff").mkdir(exist_ok=True)

class TestFF:
    def test_atom_type(self):
        re = subprocess.run(f"craton ff atom_type -i {sdf_file} -o ./runjob/ff/atom_type",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton ff atom_type -i {ligand_file} -o ./runjob/ff/atom_type",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_assign_ff(self):
        re = subprocess.run(f"craton ff assign_ff -i {sdf_file} -o ./runjob/ff/assign_ff",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton ff assign_ff -i {ligand_file} -o ./runjob/ff/assign_ff",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        
