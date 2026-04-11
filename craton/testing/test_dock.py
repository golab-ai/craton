import os
import subprocess
from pathlib import Path

CRATON_DIR = Path(__file__).parent
#ROOT_DIR = CRATON_DIR.parent
ROOT_DIR = CRATON_DIR

protein_file = f"{ROOT_DIR}/data/tyk2_protein.pdb"
ligand_file = f"{ROOT_DIR}/data/tyk2_ligands.sdf"
sdf_file = f"{ROOT_DIR}/data/ligand1prepped.sdf"

Path("./runjob").mkdir(exist_ok=True)
Path("./runjob/dock").mkdir(exist_ok=True)

class TestDock:
    #def test_pocket(self):
        #os.system("export KMP_DUPLICATE_LIB_OK=TRUE")
    #    re = subprocess.run(f"craton dock pocket -i {protein_file} -o ./runjob/dock/pocket",shell=True,capture_output=True,text=True)
    #    if re.returncode != 0:
    #        assert False, f"Simulation failed with error: {re.stderr}"


    def test_assign_dock(self):
        #os.system("export KMP_DUPLICATE_LIB_OK=TRUE")
        re = subprocess.run(f"craton dock dock -p {protein_file} -i {ligand_file} -center 1.508,7.143,-27.885 -box 9.600,9.600,10.200 -o ./runjob/dock/dock",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        
        
