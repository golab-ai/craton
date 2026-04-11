import os
import subprocess
from pathlib import Path

CRATON_DIR = Path(__file__).parent
#ROOT_DIR = CRATON_DIR.parent
ROOT_DIR = CRATON_DIR

protein_file = f"{ROOT_DIR}/data/tyk2_protein.pdb"
ligand_file = f"{ROOT_DIR}/data/tyk2_ligands.sdf"
hfe_file = f"{ROOT_DIR}/data/hfe_test.csv"
liquid_file = f"{ROOT_DIR}/data/liquid_test.csv"
logp_file = f"{ROOT_DIR}/data/logp_test.csv"

Path("./runjob").mkdir(exist_ok=True)
Path("./runjob/simulation").mkdir(exist_ok=True)
class TestSimulation:

    def test_rbfe(self):
        re = subprocess.run(f"craton simulation rbfe --protein {protein_file} --ligands {ligand_file} -o ./runjob/simulation/rbfe",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_abfe(self):
        re = subprocess.run(f"craton simulation abfe --protein {protein_file} --ligands {ligand_file} -o ./runjob/simulation/abfe",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_hfe(self):
        re = subprocess.run(f"craton simulation ahfe --ligands {hfe_file} -o ./runjob/simulation/hfe",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_complex(self):
        re = subprocess.run(f"craton simulation complex --protein {protein_file} --ligands {ligand_file} -o ./runjob/simulation/complex",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_liquid(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -o ./runjob/simulation/liquid",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_logp(self):
        re = subprocess.run(f"craton simulation alogp --ligands {logp_file} --solvent water,octanol -o ./runjob/simulation/logp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_protein(self):
        re = subprocess.run(f"craton simulation protein --molecules {protein_file} -o ./runjob/simulation/protein",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_solution(self):
        re = subprocess.run(f"craton simulation solution --molecules {liquid_file} -o ./runjob/simulation/solution",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_vacuum(self):
        re = subprocess.run(f"craton simulation vacuum --molecules {liquid_file} -o ./runjob/simulation/vacuum",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_ap(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p ap -o ./runjob/property_ap -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_cp(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p cp -o ./runjob/simulation/property_cp -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_dc(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p dc -o ./runjob/simulation/property_dc -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_density(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p den -o ./runjob/simulation/property_den -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_er(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p er -o ./runjob/simulation/property_er -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_hov(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p hov -o ./runjob/simulation/property_hov -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_kt(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p kt -o ./runjob/simulation/property_kt -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_st(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p st -o ./runjob/simulation/property_st -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_td(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p td -o ./runjob/simulation/property_td -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_vis(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p vis -o ./runjob/simulation/property_vis -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_propety_all(self):
        re = subprocess.run(f"craton simulation liquid --molecules {liquid_file} -p den:hov:cp:kt:ap:vis:td:st -o ./runjob/simulation/property_all -eng lmp",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        
if __name__ == "__main__":
    print(protein_file,logp_file)