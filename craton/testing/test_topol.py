import os
import subprocess
from pathlib import Path

CRATON_DIR = Path(__file__).parent
#ROOT_DIR = CRATON_DIR.parent
ROOT_DIR = CRATON_DIR

#protein_file = f"{ROOT_DIR}/data/fep/r_group/tyk2_protein.pdb"
ligand_file = f"{ROOT_DIR}/data/liquid_test.csv"
sdf_file = f"{ROOT_DIR}/data/ligand1prepped.sdf"
smiles = "CCc1cc2c(SCC(=O)N3CCNC3=O)ncnc2s1"

Path("./runjob").mkdir(exist_ok=True)
Path("./runjob/topol").mkdir(exist_ok=True)

class TestMM:
    def test_angle_belding(self):
        re = subprocess.run(f"craton stru vary -i '{smiles}' -a 5-6-7 -v 130.0 -o ./runjob/topol/angle_bending",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_angle_degree(self):
        re = subprocess.run(f"craton stru measure -i '{smiles}' -a 5-6-7 -o ./runjob/topol/angle_degree",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_bond_distance(self):
        re = subprocess.run(f"craton stru measure -i '{smiles}' -a 5-6 -o ./runjob/topol/bond_distance",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_bond_stretching(self):
        re = subprocess.run(f"craton stru vary -i '{smiles}' -a 5-6 -v 2.0 -o ./runjob/topol/bond_stretching",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_dihedral_degree(self):
        re = subprocess.run(f"craton stru measure -i '{smiles}' -a 4-5-6-7 -o ./runjob/topol/dihedral_degree",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_bond_stretching(self):
        re = subprocess.run(f"craton stru vary -i '{smiles}' -a 4-5-6-7 -v 60.0 -o ./runjob/topol/bond_stretching",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_vary(self):
        re = subprocess.run(f"craton stru vary -i '{smiles}' -a 4-5-6-7 -v 60.0 -o ./runjob/topol/vary",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_measure(self):
        re = subprocess.run(f"craton stru measure -i '{smiles}' -a 4-5-6-7 -o ./runjob/topol/measure",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_chiral(self):
        re = subprocess.run(f"craton stru topol chiral -i {ligand_file} -o ./runjob/topol/chiral",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_croase_grain(self):
        re = subprocess.run(f"craton stru topol cg -i {ligand_file} -o ./runjob/topol/cg",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_fragmentation(self):
        re = subprocess.run(f"craton stru topol frag -i {ligand_file} -o ./runjob/topol/frag",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_function_group(self):
        re = subprocess.run(f"craton stru topol fg -i {ligand_file} -o ./runjob/topol/fg",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_hybrid(self):
        re = subprocess.run(f"craton stru topol hybrid -i {ligand_file} -o ./runjob/topol/hybride",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_interaction_site(self):
        re = subprocess.run(f"craton stru topol interaction_site -i {ligand_file} -o ./runjob/topol/interaction_site",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_ring(self):
        re = subprocess.run(f"craton stru topol ring -i {ligand_file} -o ./runjob/topol/ring",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_torsion(self):
        re = subprocess.run(f"craton stru topol torsion -i {ligand_file} -o ./runjob/topol/torsion",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    