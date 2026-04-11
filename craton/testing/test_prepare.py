import os
import subprocess
from pathlib import Path

CRATON_DIR = Path(__file__).parent
#ROOT_DIR = CRATON_DIR.parent
ROOT_DIR = CRATON_DIR

#protein_file = f"{ROOT_DIR}/data/fep/r_group/tyk2_protein.pdb"
ligand_file = f"{ROOT_DIR}/data/ligands.csv"
sdf_file = f"{ROOT_DIR}/data/ligand1prepped.sdf"
smiles = "O(C(OCC)C)CC"
IUPAC_name = "acetal"
CAS_number = "105-57-7"


Path("./runjob").mkdir(exist_ok=True)
Path("./runjob/prepare").mkdir(exist_ok=True)
class TestPrepare:

    def test_uniprot(self):
        re = subprocess.run(f"craton prepare uniprot -t tyk2 -o ./runjob/prepare/TYK2 < uniprot.in ",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_pdb(self):
        re = subprocess.run(f"craton prepare pdb -i ./runjob/prepare/TYK2/P29597_info.txt -o ./runjob/prepare/PDB",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
    
    def test_prepare_protein(self):
        re = subprocess.run(f"craton prepare protein -i ./runjob/prepare/PDB/8S99.pdb -of tyk2_prepare -o ./runjob/prepare/prepare_pdb",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
    
    def test_prepare_ligand(self):
        re = subprocess.run(f"craton prepare ligand -i {ligand_file} -of ligands_prepare.sdf -o ./runjob/prepare/prepare_sdf",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare ligand -i {sdf_file} -of ligands_prepare_2.sdf -o ./runjob/prepare/prepare_sdf",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_molecule_info(self):
        re = subprocess.run(f"craton prepare mol_info -i '{smiles}' -it smiles -o ./runjob/prepare/molecule_info",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare mol_info -i {CAS_number} -it cas -o ./runjob/prepare/molecule_info",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare mol_info -i {IUPAC_name} -it iupac -o ./runjob/prepare/molecule_info",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_mutation(self):
        re = subprocess.run(f"craton prepare mutation -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r ARG_583_A -m LEU -of tyk2_mutation_LEU -o ./runjob/prepare/mutation",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare mutation -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r ARG_583_A -m AIB -of tyk2_mutation_AIB -o ./runjob/prepare/mutation",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_modify(self):
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r SER_593_A -m pho -of tyk2_modify_ser_pho -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r SER_593_A -m n-met -of tyk2_modify_ser_n_met -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r THR_599_A -m pho -of tyk2_modify_thr_pho -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r TYR_604_A -m pho -of tyk2_modify_tyr_pho -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r TYR_604_A -m suf -of tyk2_modify_tyr_suf -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r ARG_607_A -m met -of tyk2_modify_arg_met -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r LYS_642_A -m met -of tyk2_modify_lys_met -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton prepare modify -i ./runjob/prepare/prepare_pdb/tyk2_prepare.pdb -r MET_662_A -m n-met -of tyk2_modify_met_n__met -o ./runjob/prepare/modify",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"