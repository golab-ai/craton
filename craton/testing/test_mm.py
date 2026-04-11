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
Path("./runjob/mm").mkdir(exist_ok=True)

class TestMM:
    def test_center(self):
        re = subprocess.run(f"craton mm center -i {sdf_file} -p center -o ./runjob/mm/center",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm center -i {ligand_file} -p center -o ./runjob/mm/center",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        

    def test_cob(self):
        re = subprocess.run(f"craton mm center -i {sdf_file} -p cob -o ./runjob/mm/cob",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm center -i {ligand_file} -p cob -o ./runjob/mm/cob",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        

    def test_cog(self):
        re = subprocess.run(f"craton mm center -i {sdf_file} -p cog -o ./runjob/mm/cog",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm center -i {ligand_file} -p cog -o ./runjob/mm/cog",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_com(self):
        re = subprocess.run(f"craton mm center -i {sdf_file} -p com -o ./runjob/mm/com",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm center -i {ligand_file} -p com -o ./runjob/mm/com",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_size(self):
        re = subprocess.run(f"craton mm center -i {sdf_file} -p size -o ./runjob/mm/size",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm center -i {ligand_file} -p size -o ./runjob/mm/size",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_multipole(self):
        re = subprocess.run(f"craton mm multipole -p multipole -i {sdf_file} -o ./runjob/mm/multipole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm multipole -p multipole -i {ligand_file} -o ./runjob/mm/multipole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_dipole(self):
        re = subprocess.run(f"craton mm multipole -p dipole -i {sdf_file} -o ./runjob/mm/dipole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm multipole -p dipole -i {ligand_file} -o ./runjob/mm/dipole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_quadrupole(self):
        re = subprocess.run(f"craton mm multipole -p quadrupole -i {sdf_file} -o ./runjob/mm/quadrupole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm multipole -p quadrupole -i {ligand_file} -o ./runjob/mm/quadrupole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
    
    def test_octupole(self):
        re = subprocess.run(f"craton mm multipole -p octupole -i {sdf_file} -o ./runjob/mm/octupole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm multipole -p octupole -i {ligand_file} -o ./runjob/mm/octupole",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_energy(self):
        re = subprocess.run(f"craton mm calculate energy -i {sdf_file} -o ./runjob/mm/energy",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm calculate energy -i {ligand_file} -o ./runjob/mm/energy",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_force(self):
        re = subprocess.run(f"craton mm calculate force -i {sdf_file} -o ./runjob/mm/force",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm calculate force -i {ligand_file} -o ./runjob/mm/force",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_freq(self):
        #re = subprocess.run(f"craton mm calculate freq -i {sdf_file} -o ./runjob/mm/freq",shell=True,capture_output=True,text=True)
        #if re.returncode != 0:
        #    assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm calculate freq -i {ligand_file} -o ./runjob/mm/freq",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_hessian(self):
        re = subprocess.run(f"craton mm calculate hessian -i {sdf_file} -o ./runjob/mm/hessian",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm calculate hessian -i {ligand_file} -o ./runjob/mm/hessian",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_inertia(self):
        re = subprocess.run(f"craton mm inertia -i {sdf_file} -o ./runjob/mm/inertia",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm inertia -i {ligand_file} -o ./runjob/mm/inertia",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_optimized(self):
        re = subprocess.run(f"craton mm opt -i {sdf_file} -o ./runjob/mm/opt",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm opt -i {ligand_file} -o ./runjob/mm/opt",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
    
    def test_scan(self):
        re = subprocess.run(f"craton mm scan -i {sdf_file} -o ./runjob/mm/scan",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm scan -i {ligand_file} -o ./runjob/mm/scan",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_surface(self):
        re = subprocess.run(f"craton mm surface -i {sdf_file} -o ./runjob/mm/surface",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm surface -i {sdf_file} -o ./runjob/mm/surface",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    def test_volume(self):
        re = subprocess.run(f"craton mm volume -i {sdf_file} -o ./runjob/mm/volume",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"
        re = subprocess.run(f"craton mm volume -i {sdf_file} -o ./runjob/mm/volume",shell=True,capture_output=True,text=True)
        if re.returncode != 0:
            assert False, f"Simulation failed with error: {re.stderr}"

    