import os
import shutil
from typing import Iterable
from ..utils import logger
#from ..chem import FormatMolecule as FM
from .. import CRATON_CONFIGURE

import subprocess

tmp_path = CRATON_CONFIGURE["path"]["tmp"]

def get_am1bcc(molecule):
    from ..chem import FormatMolecule as FM
    antechamber = shutil.which("antechamber")
    if not antechamber:
        logger.error(f"Please make sure AmberTools is installed and the path is added to environment variable.")
        raise Exception("Getting AM1BCC charge failed")
    molfilename = molecule.mole_name.split("/")[-1]
    FM._convert([molecule],otype="mol2",ofilename=f"{molfilename}",opath=tmp_path)
    FM._convert([molecule],otype="mol",ofilename=f"{molfilename}",opath=tmp_path)


    am1bcc_mol2_path = f"{tmp_path}/{molfilename}_am1bcc.mol2"

    if not os.path.isfile(am1bcc_mol2_path):
        logger.info(f"Run am1bcc calculation.")
        net_charge = sum(molecule.formal_charge)
        try:
            # run antechamber with mol2 as input
            logger.info(f"Run antechamber with mol2 file as input.")
            p = subprocess.run(
                [
                    antechamber,
                    "-at",
                    "gaff2",  # atom type
                    "-i",
                    f"{tmp_path}/{molfilename}.mol2",  # input file name
                    "-fi",
                    "mol2",  # input file format
                    "-o",
                    am1bcc_mol2_path,  # output file name
                    "-fo",
                    "mol2",  # output file format
                    "-c",
                    "bcc",  # charge method
                    "-s",
                    "2",  # status information
                    "-eq",
                    "2",  # equalizing charge
                    "-nc",
                    str(net_charge),  # molecule net charge
                    "-pf",
                    "y",  # remove intermediate file
                    "-ek",
                    "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, maxcyc=0",
                    # am1 calculation setting
                ]
            )
            assert p.returncode == 0
        except AssertionError:
            try:
                # run antechamber with mol as input
                logger.info(f"Antechamber with mol2 file failed. Run antechamber with mol file as input.")
                p = subprocess.run(
                    [
                        antechamber,
                        "-at",
                        "gaff2",  # atom type
                        "-i",
                        f"{tmp_path}/{molfilename}.mol",  # input file name
                        "-fi",
                        "mdl",  # input file format
                        "-o",
                        am1bcc_mol2_path,  # output file name
                        "-fo",
                        "mol2",  # output file format
                        "-c",
                        "bcc",  # charge method
                        "-s",
                        "2",  # status information
                        "-eq",
                        "2",  # equalizing charge
                        "-nc",
                        str(net_charge),  # molecule net charge
                        "-pf",
                        "y",  # remove intermediate file
                        "-ek",
                        "qm_theory='AM1', grms_tol=0.0005, scfconv=1.d-10, ndiis_attempts=700, maxcyc=0",
                    ]
                )
                assert p.returncode == 0
            except AssertionError:
                logger.warning(f"Am1bcc calculation is failed. Charge is missing for molecule: {molfilename}.")
    # reload molecule object from updated mol2 file with am1bcc charge
    if os.path.isfile(am1bcc_mol2_path):
        with open(am1bcc_mol2_path) as inf:
            lines = inf.readlines()
            arr = lines[lines.index("@<TRIPOS>ATOM\n") + 1 : lines.index("@<TRIPOS>BOND\n")]
        for ii,atom in enumerate(molecule.Atoms):
            atom.ff_charge = float(arr[ii].split()[-1])
    else:
        logger.error("Mol2 with am1bcc information not found")
        raise Exception("Calc am1bcc charge failed")
