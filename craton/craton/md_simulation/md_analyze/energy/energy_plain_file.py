from pathlib import Path

import gromacs
import pandas as pd
from gromacs.cbook import Transformer, edit_mdp, trj_fitandcenter
from gromacs.fileformats.xpm import XPM
from gromacs.fileformats.xvg import XVG
from gromacs.setup import energy_minimize
from gromacs.utilities import in_dir

from ....chemistry.chemsystem import System
######from ...chemistry.share.setting import vacuum_mdp
from ...software.gromacs import GroInputFile
from ...run.md_force_field import BuildMoleculeTopologyForMD
from ...run.md_load_file import MDInputFiles, ParseMDInputStructure
from ...run.md_setting import ForceFieldSetting

energy_item = [
    "Bond",
    "Angle",
    "Proper-Dih.",
    "LJ-14",
    "Coulomb-14",
    "LJ-(SR)",
    "Disper.-corr.",
    "Coulomb-(SR)",
    "Potential",
]
pairwise_energy_name_prefix = ["Coul-SR", "Coul-14", "LJ-SR", "LJ-14"]

pd.set_option("display.max_columns", None)


class EnergyFromRawData:
    def __init__(self, moles, protein, logs="logs"):
        self.moles = moles
        self.protein = protein
        self.logs = Path(logs)
        if self.moles:
            BuildMoleculeTopologyForMD.get_force_field_and_assign_for_molecules(
                self.moles, ForceFieldSetting(charge_method="am1bcc"), output_directory=self.logs
            )

    @classmethod
    def mols_from_file(cls, protein_file=None, ligand_file=None, smiles_file=None, smiles=None, output="logs"):
        with in_dir(output):
            mdinput = ParseMDInputStructure(
                MDInputFiles(protein_file, ligand_file, smiles_file, smiles),
                output_directory=output,
                ignore_residue_missing=True,
            )
            moles, protein = None, None
            if ligand_file:
                moles = mdinput.moles
            if protein_file:
                protein = mdinput.protein
            return cls(moles, protein, logs=output)

    def write_prepare_file(self, mol):
        with in_dir(self.logs):
            name = mol.name
            system = System()
            system.mole_number = [1]
            system.mole = [mol]
            system.name = mol.name
            system.coor = mol.coordinates

            if self.protein:
                system.mole_number.append(1)
                system.mole.append(self.protein)
                system.coor = mol.coordinates + self.protein.coordinates
                system.protein_force_field = "amber99sb"

            system.lattics = [9999, 9999, 9999]
            top_file = f"{mol.name}_topol.top"
            conf_file = f"{mol.name}_conf.gro"

            gro = GroInputFile()
            gro.import_systemobj(system)
            gro.write_small_molecule_system(top_file, conf_file)

    def single_point_energy(self, mol, optimize=True, center=True, mdp=vacuum_mdp):
        kwargs = {}
        if not optimize:
            kwargs["nsteps"] = 0
        with in_dir(self.logs):

            gromacs.grompp(f=mdp, o=f"{mol.name}.tpr", c=f"{mol.name}_conf.gro", p=f"{mol.name}_topol.top", maxwarn="1")
            gromacs.mdrun(deffnm=f"{mol.name}")
            # energy_minimize(
            #     dirname=".",
            #     mdp=mdp,
            #     struct=f'{mol.name}_conf.gro',
            #     top=f'{mol.name}_topol.top',
            #     deffnm=f'{mol.name}',
            #     output=f'{mol.name}.pdb',
            #     maxwarn=1,
            #     emtol=100,
            #     **kwargs)
            gromacs.energy(f=f"{mol.name}.edr", o=f"{mol.name}.xvg", input=energy_item)
            if center:
                trj_fitandcenter(
                    f=f"{mol.name}.gro",
                    s=[f"{mol.name}.tpr", f"{mol.name}_conf.gro"],
                    o=f"{mol.name}_optimize.pdb",
                    input=("system", "system", "system"),
                )
        return XVG(self.logs / f"{mol.name}.xvg").to_df()

    def pairwise_energy(self, mol, groups1, groups2, optimize=False, mdp=vacuum_mdp, **kwargs):
        groups1_name = "_".join(groups1.strip().split()) + "_"
        groups2_name = "_".join(groups2.strip().split()) + "_"
        extend_energy_item = [f"Coul-SR:{groups1_name}-{groups2_name}", f"LJ-SR:{groups1_name}-{groups2_name}"]
        if not optimize:
            kwargs["nsteps"] = 0
        with in_dir(self.logs):
            edit_mdp(mdp, new_mdp="vacuum.mdp", energygrps=[groups1_name, groups2_name])
            gromacs.select(s=f"{mol.name}_conf.gro", on="sele.ndx", input=[groups1, groups2])
            energy_minimize(
                dirname=".",
                mdp="vacuum.mdp",
                struct=f"{mol.name}_conf.gro",
                top=f"{mol.name}_topol.top",
                deffnm=mol.name,
                output=f"{mol.name}.pdb",
                maxwarn=1,
                emtol=100,
                n="sele.ndx",
                **kwargs,
            )
            gromacs.energy(f=f"{mol.name}.edr", o=f"{mol.name}.xvg", input=energy_item + extend_energy_item)
            if optimize:
                Transformer(s=f"{mol.name}.tpr", f=f"{mol.name}.pdb", force=True).center_fit(
                    o=f"{mol.name}_optimize.pdb"
                )
        return XVG(self.logs / f"{mol.name}.xvg").to_df()

    def matrix_energy(self, mol, *args, optimize=False, mdp=vacuum_mdp, **kwargs):
        group_file = "groups.dat"
        group_name = []
        for group in args:
            group_name.append("_".join(group.strip().split()) + "_")
        if not optimize:
            kwargs["nsteps"] = 0
        with in_dir(self.logs):
            edit_mdp(mdp, new_mdp="vacuum.mdp", energygrps=group_name)
            with open(group_file, "w") as f:
                f.write(f"{len(group_name)}\n")
                for group in group_name:
                    f.write(group + "\n")
            gromacs.select(s=f"{mol.name}_conf.gro", on="sele.ndx", input=args)
            energy_minimize(
                dirname=".",
                mdp="vacuum.mdp",
                struct=f"{mol.name}_conf.gro",
                top=f"{mol.name}_topol.top",
                deffnm=mol.name,
                output=f"{mol.name}.pdb",
                maxwarn=1,
                emtol=100,
                n="sele.ndx",
                **kwargs,
            )
            gromacs.enemat(f=f"{mol.name}.edr", groups=group_file, free="no", emat="matrix_energy.xpm")
            # gromacs.energy(f=f'{mol.name}.edr', o=f'{mol.name}.xvg', input=energy_item+extend_energy_item)
            if optimize:
                Transformer(s=f"{mol.name}.tpr", f=f"{mol.name}.pdb", force=True).center_fit(
                    o=f"{mol.name}_optimize.gro"
                )
        return (
            XPM(self.logs / "totalmatrix_energy.xpm").array,
            XPM(self.logs / "LJ-SRmatrix_energy.xpm").array,
            XPM(self.logs / "Coul-SRmatrix_energy.xpm").array,
        )


if __name__ == "__main__":
    energy_calc = EnergyFromRawData.mols_from_file(
        ligand_file="/Users/haomiao/Documents/cadd/fep/dataset/Thrombin_ligands.sdf",
        protein_file="/Users/haomiao/Documents/cadd/fep/dataset/Thrombin_protein.pdb",
    )
    energy_calc.write_prepare_file(energy_calc.moles[0])
    matrix_e = energy_calc.matrix_energy(energy_calc.moles[0], "resname LIG", "resid 100", "resid 101", "resid 102")
    pairwise_e = energy_calc.pairwise_energy(energy_calc.moles[0], "resname LIG", "resid 100")

    print(matrix_e)
    print(pairwise_e)
