import multiprocessing
from collections import defaultdict
from pathlib import Path

import gromacs
import MDAnalysis
import numpy as np
import pandas as pd
#from gromacs.cbook import edit_mdp
#from gromacs.fileformats.xvg import XVG
#from gromacs.utilities import in_dir
from MDAnalysis.selections.gromacs import SelectionWriter

from ...utils import logger

energy_item = [
    "Bond",
    "Angle",
    "Proper-Dih.",
    "Improper-Dih.",
    "LJ-14",
    "Coulomb-14",
    "LJ-(SR)",
    "Disper.-corr.",
    "Coulomb-(SR)",
    "Coul.-recip.",
    "Potential",
]

CPU_NUMBER = multiprocessing.cpu_count()

if CPU_NUMBER > 64:
    CPU_NUMBER = 64


def add_energygrps(mdp_file):
    found_energygrps = False
    with open(mdp_file, "r") as f:
        for line in f:
            if line.strip().startswith("energygrps"):
                found_energygrps = True
                break
    if not found_energygrps:
        with open(mdp_file, "a") as f:
            f.write("energygrps=")


def modify_itp_use_topology_b(raw_itp, raw_top):
    new_itp = raw_itp.parent / "top_b.itp"
    new_top = raw_top.parent / "top_b.top"
    with open(new_itp, "w") as f:
        with open(raw_itp, "r") as fin:
            start_parse = False
            for line in fin:
                if start_parse:
                    if line == "[ bonds ]":
                        start_parse = False
                    token = line.split()
                    if len(token) == 11:
                        token[1] = token[8]
                        token[6] = token[9]
                        token[7] = token[10]
                        f.write("      ".join(token) + "\n")
                    else:
                        f.write(line)
                elif "[ atoms ]" in line:
                    start_parse = True
                    f.write(line)
                else:
                    f.write(line)

    raw_itp_name = raw_itp.name
    with open(new_top, "w") as f_w:
        with open(raw_top, "r") as f_r:
            for line in f_r:
                if raw_itp_name in line:
                    line = '#include "top_b.itp"\n'
                f_w.write(line)
    return new_top


def _extract_residue_from_interaction_df(df):
    selection_str = ""
    for res in df.index:
        resindex, _ = res.split("_")
        selection_str += f"resid {resindex} or "
    return selection_str[:-3]


def _get_statistical_data_from_energy_df(df):
    columns = df.columns
    residue_energy_stat = defaultdict(dict)
    for i in range(0, len(columns), 2):
        interaction_type1, interaction_residue = columns[i].split(":")
        interaction_type2, interaction_residue = columns[i + 1].split(":")
        bs_residue = interaction_residue.split("-")[1]
        total_energy = df[columns[i]] + df[columns[i + 1]]
        residue_energy_stat[bs_residue]["total_mean"] = total_energy.mean()
        residue_energy_stat[bs_residue]["total_std"] = total_energy.std()
        residue_energy_stat[bs_residue][f"{interaction_type1}_mean"] = df[columns[i]].mean()
        residue_energy_stat[bs_residue][f"{interaction_type1}_std"] = df[columns[i]].std()
        residue_energy_stat[bs_residue][f"{interaction_type2}_mean"] = df[columns[i + 1]].mean()
        residue_energy_stat[bs_residue][f"{interaction_type2}_std"] = df[columns[i + 1]].std()

    df = pd.DataFrame.from_dict(residue_energy_stat).T
    df.sort_index(key=lambda col: np.array([int(item.split("_")[0]) for item in col]), inplace=True, ascending=False)
    return df


class EnergyFromGmxData:
    def __init__(
        self,
        md_path=".",
        gro_file="prod_npt.gro",
        xtc_file="prod_npt.xtc",
        mdp_file="_prod_npt.mdp",
        top_file="topol.top",
        logs="logs",
        output_file="lie.csv",
        md_type="fep",
    ):
        md_path = Path(md_path).resolve()
        self.gro = md_path / gro_file
        self.xtc = md_path / xtc_file
        self.mdp = md_path / mdp_file
        self.top = md_path / top_file
        if md_type == "fep":
            self.itp = next(md_path.glob("*_to_*.itp"))
        if not self.top.exists():  # for two stages
            self.top = md_path.parent / top_file
        self.logs = Path(logs)
        self.output_file = output_file
        if not self.gro.exists():
            raise RuntimeError(f"cannot find gro file: {gro_file}")
        if not self.xtc.exists():
            raise RuntimeError(f"cannot find xtc file: {xtc_file}")
        if not self.mdp.exists():
            raise RuntimeError(f"cannot find mdp file: {mdp_file}")
        if not self.top.exists():
            raise RuntimeError(f"cannot find top file: {mdp_file}")
        if md_type == "fep":
            if not self.itp.exists():
                raise RuntimeError(f"cannot find itp file")

    def pairwise_energy(self, groups1, *other_group):
        if groups1 in ["LIG", "Protein"]:
            group1_name = groups1
        else:
            group1_name = "_".join(groups1.strip().split()) + "_"
        group2_names = []
        extend_energy_item = []
        for group_2 in other_group:
            if group_2 in ["LIG", "Protein"]:
                group2_name = group_2
            else:
                group2_name = "_".join(group_2.strip().split()) + "_"
            group2_names.append(group2_name)
            extend_energy_item.append(f"Coul-SR:{group1_name}-{group2_name}")
            extend_energy_item.append(f"LJ-SR:{group1_name}-{group2_name}")
        add_energygrps(self.mdp)
        with in_dir(self.logs):
            edit_mdp(self.mdp, new_mdp="new.mdp", energygrps=[group1_name] + group2_names)
            gromacs.select(s=self.gro, on="sele.ndx", input=["System", groups1, *other_group])
            gromacs.grompp(f="new.mdp", p=self.top, c=self.gro, maxwarn="3", o="rerun", n="sele.ndx")
            gromacs.mdrun(rerun=self.xtc, s="rerun.tpr", e="rerun.edr")
            gromacs.energy(f="rerun.edr", o="rerun.xvg", input=energy_item + extend_energy_item)
        return XVG(self.logs / "rerun.xvg").to_df().to_csv(self.output_file, float_format="%.3f")

    def binding_site_energy(self, cutoff=4, interaction_df=None, use_top_b=False):
        u = MDAnalysis.Universe(self.gro, str(self.xtc))
        logger.debug(interaction_df)
        if interaction_df is not None:
            selection_string = _extract_residue_from_interaction_df(interaction_df)
            binding_atoms = u.select_atoms(selection_string)
        else:
            binding_atoms = u.select_atoms(f"(around {cutoff} resname LIG) and protein")
        if use_top_b:
            self.top = modify_itp_use_topology_b(self.itp, self.top)
        group_names, extend_energy_items = [], []
        if not binding_atoms:
            raise RuntimeError("selection error ")
        with in_dir(self.logs):
            sw = SelectionWriter("sele.ndx")
            sw.write(u, name="system")
            sw.write(u.select_atoms("resname LIG"), name="LIG")
            group_names.append("LIG")
            for residue in {atom.residue for atom in binding_atoms}:
                group_name = f"{residue.resid}_{residue.resname}"
                sw.write(residue, name=group_name)
                group_names.append(group_name)
                extend_energy_items.append(f"Coul-SR:LIG-{group_name}")
                extend_energy_items.append(f"LJ-SR:LIG-{group_name}")
            sw.close()
            add_energygrps(self.mdp)
            edit_mdp(self.mdp, new_mdp="new.mdp", energygrps=group_names, free_energy="no")
            gromacs.grompp(f="new.mdp", p=self.top, c=self.gro, maxwarn="3", o="rerun", n="sele.ndx")
            gromacs.mdrun(rerun=self.xtc, s="rerun.tpr", e="rerun.edr")
            gromacs.energy(f="rerun.edr", o="rerun.xvg", input=energy_item + extend_energy_items)
        return _get_statistical_data_from_energy_df(XVG(self.logs / "rerun.xvg").to_df()[extend_energy_items])


class EnergyDifferenceFromFEP:
    def __init__(
        self,
        fep_path=".",
        xtc_file="prod_npt.xtc",
        tpr_file="prod_npt.tpr",
        logs="logs",
        mdp_file="_prod_npt.mdp",
        top_file="topol.top",
    ):
        self.fep_path = Path(fep_path).resolve()
        self.xtc_file = xtc_file
        self.tpr_file = tpr_file
        self.mdp_file = mdp_file
        self.logs = Path("logs")
        self.top_file = top_file

    def energy_diff(self, left, right, extend_energy_item=None, tpr_file=None, output="energy_diff.csv"):
        if extend_energy_item is None:
            extend_energy_item = []
        if tpr_file is None:
            tpr_file = self.tpr_file
        xtc_file = self.fep_path / left / self.xtc_file
        df_list = []
        with in_dir(self.logs):
            for item in [left, right]:
                gromacs.mdrun(rerun=xtc_file, s=self.fep_path / item / tpr_file, e=f"{left}{item}.edr")
                gromacs.energy(f=f"{left}{item}.edr", o=f"{left}{item}.xvg", input=energy_item + extend_energy_item)
                df_list.append(XVG(f"{left}{item}.xvg").to_df())
            (df_list[1] - df_list[0]).to_csv(output, float_format="%.3f")

    def energy_diff_group(self, left, right, groups1, groups2, output="energy_diff.csv"):
        if groups1 in ["LIG", "Protein"]:
            groups1_name = groups1
        else:
            groups1_name = "_".join(groups1.strip().split()) + "_"
        if groups2 in ["LIG", "Protein"]:
            groups2_name = groups2
        else:
            groups2_name = "_".join(groups2.strip().split()) + "_"
        extend_energy_item = [
            f"Coul-SR:{groups1_name}-{groups2_name}",
            f"Coul-SR:{groups1_name}-{groups1_name}",
            f"Coul-SR:{groups1_name}-rest",
            f"Coul-SR:{groups2_name}-rest",
            f"LJ-SR:{groups1_name}-{groups1_name}",
            f"LJ-SR:{groups1_name}-{groups2_name}",
            f"LJ-SR:{groups1_name}-rest",
            f"LJ-SR:{groups2_name}-rest",
            f"LJ-14:{groups1_name}-{groups1_name}",
            f"LJ-14:{groups1_name}-{groups2_name}",
            f"Coul-14:{groups1_name}-{groups1_name}",
            f"Coul-14:{groups1_name}-{groups2_name}",
        ]

        for item in [left, right]:
            tpr_file = self.fep_path / item / self.tpr_file
            mdp_file = self.fep_path / item / self.mdp_file
            add_energygrps(mdp_file)
            with in_dir(self.fep_path / item):
                edit_mdp(mdp_file, new_mdp="new.mdp", energygrps=[groups1_name, groups2_name])
                gromacs.select(s=tpr_file, on="sele.ndx", input=["System", groups1, groups2])
                gromacs.grompp(f="new.mdp", p=self.top_file, c=self.tpr_file, maxwarn="3", o="rerun", n="sele.ndx")
        self.energy_diff(left, right, extend_energy_item=extend_energy_item, tpr_file="rerun.tpr", output=output)
