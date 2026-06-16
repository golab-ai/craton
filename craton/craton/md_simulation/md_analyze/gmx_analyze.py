from collections import defaultdict
from pathlib import Path
import os
import random
from copy import deepcopy
import pandas as pd
import numpy as np
from scipy.constants import R  # 8.314472  J mol^-1 K^-1
import re
import simplejson
from typing import Dict
import shutil
import sys
#import pandas as pd


from ..fep.bar import BarForDDG, BETA
from ...utils.commons import parallel_run
from ...utils import logger
from ...utils.figure import DrawFigure as DF
from ...mm_calculator import Calculator as Calc 
from ...utils.geometry import calc_stru_para
from ...chemkit import InteractionModel as IM #get_interaction_model
from ...chemkit.interaction.report import ReportMultipleInteraction, merge_interaction_dfs
from ...chemkit.interaction.visualize_by_plotly import fep_bar_plot, fep_bar_plot_with_energy,normal_md_bar_plot,normal_md_bar_plot_with_energy
from ..fep.cycle_closure import CycleClosure
from ..mapping.graph import Graph

from .analyze_visualization import (generate_torsion_figure, 
    generate_torsion_figure_4col,
    plot_exchange_rate,
    plot_convergence_block,
    plot_convergence,
    plot_correlation)

#from .pl_interaction.fep_md_traj import FEPTrajInteraction



KCAL_MOL_2_KJ_MOL = 4.184
TEMPERATURE = 310
KB = R / 1000
BETA = 1 / (KB * TEMPERATURE)  # mol / KJ
if shutil.which("gmx_mpi") is not None:
    gmx = "gmx_mpi"
else:
    if shutil.which("gmx") is not None:
        gmx = "gmx"
    else:
        sys.exit("Gromacs未安装")


class GmxAnalyze:
    def __init__(self,temperature=310.0):
        self.temperature = temperature
        self.dump_data = {}
        self.group_idx = None
        self.u_nks = None
        
    def parse_xvg_file(self,file_name):
        with open(file_name) as inf:
            lines = inf.readlines()
        datas = []
        units = []
        labels = []
        for line in lines:
            if line[0] not in ["#","@"]:
                datas.append([float(v) for v in line.strip().split() ])
            elif line.startswith("@    xaxis  label"):
                ss = line.strip().split()
                labels.append(ss[3].strip('"'))
                units.append(ss[4].strip('"'))
            elif line.startswith("@    yaxis  label"):
                ss = line.strip().split()
                for s in ss[3:]:
                    units.append(s.strip('"').strip(","))
            elif line.startswith("@ s"):
                ss = line.strip().split()
                labels.append(ss[3].strip('"'))
        return labels, units, datas

    def parse_xvg_dh_file(self,file_name):
        with open(file_name) as inf:
            lines = inf.readlines()
        datas = []
        temperature = TEMPERATURE
        for line in lines:  
            if line[0] not in ["#","@"]:
                datas.append([float(v) for v in line.strip().split() ])
            elif line.startswith("@ subtitle"):
                ss = line.strip().split()
                temperature = float(ss[4])
                ndx = int(ss[8].strip(":"))
        #temperature = TEMPERATURE
        #BETA  = 1 / (KB * temperature)
        nn = len(datas[0]) - 4
        u_k = {ii:[] for ii in range(nn)}
        for data in datas:
            for ii in range(nn):
                u_k[ii].append( (data[ii+3] + data[-1]) *  BETA )
                #u_k[ii].append(data[ii+3] * BETA)
        u_k = pd.DataFrame(u_k, columns=[ii for ii in range(nn)])
        u_k["lambda"] = f"lambda_{ndx}"
        u_k["window"] = ndx
        return u_k.set_index(["lambda","window"])

    def parse_csv_file(self,csvfile,temperature=310.0):

        BETA = 1 / (KB * temperature)
        df = pd.read_csv(csvfile, na_filter=True, memory_map=True, sep=r"\s+")
        # remove Time, U, pV column, only save the dU column
        u_k = df.loc[:, ~df.columns.isin(["Time", "U", "pV"])]
        # add pV and convert to reduced potential
        u_k = u_k.add(df.pV, axis="rows") * BETA
        # rename window list
        u_k.columns = [f"{i}" for i in range(len(u_k.columns))]
        # add a column to df for using groupby
        u_k["lambda"] = f"lambda_{csvfile.parent.name}"
        u_k["window"] = f"{int(csvfile.parent.name)}"  # for csv
        # set lambda index for later groupby
        return u_k.set_index(["lambda", "window"])

    def parse_gro_file(self,file_name):
        with open(file_name) as inf:
            lines = inf.readlines()
        ans =int(lines[1].strip())
        total_coors = []
        for ii in range(0,len(lines),ans+3):
            coors = []
            for jj in range(2,ans+2):
                line = lines[ii+jj]
                coors.append([float(line[20:28])*10.0,float(line[28:36])*10.0,float(line[36:].strip())*10.0])
            total_coors.append(coors)
        del lines
        return total_coors

    def file_figure(self,output_dir,labels,units,datas,pre_fn=""):
        X = [v[0] for v in datas]
        for ii in range(1,len(labels)):
            Y = [v[ii] for v in datas]
            DF.line_draw([X,Y], name=f"{pre_fn}{labels[ii]}", labels=labels[ii], xylabels=[f"{labels[0]} {units[0]}",f"{labels[ii]} {units[ii]}"], save_path=output_dir)

    def _parse_args(self,prop,args):
        default = {
                    "energy":[["job_dir","."],["output_dir","."],["attrs",["Potential","Temperature","Pressure"]],["pre_fn",""]],
                    "rmsd":[["job_dir","."],["output_dir","."],["attrs",["Protein","LIG"]],["pre_fn",""]],
                    "interaction":[["job_dir","."],["output_dir","."],["molecules",None],["molecule_dir",None],["fn",""], ["save_fig",True],["dimer_type","wat_com"]],
                    "pair_interaction":[["job_dir","."],["output_dir","."],["molecules",None],["molecule_dir",None]],
                    "bar":[["job_dir","."],["output_dir","."],["pre_fn",""]],
                    "ddg":[["job_dir","."],["output_dir","."]],
                    "exchange":[["job_dir","."],["output_dir","."],["pre_fn",""]],
                    "block_ddg":[["job_dir","."],["output_dir","."],["pre_fn",""],["block_num",5]],
                    "accum_ddg":[["job_dir","."],["output_dir","."],["pre_fn",""],["block_num",5]],
                    "rbfe":[["job_dir","."],["output_dir","."],["attrs",None],["molecule_dir",None],["molecules",None]],
                    "all_rbfe":[["job_dir","."],["output_dir","."],["attrs",None],["expt_file",None],["pka_file",None],["parallel",True]],
                    "all_dimer":[["job_dir","."],["output_dir","."],["molecules",None],["molecule_dir",None],["fn",""], ["save_fig",True],["dimer_type","wat_com"],["parallel",True]],
                    
                    }
        values = []
        for rr in default[prop]:
            if rr[0] in args:
                if args[rr[0]] is not None:
                    values.append(args[rr[0]])
                else:
                    values.append(rr[1])
            else:
                values.append(rr[1])
        return values

    def get_index_number(self,file_name):
        with open(file_name) as inf:
            lines = inf.readlines()
        group_idx = {}
        nn = 0
        for line in lines:
            if line.startswith("["):
                #group = line.strip("[")
                #group = group.strip("]")
                #group = group.strip()
                group_idx[line[2:-3]] = nn
                nn += 1
        return group_idx

    def generate_group(self,job_dir,args=None):
        
        with open(f"{job_dir}/group.in",'w') as outf:
            outf.write("0\n\n")
            outf.write("q\n\n\n")
        os.system(f"{gmx} make_ndx -f {job_dir}/prod_npt.tpr -o {job_dir}/index.ndx < {job_dir}/group.in")
        return self.get_index_number(f"{job_dir}/index.ndx")

    def trj_pbc(self,job_dir,group_idx,attrs = None,args=None):
        if attrs is None:
            attrs = ["System"]
        with open(f"{job_dir}/trj_pbc.in",'w') as outf:
            for group in attrs:
                outf.write("%d "%group_idx[group])
            outf.write("\n\n\n")
        os.system(f"{gmx} trjconv -f {job_dir}/prod_npt.xtc -s {job_dir}/prod_npt.tpr -pbc whole -o {job_dir}/prod_npt-pbc.xtc -n {job_dir}/index.ndx < {job_dir}/trj_pbc.in")

    #def get_energy(self,job_dir,output_dir,attrs=None,pre_fn=""):
    def get_energy(self,args):
        job_dir,output_dir,attrs,pre_fn = self._parse_args("energy",args)
        with open(f"{job_dir}/energy.in",'w') as outf:
            for attr in attrs:
                outf.write("%s\n"%attr)
            outf.write("\n\n")
        os.system(f"{gmx} energy -f {job_dir}/prod_npt.edr -o {job_dir}/potention_energy.xvg < {job_dir}/energy.in")
        labels, units,datas = self.parse_xvg_file(f"{job_dir}/potention_energy.xvg")
        self.file_figure(output_dir,labels,units,datas,pre_fn=pre_fn)
        
    def get_rmsd(self,args):
        job_dir,output_dir,attrs,pre_fn = self._parse_args("rmsd",args)
        if not Path(job_dir,"prod_npt-pbc.xtc").exists():
            ndx = self.generate_group(job_dir)
            self.trj_pbc(job_dir,ndx)
        for group in attrs:
            with open(f"{job_dir}/{group}_rmsd.in",'w') as outf:
                outf.write("%s\n"%group)
                outf.write("%s\n\n\n"%group)
            os.system(f"{gmx} rms -f {job_dir}/prod_npt-pbc.xtc -s {job_dir}/prod_npt.tpr -o {job_dir}/rmsd-{group}.xvg -n {job_dir}/index.ndx < {job_dir}/{group}_rmsd.in")
            labels, units,datas = self.parse_xvg_file(f"{job_dir}/rmsd-{group}.xvg")
            labels[1] = group
            self.file_figure(output_dir,labels,units,datas,pre_fn=pre_fn)

    def get_cluster(self,job_dir,args=None):
        job_dir,output_dir,group_idx,cluster_group = self._parse_args(job_dir,args,attrs=["output_dir","group_idx","cluster_group"])
        os.system(f"{gmx} cluster  -f prod_npt-pbc.xtc -s prod_npt.tpr -n index.ndx -cl cluster.pdb -cuttoff X")

    def get_lig_torsion_statistical(self,job_dir,output_dir,molecule,torsion_energy,xvg_flag=True,save_figure=True):
        
        torsion_distribute = {term:[] for term in torsion_energy.keys()}
        torsion_term = {term:[int(an) for an in term.split("-")] for term in torsion_distribute.keys()}
        
        if not Path(job_dir,"prod_npt-all.gro").exists:
            with open(f"{job_dir}/trj.in",'w') as outf:
                outf.write("System\n")
                outf.write("\n\n\n\n")
            os.system(f"{gmx} trjconv -f {job_dir}/prod_npt-pbc.xtc -s {job_dir}/prod_npt.tpr -o {job_dir}/prod_npt-all.gro < {job_dir}/trj.in")
        coors = self.parse_gro_file(f"{job_dir}/prod_npt-lig.gro")
        
        for coor in coors:
            for term,ans in torsion_term.items():
                torsion_distribute[term].append(calc_stru_para([coor[ans[0]],coor[ans[1]],coor[ans[2]],coor[ans[3]]]))
        
        return torsion_distribute

    def get_pair_lig_torsion_statistical(self,job_dir,output_dir,molecules,torsion_energys):
        from craton import molxpert as MX
        molecules = MX.molecule_structure(molecules)
        molecules = MX.molecule_torsion(molecules)
        torsion_energys = Calc._torsion_scan(molecules)
        
        
        bfepath = Path(job_dir, "rbfe")
        hfepath = Path(job_dir,"rhfe")
        pair_name = job_dir.parent.name 
        
        dds = [int(dd) for dd in os.listdir(bfepath) if Path(bfepath,dd).is_dir() if dd not in ["job_info"]]
        adir = "0"
        bdir = str(max(dds))
        
        a_bfe_distribution = self.get_lig_torsion_statistical(Path(bfepath,adir),output_dir,molecules[0],torsion_energys[0],xvg_flag=True,save_figure=False)
        b_bfe_distribution = self.get_lig_torsion_statistical(Path(bfepath,bdir),output_dir,molecules[1],torsion_energys[1],xvg_flag=True,save_figure=False)
        
        a_hfe_distribution = self.get_lig_torsion_statistical(Path(hfepath,adir),output_dir,molecules[0],torsion_energys[0],xvg_flag=False,save_figure=False)
        b_hfe_distribution = self.get_lig_torsion_statistical(Path(hfepath,bdir),output_dir,molecules[1],torsion_energys[1],xvg_flag=False,save_figure=False)
    
    def _set_molecule_dir(self,job_dir):
        if Path(job_dir,"job_info").exists():
            return Path(job_dir,"job_info"), Path(job_dir,"job_info","md_setting.json")
        parent_1 = job_dir.parent
        if Path(parent_1,"job_info").exists():
            return Path(parent_1,"job_info"), Path(parent_1,"job_info","md_setting.json")
        
    def _set_molecule_dir_old(self,job_dir):
        protein_file = None
        ligand_file = None
        if job_dir.name.isdigit():
            if job_dir.name == "0":
                ligand_file_name = "left.mol"
            else:
                ligand_file_name = "right.mol"
        else:
            ligand_file_name = f"{job_dir.name}.mol"
            
        if Path(job_dir,"protein.pdb").exists and Path(job_dir,ligand_file_name).exists():
            return Path(job_dir,"protein.pdb"),Path(job_dir,ligand_file_name)

        parent_1 = job_dir.parent
        if Path(parent_1,"protein.pdb").exists() and Path(parent_1,ligand_file_name).exists():
            return Path(parent_1,"protein.pdb"),Path(parent_1,ligand_file_name)
        else:
            if Path(parent_1,"intermedital").exists():
                if Path(parent_1,"intermedital","molecules","protein.pdb").exists():
                    protein_file = Path(parent_1,"intermedital","molecules","protein.pdb")
                if Path(parent_1,"intermedital","molecules",ligand_file_name).exists():
                    ligand_file = Path(parent_1,"intermedital","molecules",ligand_file_name)
                else:
                    if Path(parent_1,"intermedital","molecules",job_dir.name,ligand_file_name).exists():
                        ligand_file = Path(parent_1,"intermedital","molecules",job_dir.name,ligand_file_name)
                    #else:
                    #    if Path(parent_1,"intermedital","molecules",job_dir.parent.name,ligand_file).exists():
                    #        ligand_file = Path(parent_1,"intermedital","molecules",job_dir.parent.name,ligand_file)
                return protein_file,ligand_file
            
        parent_2 = job_dir.parent.parent
        if Path(parent_2,"protein.pdb").exists() and Path(parent_2,ligand_file_name).exists():
            return Path(parent_2,"protein.pdb"),Path(parent_2,ligand_file_name)
        else:
            if Path(parent_2,"intermedital").exists():
                if Path(parent_2,"intermedital","molecules","protein.pdb").exists():
                    protein_file = Path(parent_2,"intermedital","molecules","protein.pdb")
                if Path(parent_2,"intermedital","molecules",ligand_file_name).exists():
                    ligand_file = Path(parent_2,"intermedital","molecules",ligand_file_name)
                else:
                    if Path(parent_2,"intermedital","molecules",job_dir.name,ligand_file_name).exists():
                        ligand_file = Path(parent_2,"intermedital","molecules",job_dir.name,ligand_file_name)
                    else:
                        if Path(parent_2,"intermedital","molecules",job_dir.parent.name,ligand_file_name).exists():
                            ligand_file = Path(parent_2,"intermedital","molecules",job_dir.parent.name,ligand_file_name)
                        #else:
                        #    if Path(parent_2,"intermedital","molecules",job_dir.parent.parent.name,ligand_file).exists():
                        #        ligand_file = Path(parent_2,"intermedital","molecules",job_dir.parent.parent.name,ligand_file)
                return protein_file,ligand_file
            
        parent_3 = job_dir.parent.parent.parent
        
        if Path(parent_3,"intermedital").exists():
            if Path(parent_3,"intermedital","molecules","protein.pdb").exists():
                protein_file = Path(parent_3,"intermedital","molecules","protein.pdb")
            if Path(parent_3,"intermedital","molecules",ligand_file_name).exists():
                ligand_file = Path(parent_3,"intermedital","molecules",ligand_file_name)
            else:
                if Path(parent_3,"intermedital","molecules",job_dir.name,ligand_file_name).exists():
                    ligand_file = Path(parent_3,"intermedital","molecules",job_dir.name,ligand_file_name)
                else:
                    if Path(parent_3,"intermedital","molecules",job_dir.parent.name,ligand_file_name).exists():
                        ligand_file = Path(parent_3,"intermedital","molecules",job_dir.parent.name,ligand_file_name)
                    else:
                        if Path(parent_3,"intermedital","molecules",job_dir.parent.parent.name,ligand_file_name).exists():
                            ligand_file = Path(parent_3,"intermedital","molecules",job_dir.parent.parent.name,ligand_file_name)
                        #else:
                        #    if Path(parent_3,"intermedital","molecules",job_dir.parent.parent.parent.name,ligand_file).exists():
                        #        ligand_file = Path(parent_3,"intermedital","molecules",job_dir.parent.parent.parent.name,ligand_file)
            return protein_file,ligand_file
            
        return protein_file,ligand_file
                
    def get_dimer(self,args,idx=None):
        def get_interaction_script(acts):
            ss = ""
            for act in acts.__interaction__:
                ss += f"{act.acceptor.site_script}&{act.donor.site_script}$"
            return ss[:-1]
        
        random_n = 10
        from ...chem.chemsystem import System
        job_dir,output_dir,molecules,molecule_dir,fn, save_fig, dimer_type = self._parse_args("interaction",args)
        output_dir=Path(output_dir)
        job_dir = Path(job_dir).resolve()
        if output_dir != "db":
            Path(output_dir).mkdir(exist_ok=True)
        if molecules is None:
            if molecule_dir is None:
                molecule_path,md_setting_file = self._set_molecule_dir(job_dir)

            from craton import molxpert as MX
            import json
            md_setting = json.loads(open(md_setting_file).read())
            molecule1 = MX.molecule_create(f"{molecule_path}/{md_setting['molecules'][0]}.mtx",parallel=False)[0]
            #molecule1 = MX.molecule_structure(molecule1)[0]
            

            molecule2  = MX.molecule_create(f"{molecule_path}/{md_setting['molecules'][1]}.mtx",parallel=False)[0]
            molecule2 = MX.protein_structure(molecule2)
            molecules = [molecule1,molecule2]
        _dimer_types_ = dimer_type.split("_")
        
        molecule1 = molecules[0]
        molecule2 = molecules[1]
        for atom in molecule1.Atoms:
            atom.residue = _dimer_types_[0]
        for atom in molecule2.Atoms:
            atom.residue = _dimer_types_[1]
        
        molecule1.shift_an = 0
        molecule2.shift_an = len(molecule1.Atoms)
        
        try:
            if not Path(job_dir,"prod_npt-all.gro").exists():
                if not Path(job_dir,"prod_npt-pbc.xtc").exists():
                    ndx = self.generate_group(job_dir)
                    self.trj_pbc(job_dir,ndx)
                with open(f"{job_dir}/trj.in",'w') as outf:
                    outf.write("System\n")
                    outf.write("\n\n\n\n")
                os.system(f"{gmx} trjconv -f {job_dir}/prod_npt-pbc.xtc -s {job_dir}/prod_npt.tpr -pbc whole -o {job_dir}/prod_npt-all.gro < {job_dir}/trj.in")
            total_coors = self.parse_gro_file(f"{job_dir}/prod_npt-all.gro")
            total_interactions = IM.get_interaction_model(molecule1,molecule2,coordinates=total_coors,probe_flag=False,parallel=False)
        except:
            logger.error(f"{job_dir} is error")
            if idx is None:
                return []
            else:
                return [],idx
        datas = []
        
        tmp_dicts_inter = {}
        none_inter = []
        for ii,interaction in enumerate(total_interactions):
            if len(interaction.__interaction__) > 0:# is not None:
                ss = get_interaction_script(interaction)
                if ss not in tmp_dicts_inter:
                    tmp_dicts_inter[ss] = []
                tmp_dicts_inter[ss].append([ii,interaction.__interaction__[0],interaction.__interaction__[0].type,interaction.__interaction__[0].subtype,
                                        ss,interaction.__interaction__[0].distance])
            else:
                none_inter.append([ii,None,None,None,None,None])
                
        act_ndx_info = []  
        for kk,vv in tmp_dicts_inter.items():
            vvv = sorted(vv,key=lambda x:x[5])
            act_ndx_info.append(vvv[0])
        if random_n > 0:
            if len(none_inter) <= random_n:
                act_ndx_info.extend(none_inter)
            else:
                rad_ndx =random.sample([ii for ii in range(len(none_inter))],random_n)
                for ii in rad_ndx:
                    act_ndx_info.append(none_inter[ii]) 
        
        for nn,rr in enumerate(act_ndx_info):
            ii = rr[0]
            inter = rr[1]
            interaction_type = rr[2]
            interaction_subtype = rr[3]
            interaction_script = rr[4]
            
            if inter is not None:
                interaction = {"type":inter.type,"subtype":inter.subtype,"distance":inter.distance,"angle":inter.angle,"offset":inter.offset,"auxi_angle":inter.auxi_angle,"ring_in_protein":inter.ring_in_protein}
                if inter.acceptor.site_script[:3] == molecule1.Atoms[0].residue:
                    interaction["site1"] = inter.acceptor.site_script
                    interaction["site2"] = inter.donor.site_script
                else:
                    interaction["site2"] = inter.acceptor.site_script
                    interaction["site1"] = inter.donor.site_script
            else:
                interaction = None
            
            this_system = System("dimer")
            this_system.name = f"{molecule1.mole_name}_{molecule2.mole_name}_{interaction_type}_{ii}"
            this_molecule1 = deepcopy(molecule1)
            this_molecule2 = deepcopy(molecule2)
            coor = total_coors[ii]
            shift_an = this_molecule1.shift_an
            for atom in this_molecule1.Atoms:
                atom.coor = coor[atom.ID + shift_an]
            shift_an = this_molecule2.shift_an
            for atom in this_molecule2.Atoms:
                atom.coor = coor[atom.ID + shift_an]
            this_system.molecules = [this_molecule1,this_molecule2]
            this_system.molecule_number = [1,1]
            this_system.molecule_types = [this_molecule1.mole_name,this_molecule2.mole_name]
            this_system.dimer_type = dimer_type
            this_system.interaction = interaction_type
            this_system.interaction_subtype = interaction_subtype
            this_system.interaction_script = interaction_script
            this_system.source = "dimer_simulation"
            this_system.interaction_detail = interaction
            datas.append(this_system)
        if output_dir.name == "db":
            from ...chem.database.mongodb import CompoundDB
            this_compounddb = CompoundDB(configure={"datasearch":{"compound_style":"dimer"}})
            this_compounddb.insert_to_db(datas)
                
        
        #model_report = ReportMultipleInteraction(total_interactions)
        #df = model_report.info_df(show_figure=False)
        #model_report.info_detail(f"{output_dir}/{fn}interaction_detail.txt")
        
        if idx is None:
            return total_interactions
        else:
            logger.info(f"##########################Done:  {idx} {job_dir} #######################")
            return total_interactions,idx
    
    def get_pl_interaction(self,args):
        job_dir,output_dir,molecules,molecule_dir,fn, save_fig, __ = self._parse_args("interaction",args)
        job_dir = Path(job_dir).resolve()
        Path(output_dir).mkdir(exist_ok=True)
        if molecules is None:
            if molecule_dir is None:
                molecule_path,md_setting_file = self._set_molecule_dir(job_dir)

            from craton import molxpert as MX
            import json
            md_setting = json.loads(open(md_setting_file).read())
            try:
                protein = MX.molecule_create(f"{molecule_path}/{md_setting['protein']}.mtx", parallel=False)[0]
                protein = MX.protein_structure(protein)
                ligand  = MX.molecule_create(f"{molecule_path}/{md_setting['ligand']}.mtx", parallel=False)[0]
            except:
                protein = MX.molecule_create(f"{molecule_path}/{md_setting['molecules'][1]}.mtx", parallel=False)[0]
                protein = MX.protein_structure(protein)
                ligand  = MX.molecule_create(f"{molecule_path}/{md_setting['molecules'][0]}.mtx", parallel=False)[0]
            molecules = [protein,ligand]
            #protein = MX.protein_create(protein_molecule_dir)
            #protein = MX.protein_structure(protein)
            #ligands = MX.molecule_create(str(ligand_molecule_dir),extra_var={"smiles_flag":False},show_figure=False)
        
            
        protein = molecules[0]
        ligand = molecules[1]
        
        ligand.shift_an = 0
        protein.shift_an = len(ligand.Atoms)
        if not Path(job_dir,"prod_npt-all.gro").exists():
            if not Path(job_dir,"prod_npt-pbc.xtc").exists():
                ndx = self.generate_group(job_dir)
                self.trj_pbc(job_dir,ndx)
            with open(f"{job_dir}/trj.in",'w') as outf:
                outf.write("System\n")
                outf.write("\n\n\n\n")
            os.system(f"{gmx} trjconv -f {job_dir}/prod_npt-pbc.xtc -s {job_dir}/prod_npt.tpr -pbc whole -o {job_dir}/prod_npt-all.gro < {job_dir}/trj.in")
        total_coors = self.parse_gro_file(f"{job_dir}/prod_npt-all.gro")
        total_interactions = IM.get_interaction_model(protein,ligand,coordinates=total_coors,parallel=False)
        
        
        model_report = ReportMultipleInteraction(total_interactions)
        df = model_report.info_df(show_figure=False)
        model_report.info_detail(f"{output_dir}/{fn}interaction_detail.txt")
        if save_fig:
            normal_md_bar_plot(df,output_file=f"{output_dir}/{fn}interaction.html")
        
        return total_interactions, df
    
    def get_pair_pl_interaction(self,args):
        job_dir,output_dir,molecules,molecule_dir = self._parse_args("pair_interaction",args)
        job_dir = Path(job_dir).resolve()
        Path(output_dir).mkdir(exist_ok=True)
        
        bfepath = Path(job_dir, "bfe")
        hfepath = Path(job_dir,"hfe")
        pair_name = job_dir.parent.name    
    
        dds = [int(dd) for dd in os.listdir(bfepath) if Path(bfepath,dd).is_dir() if dd not in ["job_info"]]
        adir = "0"
        bdir = str(max(dds))
            
        #a_bfe_interaction, a_bfe_df = self.get_pl_interaction(Path(bfepath,adir),output_dir,molecules=[deepcopy(molecules[0]),molecules[1]],fn="a_bfe_")
        a_bfe_interaction, a_bfe_df = self.get_pl_interaction({"job_dir":Path(bfepath,adir),"output_dir":output_dir,"fn":"a_bfe_","save_fig":False})
        #b_bfe_interaction, b_bfe_df = self.get_pl_interaction(Path(bfepath,bdir),output_dir,molecules=[deepcopy(molecules[0]),molecules[2]],fn="b_bfe_")
        b_bfe_interaction, b_bfe_df = self.get_pl_interaction({"job_dir":Path(bfepath,bdir),"output_dir":output_dir,"fn":"b_bfe_","save_fig":False})
        
        merge_a_bfe_df, merge_b_bfe_df = merge_interaction_dfs(a_bfe_df, b_bfe_df)
        
        
        fep_bar_plot(
                merge_a_bfe_df, merge_b_bfe_df,
                pair_name, output_dir / "interaction_result.html"
            )
    
        return merge_a_bfe_df, merge_b_bfe_df
        
    def get_unks(self,job_dir):
        dds = [dd for dd in os.listdir(job_dir) if os.path.isdir(f"{job_dir}/{dd}") and dd not in ["job_info"]]
        try:
            csvfiles = [Path(f"{job_dir}/{dd}/prod_npt.csv") for dd in dds]
            u_nks = [self.parse_csv_file(subfile, temperature=self.temperature) for subfile in csvfiles]
        except:
            xvgfiles = [Path(f"{job_dir}/{dd}/prod_npt.xvg") for dd in dds]
            u_nks = [self.parse_xvg_dh_file(subfile) for subfile in xvgfiles]
        return u_nks

    def get_bar(self,args):
        job_dir,output_dir,pre_fn = self._parse_args("bar",args)
        #job_dir,output_dir = self._parse_args(job_dir,args,attrs=["output_dir"])
        u_nks = self.get_unks(job_dir)
        bar_data = {}
        bar_estimator = BarForDDG(u_nks)
        dFs, df_error, df_bootstrapping_error, sa, sb = bar_estimator.run()
        bar_data["ddg"] = dFs
        bar_data["ddg_error"] = df_error
        bar_data["ddg_bootstrapping_error"] = df_bootstrapping_error
        bar_data["sa"] = sa
        bar_data["sb"] = sb
        bar_data["window"] = [f"{i} -> {i+1}" for i in range(len(dFs))]

        _final = {
            "ddg": bar_data["ddg"].sum(),
            "ddg_error": np.linalg.norm(bar_data["ddg_error"]),
            "ddg_bootstrapping_error": np.linalg.norm(bar_data["ddg_bootstrapping_error"]),
            "window": "total",
        }
        
        bar_data = pd.concat((pd.DataFrame.from_records(bar_data), pd.DataFrame.from_records([_final])), ignore_index=True)
        bar_data.set_index("window", inplace=True)

        bar_data.to_csv(Path(output_dir, f"{pre_fn}bar.csv"))
        sub_df = bar_data[["ddg", "ddg_error"]]
        return sub_df.to_dict(orient="index")
        self.dump_data["energy"] = sub_df.to_dict(orient="index")
                
    def get_ddg(self,args):
        job_dir,output_dir = self._parse_args("ddg",args)
        bfepath = Path(job_dir, "bfe")
        hfepath = Path(job_dir,"hfe")

        bfe = self.get_bar({"job_dir":bfepath,"output_dir":output_dir,"pre_fn":"bfe_"})
        hfe = self.get_bar({"job_dir":hfepath,"output_dir":output_dir,"pre_fn":"hfe_"})
        
        ddg = bfe["total"]["ddg"] - hfe["total"]["ddg"]
        ddg_error = max([bfe["total"]["ddg_error"] , hfe["total"]["ddg_error"]])
        datas = { "ddg":ddg,
                 "ddg_error":ddg_error,
                 "bfe_ddg":bfe["total"]["ddg"],
                 "bfe_ddg_error":bfe["total"]["ddg_error"],
                 "hfe_ddg":hfe["total"]["ddg"],
                 "hfe_ddg_error":hfe["total"]["ddg_error"]
                 }
        return datas
        
    def get_exchange(self,args):
        job_dir,output_dir,pre_fn = self._parse_args("exchange",args)
        def _count_array(count_dict, lambda_array):
            for i, item in enumerate(lambda_array):
                count_dict[i][item] += 1
        
        analyze_data = defaultdict(lambda: defaultdict(lambda: 0))
        #job_dir,output_dir = self._parse_args(job_dir,args,attrs=["output_dir"])
        log_file = "0/prod_npt.log"
        with open(Path(f"{job_dir}/{log_file}")) as f:
            lambda_array = None
            count = 0
            for line in f:
                if line.startswith("Replica exchange in lambda"):
                    lambda_array = [int(item) for item in next(f).split()]
                    lambda_array = [f"replica_{item:02d}" for item in lambda_array]
                    _count_array(analyze_data, lambda_array)
                if lambda_array and line.startswith("Repl ex"):
                    count += 1
                    token = line.split()[2:]
                    for i in range(len(token)):
                        if token[i] == "x":
                            left = int(token[i - 1])
                            right = int(token[i + 1])
                            lambda_array[left], lambda_array[right] = lambda_array[right], lambda_array[left]
                    _count_array(analyze_data, lambda_array)

        frames = []
        for _, value in analyze_data.items():
            frames.append(pd.Series(value))
        df = pd.concat(frames, axis=1)
        df = df.sort_index().T
        ####self.dump_data = df.to_dict(orient="index")
        plot_exchange_rate(df, Path(output_dir, f"{pre_fn}replica_exchange.html"))
        return df.to_dict(orient="index")

    def get_block_ddg(self,args):
        job_dir,output_dir,pre_fn,block_num = self._parse_args("block_ddg",args)
        from alchemlyb.estimators import MBAR as _MBAR
        #job_dir,output_dir,u_nks,num_points = self._parse_args(job_dir,args,attrs=["output_dir","u_nks","block_num"])
        u_nks = self.get_unks(job_dir=job_dir)
        num_points = block_num
        
        analyze_data = pd.DataFrame()
        analyze_data.attrs = {"temperature": self.temperature, "energy_unit": "kcal/mol"}
        
        forward, forward_error = [], []
        for i in range(1, num_points + 1):
            slice_right = int(len(u_nks[0]) / num_points * i)
            slice_left = int(len(u_nks[0]) / num_points * (i - 1))
            u_nk_current = pd.concat([data[slice_left:slice_right] for data in u_nks])
            estimate = _MBAR(method="L-BFGS-B").fit(u_nk_current)
            forward.append(estimate.delta_f_.iloc[0, -1])
            forward_error.append(estimate.d_delta_f_.iloc[0, -1])

        analyze_data["ddg"] = forward
        analyze_data["ddg_error"] = forward_error
        analyze_data /= BETA * KCAL_MOL_2_KJ_MOL

        ax = plot_convergence_block(analyze_data["ddg"], analyze_data["ddg_error"])
        ax.figure.savefig(Path(output_dir, f"{pre_fn}convergence_accumulate_block.png"))
        return analyze_data.to_dict(orient="index")

    def get_accumulate_ddg(self,args):
        job_dir,output_dir,pre_fn,block_num = self._parse_args("block_ddg",args)
        from alchemlyb.estimators import MBAR as _MBAR
        ###job_dir,output_dir,u_nks,num_points = self._parse_args(job_dir,args,attrs=["output_dir","u_nks","block_num"])
        u_nks = self.get_unks(job_dir=job_dir)
        num_points = block_num
        
        
        analyze_data = pd.DataFrame()
        analyze_data.attrs = {"temperature": self.temperature, "energy_unit": "kcal/mol"}
        ddg_final = False
        #num_points = self._kwargs.get("num_points", 5)
        #ddg_final, ddg_error_final = self._kwargs.get("ddg"), self._kwargs.get("ddg_error")
        forward, forward_error, backward, backward_error = [], [], [], []
        for i in range(1, num_points + 1):
            slice_right = int(len(u_nks[0]) / num_points * i)
            u_nk_current = pd.concat([data[:slice_right] for data in u_nks])
            estimate = _MBAR(method="L-BFGS-B").fit(u_nk_current)
            forward.append(estimate.delta_f_.iloc[0, -1])
            forward_error.append(estimate.d_delta_f_.iloc[0, -1])

            u_nk_current = pd.concat([data[-slice_right:] for data in u_nks])
            estimate = _MBAR(method="L-BFGS-B").fit(u_nk_current)
            backward.append(estimate.delta_f_.iloc[0, -1])
            backward_error.append(estimate.d_delta_f_.iloc[0, -1])

        if ddg_final:
            forward.append(ddg_final)
            forward_error.append(ddg_error_final)
            backward.append(ddg_final)
            backward_error.append(ddg_error_final)

        analyze_data["forward"] = forward
        analyze_data["forward_error"] = forward_error
        analyze_data["backward"] = backward
        analyze_data["backward_error"] = backward_error
        analyze_data /= BETA * KCAL_MOL_2_KJ_MOL

        ax = plot_convergence(
            **analyze_data,
            units="kcal/mol",
        )
        ax.figure.savefig(Path(output_dir, f"{pre_fn}convergence_accumulate_plot.png"))
        return analyze_data.to_dict(orient="index")

    #def get_rbfe(self,job_dir,attrs=None,output_dir=Path("."),molecule_dir=None,idx = None):
    def get_rbfe(self,args,idx=None):
        job_dir,output_dir,attrs,molecule_dir,molecules= self._parse_args("rbfe",args)
        Path(output_dir).mkdir(exist_ok=True)
        datas = {}
        if attrs == "all":
            attrs = ["interaction","energy","rmsd","exchange"]
        if attrs is not None:
            for run_dir in ["bfe","hfe"]:
                this_job_dir = Path(job_dir,run_dir)
                dds = [int(dd) for dd in os.listdir(Path(this_job_dir)) if Path(this_job_dir,dd).is_dir() and dd not in ["job_info"]]
            
                adir = Path(this_job_dir,"0")
                bdir = Path(this_job_dir,str(max(dds)))
                ndxa = self.generate_group(adir)
                self.trj_pbc(adir,ndxa)
                ndxb = self.generate_group(bdir)
                self.trj_pbc(bdir,ndxb)
            
            
                if run_dir == "hfe":
                    rmsd_attrs = ["LIG"]
                else:
                    rmsd_attrs = ["Protein","LIG"]
                if "energy" in attrs:
                    self.get_energy({"job_dir":adir,"output_dir":output_dir,"pre_fn":f"{run_dir}_a_"})
                    self.get_energy({"job_dir":bdir,"output_dir":output_dir,"pre_fn":f"{run_dir}_b_"})
                if "rmsd" in attrs:
                    if run_dir == "bfe":
                        self.get_rmsd({"job_dir":adir,"output_dir":output_dir,"attrs":rmsd_attrs,"pre_fn":f"{run_dir}_a_"})
                        self.get_rmsd({"job_dir":bdir,"output_dir":output_dir,"attrs":rmsd_attrs,"pre_fn":f"{run_dir}_b_"})
            
                if "exchange" in attrs:
                    datas[f"{run_dir}_exchanged"] = self.get_exchange({"job_dir":this_job_dir,"output_dir":output_dir,"pre_fn":f"{run_dir}_"})
                if "accum_ddg" in attrs:
                    datas[f"{run_dir}_accumul_ddg"] = self.get_accumulate_ddg({"job_dir":this_job_dir,"output_dir":output_dir,"pre_fn":f"{run_dir}_"})
                if "block_ddg" in attrs:
                    datas[f"{run_dir}_block_ddg"] = self.get_block_ddg({"job_dir":this_job_dir,"output_dir":output_dir,"pre_fn":f"{run_dir}_"})
            
            if "interaction" in attrs:
                left_interactions, right_interactions = self.get_pair_pl_interaction({"job_dir":Path(job_dir),
                                                                                          "output_dir":Path(output_dir),
                                                                                          "molecule_dir":molecule_dir,
                                                                                          "molecules":molecules,})
                datas["left_interactions"] = left_interactions
                datas["right_interactions"] = right_interactions
                
        os.system(f"cp {job_dir}/bfe/job_info/*_atom_mapping.png {output_dir}")
        os.system(f"cp {job_dir}/bfe/job_info/*_atom_mapping_nonH.png {output_dir}")
        ddg_datas = self.get_ddg({"job_dir":job_dir,"output_dir":output_dir})
        datas["ddg"] = ddg_datas
        
        if idx == None:
            return datas
        else:
            return datas,idx

    def recalculate_cycle_closure(self,ddg_df,dg_df,input=None,pka_file=None):
        if input is not None:
            with open(input, "r") as f:
                cc_data = simplejson.load(f)
            ddg_df = pd.DataFrame.from_dict(cc_data["edge"], orient="index")
            dg_df = pd.DataFrame.from_dict(cc_data["node"], orient="index")
            ddg_df["name"] = ddg_df.index
            dg_df["name"] = dg_df.index
        cc_ddg_df, cc_dg_df = CycleClosureFreeEnergy.run_cycle_closure(ddg_df, dg_df,pka_file=pka_file)
        return cc_ddg_df,cc_dg_df

    #def get_all_rbfe(self,job_dir,output_dir,attrs=None,expt_file=None,pka_file=None,parallel=False):
    def get_all_rbfe(self,args):
        job_dir,output_dir,attrs,expt_file,pka_file,parallel = self._parse_args("all_rbfe",args)
        Path(output_dir).mkdir(exist_ok=True)
        dds = [dd for dd in os.listdir(job_dir) if dd.find("_to_") != -1]
        output_dirs = [Path(output_dir,dd) for dd in dds]
        input_dirs = [Path(job_dir,dd) for dd in dds]
        
        if parallel:
            _args = [{"job_dir":this_job_dir,"output_dir":output_dirs[ii],"attrs":attrs,"idx":ii} for ii,this_job_dir in enumerate(input_dirs) ]
            _arr = parallel_run(self.get_rbfe,_args)
        else:
            _arr = []
            for ii,this_job_dir in enumerate(input_dirs):
                _arr.append(self.get_rbfe({"job_dir":this_job_dir,"output_dir":output_dirs[ii],"attrs":attrs}))
        datas = {d_name:_arr[ii] for ii,d_name in enumerate(dds)}
        
        
        all_pair_bar_result = defaultdict(list)
        for key,value in datas.items():
            all_pair_bar_result["name"].append(key)
            all_pair_bar_result["ddg"].append(value["ddg"]["ddg"])
            all_pair_bar_result["ddg_error"].append(value["ddg"]["ddg_error"])
            
            #for sim_type in ["bfe", "hfe"]:
            #    all_pair_bar_result[sim_type].append(value["ddg"][f"{sim_type}"])
            #    all_pair_bar_result[f"{sim_type}_error"].append(value["bar"][sim_type]["energy"]["total"]["ddg_error"])
        ddg_df = pd.DataFrame(all_pair_bar_result)
        dg_df = pd.DataFrame(columns=["name", "dg"])
        if expt_file:
            dg_df = CycleClosureFreeEnergy._get_exp_dg(expt_file)
            CycleClosureFreeEnergy._add_exp_ddg_column(ddg_df, dg_df)
        cc_ddg_df,cc_dg_df = self.recalculate_cycle_closure(ddg_df,dg_df,pka_file=pka_file)
        WriteExcel.write_rbfe_excel(cc_ddg_df,cc_dg_df,output_dir,exp_file=expt_file,two_stages=False,pic_dir=job_dir)

        image_dir = Path(output_dir, "img")
        image_dir.mkdir(exist_ok=True)
        molecules_copied = set()
        for d_name in dds:
            lig1, lig2 = d_name.split("_to_")
            for lig in [lig1, lig2]:
                if lig not in molecules_copied:
                    src = Path(job_dir, d_name, "bfe", "job_info", f"{lig}.png")
                    if src.exists():
                        shutil.copy(src, image_dir)
                    molecules_copied.add(lig)

        return datas
        
    def get_abfe(self,args, idx=None): #job_dir,output_dir=Path("."),molecule_dir=None,idx = None):
        job_dir,output_dir,attrs,molecule_dir,molecules = self._parse_args("rbfe",args)
        Path(output_dir).mkdir(exist_ok=True)
        datas = {}
        if attrs == "all":
            attrs = ["interaction","energy","rmsd"]
        if attrs is not None:
            for run_dir in ["bfe","hfe"]:
                this_job_dir = Path(job_dir,run_dir)
                #dds = [int(dd) for dd in os.listdir(Path(this_job_dir)) if Path(this_job_dir,dd).is_dir() and dd not in ["job_info"]]
            
                adir = Path(this_job_dir,"0")
                #bdir = Path(this_job_dir,str(max(dds)))
                ndxa = self.generate_group(adir)
                self.trj_pbc(adir,ndxa)
                #ndxb = self.generate_group(bdir)
                #self.trj_pbc(bdir,ndxb)
                if run_dir == "hfe":
                    rmsd_attrs = ["LIG"]
                else:
                    rmsd_attrs = ["Protein","LIG"]
                if "energy" in attrs:
                    self.get_energy({"job_dir":adir,"output_dir":output_dir,"pre_fn":f"{run_dir}_a_"})
                if "rmsd" in attrs:
                    if run_dir == "bfe":
                        self.get_rmsd({"job_dir":adir,"output_dir":output_dir,"attrs":rmsd_attrs,"pre_fn":f"{run_dir}_a_"})
            
                if "exchange" in attrs:
                    datas[f"{run_dir}_exchanged"] = self.get_exchange({"job_dir":this_job_dir,"output_dir":output_dir,"pre_fn":f"{run_dir}_"})
                if "accum_ddg" in attrs:
                    datas[f"{run_dir}_accumul_ddg"] = self.get_accumulate_ddg({"job_dir":this_job_dir,"output_dir":output_dir,"pre_fn":f"{run_dir}_"})
                if "block_ddg" in attrs:
                    datas[f"{run_dir}_block_ddg"] = self.get_block_ddg({"job_dir":this_job_dir,"output_dir":output_dir,"pre_fn":f"{run_dir}_"})
            
            if "interaction" in attrs:
                interactions_models,interactions_df = self.get_pl_interaction({"job_dir":Path(job_dir,"bfe","0"),
                                                                                          "output_dir":Path(output_dir),
                                                                                          "molecule_dir":molecule_dir,
                                                                                          "molecules":molecules,})
                datas["interactions"] = interactions_models
        
        
        ddg_datas = self.get_ddg({"job_dir":job_dir,"output_dir":output_dir})
        datas["ddg"] = ddg_datas
        
        #interactions = self.get_pl_interaction(Path(job_dir),Path(output_dir),molecules=[protein,ligand],save_fig=True)
        if idx == None:
            return datas
        else:
            return datas,idx
    
    def get_all_dimer(self,args):
        job_dir,output_dir,molecules,molecule_dir,fn, save_fig, dimer_type,parallel = self._parse_args("all_dimer",args)
        dds = [dd for dd in os.listdir(job_dir) if os.path.isdir(f"{job_dir}/{dd}")]
        
        if output_dir.name != "db":
            Path(output_dir).mkdir(exist_ok=True)
            output_dirs = [{"output_dir":Path(output_dir,dd),"dimer_type":dimer_type,"save_fig":save_fig} for dd in dds]
        else:
            output_dirs = [{"output_dir":output_dir,"dimer_type":dimer_type,"save_fig":save_fig} for dd in dds]
        input_dirs = [Path(job_dir,dd) for dd in dds]
        if parallel:
            for ii, job_dir in enumerate(input_dirs):
                output_dirs[ii]["job_dir"] = job_dir
            _arr = parallel_run(self.get_dimer,output_dirs)
        else:
            _arr = []
            for ii,this_job_dir in enumerate(input_dirs):
                print("##########################runing to :  ",ii," #######################")
                try:
                    this_arg = output_dirs[ii]
                
                    this_arg["job_dir"] =this_job_dir
                    _arr.append(self.get_dimer(this_arg))
                except:
                    logger.error(f"error md simulation: {this_job_dir}")
    
    def get_total_dimer(self,args):
        job_dir,output_dir,molecules,molecule_dir,fn, save_fig, dimer_type,parallel = self._parse_args("all_dimer",args)
        ds = [dd for dd in os.listdir(job_dir) if os.path.isdir(f"{job_dir}/{dd}")]
        for dd in ds:
            self.get_all_dimer({"job_dir":Path(job_dir,dd),"output_dir":output_dir,"dimer_type":dimer_type,"save_fig":save_fig})               
    
    def get_all_abfe(self,args):
        job_dir,output_dir,attrs,expt_file,pka_file,parallel = self._parse_args("all_rbfe",args)
        Path(output_dir).mkdir(exist_ok=True)
        dds = [Path(job_dir,dd) for dd in os.listdir(job_dir)] 
        dds = [dd for dd in dds if dd.is_dir()]
        d_names = [dd.name for dd in dds]
        
        output_dirs = [Path(output_dir,dd) for dd in d_names]
        if parallel:
            _args = [{"job_dir":this_job_dir,"output_dir":output_dirs[ii],"attrs":attrs,"idx":ii} for ii,this_job_dir in enumerate(dds) ]
            _arr = parallel_run(self.get_abfe,_args)
        else:
            _arr = []
            for ii,this_job_dir in enumerate(dds):
                _arr.append(self.get_abfe({"job_dir":this_job_dir,"output_dir":output_dirs[ii],"attrs":attrs}))
        
        
        _datas = {d_name:_arr[ii] for ii,d_name in enumerate(d_names)}
        datas = {}
        _label = {"ddg":"cc_dg","ddg_error":"dg_error","bfe_ddg":"bfe_dg","bfe_ddg_error":"bfe_dg_error","hfe_ddg":"hfe_dg","hfe_ddg_error":"hfe_dg_error"}
              
        datas = {lig_name:{attr_t:data["ddg"][attr_r] for attr_r, attr_t in _label.items()} for lig_name,data in _datas.items()}

        if expt_file is not None:
            expt_datas = {}
            with open(expt_file) as inf:
                lines = inf.readlines()
            for line in lines[1:]:
                ss = line.split(",")
                expt_datas[ss[0]] = float(ss[1])
            for lig,data in datas.items():
                if lig in expt_datas:
                    data["dg"] = expt_datas[lig]
                    data["ue"] = abs(data["cc_dg"] - expt_datas[lig])
                else:
                    data["dg"] = "nan"
                    data["ue"] = -1
            attrs = ["ligand","dg","cc_dg","ue","dg_error","bfe_dg","bfe_dg_error","hfe_dg","hfe_dg_error"]
        else:
            attrs = ["ligand","cc_dg","dg_error","bfe_dg","bfe_dg_error","hfe_dg","hfe_dg_error"]
        
        with open(f"{output_dir}/abfe_result.csv",'w') as outf:
            outf.write(",".join(attrs))
            outf.write("\n")
            for lig,data in datas.items():
                outf.write(f"{lig}, ")
                for attr in attrs[1:]:
                    outf.write("%.3f, " %data[attr])
                outf.write("\n")
        WriteExcel.write_abfe_excel(datas,output_dir,pic_dir=job_dir)
        return datas

    def get_complex(self,job_dir,args=None):
        job_dir,output_dir = self._parse_args(job_dir,args,attrs=["output_dir"])
        group_idx = self.generate_group(job_dir,args=None)
        if args is None:
            args = {}

        args.update({"output_dir":output_dir,
                        "group_idx":group_idx,
                        "pbc_groups":["LIG","Protein"],
                        "rmsd_group":["Protein","LIG"],
                        "energy_attrs":["Potential","Temperature","Pressure"]}
        )
        self.trj_pbc(job_dir,args=args)
        self.get_rmsd({"job_dir":job_dir,"output_dir":output_dir})
        #self.get_energy(job_dir,args=args)
        self.get_energy({"job_dir":job_dir,"output_dir":output_dir})

    def get_protein(self,job_dir,args=None):
        job_dir,output_dir = self._parse_args(job_dir,args,attrs=["output_dir"])
        group_idx = self.generate_group(job_dir,args=None)
        if args is None:
            args = {}

        args.update({"output_dir":output_dir,
                        "group_idx":group_idx,
                        "pbc_groups":["Protein"],
                        "rmsd_group":["Protein"],
                        "energy_attrs":["Potential","Temperature","Pressure"]}
        )
        self.trj_pbc(job_dir,args=args)
        self.get_rmsd({"job_dir":job_dir,"output_dir":output_dir,"attrs":["Protein"]})
        self.get_energy({"job_dir":job_dir,"output_dir":output_dir})

    def get_normal(self,job_dir,args=None):
        job_dir,output_dir = self._parse_args(job_dir,args,attrs=["output_dir"])
        if args is None:
            args = {}

        args.update({"output_dir":output_dir,"energy_attrs":["Potential","Temperature","Pressure"]})
        self.get_energy(job_dir,args=args)

    def run(self,ptype,args):
        __FUNC = {
            "trj_pbc":self.trj_pbc,
            "trj-pbc":self.trj_pbc,
            "bar": self.get_bar,
            "ddg": self.get_ddg,
            "energy": self.get_energy,
            "rmsd": self.get_rmsd,
            "torsion": self.get_lig_torsion_statistical,
            "interaction": self.get_pl_interaction,
            "pair-interaction": self.get_pair_pl_interaction,
            "fep-exchange": self.get_exchange,
            "block-ddg":self.get_block_ddg,
            "accum-ddg":self.get_accumulate_ddg,
            "pair_interaction": self.get_pair_pl_interaction,
            "fep_exchange": self.get_exchange,
            "block_ddg":self.get_block_ddg,
            "accum_ddg":self.get_accumulate_ddg,
            "rbfe": self.get_rbfe,
            "abfe": self.get_abfe,
            "complex":self.get_complex,
            "protein": self.get_protein,
            "normal": self.get_normal,
            "all-rbfe": self.get_all_rbfe,
            "all-abfe": self.get_all_abfe,
            "all_rbfe": self.get_all_rbfe,
            "all_abfe": self.get_all_abfe,
            "dimer": self.get_dimer,
            "all_dimer":self.get_all_dimer,
            "all-dimer":self.get_all_dimer,
            "total_dimer": self.get_total_dimer,
            "total-dimer":self.get_total_dimer,
        }
        
        return __FUNC[ptype](args)

class CycleClosureFreeEnergy:
    def __init__(self,ddg_df,expt_file=None,pka_file=None):
        self.ddg_df = ddg_df
        self.expt_file = expt_file
        self.pka_file = pka_file
    @staticmethod
    def _get_exp_dg(file):
        df = pd.DataFrame(columns=["name", "dg"])
        if Path(file).suffix == ".csv":
            df = pd.read_csv(file, converters={"name": str})
            df.columns = ["name", "dg"]
            df = df.astype({"name": str})
        elif Path(file).suffix == ".gpickle":
            exp_dict: Dict[str, float] = {}
            graph = Graph.load(file)
            for node, data in graph.nodes_iter(data=True):
                if dg := data.get("dg"):
                    exp_dict[node.name] = float(dg)
            df = pd.DataFrame(list(exp_dict.items()), columns=["name", "dg"])
            df = df.astype({"name": str})
        if df is None:
            raise RuntimeError(f"{file} cannot extract dg")
        return df

    @staticmethod
    def _add_exp_ddg_column(ddg_df: pd.DataFrame, dg_df: pd.DataFrame):
        dg = dict(zip(dg_df["name"], dg_df["dg"]))
        exp_ddg = []
        for _, row in ddg_df.iterrows():
            name = row["name"]
            n1, n2 = str(name).split("_to_")
            if n1 in dg and n2 in dg:
                exp_ddg.append(dg[n2] - dg[n1])
            else:
                exp_ddg.append(None)
        ddg_df["exp_ddg"] = exp_ddg
        ddg_df["ue"] = abs(ddg_df["exp_ddg"] - ddg_df["ddg"])

    def cycleclosure(self):

        dg_df = pd.DataFrame(columns=["name", "dg"])
        if self.exp_file:
            dg_df = self._get_exp_dg(self._exp_file)
            self._add_exp_ddg_column(ddg_df, dg_df)
        ddg_df, dg_df = self.run_cycle_closure(self.ddg_df, dg_df, pka_file=self._pka_file)

        if "ue" in ddg_df.columns:
            ddg_df = ddg_df.sort_values("ue")
        with pd.ExcelWriter(Path(self.output_dir) / "all_pair_result.xlsx") as writer:
            workbook = writer.book
            format = workbook.add_format()
            format.set_font_size(40)
            ddg_df.to_excel(writer, sheet_name="ddg")
            dg_df.to_excel(writer, sheet_name="dg")

            if self._exp_file:
                self.write_statistical_result(dg_df, writer)

            # add pictures to ddg_df
            if not self._two_stages:
                self.add_rdkit_pic_to_sheet(writer, ddg_df, dg_df)
        if (p := Path("correlation_figure.png")).exists():
            p.unlink()

        # for ui / output json
        basic_result = {}
        basic_result["edge"] = ddg_df.to_dict(orient="index")
        basic_result["node"] = dg_df.to_dict(orient="index")
        with open(Path(self.output_dir) / "analyze_basic.json", "w") as f:
            simplejson.dump(basic_result, f, ignore_nan=True)

    @staticmethod
    def run_cycle_closure(ddg_df, dg_df, pka_file=None):
        cc_graph = CycleClosure.ccgraph_from_csv_or_dataframe(ddg_df, dg_df)
        cc_graph.run()

        if pka_file:
            cc_graph.pka_correction(pka_file)

        cc_ddg, cc_ddg_error = {}, {}
        for edge, data in cc_graph.graph.edges_iter(data=True):
            cc_ddg[edge.name] = data.get("cc_ddg")
            cc_ddg_error[edge.name] = data.get("cc_ddg_error")

        ddg_df = ddg_df.merge(pd.DataFrame(list(cc_ddg.items()), columns=["name", "cc_ddg"]), on="name").merge(
            pd.DataFrame(list(cc_ddg_error.items()), columns=["name", "cc_ddg_error"])
        )

        cc_dg, cc_dg_error = {}, {}
        for node, data in cc_graph.graph.nodes_iter(data=True):
            cc_dg[node.name] = data.get("cc_dg")
            cc_dg_error[node.name] = data.get("cc_dg_error")

        dg_df = dg_df.merge(pd.DataFrame(list(cc_dg.items()), columns=["name", "cc_dg"]), on="name", how="right").merge(
            pd.DataFrame(list(cc_dg_error.items()), columns=["name", "cc_dg_error"]), on="name"
        )
        ddg_df.set_index("name", inplace=True)
        dg_df.set_index("name", inplace=True)
        return ddg_df, dg_df

    def get_rbfe_ddg(self,job_dir=None,args=None):
        all_pair_bar_result = defaultdict(list)
        dds = [dd for dd in os.listdir(self.job_dir) if dd.find("_to_") != -1]
        dds = [f"{self.output_dir}/{dd}" for dd in dds]
        for dd in dds:
            all_pair_bar_result["name"].append(Path(dd).name)
            with open(f"{dd}/rbfe/bar.csv") as inf:
                lines = inf.readlines()
                ss = lines.split(",")
                all_pair_bar_result["rbfe"].append(round(float(ss[1]),3))
                all_pair_bar_result["rbfe_error"].append(round(float(ss[3]),4))
            with open(f"{dd}/rhfe/bar.csv") as inf:
                lines = inf.readlines()
                ss = lines.split(",")
                all_pair_bar_result["rhfe"].append(round(float(ss[1]),3))
                all_pair_bar_result["rhfe_error"].append(round(float(ss[3]),4))
        df = pd.DataFrame(all_pair_bar_result)
        df["ddg"] = df["rbfe"] - df["rhfe"]
        df["ddg_error"] = (df["rbfe_error"] ** 2 + df["rhfe_error"] ** 2) ** 0.5

        dg_df = pd.DataFrame(columns=["name", "dg"])
        if self.expt_file:
            dg_df = self._get_exp_dg(self.expt_file)
            self._add_exp_ddg_column(ddg_df, dg_df)
        ddg_df, dg_df = self.run_cycle_closure(ddg_df, dg_df, pka_file=self.pka_file)

class WriteExcel:
    def __init__(self):
        pass

    @staticmethod
    def write_rbfe_excel(ddg_df,dg_df,output_dir,exp_file=None,two_stages=False,pic_dir=None):
        if pic_dir is None:
            pic_dir = Path(output_dir.parent,"intermedital","png")
        if "ue" in ddg_df.columns:
            ddg_df = ddg_df.sort_values("ue", ascending=False)
        with pd.ExcelWriter(Path(output_dir) / "all_pair_result.xlsx") as writer:
            workbook = writer.book
            format = workbook.add_format()
            format.set_font_size(40)
            ddg_df.to_excel(writer, sheet_name="ddg")
            dg_df.to_excel(writer, sheet_name="dg")

            if exp_file:
                WriteExcel.write_statistical_result(dg_df, writer)

            # add pictures to ddg_df
            if not two_stages:
                WriteExcel.add_pic_to_sheet(writer, ddg_df, dg_df,pic_dir)
        if (p := Path("correlation_figure.png")).exists():
            p.unlink()

        # for ui / output json
        basic_result = {}
        basic_result["edge"] = ddg_df.to_dict(orient="index")
        basic_result["node"] = dg_df.to_dict(orient="index")
        with open(Path(output_dir) / "analyze_basic.json", "w") as f:
            simplejson.dump(basic_result, f, ignore_nan=True)

    @staticmethod
    def write_abfe_excel(dg,output_dir,pic_dir=None):
        pf = pd.DataFrame.from_dict(dg,orient='index')
        if "ue" in pf.columns:
            _dg = {lig_name:dd for lig_name,dd in dg.items() if dd["ue"] != "nan"}
            _pf = pd.DataFrame.from_dict(_dg,orient='index')
            _pf = _pf.sort_values("ue",ascending=False)
        if pic_dir is None:
            pic_dir = Path(output_dir.parent,"intermedital","png")
        
        
        with pd.ExcelWriter(Path(output_dir) / "abfe_all_result.xlsx") as writer:
            workbook = writer.book
            format = workbook.add_format()
            format.set_font_size(40)
            pf.to_excel(writer, sheet_name="dg")

            if "ue" in pf.columns:
                WriteExcel.write_statistical_result(_pf, writer)

            # add pictures to ddg_df
            WriteExcel.add_pic_to_sheet_abfe(writer, pf, pic_dir)
        if (p := Path("correlation_figure.png")).exists():
            p.unlink()


    @staticmethod
    def write_statistical_result(dg_df, writer):
        from itertools import combinations

        cc_dg = dg_df["cc_dg"].replace({np.nan: None}).to_dict()
        dg = dg_df["dg"].replace({np.nan: None}).to_dict()

        cnt, rmse, less_than_1, between_1_2, greater_than_2, opposite = 0, 0, 0, 0, 0, 0
        for n1, n2 in combinations(cc_dg.keys(), 2):
            if cc_dg.get(n1) is None or cc_dg.get(n2) is None or dg.get(n1) is None or dg.get(n2) is None:
                continue
            cnt += 1
            calc_dg = cc_dg.get(n1) - cc_dg.get(n2)
            exp_dg = dg.get(n1) - dg.get(n2)
            if calc_dg * exp_dg < 0:
                opposite += 1
            ue = abs(calc_dg - exp_dg)
            if ue < 1:
                less_than_1 += 1
            elif 1 <= ue <= 2:
                between_1_2 += 1
            elif ue > 2:
                greater_than_2 += 1
            rmse += ue**2

        if cnt == 0:
            return

        rmse = np.sqrt(rmse / cnt)
        less_than_1 /= cnt
        between_1_2 /= cnt
        greater_than_2 /= cnt
        opposite /= cnt
        correlation = dg_df[["cc_dg", "dg"]].corr().iloc[0, 1]

        stat_res = {
            "all_rmse": rmse,
            "<1": less_than_1,
            "[1,2]": between_1_2,
            ">2": greater_than_2,
            "opposite": opposite,
            "correlation": correlation,
        }

        workbook = writer.book
        worksheet = workbook.add_worksheet("statistical")

        for col, (key, value) in enumerate(stat_res.items()):
            worksheet.write(0, col, key)
            worksheet.write(1, col, f"{value:.3f}")

        ax = plot_correlation(dg_df["cc_dg"], dg_df["dg"])
        ax.figure.savefig("correlation_figure.png")
        worksheet.insert_image("A3", "correlation_figure.png", {"x_offset": 2})

    @staticmethod
    def add_pic_to_sheet(writer, ddg_df, dg_df,pic_dir):
        n = len(ddg_df)
        worksheet = writer.sheets["ddg"]
        worksheet.set_column_pixels(0, 0, 350 * 2)
        mapper = {}
        for i, edge_name in enumerate(ddg_df.index):
            this_pic_dir = f"{pic_dir}/{edge_name}/bfe/job_info"
            ligand1, ligand2 = edge_name.split("_to_")
            ligand1_picture = Path(this_pic_dir)  / f"{ligand1}_atom_mapping_nonH.png"
            ligand2_picture = Path(this_pic_dir) /  f"{ligand2}_atom_mapping_nonH.png"
            ligand1_raw_picture = Path(this_pic_dir)  / f"{ligand1}.png"
            ligand2_raw_picture = Path(this_pic_dir)  / f"{ligand2}.png"
            
            worksheet.set_row_pixels(i + 1, 350)
            worksheet.insert_image(i + 1, 0, ligand1_picture, {"x_offset": 0, "y_offset": 25})
            worksheet.insert_image(i + 1, 0, ligand2_picture, {"x_offset": 350, "y_offset": 25})
            mapper[ligand1] = ligand1_raw_picture
            mapper[ligand2] = ligand2_raw_picture

        worksheet = writer.sheets["dg"]
        worksheet.set_column_pixels(0, 0, 350)
        for i, name in enumerate(dg_df.index):
            worksheet.set_row_pixels(i + 1, 350)
            worksheet.insert_image(i + 1, 0, mapper[name], {"x_offset": 0, "y_offset": 25})

    @staticmethod
    def add_pic_to_sheet_abfe(writer, dg_df,pic_dir):
        n = len(dg_df)

        worksheet = writer.sheets["dg"]
        worksheet.set_column_pixels(0, 0, 350)
        for i, name in enumerate(dg_df.index):
            ligand_picture = f"{pic_dir}/{name}/bfe/job_info/{name}.png"
            worksheet.set_row_pixels(i + 1, 350)
            worksheet.insert_image(i + 1, 0, ligand_picture, {"x_offset": 0, "y_offset": 25})
