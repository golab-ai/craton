from copy import deepcopy
import json
from itertools import chain
from pathlib import Path

from ..chem import FormatMolecule as FM
from ..chem import _show_molecule

from ..utils.commons import parallel_run

from .gmxmd import GmxEngine, GmxSlurm
from .lmpmd import LmpEngine, LmpSlurm


def run_write_system(system):
    
    Path(system.output_dir).mkdir(parents=True, exist_ok=True)
    if system.mdengine in ["gmx","gromacs"]:
        GmxEngine(system)
    elif system.mdengine in ["lammps","lmp"]:
        LmpEngine(system)


def write_system(systems,parallel=True):
    if parallel:
        parallel_run(run_write_system,systems,keep_order=False,return_result=False)
    else:
        for system in systems:
            run_write_system(system)

def run_write_bash(system):
    Path(system.output_dir).mkdir(parents=True, exist_ok=True)
    if system.mdengine in ["gmx","gromacs"]:
        GmxSlurm(system)
    elif system.mdengine in ["lammps","lmp"]:
        LmpSlurm(system)

def write_bash(systems,parallel=True):
    if parallel:
        parallel_run(run_write_bash,systems,keep_order=False,return_result=False)
    else:
        for system in systems:
            run_write_bash(system)

def run_rfe_info(sm):
    parent_dir = sm.output_dir
    Path(parent_dir).mkdir(exist_ok=True)
    this_dir = f"{parent_dir}/job_info"
    Path(this_dir).mkdir(exist_ok=True)

    job_json = {"simulation_type": sm.simulation_type,
                "mdengine": sm.mdengine,
                "ligand": sm.molecules[0].mole_name,
                "protein":sm.molecules[1].mole_name,
                "left": sm.atom_mapping[0].mole_name,
                "right": sm.atom_mapping[1].mole_name,
                "total_atom_mapping": sm.atom_mapping[2],
                "atom_mapping": sm.atom_mapping[3],
                "atom_mapping_nonH": sm.atom_mapping[4],
                "md_setting": {kk:vv for kk,vv in sm.md_setting.items() if kk != "free_energy_auixed"},
                }
    with open(f"{this_dir}/md_setting.json",'w') as outf:
        outf.write(json.dumps(job_json))
    for molecule in sm.molecules:
        FM._convert(molecule,otype="mtx",ofilename=molecule.mole_name,opath=this_dir,extra_var="all")
    for ii,molecule in enumerate(sm.atom_mapping[:2]):
        mole_name = "left" if ii == 0 else "right"
        highlights = [nn for nn in sm.atom_mapping[3].keys()] if ii == 0 else [nn for nn in sm.atom_mapping[3].values()]
        highlights_nonH = [nn for nn in sm.atom_mapping[4].keys()] if ii == 0 else [nn for nn in sm.atom_mapping[4].values()]
        FM._convert(molecule,otype="mtx",ofilename=mole_name,opath=this_dir,extra_var="all")
        FM._convert(molecule,otype="mol",ofilename=mole_name,opath=this_dir)
        FM._convert(molecule,otype="png",ofilename=mole_name,opath=this_dir)
        _show_molecule(molecule,attrs="normal",extra={"highlights":highlights,"type":"atoms","fname_pre":"atom_mapping"},save_file=True,opath=this_dir,TD_flag=False,remove_H_flag=False)
        _show_molecule(molecule,attrs="normal",extra={"highlights":highlights_nonH,"type":"atoms","fname_pre":"atom_mapping_nonH"},save_file=True,opath=this_dir,TD_flag=False)

def run_normal_info(sm):
    ST = sm.simulation_type
    parent_dir = sm.output_dir
    Path(parent_dir).mkdir(exist_ok=True)
    this_dir = f"{parent_dir}/job_info"
    Path(this_dir).mkdir(exist_ok=True)

    molecule_smiles = []
    molecule_inchi_key = []


    for kk,molecule in enumerate(sm.molecules):
        FM._convert(molecule,otype="mtx",ofilename=molecule.mole_name,opath=this_dir,extra_var="all")
        if ST in ["rbfe","rhfe"] and kk == 0:
            molecule_smiles.append(molecule.mole_name)
            molecule_inchi_key.append(molecule.mole_name)
        else:
            if molecule.style not in ["pdb","protein","template","dna","rna","DNA","RNA","Protein"]:
                if len(molecule.Atoms) <= 200:
                    FM._convert(molecule,otype="mol",ofilename=molecule.mole_name,opath=this_dir)
                    FM._convert(molecule,otype="png",ofilename=molecule.mole_name,opath=this_dir)
                    molecule_smiles.append(molecule.smiles)
                    molecule_inchi_key.append(molecule.inchi_key)
                else:
                    molecule_smiles.append(molecule.mole_name)
                    molecule_inchi_key.append(molecule.mole_name)
            else:
                molecule_smiles.append(molecule.mole_name)
                molecule_inchi_key.append(molecule.mole_name)
    if ST in ["rbfe","rhfe"]:
        for ii,molecule in enumerate(sm.atom_mapping[:2]):
            mole_name = "left" if ii == 0 else "right"
            highlights = [nn for nn in sm.atom_mapping[3].keys()] if ii == 0 else [nn for nn in sm.atom_mapping[3].values()]
            highlights_nonH = [nn for nn in sm.atom_mapping[4].keys()] if ii == 0 else [nn for nn in sm.atom_mapping[4].values()]
            FM._convert(molecule,otype="mtx",ofilename=mole_name,opath=this_dir,extra_var="all")
            FM._convert(molecule,otype="mol",ofilename=mole_name,opath=this_dir)
            FM._convert(molecule,otype="png",ofilename=mole_name,opath=this_dir)
            _show_molecule(molecule,attrs="normal",extra={"highlights":highlights,"type":"atoms","fname_pre":"atom_mapping"},save_file=True,opath=this_dir,TD_flag=False,remove_H_flag=False)
            _show_molecule(molecule,attrs="normal",extra={"highlights":highlights_nonH,"type":"atoms","fname_pre":"atom_mapping_nonH"},save_file=True,opath=this_dir,TD_flag=False)


        
    md_setting = {kk:vv for kk,vv in sm.md_setting.items() if kk != "free_energy_auixed"}
    md_info_json = {
                    "simulation_type":sm.simulation_type,
                    "molecule_name":[molecule.mole_name for molecule in sm.molecules],
                    "molecule_number":sm.molecule_number,
                    "molecule_role":[molecule.molecule_role for molecule in sm.molecules],
                    "molecule_smiles":molecule_smiles,
                    "molecule_inchi_key":molecule_inchi_key,
                    "mdengine":sm.mdengine,
                    "output_dir":sm.output_dir,
                    "info_dir":this_dir,
                    "md_setting":md_setting,
                    "name":sm.name,
                    }
    if ST in ["rbfe","rhfe"]:
        md_info_json["left"] = sm.atom_mapping[0].mole_name
        md_info_json["right"] = sm.atom_mapping[1].mole_name
        md_info_json["left_smiles"] = sm.atom_mapping[0].smiles
        md_info_json["right_smiles"] = sm.atom_mapping[1].smiles
        md_info_json["left_inchi_key"] = sm.atom_mapping[0].inchi_key
        md_info_json["right_inchi_key"] = sm.atom_mapping[1].inchi_key
        md_info_json["total_atom_mapping"] = sm.atom_mapping[2]
        md_info_json["atom_mapping"] = sm.atom_mapping[3]
        md_info_json["atom_mapping_nonH"] = sm.atom_mapping[4]
    with open(f"{this_dir}/md_setting.json",'w') as outf:
        outf.write(json.dumps(md_info_json))
    
def run_info(sm):
    #if sm.simulation_type in ["rbfe","rhfe"]:
    #    run_rfe_info(sm)
    #else:
    run_normal_info(sm)
            
def write_job_info(systems,parallel=True):
    _systems = [sm for sm in systems]
    if parallel:
        parallel_run(run_info,_systems,keep_order=False,return_result=False)
    else:
        for system in _systems:
            run_info(system)