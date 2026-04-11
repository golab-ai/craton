import gzip
import json
import multiprocessing as mp
from pathlib import Path
from copy import deepcopy
from typing import *

import munch
import numpy as np

from ...utils import logger
#from ...mm_calculator import Calculator as Calc


#from ..conformation.conformation import ConformType
from ...chemkit.conformation.conformation import ConformType
#from ..force_field import  prepare_force_field

from .intra_pairs_analyze import IntraAnalyze as IA
#from ..force_field import _special_dihedrals_infos

def get_parameters_results(molecules, ts_molecules, this_ff):
    values = {}
    _label = {"bondterm":"Bonds","angleterm":"Angles","dihedralterm":"Dihedrals","improperterm":"Impropers"}
    for term,params in this_ff.items():
        if term in ["bondterm","angleterm","dihedralterm","improperterm",]:
            values[_label[term]] = {}
            for para, v1 in params.items():
                values[_label[term]][para] = [[],[],1] if v1["ptag"] in ["isfitting","pseudofitting"] else [[],[],0]

    _val_type = {"Bonds":"value","Angles":"value","Dihedrals":"value","Impropers":"value"}
    for ii,molecule in enumerate(molecules):
        ts_molecule = ts_molecules[ii]
        for term in ["Bonds","Angles","Dihedrals","Impropers"]:
            if hasattr(molecule,term):
                molecule_topols = getattr(molecule,term)
                ts_molecule_topols = getattr(ts_molecule,term)
                for jj, topol in enumerate(molecule_topols):
                    ts_topol = ts_molecule_topols[jj]
                    topol_name = topol.atom_type_used_name
                    values[term][topol_name][0].append(round(getattr(ts_topol,_val_type[term]),4))
                    values[term][topol_name][1].append(round(getattr(topol,_val_type[term]),4))

    return values

def get_molecule_info(molecules):
    datas = {
        "total": {
            "pes_number": [
                0,
                0,
            ],
            "conformations": 0,
            "miniminum_number": 0,
            "molecules_number": 0,
        },
    }
    inchi_key_dict = {}
    for mol in molecules:
        if mol.inchi_key not in inchi_key_dict:
            datas["total"]["molecules_number"] += 1
            datas[mol.inchi_key] = {"conformations": 0, "pes_number": [0, 0], "miniminum_number": 0}
            inchi_key_dict[mol.inchi_key] = []
        datas[mol.inchi_key]["conformations"] += 1
        datas["total"]["conformations"] += 1
        if hasattr(mol, "constrain"):
            #scan_term_name = "-".join([str(n) for n in mol.scan_term[0]])
            scan_term_name = mol.constrain[0].name
            if scan_term_name not in inchi_key_dict[mol.inchi_key]:
                datas[mol.inchi_key]["pes_number"][0] += 1
                datas[mol.inchi_key]["pes_number"][1] += 1
                datas["total"]["pes_number"][0] += 1
                datas["total"]["pes_number"][1] += 1
                inchi_key_dict[mol.inchi_key].append(scan_term_name)
            else:
                datas[mol.inchi_key]["pes_number"][1] += 1
                datas["total"]["pes_number"][1] += 1
        else:
            datas[mol.inchi_key]["miniminum_number"] += 1
            datas["total"]["miniminum_number"] += 1
    return datas

def get_parameter_info(this_ff) :
    __tt = {"atomtype", "bondterm", "angleterm", "dihedralterm", "improperterm"}
    datas = dict()
    fix_datas = dict()
    for key, value in this_ff.items():
        if key in __tt:
            for term, para in value.items():
                if "isfitting" in para:
                    if key not in datas:
                        datas[key] = {}
                    datas[key][term] = [round(v, 4) for v in para["parameter"]]
                else:
                    if key not in fix_datas:
                        fix_datas[key] = {}
                    fix_datas[key][term] = [round(v, 4) for v in para["parameter"]]
    datas["fix_datas"] = fix_datas
    return datas

def qm_mm_analyze(
        ts_molecules,
        results_path,
        force_field=None,
        atom_type_file=None,
        optimizer="openmm",
        done_fitting=None,
        init_this_ff=None,
        validation_terms=None,
        optimize_flag=True
    ):
    from ...mm_calculator import Calculator as Calc
    if done_fitting is None:
        ts_molecules,this_ff = prepare_force_field(ts_molecules,atom_type_file=atom_type_file,force_field_file=force_field)

    else:
        this_ff = force_field
        if init_this_ff is None:
            init_this_ff = this_ff
    #ts_molecules,this_ff = prepare_force_field(ts_molecules,atom_type_file=atom_type_file,force_field_file=force_field)
    
    #if init_this_ff is None:
    #    init_this_ff = this_ff

    if validation_terms is None:
        validation_terms = ["energy", "pes", "esp_charge", "Bonds", "Angles", "Dihedrals", "rmsd", "Pair1n","hessian","freq","force"]
    if not optimize_flag:
        validation_terms = list(set(validation_terms)-set(["Bonds","Angles","Dihedrals","rmsd","Pair1n"]))

    molecule_info = get_molecule_info(ts_molecules)
    
    parameter_info = get_parameter_info(init_this_ff)

    ts_molecules_opt = [molecule for molecule in ts_molecules if molecule.conform_type in ConformType.TORSION_SCAN_TYPES + [ConformType.LOCAL_MINIMUM]]
    ts_molecules_stretch = [molecule for molecule in ts_molecules if molecule.conform_type == ConformType.STRETCH]
    molecules_opt = [deepcopy(mol) for mol in ts_molecules_opt]
    molecules_stretch = [deepcopy(mol) for mol in ts_molecules_stretch]

    #### validation force
    force_flag = False
    if "force" in validation_terms:
        molecules_force = Calc._mix_energy(molecules_opt+molecules_stretch,terms=["force"])
        data_f,stat_f = IA.intra_molecule_pairs(molecules_force,ts_molecules_opt+ts_molecules_stretch,terms=["force"])
        validation_terms = list(set(validation_terms) - set(["force"]))
        force_flag = True

    #data_f, stat_f, _ = pairs_analyze(molecules_opt+molecules_stretch,ts_molecules_opt+ts_molecules_stretch,terms=["force"],)

    ### optimize the comformers of opt comformers
    if optimize_flag:
        logger.info(f"{len(molecules_opt)} conformers optimizing in intra_validation ......")
        molecules_opt = Calc._optimize(molecules_opt, optimizer=optimizer)
        logger.info(f"{len(molecules_opt)} conformers optimizing in intra_validation DONE")

        #terms = ["energy", "pes", "esp_charge", "Bonds", "Angles", "Dihedrals", "rmsd", "Pair1n","hessian","freq"]

    ### validation other attributions
    logger.info(f"Calculate energy force hessian frequency ...... ")
    molecules_opt = Calc._mix_energy(molecules_opt,terms=validation_terms)
    logger.info(f"Calculate energy force hessian frequency DONE ")
    logger.info(f"Analyze data ...... ")


    data,statics = IA.intra_molecule_pairs(molecules_opt,ts_molecules_opt,terms=validation_terms)
    if "torsion_infos" not in this_ff:
        _special_dihedrals_infos(ts_molecules_opt,this_ff)
    for inchi_key,vv in data["pes"].items():
        for torsion,dd in vv.items():
            dd.append(this_ff["torsion_infos"][inchi_key][torsion])
    logger.info(f"Analyze data DONE ")
    #data, statics = pairs_analyze(molecules_opt,ts_molecules_opt,terms=terms,)

    # get parameters results
    parameter_results = get_parameters_results(molecules_opt, ts_molecules_opt, this_ff)

    if force_flag:
        data["force"] = data_f["force"]
        statics["force"] = stat_f["force"]

    if len(data["energy"]["check_conformer"]) > 0:
        from ...chem import FormatMolecule as FM
        Path(f"{results_path}/check_ts_conformer").mkdir(exist_ok=True)
        Path(f"{results_path}/check_conformer").mkdir(exist_ok=True)
        check_ts_molecule = [ts_molecules_opt[ii] for ii in data["energy"]["check_conformer"]]
        check_molecule = [molecules_opt[ii] for ii in data["energy"]["check_conformer"]]
        FM._convert(check_ts_molecule,otype="gjf",opath= f"{results_path}/check_ts_conformer",)
        FM._convert(check_molecule,otype="gjf",opath=f"{results_path}/check_conformer")

    results_datas = {"initi_info": {"init_moles": molecule_info,"init_para": parameter_info,},
                     "results": {"properties": data,"property_statics":statics,"parameters": parameter_results},
                    }
    if done_fitting is not None:
         results_datas["process"] = done_fitting,
    # output as gzip file
    with gzip.open(f"{results_path}/intra_validation_results.json.gz", "wt") as outf:
        outf.write(json.dumps(results_datas, indent=2))

    return munch.munchify(results_datas)