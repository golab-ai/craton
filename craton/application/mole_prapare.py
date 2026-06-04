from craton import molxpert as MX
import csv,os,sys
import json
from pathlib import Path
from copy import deepcopy    

class MolPrepare:
    def __init__(self):
        pass

    @staticmethod
    def uniprot(target=None,uniprot_id=None,output_directory="."):
        MX.uniport_info(target=target,uniprot_id=uniprot_id,output_directory=output_directory)

    @staticmethod   
    def info(inps,in_type,print_flag=True,output_directory="."):
        molecules = []
        if in_type not in ["name","cas-no"]:
            molecules = MX.molecule_create(inps)

        if len(molecules) > 0:
            MX.pubchem_info([molecule.smiles for molecule in molecules],"smiles",print_flag=print_flag,opath=output_directory)
        else:
            MX.pubchem_info(inps,in_type,print_flag=print_flag,opath=output_directory)

    @staticmethod
    def get_pdb(info_file=None,pdb_id=None,output_dir=".",output_format="pdb"):
        MX.pdb_file(info_file=info_file,pdb_id=pdb_id,output_dir=output_dir,output_format=output_format)

    @staticmethod
    def protein_prepare(protein_file,ofilename=None,opath="."):
        if isinstance(protein_file, str):
            protein = MX.molecule_create(protein_file)[0]
        else:
            protein = protein_file
        protein = MX.protein_prepare(protein)
        if ofilename is None:
            ofilename = protein.mole_name
        else:
            if ofilename.endswith(".pdb"):
                ofilename = ofilename[:-4]
        MX.format_convert(protein,otype="pdb",ofilename=ofilename,opath=opath)

    @staticmethod
    def protein_mutation(protein_file,residue,mutation,ofilename=None,opath="."):
        change = [[[[residue,mutation]],"mutation"]]
        if isinstance(protein_file, str):
            protein = MX.molecule_create(protein_file)[0]
        else:
            protein = protein_file

        protein = MX.protein_prepare(protein)
        from craton.craton.chemkit.structure.structure import protein_ring_and_charge_group
        protein = protein_ring_and_charge_group(protein)
        
        protein2 = MX.protein_process(protein,change)
        if ofilename is None:
            ofilename = f"{protein2.mole_name}_{residue}_{mutation}"
        else:
            if ofilename.endswith(".pdb"):
                ofilename = ofilename[:-4]
        MX.format_convert(protein2,otype="pdb",ofilename=ofilename,opath=opath)

    @staticmethod
    def protein_modify(protein_file,residue,modify_fg,ofilename=None,opath="."):
        if modify_fg not in ["pho","suf","met","n-met"]:
            print(f"the modify target set is error, the modify type must be in pho suf, met, n-met, now is {modify_fg}")
            sys.exit(1)
        change = [[[[residue,modify_fg]],"modify"]]
        if isinstance(protein_file, str):
            protein = MX.molecule_create(protein_file)[0]
        else:
            protein = protein_file
        protein = MX.protein_prepare(protein)
        from craton.craton.chemkit.structure.structure import protein_ring_and_charge_group
        protein = protein_ring_and_charge_group(protein)

        protein2 = MX.protein_process(protein,change)
        if ofilename is not None:
            if ofilename.endswith(".pdb"):
                ofilename = ofilename[:-4]
        else:
            ofilename = f"{protein2.mole_name}_{residue}_{modify_fg}"
        MX.format_convert(protein2,otype="pdb",ofilename=ofilename,opath=opath)

    @staticmethod
    def ligand_prepare(inputs,ph_min=7.4,ph_max=7.4,ofilename=None,output_directory="."):
        molecules = MX.molecule_create(inputs)
        molecules2 = MX.molecule_prepare(molecules,ph_min=ph_min,ph_max=ph_max)
        for ii,molecule in enumerate(molecules2):
            molecule.mole_name = molecules[ii].mole_name
            if molecule.inchi_key != molecules[ii].inchi_key:
                if molecules[ii].mole_name != molecules[ii].inchi_key:
                    molecule.mole_name = f"{molecules[ii].mole_name}_charged"
            molecule.associated_data = {}
            molecule.associated_data["smiles"] = molecule.smiles
            molecule.associated_data["origin_smiles"] = molecules[ii].smiles
            molecule.associated_data["inchi_key"] = molecule.inchi_key

        MX.format_convert(molecules2,otype="sdf",ofilename=ofilename,opath=output_directory)
