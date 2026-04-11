import textwrap
from pathlib import Path

import click
import yaml
#import simplejson

from craton.craton.utils import logger
from craton import molxpert as MX
from craton.application.chem_info import chem_info, molecule_descript
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])



@click.group(name="tool",context_settings=CONTEXT_SETTINGS)
def tool():
    """Analyzing the fep molecule_dynamics molecule_dynamics
    The analyzing sub-group containing following commands:

    \b
    * md, including bar, mbar, exchange, rmsd, torsion, interaction et al
    * rbfe, all pair directorys
    * all_extra, all pair directorys
    """

@tool.command("save")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def save_info(inputs,output_dir):
    """
    create descript of a molecule, and save as a file

    \b
    inputs: molecules
    output_dir: the path for save file
    """

    molecule_descript(inputs,output_dir)


@tool.command("torsion")
@click.option("-i","--input_directory",help="The directory of gaussian log files",default="./",show_default=True,)
@click.option("-o", "--output", help="The directory of output", default="./", show_default=True)
#@click.option("-exp", "--exp_file", help="Experiemnt dG file, this could be csv file or gpickle file")
#@click.option("--two_stages/--no-two_stages", default=False, show_default=True, help="If two stages are used?")
#@click.option("-pka", "--pka-file", help="pka file", default=None, show_default=True)
def torsion_scan_analyze(input_directory, output):
    """
    torsion  result analyze base on gaussian 16 output files:
    "-i": The directory of gaussian log files
    "-o": The directory of output
    """
    i_path = input_directory
    o_path = output
    molecules = MX.molecule_create(i_path,parallel=False)
    smiles_arr = [m.smiles for m in molecules]
    model_molecules = MX.molecule_create(smiles_arr)
    model_molecules = MX.molecule_structure(model_molecules)
    model_molecules = MX.molecule_torsion(model_molecules)

    MX.molecule_topolgy_update(molecules,model_molecules)

    scan_curves = MX.scan_curve(molecules)
    for inchi_key,torsions in scan_curves.items():
            for torsion_name,mm in torsions.items():
                    mm = sorted(mm,key=lambda m:m.constrain[0].fix_value)
                    fname = f"{inchi_key}_{torsion_name}"
                    MX.format_convert(mm,otype="sdf",ofilename=fname,opath=o_path)
                    atoms = [int(an) for an in torsion_name.split("-")]
                    MX.molecule_show(mm[0],attrs=["torsion"],opath=o_path,extra=atoms)


    datas = MX.scan_curve_data(scan_curves)
    for name,XX in datas.items():
        args = {"name":f"pes_{name}","xylabels":["trosion","energy(kcal/mol)"],"save_path":o_path}
        MX.figure_show(XX,"pes",args=args)


@tool.command("chem-info")
@click.option("-i","--input_file",help="The directory of gaussian log files",default="./",show_default=True,)
#@click.option("-o", "--output_file", help="The directory of output", default="./", show_default=True)
def get_chem_info(input_file):
    chem_info(input_file)
    

@tool.command("pubchem")
@click.option("-i","--input",help="input similes or name or file",default=None,show_default=True,)
@click.option("-it","--input_type",help="the type of input",default="smiles",show_default=True,)
@click.option("-ot","--output_type",help="the type of output",default="name",show_default=True,)
#@click.option("-o", "--output_file", help="The directory of output", default="./", show_default=True)
def get_pubchem_info(input,input_type,output_type):
    """
    type: smiles, inchikey,name, nickname,formula,cas_no,weight,
    """
    if Path(input).exists():
        results = MX.get_pubchem(input,input_type=input_type,output_type=output_type,file_flag=True)
    else:
        results = MX.get_pubchem(input,input_type=input_type,output_type=output_type)
        print(results[0],results[1])
    
@tool.command("cg")
@click.option("-i","--input_molecules",help="The directory of gaussian log files",default="./",show_default=True,)
@click.option("-o", "--output", help="The directory of output", default="./", show_default=True)
def croase_grain(input_molecules,output):
    molecules = MX.molecule_create(input_molecules)
    molecules = MX.molecule_structure(molecules)
    molecules = MX.atom_cluster(molecules)
    cg_molecules = MX.cg_molecule(molecules)

    MX.format_convert(cg_molecules,otype="mtx",opath=output,extra_var="all")

@tool.command("property")
@click.option("-i","--inputs",help="The directory of gaussian log files",default="514-10-3",show_default=True,)
@click.option("-p","--props",help="The properties",default="density",show_default=True)
@click.option("-t","--temperature",help="The temperatures",default=300.0,show_default=True)
@click.option("-f", "--yaml_file", type=click.File("r"), help="if need more fined control, please specify a yaml file")
#@click.option("-o", "--output", help="The directory of output", default="./", show_default=True)
def thermo_dyna_property(inputs,props,temperature,yaml_file):
    molecules = inputs
    property = props
    temperatures = temperature
    config = {"molecules":inputs,"property":props,"temperatures":temperature,"molecule_type":"CAS_number"}
    if yaml_file is not None:
        y_config = yaml.safe_load(yaml_file.read())
        for attr in ["molecules","property","temperatures","moleucle_type"]:
            if attr in y_config:
                config[attr] = y_config[attr]
    res = MX.get_property(config["molecules"],config["property"],molecule_type=config["molecule_type"],temperatures=config["temperatures"])
    print(res)
    
@tool.command("peptide")
@click.option("-n","--residue_number",help="the residue number for generation",default=1,show_default=True,)
@click.option("-r","--c_terminal",help="the C terminal cap",default=None,show_default=True,)
@click.option("-l","--n_terminal",help="the N terminal cap",default=None,show_default=True,)
@click.option("-t","--template",help="the template",default=None,show_default=True,)
@click.option("-o", "--output", help="The directory of output", default="./", show_default=True)
####n,left_cap="ACE",right_cap="NME",terminal_flag=True,templates=None
def generate_peptide(residue_number,c_terminal,n_terminal,template,output):
    """
    * c_terminal: none or ACE-OH0
    * n_terminal: none or NME-H00
    * template: none or 0:ALA-GLY_1:ASP-HIS-SER_3:VAL-GLU
    """
    if c_terminal == "none" or c_terminal is None:
        c_terminal = "NME"
    else:
        c_terminal = [rr if rr != "none" else None for rr in c_terminal.split("-")]
    if n_terminal == "none" or n_terminal is None:
        n_terminal = "ACE"
    else:
        n_terminal = [rr if rr != "none" else None for rr in n_terminal.split("-")]
    if template == "none" or template is None:
        template = None
    else:
        template=template.split("-")
        #template = {}
        #for rr in sss:
        #    ss = rr.split(":")
        #    template[int(ss[0])] = ss[1].split("-")
    arrs = MX.create_peptide(residue_number,left_cap=n_terminal,right_cap=c_terminal,templates=template)

@tool.command("rna")
@click.option("-n","--residue_number",help="the residue number for generation",default=1,show_default=True,)
@click.option("-t","--template",help="the template",default="rna",show_default=True,)
@click.option("-o", "--output", help="The directory of output", default="./", show_default=True)
####n,left_cap="ACE",right_cap="NME",terminal_flag=True,templates=None
def generate_peptide(residue_number,template,output):
    """
    * template: none or 0:ALA-GLY_1:ASP-HIS-SER_3:VAL-GLU
    """
    if template == "none" or template is None:
        template = None
    elif template in ["rna","dna"]:
        template = template
    else:
        template=template.split("-")
        #template = {}
        #for rr in sss:
        #    ss = rr.split(":")
        #    template[int(ss[0])] = ss[1].split("-")
    arrs = MX.create_dnarna(residue_number,templates=template)
    print(arrs)

@tool.command("dna")
@click.option("-n","--residue_number",help="the residue number for generation",default=1,show_default=True,)
@click.option("-t","--template",help="the template",default="dna",show_default=True,)
@click.option("-o", "--output", help="The directory of output", default="./", show_default=True)
####n,left_cap="ACE",right_cap="NME",terminal_flag=True,templates=None
def generate_peptide(residue_number,template,output):
    """
    * template: none or 0:ALA-GLY_1:ASP-HIS-SER_3:VAL-GLU
    """
    if template == "none" or template is None:
        template = None
    elif template in ["rna","dna"]:
        template = template
    else:
        template=template.split("-")
        #template = {}
        #for rr in sss:
        #    ss = rr.split(":")
        #    template[int(ss[0])] = ss[1].split("-")
    arrs = MX.create_dnarna(residue_number,templates=template)
    print(arrs)

@tool.command("atf2json")
@click.option("-f","--file_name",help="the atom type define txt file name",show_default=True,)
def atf_to_json(file_name):
    MX._convert_atom_type_file_to_json_(file_name)

@tool.command("fff2json")
@click.option("-f","--file_name",help="the force field parameter txt file name",show_default=True,)
def fff_to_json(file_name):
    MX._convert_force_field_file_to_json_(file_name)

@tool.command("aminoacid")
@click.option("-f","--file_name",help="the rtp file name",show_default=True,)
def convert_aminoacid_template(file_name):
    MX._aminoacid_json(file_name)
    
@tool.command("amberff")
@click.option("-a","--atf",help="the atom type define file",show_default=True,)
@click.option("-v","--nonbf",help="the nonb force field file",show_default=True,)
@click.option("-b","--bondf",help="the bonded force field file",show_default=True,)
def convert_amber_ff(atf,nonbf,bondf):
    MX._amberff_to_ff(atf,nonbf,bondf)