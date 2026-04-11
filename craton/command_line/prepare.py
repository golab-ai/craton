import textwrap
import sys
from pathlib import Path

import click
import yaml
#import simplejson

from craton.craton.utils import logger
from craton import molxpert as MX
from craton.application.mole_prapare import MolPrepare as MP
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(name="prepare",context_settings=CONTEXT_SETTINGS)
def prepare():
    """prepare the biomacromolecule and small molecule

    \b
    * mol_info: get the information of moleucles including iupac name, cas-no, smiles et al
    * ligand_prepare: small molecule prepare, including ionization et al
    * uniprot: get uniprot information
    * pdb: get pdb file from https://www.rcsb.org/
    * protein_prepare: prepare pdb file, including ionization, add terminal of N and C, and add missing hydrogen atom et al
    * mutation: residue mutation of protein
    """

@prepare.command("mol_info")
@click.option("-i","--inputs",help="the molecule for get infomation",default=".",show_default=True,)
@click.option("-it","--in_type",help="the type of inputs",default="smiles",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def molecule_info(inputs,in_type,output_dir):
    """
    get the information of moleucles including iupac name, cas-no, smiles et al

    \b
    inputs: inputs moleucle 
    in_type: type of input, such as name, smiles, cas-no
    output_dir: output directory
    """
    MP.info(inputs,in_type,print_flag=True,output_directory=output_dir)

@prepare.command("ligand")
@click.option("-i","--inputs",help="input files",default=None,show_default=True,)
@click.option("-pi","--ph_min",help="minize value of ph",default=7.4,show_default=True,)
@click.option("-pa","--ph_max",help="maxmum value of ph",default=7.4,show_default=True,)
@click.option("-of","--output_file",help="output pdb file",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def ligand_prepare(inputs,ph_min,ph_max,output_file,output_dir):
    """
    small molecule prepare, including ionization et al
    """
    MP.ligand_prepare(inputs,ph_min=ph_min,ph_max=ph_max,ofilename=output_file,output_directory=output_dir)

@prepare.command("uniprot")
@click.option("-t","--target",help="the name of target protein",default=None,show_default=True,)
@click.option("-i","--uniprot_id",help="uniprot id",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def uniprot_info(target,uniprot_id,output_dir):
    """
    uniprot information including sequence, pdb file, and fasta file

    \b
    target: target protein 
    uniprot_id: uniprot id
    output_dir: output directory
    """
    MP.uniprot(target=target,uniprot_id=uniprot_id,output_directory=output_dir)

@prepare.command("pdb")
@click.option("-i","--info_file",help="the file including pdb ids",default=None,show_default=True,)
@click.option("-id","--pdb_id",help="pdb id",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def pdb_info(info_file,pdb_id,output_dir):
    """
    get pdb file from https://www.rcsb.org/
    """
    MP.get_pdb(info_file=info_file,pdb_id=pdb_id,output_dir=output_dir)

@prepare.command("protein")
@click.option("-i","--protein_file",help="input dpb ifle",default=None,show_default=True,)
@click.option("-of","--output_file",help="output pdb file",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def protein_prepare(protein_file,output_file,output_dir):
    """
    prepare pdb file, including ionization, add terminal of N and C, and add missing hydrogen atom et al
    """
    MP.protein_prepare(protein_file,ofilename=output_file,opath=output_dir)

@prepare.command("mutation")
@click.option("-i","--protein_file",help="input dpb file",default=None,show_default=True,)
@click.option("-r","--residue",help="target residue,such as 127-ARG-H",show_default=True,)
@click.option("-m","--mutation",help="mutation residue,such as LEU",show_default=True,)
@click.option("-of","--output_file",help="output pdb file",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def protein_mutation(protein_file,residue,mutation,output_file,output_dir):
    """
    residue mutation of protein
    """
    MP.protein_mutation(protein_file,residue,mutation,ofilename=output_file,opath=output_dir)

@prepare.command("modify")
@click.option("-i","--protein_file",help="input dpb file",default=None,show_default=True,)
@click.option("-r","--residue",help="target residue,such as 127-ARG-H",show_default=True,)
@click.option("-m","--modify_fg",help="modify type,such as pho",show_default=True,)
@click.option("-of","--output_file",help="output pdb file",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def protein_mutation(protein_file,residue,modify_fg,output_file,output_dir):
    """
    residue mutation of protein
    """
    if modify_fg not in ["pho","suf","met","n-met"]:
        print("the modify target set is error, the modify type must be in pho suf, met, n-met, now is {_tmp}")
        sys.exit(1)

    MP.protein_modify(protein_file,residue,modify_fg,ofilename=output_file,opath=output_dir)
