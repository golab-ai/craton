import textwrap
import sys
from pathlib import Path

import click
import yaml
#import simplejson

from craton.craton.utils import logger
from craton import molxpert as MX
from craton.application.mole_dock import MolDock
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(name="dock",context_settings=CONTEXT_SETTINGS)
def dock():
    """run dock

    \b
    * pocket: get pocket information of a protein
    * dock: run docking 
    """

@dock.command("pocket")
@click.option("-i","--inputs",help="protein file",default=".",show_default=True,)
@click.option("-o","--output_directory",help="output directory",default=".",show_default=True,)
def pocket_info(inputs,output_directory):
    """
    get the pocket information of protein

    \b

    """
    results = MolDock.pocket(inputs,output_directory=output_directory)

@dock.command("dock")
@click.option("-p","--protein",help="protein file",default=".",show_default=True,)
@click.option("-i","--ligands",help="ligands file",default=".",show_default=True,)
@click.option("-center","--center",help="the center coordinates of pocket",default=None,show_default=True,)
@click.option("-box","--box_size",help="the box size of pocket",default=None,show_default=True,)
@click.option("-charge","--charge_method",help="charge method of ligands",default="binc",show_default=True,)
@click.option("-o","--output_dir",help="output directory",default=".",show_default=True,)
def docking_runing(protein,ligands,center,box_size,charge_method,output_dir):
    """
    get the pocket information of protein

    \b

    """
    if center is None:
        center = None
    else:
        center = [float(x) for x in center.split(",")]
    if box_size is None:
        box_size = None
    else:
        box_size = [float(x) for x in box_size.split(",")]  
    DOCK = MolDock(protein,ligands,center=center,box_size=box_size,charge_method=charge_method,output_directory=output_dir)
    DOCK.run_vina()

