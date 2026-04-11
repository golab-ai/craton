import textwrap
from pathlib import Path

import click
import yaml
#import simplejson

from craton.craton.utils import logger
from craton import molxpert as MX
from craton.application.mm_calculation import MMCalculation as MMC
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(name="mm",context_settings=CONTEXT_SETTINGS)
def mm():
    """calculation based on force field

    \b
    * calculate, calculate energy, force, hessian, frequency et al
    * opt, optimize structrue
    * scan, torsion scan
    """
@mm.command("calculate")
@click.argument(
    "prop",
    type=click.Choice(
        ["energy","force","hessian","freq","frequency"]
    ),
)
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def calculate_energy(prop,inputs,output_dir):
    """
    calculate energy, force, hessian and frequency

    \b
    prop: energy, force, hessian,freq
    inputs: molecules
    output_dir: the path for save file
    """
    __label = {
        "energy":"energy",
        "force":"force",
        "hessian":"hessian",
        "freq": "freq",
        "frequency":"frequency"
    }
    mm_calculator = MMC(inputs,prop=__label[prop],output_directory=output_dir)
    mm_calculator.energy_calculate()

@mm.command("opt")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def save_info(inputs,output_dir):
    """
    create descript of a molecule, and save as a file

    \b
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,output_directory=output_dir)
    mm_calculator.molecule_optimize()

@mm.command("scan")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def save_info(inputs,output_dir):
    """
    create descript of a molecule, and save as a file

    \b
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,output_directory=output_dir)
    mm_calculator.torsion_scan()

@mm.command("multipole")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-p","--prop",help="dipole, quadrupole, octupole, multipole",default="energy",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def multipole_info(prop,inputs,output_dir):
    """
    calculate multipole, such as dipole, quadrupole, octupole

    \b
    prop: dipole, quadrupole, octupole
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,prop=prop,output_directory=output_dir)
    mm_calculator.multipole_moment()

@mm.command("volume")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def volume_info(inputs,output_dir):
    """
    calculate multipole, such as dipole, quadrupole, octupole

    \b
    prop: dipole, quadrupole, octupole
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,output_directory=output_dir)
    mm_calculator.volume_surface()

@mm.command("surface")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def surface_info(inputs,output_dir):
    """
    calculate multipole, such as dipole, quadrupole, octupole

    \b
    prop: dipole, quadrupole, octupole
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,output_directory=output_dir)
    mm_calculator.volume_surface()

@mm.command("center")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-p","--prop",help="center, cog, com, cob, size",default="cog",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def center_info(inputs,prop,output_dir):
    """
    calculate multipole, such as dipole, quadrupole, octupole

    \b
    prop: dipole, quadrupole, octupole
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,prop=prop,output_directory=output_dir)
    mm_calculator.center()

@mm.command("inertia")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def inertia_info(inputs,output_dir):
    """
    calculate multipole, such as dipole, quadrupole, octupole

    \b
    prop: dipole, quadrupole, octupole
    inputs: molecules
    output_dir: the path for save file
    """
    mm_calculator = MMC(inputs,output_directory=output_dir)
    mm_calculator.inertia_moment()