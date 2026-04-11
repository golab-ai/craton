import textwrap
from pathlib import Path

import click
import yaml
#import simplejson

from craton.craton.utils import logger
from craton import molxpert as MX
from craton.application.chem_info import ChemInfo as CI
from craton.application.chem_info import MoleEdit as ME
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(name="stru",context_settings=CONTEXT_SETTINGS)
def stru():
    """Analyzing topol info of moleucle, calculate structure parameter and vary molecule structure:

    \b
    * topol: analyze topol structure of molecules
    * measure: measure the structure parameter of molecule
    * vary: bond strecching, angle blending and dihderal rotating
    """

@stru.command("topol")
@click.argument(
    "topol_typ",
    type=click.Choice(
        ["ring","chiral","torsion","hybrid","interaction_site","frag","cg","fg","image"]
    ),
)
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def topol_analyze(topol_typ,inputs,output_dir):
    """
    analyze topol structure of molecules

    \b
    topol_typ: ring, chiral, torsion, hybrid, interaction_site, frag(fragmentation), cg(corase grain bead), fg(function group), image
    inputs: molecules
    output_dir: the path for save file
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    chemi = CI(inputs,output_directory=output_dir)

    _label = {
        "ring": "ring",
        "chiral": "chiral_atom",
        "torsion": "torsion",
        "hybrid": "hybrid",
        "interaction_site": "interaction_site",
        "frag": "fragmentation",
        "cg": "atom_cluster",
        "fg": "function_group",
        "image": "molecule_image",
    }
    action_name = _label[topol_typ]
    action = getattr(chemi, action_name, None)
    if action is None:
        raise click.ClickException(f"ChemInfo has no action: {action_name}")
    result = action() if callable(action) else action
    if callable(result):
        result()

    # Fallback: if topol action does not write files explicitly, dump molecule info.
    out_path = Path(output_dir)
    if not any(out_path.iterdir()) and hasattr(chemi, "molecules"):
        try:
            MX.save_molecule_info(chemi.molecules, fpath=output_dir, parallel=False)
        except Exception as e:
            logger.warning(f"topol fallback output failed: {e}")

@stru.command("measure")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-a","--atoms",help="atom ID, such as: 2-3, 4-5-6,7-8-9-10",default=None,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def stru_parameter(inputs,atoms,output_dir):
    """
    measure structure parameter of bond, angle, and dihedral

    \b
    measure_type: bond, angle or dihedral
    inputs: molecules
    atoms: atom ID
    output_dir: the path for save file
    """
    if atoms is None:
        raise click.BadParameter("option '--atoms' is required, e.g. 2-3 or 4-5-6")
    atoms = [int(an) for an in atoms.split("-")]

    mole_edit = ME(inputs,output_directory=output_dir)
    mole_edit.calculate_structure_parameter(atoms)

@stru.command("vary")
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-a","--atoms",help="atom ID, such as: 3-7,4-5-6,11-0-7-15",default=None,show_default=True,)
@click.option("-v","--value",help="vary value",default=None,show_default=True,)
@click.option("-dv","--del_value",help="是否del值",default=False,show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def change_stur_parameter(inputs,atoms,value,del_value,output_dir):
    """
    create descript of a molecule, and save as a file

    \b
    inputs: molecules
    output_dir: the path for save file
    """
    if atoms is None:
        raise click.BadParameter("option '--atoms' is required, e.g. 3-7 or 11-0-7-15")
    atoms = [int(an) for an in atoms.split("-")]

    mole_edit = ME(inputs,output_directory=output_dir)
    mole_edit.change_structure_parameter(atoms,value,del_value=del_value)

@stru.command("rmsd")
@click.option("-i","--inputs",help="reference molecule",default=".",show_default=True,)
@click.option("-t","--inputs2",help="target molecule",default=".",show_default=True,)
@click.option("-o","--output_dir",help="save path",default=".",show_default=True)
def molecule_rmsd(inputs,inputs2,output_dir):
    """
    calculate the rmsd of two molecule

    \b
    inputs: reference molecule
    inputs2: target molecule
    output_dir: the path for save file
    """

    mole_edit = ME(inputs,output_directory=output_dir)
    mole_edit.get_RMSD(inputs2)
