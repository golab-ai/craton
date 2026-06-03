import click
import yaml
from copy import deepcopy

from craton.application.md_simulation import Simulation_run
from craton.application.qm_conformation_search import QMConformationSearch as QMC
from craton.application.qm_conformation_search import meta_qm_calc

MD_SIMULATION_TYPES = [
    "vacuum","solution","liquid","complex","protein",
    "rbfe","rhfe","rlogp","rlogs","mem-rbfe","cov-rbfe",
    "mutation","pep-rbfe","rna-rbfe",
    "abfe","ahfe","hfe","alogp",
    "bilayer","biomembrane","mem-protein"
    ]
QM_SIMULATION_TYPES = ["Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q100",]

@click.command("simulation")
@click.argument(
    "simulation_type",
    type=click.Choice(
        ["rbfe", "abfe", "rhfe","ahfe","hfe",
         "vacuum", "solution", "complex", "protein", "alogp","rlogp","alogs","rlogs", "liquid",
        "Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q100",
        "yaml","autoqm"]
    ),
)
@click.option("--ligands", "ligands", help="ligands file, support the format .sdf")
@click.option("--protein", "protein", help="protein file, support the format .pdb")
@click.option("--molecules", "molecules", help="molecules, support the format .mol, smiles, ......")
@click.option("--coligands", "coligands", help="co-ligands file, suport the format .sdf")
@click.option("--repeat","repeat", default=1,help="the repeat number of calculation job")
@click.option("-o", "--output", "output_directory", default="output", help="output directory")
@click.option("-t", "--simulation_time", default=None, help="Simulation steps (2fs per step)", show_default=True)
@click.option("-c", "--charge_method", default=None, help="charge method")
@click.option("-n", "--molecule_number", default=None, help="the molecule number of auto-qmcalc")
@click.option("-p", "--property", default="normal", help="property, such as density")
@click.option("-eng", "--mdengine", default="gmx", help="md engine, such as gmx, lmp")
@click.option("--ncpu", "ncpu", default=None, type=int, help="number of CPUs (EnvironmentSetting.ncpu, e.g. SLURM cpus-per-task)")
@click.option("-f", "--yaml_file", type=click.File("r"), help="if need more fined control, please specify a yaml file")
def run_simulaiton(
    simulation_type, ligands, protein, molecules, coligands, repeat, output_directory, simulation_time, charge_method, molecule_number,property,mdengine,ncpu,yaml_file
):
    """Preparing Gromacs MD molecule_dynamics files, most of setting are defined in the 'default_settings.py' file.
    If need more fined control, you can specify a yaml file, the template of which could
    be found in the 'tests' directory. For example "inf_md molecule_dynamics yaml -f test.yaml'

     \b
    Different simulations should provide different files, as showing by the following rules:
        * pure-liquid, smiles or smiles_file
        * vacuum, smiles or smiles_file
        * solution, smiles or smiles_file
        * biomolecule, protein_file
        * complex,  protein_file, ligand_file
        * rbfe(relative binding free energy), protein_file, ligand_file
        * hfe(hydration free energy), smiles or smiles_file
        * abfe(absolute binding free neergy), protein_file, ligand_file
        * yaml, provide yaml file for more settings.
    """
    total_config = {}
    if simulation_type == "autoqm":
        meta_qm_calc(nn=molecule_number)
        return
    elif simulation_type == "yaml":
        if yaml_file is None:
            raise RuntimeError("When choose yaml, the yaml file must be supported followed by -f ")
        config = yaml.safe_load(yaml_file.read())
        multitask = False
        for aa,bb in config.items():
            if aa.find("task") != -1:
                multitask = True
                break

        if multitask:
            _config = deepcopy(config)

        else:
            _config = {"task0":config}

    else:
        _config = {"task0": {
            "simulation_type": simulation_type,
            "ligands": ligands,
            "protein": protein,
            "molecules": molecules,
            "repeat":repeat,
            "simulation_time":simulation_time,
            "output_directory": output_directory,
            "coligands": coligands,
            "property": property.split(":"),
            "charge_method": charge_method,
            "mdengine": mdengine,
        }}

    total_config = create_repeat_config(_config)
    if ncpu is not None:
        for bb in total_config.values():
            bb["ncpu"] = ncpu

    for aa,bb in total_config.items():
        if bb["simulation_type"] in MD_SIMULATION_TYPES:
            run_md_simulaiton(bb)
        elif bb["simulation_type"] in QM_SIMULATION_TYPES:
            bb["stage"] = bb["simulation_type"]
            run_qm_simulation(bb)
    
def run_md_simulaiton(config):
    click.echo(config)
    Simulation_run(config)

def run_qm_simulation(config):
    click.echo(config)
    QMC(config)

def create_repeat_config(config):
    total_config = {}
    for task,this_config in config.items():
        if "repeat" in this_config and this_config["repeat"] > 1:
            if "output_directory" not in this_config:
                output_directory = "./output"
            else:
                output_directory = this_config["output_directory"]
            for ii in range(this_config["repeat"]):
                _task = f"{task}_{ii}"
                _output_directory = f"{output_directory}_{ii+1}"
                _this_config = deepcopy(this_config)
                _this_config["output_directory"] = _output_directory
                total_config[_task] = _this_config
        else:
            total_config[task] = this_config 
    return total_config
