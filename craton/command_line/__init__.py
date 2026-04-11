import click

#from .fep_analyze import fep_analyze
#from .normal_analyze import normal_analyze
from .simulation import run_simulaiton
from .fep_analyze import analyze
from .cratondb import data
from .force_field import ff
from .tool import tool
from .mm import mm
from .stru import stru
from .prepare import prepare
from .dock import dock
##from ..tools import tools_app

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def craton_app():
    """Simulation application, it contains the 'molecule_dynamics', 'analyze' and 'extend' sub-groups."""


craton_app.add_command(run_simulaiton)
# craton_app.add_command(md_check)
craton_app.add_command(analyze)
craton_app.add_command(data)
craton_app.add_command(ff)
craton_app.add_command(dock)
craton_app.add_command(mm)
craton_app.add_command(stru)
craton_app.add_command(tool)
craton_app.add_command(prepare)

###simulation_app.add_command(tools_app)
#simulation_app.add_command(normal_analyze)
# simulation_app.add_command(extend)
