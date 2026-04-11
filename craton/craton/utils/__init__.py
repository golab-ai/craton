from .common._logger import _get_logger
logger = _get_logger("craton.log", log_directory=".")

from .figure import DrawFigure as DF
from .pubchem_query import run
from .fetch_pdb import pdb
from .fetch_uniprot import uniprot


def _show_figure(XX,figure_type,args=None):
    __func = {
                "diagonal": DF.diagonal_draw,
                "pes":      DF.pes1d_draw,
                "violin":   DF.violin,
                "bar":      DF.bar,
                "pie":      DF.pie,
                "line":     DF.line_draw,
              }
    if args is not None:
        file_name = __func[figure_type](XX,**args)
    else:
        file_name = __func[figure_type](XX)
    return file_name

def _pubchem_info(strs,typ,print_flag=True,opath="."):
    if not isinstance(strs,list):
        strs = [strs]
    infos = []
    for str in strs:
        infos.append(run(str,typ,print_flag=print_flag,opath=opath))
    return infos

def _get_pdb_file(info_file=None,pdb_id=None,output_dir=".",output_format="pdb"):
    return pdb(info_file=info_file,pdb_id=pdb_id,output_dir=output_dir,output_format=output_format)

def _get_uniport(target=None,uniprot_id=None,output_directory="."):
   return uniprot(target=target,uniprot_id=uniprot_id,output_directory=output_directory)