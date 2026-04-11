from .stero_change import calculate_structure, change_structure
from .update_topolgy import template_molecule_topolgy, expand_conformers

class MolEdit:
    def __init__(self) -> None:
        pass

    @staticmethod
    def _structure_calculate(molecule,patoms):
        return calculate_structure(molecule,patoms)
    
    @staticmethod
    def _structure_change(molecule,patoms,value,del_value=False,improper_flag=False):
        return change_structure(molecule,patoms,value,del_value=del_value,improper_flag=improper_flag)

    @staticmethod
    def _molecule_topolgy_update(molecules1,molecules2,match_key=None):
        template_molecule_topolgy(molecules1,molecules2,match_key=match_key)

    @staticmethod
    def _conformer_expand(molecules1,molecules2,attrs=["coordinates"],match_key="inchi_key"):
        return expand_conformers(molecules1,molecules2,attrs=attrs,match_key=match_key)