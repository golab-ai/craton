from .structure import get_2d_connectivity
from .structure import remove_hydrogen_atoms
from .ring import find_cyclo
from .ring import cyclo_property
from .ring import cyclo_blocks
from .structure import assign_structure_info
from .structure import assign_local_info
from .structure import assign_connect_type
from .structure import assign_conjugate_info
from .structure import divide_charge 
from .structure import assign_hybrid
from .structure import flexible_torsion
from .structure import scan_torsion
from .structure import get_fragment_inchi_key
from .structure import check_chiral_atom
from .chain import get_chain,search_chain
from .structure import find_mole,get_2d_connectivity,remove_hydrogen_atoms
from .structure import _update_molecule_topol_value

from .connectivity import molecule_coordinate_to_bond_type
from ...utils.commons import parallel_run

from .structure import protein_ring_and_charge_group

from .model import Model


#####
class Structure:
    def __init__(self):
        pass
    
    @staticmethod
    def _update_topol_values(molecules,parallel=True):
        if not isinstance(molecules,list):
            molecules = [molecules]
        if parallel:
            return parallel_run(_update_molecule_topol_value,molecules)
        else:
            return [_update_molecule_topol_value(molecule) for molecule in molecules]


    @staticmethod
    def _find_molecule(connects):
        return find_mole(connects)

    @staticmethod
    def _single_basic_structure_analyze(molecule,ignore_existing=False,idx=None):
        from ..biomacromolecule import Protein #_pdb_prepare
        if not ignore_existing and "structure" in molecule.steps:
            pass
        else:
            if molecule.style not in ["pdb","protein","template","dna","rna","DNA","RNA","Protein"]:
                element_number, all_connect = get_2d_connectivity(molecule)
                reduce_connect, __, ___ = remove_hydrogen_atoms(element_number, all_connect)
                cyclos = find_cyclo(reduce_connect)
                cyclos = cyclo_property(molecule, cyclos)
                rings,ring_blocks,ring_block_components = cyclo_blocks(cyclos)
                assign_structure_info(molecule,rings,ring_blocks,ring_block_components)
                assign_local_info(molecule,reduce_connect)
                assign_connect_type(molecule)
                assign_conjugate_info(molecule)
                if len(molecule.Atoms) > 1:
                    divide_charge(molecule)
                molecule.steps.append("structure")
            else:
                molecule = Protein._pdb_prepare(molecule)
                molecule = Structure._protein_structure(molecule)
                molecule.steps.append("atom type")
                molecule.steps.append("force field")
        if idx is not None:
            return molecule,idx
        else:
            return molecule

    @staticmethod
    def _basic_structure_analyze(molecules,ignore_existing=False,parallel=True):
        new_molecules = []
        if not isinstance(molecules,list):
            molecules=[molecules]
        if parallel:
            new_molecules = parallel_run(Structure._single_basic_structure_analyze,molecules, kwds=[{"ignore_existing":ignore_existing} for molecule in molecules])
        else:
            for molecule in molecules:
                new_molecules.append(Structure._single_basic_structure_analyze(molecule,ignore_existing=ignore_existing))
        return new_molecules

    @staticmethod
    def _assign_hybrid(molecules,parallel=True):
        if not isinstance(molecules,list):
            molecules=[molecules]

        if parallel:
            molecules = parallel_run(assign_hybrid,molecules)
        else:
            for molecule in molecules:
                assign_hybrid(molecule)
        for molecule in molecules:
            molecule.steps.append("hybrid")
        return molecules

    @staticmethod
    def _protein_structure(protein):
        return protein_ring_and_charge_group(protein)

    @staticmethod
    def _get_flexible_torsion(molecules,parallel=True):
        if not isinstance(molecules,list):
            molecules=[molecules]

        if parallel:
            molecules = parallel_run(flexible_torsion,molecules)
            molecules = parallel_run(scan_torsion,molecules)
        else:
            for molecule in molecules:
                flexible_torsion(molecule)
                scan_torsion(molecule)
        for molecule in molecules:
            molecule.steps.append("torsion")
        return molecules

    @staticmethod
    def _get_chain(molecules,parallel=True):
        if not isinstance(molecules,list):
            molecules=[molecules]

        if parallel:
            tmp_ = parallel_run(get_chain,molecules)
            datas = [rr for rr in tmp_]
        else:
            datas = []
            for molecule in molecules:
                datas.append(get_chain(molecule))
        return datas
    
    @staticmethod
    def _search_chain(terminal_atom_arr,connect_dict,additional_conditions=None):
        
        return search_chain(terminal_atom_arr,connect_dict,additional_conditions=additional_conditions)

    @staticmethod
    def _get_connectivity_bond_type(molecule):

        ####TODO
        molecule_coordinate_to_bond_type(molecule)

    @staticmethod
    def _run_model_atom(molecule,idx=None):
        interaction_model = Model(molecule).run()
        molecule._interaction_model = interaction_model
        for kk,vv in interaction_model.items():
            for term in vv:
                if term is not None:
                    if term.atom is not None:
                        IDs = [term.atom.ID]
                    else:
                        IDs = [atom.ID for atom in term.atoms]
                    for ID in IDs:
                        if not hasattr(molecule.Atoms[ID],"_interaction_model"):
                            molecule.Atoms[ID]._interaction_model = [term.type]
                        else:
                            molecule.Atoms[ID]._interaction_model.append(term.type)
        if idx is not None:
            return molecule,idx
        else:
            return molecule            
            
    @staticmethod
    def _model_atom(molecules,parallel=False):
        if not isinstance(molecules,list):
            molecules=[molecules]

        if parallel:
            new_molecules = parallel_run(Structure._run_model_atom,molecules)
        else:
            new_molecules = []
            for molecule in molecules:
                new_molecules.append(Structure._run_model_atom(molecule))
        return new_molecules
    
    @staticmethod
    def _get_chiral_atom(molecules,parallel=True):
        if not isinstance(molecules,list):
            molecules = [molecules]
        if parallel:
            new_moleucles = parallel_run(check_chiral_atom,molecules)
        else:
            new_moleucles = []
            for molecule in molecules:
                new_moleucles.append(check_chiral_atom(molecule))
        return new_moleucles
        