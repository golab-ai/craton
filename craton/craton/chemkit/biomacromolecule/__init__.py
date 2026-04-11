from copy import deepcopy
from .protein_prepare import ProteinPrepare
from .protein_process import ProteinProcess
from .protein_utils import protein_atom_mapping,create_mutation
from .protein_interaction import ProteinIntraInteraction
from .register_amino_acid import assign_atom_name_for_amino_acid, create_amino_acid_template, register_amino_acid

class Protein:
    def __init__(self) -> None:
        pass

    @staticmethod
    def _pdb_prepare(protein):
        PP = ProteinPrepare(protein)
        PP.run()
        return PP.protein

    @staticmethod
    def _pdb_run_process(protein,arg):
        ProteinP = ProteinProcess(protein)
        n_terminal = None
        c_terminal = None
        terminal = None
        if len(arg) >= 4:
            n_terminal = arg[2]
            c_terminal = arg[3]
            if len(arg) == 5:
                terminal = arg[4]
        if arg[1] == "mutation":
            protein2 = ProteinP.mutation_residue(arg[0],n_terminal=n_terminal,c_terminal=c_terminal,terminal=terminal)
        elif arg[1] == "modify":
            protein2 = ProteinP.modify_residue(arg[0],n_terminal=n_terminal,c_terminal=c_terminal,terminal=terminal)
        elif arg[1] == "delete":
            protein2 = ProteinP.delete_residue(arg[0],n_terminal=n_terminal,c_terminal=c_terminal,terminal=terminal)
        elif arg[1] == "add":
            protein2 = ProteinP.add_residue(arg[0],n_terminal=n_terminal,c_terminal=c_terminal,terminal=terminal)
        return protein2

    @staticmethod
    def _pdb_process(protein,args):
        this_protein = protein
        for arg in args:
            this_protein = Protein._pdb_run_process(this_protein,arg)
        return this_protein


    @staticmethod
    def _pdb_atom_mapping(protein1,protein2):
        return protein_atom_mapping(protein1,protein2)
    
    @staticmethod
    def _sequence_create_mutation(protein,sequences):
        return create_mutation(protein,sequences)
    
    @staticmethod
    def _assign_AA_atom_name(molecules):
        new_molecules = []
        for molecule in molecules:
            new_molecules.append(assign_atom_name_for_amino_acid(molecule))
        return new_molecules
    
    @staticmethod
    def _create_AA_template(molecules):
        return create_amino_acid_template(molecules)
    
    @staticmethod
    def _register_non_AA(molecules):
        register_amino_acid(molecules)