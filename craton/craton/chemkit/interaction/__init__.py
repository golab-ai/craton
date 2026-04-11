from .interaction_model import InteractionSite
from .interaction_detection import interaction_detection
from ...utils.commons import parallel_run


class InteractionModel:
    def __init__(self) -> None:
        pass
    
    @staticmethod
    def run_get_interaction_model(coor,molecule=None,probe=None,probe_flag=True,idx=None):
        if coor is not None:
            shift_an = molecule.shift_an
            for atom in molecule.Atoms:
                atom.coor = coor[atom.ID + shift_an]
            shift_an = probe.shift_an
            for atom in probe.Atoms:
                atom.coor = coor[atom.ID + shift_an]

        IS = InteractionSite(molecule,probe)
        if probe_flag:
            binding_residues,atoms = IS.find_binding_residue(molecule,probe)
        else:
            atoms = None
        protein_sites = IS.get_interaction_site(molecule,atoms=atoms)
        ligand_sites = IS.get_interaction_site(probe)
        inter_model = interaction_detection(ligand_sites,protein_sites)
        if idx is not None:
            return inter_model, idx
        else:
            return inter_model

    @staticmethod
    def get_interaction_model(molecule,probe,probe_flag=True, coordinates=None,parallel=True):
        if coordinates is None:
            return InteractionModel.run_get_interaction_model(None,molecule=molecule,probe=probe,probe_flag=probe_flag)
        else:
            if parallel:
                #args = [molecule for ii in range(len(coordinates))]
                kwds = [{"probe":probe,"molecule":molecule,"probe_flag":probe_flag} for coor in coordinates]
                total_inter_model = parallel_run(InteractionModel.run_get_interaction_model,coordinates,kwds=kwds,)
            else:
                total_inter_model = []
                for coor in coordinates:
                    total_inter_model.append(InteractionModel.run_get_interaction_model(coor,molecule=molecule,probe=probe,probe_flag=probe_flag))
            return total_inter_model


