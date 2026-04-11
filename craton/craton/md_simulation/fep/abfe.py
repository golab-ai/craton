from itertools import chain
from copy import deepcopy
# from ...chemkit import MolScalpel as MS
from ...chemkit import Structure as Stru
#return Stru._protein_structure(protein)
from collections import defaultdict

from ...chemkit import InteractionModel as IM #get_interaction_model

from ...utils.geometry import calc_stru_para

def get_key_interaction(interaction_models):
    def get_atom(model):
        if model.adj_atom is not None:
            return model.adj_atom.ID
        if model.atoms is not None:
            return model.atoms[0].ID
        return model.atom.ID

    model_list = ["saltbridge_lneg","saltbridge_pneg","metal_complexation",
                  "hbonds_ldon","hbonds_pdon","pi_stacking","pication_paro",
                  "pication_laro","chpi_paro","chpi_laro","halogen_bonds",
                  "hydrophobic_contacts",]
    interaction_atoms = []
    for attr in model_list:
        vvs = getattr(interaction_models,attr,[])
        for vv in vvs :
        #vv = len(interaction_models[attr])
        #if len(vv) > 0:
            acc = vv.acceptor
            acc_atom = get_atom(vv.acceptor)
            don_atom = get_atom(vv.donor)
            
            if acc.group == "LIG":
                interaction_atoms.append([don_atom,acc_atom])
                #return don_atom, acc_atom
            else:
                interaction_atoms.append([acc_atom,don_atom])
                #return acc_atom, don_atom
    return interaction_atoms
        
def get_ligand_three_atom(ligand_atom,molecule):
    mm = deepcopy(molecule)
    mm = MS._atom_cluster([mm])[0]
    #if len(mm.elem_frag) > 2:
    for ii,frag in mm.elem_frag.items():
        if ligand_atom in frag["components"]:
            break
    flag = False
    for jj in frag["connects"]:
        kks = [fn for fn in mm.elem_frag[jj]["connects"] if fn != ii]
        if len(kks) > 0:
            kk = kks[0]
            flag = True
    if flag:
        atoms = [ligand_atom]
        for an in mm.elem_frag[jj]["components"]:
            if mm.Atoms[an].atom_cluster_tag.find("center") != -1:
                atoms.append(an)
                break
        for an in mm.elem_frag[kk]["components"]:
            if mm.Atoms[an].atom_cluster_tag.find("center") != -1:
                atoms.append(an)
                break
        return atoms
    
    else:
        atoms = [ligand_atom]
        jj = [an for an in mm.Atoms[ligand_atom] if mm.Atoms[an].elem not in ["H","F","Cl","Br","I"]][0]
        kk = [an for an in mm.Atoms[jj] if an != ligand_atom and mm.Atoms[an].elem not in ["H","F","Cl","Br","I"]][0]
        return [ligand_atom,jj,kk]

def get_intermolecular_restrain(system):
    mole = deepcopy(system.molecules[0])
    mole_n = len(mole.Atoms)
    
    protein = deepcopy(system.molecules[1])
    #for atom in protein.Atoms:
    #    atom.ID = atom.ID + mole_n
    protein = Stru._protein_structure(protein)

    interaction_models = IM.get_interaction_model(protein,probe=mole)

    interaction_atoms = get_key_interaction(interaction_models)

    if len(interaction_atoms) < 3:
        pass
    else:
        P = [interaction_atoms[ii][0] + mole_n for ii in range(3)]
        L = [interaction_atoms[ii][1] for ii in range(3)]


    P_coors = [system.coordinates[an] for an in P]
    L_coors = [system.coordinates[an] for an in L]
    intermolecular_interaction = defaultdict(list)
    intermolecular_interaction["bonds"] = [
        [L[0], P[0], calc_stru_para([L_coors[0], P_coors[0]])]
    ]
    intermolecular_interaction["angles"] = [
        [
            P[1],
            P[0],
            L[0],
            calc_stru_para([P_coors[1], P_coors[0], L_coors[0]]),
        ],
        [
            P[0],
            L[0],
            L[1],
            calc_stru_para([P_coors[0], L_coors[0], L_coors[1]]),
        ],
    ]
    intermolecular_interaction["dihedrals"] = [
        [
            P[2],
            P[1],
            P[0],
            L[0],
            calc_stru_para([P_coors[2], P_coors[1], P_coors[0], L_coors[0]]),
        ],
        [
            L[1],
            L[0],
            P[0],
            P[1],
            calc_stru_para([L_coors[1], L_coors[0], P_coors[0], P_coors[1]]),
        ],
        [
            L[2],
            L[1],
            L[0],
            P[0],
            calc_stru_para([L_coors[2], L_coors[1], L_coors[0], P_coors[0]]),
        ],
    ]
    return intermolecular_interaction
