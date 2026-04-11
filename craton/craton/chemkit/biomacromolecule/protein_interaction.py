from pathlib import Path

from ...chem.molecule import Molecule

from ..structure import Structure as Stru
from ..interaction.interaction_model import InteractionSite
from ..interaction.interaction_detection import interaction_detection
#from ...chem import FormatMolecule

from .protein_prepare import ProteinPrepare

class ProteinIntraInteraction:
    def __init__(self,output_dir) -> None:
        self.path = output_dir
    
    def add_terminal(self,tmp_mole,tmp_atom,pre):
        atoms = []
        for ii,atom in enumerate(tmp_mole.Atoms):
            atoms.append(atom) 
            atoms[-1].residue_name = tmp_atom.residue
            atoms[-1].residue_ID = tmp_atom.residue_ID
            atoms[-1].chain_name = tmp_atom.chain_name
            atoms[-1].ID = f"{pre}{atoms[-1].ID}"
            atoms[-1].connectivity = [f"{pre}{an}" for an in atoms[-1].connectivity]
        return atoms

    def get_molecule(self,protein,group,group_name):
        pr_molecule = Molecule("molecule")
        pr_molecule.mole_name = group_name
        pr_molecule.Atoms = []
        for an in group.atoms:
            pr_molecule.Atoms.append(protein.Atoms[an])
        
        
        PP = ProteinPrepare(pr_molecule)
        PP.run()
        pr_molecule = PP.protein
        
        
        return pr_molecule
            
    def run_get_protein_intra_interaction(self,protein,protein_sites,group,residue_dict,ignore_atoms):
        from ...chem import FormatMolecule
        __label ={
                "saltbridge_lneg":["acceptor","atoms","donor","atoms"],
                "saltbridge_pneg":["donor","atoms","acceptor","atoms"],
                "hbonds_ldon":["acceptor","atom","donor","adj_atom"],
                "hbonds_pdon":["donor","adj_atom","acceptor","atom"],
                "pi_stacking":["acceptor","atoms","donor","atoms"],
                "pication_laro":["donor","atoms","acceptor","atoms"],
                "pication_paro":["acceptor","atoms","donor","atoms"],
                #"weakhbond":[],
                "chpi_laro":["donor","adj_atom","acceptor","atoms"],
                "chpi_paro":["acceptor","atoms","donor","adj_atom"],
                "hydrophobic_contacts":["acceptor","atom","donor","atom"],
                "halogen_bonds":["acceptor","atom","donor","atom"],
                #"waterbridges":[],
            }
    
        IS = InteractionSite(None,None)
        #protein_sites = IS.get_interaction_site(protein)
        probe_residue = Molecule("residue")
        probe_residue.Atoms = [protein.Atoms[an] for an in group.atoms]
        res_na = f"{group.group_name}-{group.group_idx}-{group.group_chain_name}"
    
        ligands_sites = {}
        for aa,bb in protein_sites.items():
            ligands_sites[aa] = []
            for bbb in bb:
                ss = f"{bbb.group}-{bbb.group_id}-{bbb.chain}"
                if ss == res_na:
                    ligands_sites[aa].append(bbb)
    
        br_10,probe_atoms = IS.find_binding_residue(protein,probe_residue)
        probe_atoms = list(set(probe_atoms).difference(ignore_atoms))
        
        protein_sites = IS.get_interaction_site(protein,atoms=probe_atoms)
        
        inter_model = interaction_detection(ligands_sites,protein_sites)
        
        
        action_rsidue = {}
        for attr in __label:
            interactions = getattr(inter_model,attr)
            if len(interactions) > 0:
                for interaction in interactions:
                    pr = getattr(interaction,__label[attr][0])
                    lr = getattr(interaction,__label[attr][2])
                    pr_name = f"{pr.group}-{pr.group_id}-{pr.chain}"
                    if f"{res_na}${pr_name}" not in self.interaction_residues and f"{pr_name}${res_na}" not in self.interaction_residues:
                        self.interaction_residues.append(f"{res_na}${pr_name}")
                        #if pr_name not in action_rsidue:
                        action_rsidue[pr_name] = [residue_dict[pr_name],[[attr,getattr(pr,__label[attr][1]),getattr(lr,__label[attr][3])]]]
                    else:
                        if pr_name in action_rsidue:
                            action_rsidue[pr_name].append([attr,getattr(pr,__label[attr][1]),getattr(lr,__label[attr][3])])
        if action_rsidue:
            pre_path = f"{self.path}/{res_na}"
            Path(pre_path).mkdir(exist_ok=True)
            this_path = f"{pre_path}/dimer"
            Path(this_path).mkdir(exist_ok=True)
            
            res_molecule = self.get_molecule(protein,group,res_na)
            res_atom_name = {ai.atom_name:ai.ID for ai in res_molecule.Atoms}
            FormatMolecule._convert(res_molecule,otype="mtx",opath=pre_path,ofilename=res_na,extra_var="all")
            FormatMolecule._convert(res_molecule,otype="sdf",opath=pre_path,ofilename=res_na,extra_var="all")

            for res,action in action_rsidue.items():
                group_name = f"{res}-{action[1][0][0]}"
                inter_molecule = self.get_molecule(protein,action[0],group_name)
                inter_atom_name = {ai.atom_name:ai.ID for ai in inter_molecule.Atoms}
                
                inter_molecule.associated_data = {"interaction_type":action[1][0][0]}
                res_interaction_atoms = action[1][0][2] if isinstance(action[1][0][2],list) else [action[1][0][2]] 
                inter_interaction_atoms = action[1][0][1] if isinstance(action[1][0][1],list) else [action[1][0][1]]
                
                inter_molecule.associated_data["residue1"] = [res_atom_name[an.atom_name] for an in res_interaction_atoms]
                inter_molecule.associated_data["residue2"] = [inter_atom_name[an.atom_name] for an in inter_interaction_atoms]
                
                FormatMolecule._convert(inter_molecule,otype="mtx",opath=this_path,ofilename=group_name,extra_var="all")
                FormatMolecule._convert(inter_molecule,otype="sdf",opath=this_path,ofilename=group_name,extra_var="all")
            
            
    def get_protein_intra_interaction(self,protein,parallel=False):
        protein = Stru._protein_structure(protein)
        residue_dict = {f"{group.group_name}-{group.group_idx}-{group.group_chain_name}":group for group in protein.Groups}
        IS = InteractionSite(None,None)
    
        protein_sites = IS.get_interaction_site(protein)
    
        nn = len(protein.Groups)
        self.interaction_residues = []
        
        for kk, group in enumerate(protein.Groups):
            ignore_atoms = []
            for jj in range(-2,3):
                if kk+jj >=0 and kk+jj < nn:
                    ignore_atoms.extend(protein.Groups[kk+jj].atoms)
            self.run_get_protein_intra_interaction(protein,protein_sites,group,residue_dict,set(ignore_atoms))
                            
