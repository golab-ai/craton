from copy import deepcopy
from ...chem.molecule import Molecule
from ...chem.atom import Atom
from ..structure.structure import protein_ring_and_charge_group
from .protein_utils import template_molecule, amino_acid
from .protein_prepare import ProteinPrepare

modify_groups = {"Pho":[{"element":"P","atom_name":"MP","atom_type_name":"P","connectivity":["MO1","MO2","MO3"],"bond_type":[2,1,1],"formal_charge":0,"plate":"no","ff_charge":0.328},
                        {"element":"O","atom_name":"MO1","atom_type_name":"O2","connectivity":["MP"],"bond_type":[2],"formal_charge":0,"plate":"no","ff_charge":-0.776},
                        {"element":"O","atom_name":"MO2","atom_type_name":"O2","connectivity":["MP"],"bond_type":[1],"formal_charge":-1,"plate":"no","ff_charge":-0.776},
                        {"element":"O","atom_name":"MO3","atom_type_name":"O2","connectivity":["MP"],"bond_type":[1],"formal_charge":-1,"plate":"no","ff_charge":-0.776}],
                 "MET":[{"element":"C","atom_name":"MC","atom_type_name":"CT","connectivity":["MH1","MH2","MH3"],"bond_type":[1,1,1],"formal_charge":0,"plate":"no","ff_charge":-0.3369},
                        {"element":"H","atom_name":"MH1","atom_type_name":"HP","connectivity":["MC"],"bond_type":[1],"formal_charge":0,"plate":"no","ff_charge":0.1123},
                        {"element":"H","atom_name":"MH2","atom_type_name":"HP","connectivity":["MC"],"bond_type":[1],"formal_charge":0,"plate":"no","ff_charge":0.1123},
                        {"element":"H","atom_name":"MH3","atom_type_name":"HP","connectivity":["MC"],"bond_type":[1],"formal_charge":0,"plate":"no","ff_charge":0.1123}],
                 "Suf":[{"element":"S","atom_name":"MS","atom_type_name":"S","connectivity":["MO1","MO2","MO3"],"bond_type":[2,2,1],"formal_charge":0,"plate":"no","ff_charge":1.8900},
                        {"element":"O","atom_name":"MO1","atom_type_name":"O2","connectivity":["MS"],"bond_type":[2],"formal_charge":0,"plate":"no","ff_charge":-0.6300},
                        {"element":"O","atom_name":"MO2","atom_type_name":"O2","connectivity":["MS"],"bond_type":[2],"formal_charge":0,"plate":"no","ff_charge":-0.6300},
                        {"element":"O","atom_name":"MO3","atom_type_name":"O2","connectivity":["MS"],"bond_type":[1],"formal_charge":-1,"plate":"no","ff_charge":-0.6300},],
                }

class ProteinProcess:
    def __init__(self,protein) -> None:
        self.protein = deepcopy(protein)
        self._template = {"pho":[{"SER":"OG","THR":"OG1","TYR":"OH"},"Pho"],
                          "met":[{"ARG":"NH2","LYS":"NZ"},"MET"],
                          "n-met":[{"all":"N",},"MET"],
                          "suf":[{"TYR":"OH"},"Suf"],
                          }

    

    def get_group_template(self,group):
        group_template = modify_groups[group]
        Atoms = [Atom(style="atom") for ii in range(len(group_template))]
        for ii,atom in enumerate(Atoms):
            for kk,vv in group_template[ii].items():
                setattr(atom,kk,vv)
        return Atoms
    
    def modify_residue(self,modifies,n_terminal=None,c_terminal=None,terminal=None,create_3D=True):
        def parse_modify(arr):
            tmp = []
            for rs in arr:
                rr = rs[0]
                typ = rs[1]
                ss = rr.split("_")
                if len(ss) == 3:
                    residue = ss[0]
                    residue_id = ss[1]
                    chain = ss[2]
                    try:
                        #anchor = f"{self._template[typ][0][residue]}_{rr}"
                        anchor = self._template[typ][0][residue]
                    except:
                        #anchor = f"{self._template[typ][0]['all']}_{rr}"
                        anchor = self._template[typ][0]['all']
                    cutpoint = None
                    residue_point = rr
                elif len(ss) == 4:
                    residue = ss[1]
                    residue_id = ss[2]
                    chain = ss[3]
                    anchor = ss[0]
                    cutpoint = None
                    residue_point = "_".join(ss[1:])
                elif len(ss) == 5:
                    residue = ss[1]
                    residue_id = ss[2]
                    chain = ss[3]
                    anchor = ss[0]
                    residue_point = "_".join(ss[1:4])
                    cutpoint = ss[4]
                tmp.append([residue_point,typ,anchor,cutpoint,residue,residue_id,chain])
            return {tt[0]:tt for tt in tmp}
            
        _tmp = set([modify[1] for modify in modifies])
        
        if len(_tmp.difference(set(self._template.keys()))) > 0:
            assert False, f"the modify target set is error, the modify type must be in pho suf, met, eth, n-met, n-eth, now is {_tmp}"
        modify_groups = {mo:self.get_group_template(self._template[mo][1]) for mo in _tmp}
        
        modify_dict = parse_modify(modifies)
        Atoms = []
        nn = len(self.protein.Atoms)
        for group in self.protein.Groups:
            if group.group_str in modify_dict:
                modify = modify_dict[group.group_str]
                modify_template = deepcopy(modify_groups[modify[1]])
                anchor_atom = [self.protein.Atoms[an] for an in group.atoms if self.protein.Atoms[an].atom_name == modify[2]][0]
                if modify[3] is None:
                    cutpoint_atom = [self.protein.Atoms[an] for an in anchor_atom.connectivity if self.protein.Atoms[an].elem == "H"][0]
                else:
                    cutpoint_atom = [self.protein.Atoms[an] for an in group.atoms if self.protein.Atoms[an].atom_name == modify[3]][0]
                    
                an = anchor_atom.connectivity.index(cutpoint_atom.ID)
                
        
                remove_ans = self.protein.find_side_componend(cutpoint_atom.ID,anchor_atom.ID)
                remove_ans.append(cutpoint_atom.ID)
                cut_charge = sum([self.protein.Atoms[an].ff_charge for an in remove_ans])
                
                
                for an in group.atoms:
                    atom0 = self.protein.Atoms[an]
                    if atom0.ID not in remove_ans:
                        Atoms.append(atom0)
                        if atom0.ID == anchor_atom.ID:
                            cc = atom0.connectivity.index(cutpoint_atom.ID)
                            atom0.connectivity = [ai for ii,ai in enumerate(anchor_atom.connectivity) if ii != cc]
                            atom0.bond_type = [ai for ii,ai in enumerate(anchor_atom.bond_type) if ii != cc]
                            atom0.connectivity.append(nn)
                            atom0.bond_type.append("1")
                            _ndxs = {atom.atom_name:kk + nn for kk,atom in enumerate(modify_template)}
                            for atom in modify_template:
                                atom.ID = _ndxs[atom.atom_name]
                                atom.connectivity = [_ndxs[an] for an in atom.connectivity]
                                atom.residu = modify[4]
                                atom.residue_ID = modify[5]
                                atom.chain_name = modify[6]
                                atom.charge_group = atom0.charge_group
                            modify_template[0].connectivity.append(atom0.ID)
                            modify_template[0].bond_type.append("1")
                            if create_3D:
                                modify_template[0].ff_charge = modify_template[0].ff_charge + cut_charge
                            Atoms.extend(modify_template)
                            nn += len(modify_template)
            else:
                for an in group.atoms:
                    Atoms.append(self.protein.Atoms[an])
                    
        _atom_ID = {atom.ID:ii for ii,atom in enumerate(Atoms)}
        for atom in Atoms:
            atom.connectivity = [_atom_ID[an] for an in atom.connectivity]
            atom.ID = _atom_ID[atom.ID]  
                  
        modify_protein = Molecule(style="protein")
        modify_protein.mole_name = self.protein.mole_name
        modify_protein.Atoms = Atoms
        modify_protein.create_topols()
        modify_protein.create_intra_nonbond_macromole()
        modify_protein.create_improper()
        
        if create_3D:
            PP = ProteinPrepare(modify_protein)
            PP.create_coor_for_loss_atoms()
            PP.assign_force_field()
            PP.optimize_position_of_loss_atom()
            modify_protein = PP.protein
            #return PP.protein  
        modify_protein = protein_ring_and_charge_group(modify_protein) 
        return modify_protein    

    def delete_residue(self,residue,n_terminal=None,c_terminal=None,terminal=None,shift_flag=True):
        Atoms = []
        terminal_del = []
        shift = [0,0]
        if terminal == "both":
            terminal_del = ["ACE","NME","NMA","MEC","MEN","NHE","NH2"]
        if terminal == "n_terminal":
            shift = [len(residue)-1,0]
            terminal_del = ["ACE","MEN",]
        if terminal == "c_terminal":
            terminal_del = ["NME","NMA","MEC","NHE","NH2"]
        
        for group in self.protein.Groups:
            if group.group_str not in residue:
                if group.group_name not in terminal_del:
                    for an in group.atoms:
                        Atoms.append(self.protein.Atoms[an])
        
        tmp_protein = Molecule("protein")
        tmp_protein.Atoms = Atoms
        tmp_protein.molecule_name = self.protein.molecule_name
        
        PP = ProteinPrepare(tmp_protein,n_terminal=n_terminal,c_terminal=c_terminal)
        PP.run()
        protein = PP.protein
        if hasattr(self.protein,"relative_shift"):
            protein.relative_shift = self.protein.relative_shift

        if not hasattr(protein,"relative_shift"):
            protein.relative_shift = shift
        else:
            
            if not shift_flag:
                protein.relative_shift = shift
            else:
                protein.relative_shift = [protein.relative_shift[0]+shift[0],protein.relative_shift[1]]
        #protein.mole_name = "_".join([self.protein,residue])
        protein = protein_ring_and_charge_group(protein)
        return protein
    
    def add_residue(self,target,residue,terminal=None,n_terminal=None,c_terminal=None):
        Atoms = []
        for group in target.Groups:
            if group.group_str not in residue:
                for an in group.atoms:
                    Atoms.append(target.Atoms[an])
        tmp_protein = Molecule("protein")
        tmp_protein.Atoms = Atoms
        tmp_protein.molecule_name = self.protein.molecule_name
        
        PP = ProteinPrepare(tmp_protein,n_terminal=n_terminal,c_terminal=c_terminal)
        PP.run()
        protein = PP.protein
        protein.mole_name = "_".join([target,mutation])
        protein = protein_ring_and_charge_group(protein)
        return protein

    def mutation_residue(self,contents,n_terminal=None,c_terminal=None,terminal=None):
        mutation_residue = {content[1]:template_molecule(amino_acid[content[1]],content[1]) for content in contents}
        targets = [content[0] for content in contents]
        Atoms = []
        for group in self.protein.Groups:
            if group.group_str in targets:
                idx = targets.index(group.group_str)
                mutation_res = deepcopy(mutation_residue[contents[idx][1]])
                save_coor_atoms = {self.protein.Atoms[an].atom_name:self.protein.Atoms[an].coordinates for an in group.atoms}
                for atom in mutation_res.Atoms:
                    atom.residue_ID = group.group_idx
                    atom.chain_name = group.group_chain_name
                    if atom.atom_name in  ["N", "H", "CA", "HA","CB","C", "O","OC1", "OC2","H1", "H2", "H3","HB","HB1","HB2"]:
                        if atom.atom_name in save_coor_atoms:
                            atom.coordinates = save_coor_atoms[atom.atom_name]
                Atoms.extend(mutation_res.Atoms)
            else:
                for an in group.atoms:
                    Atoms.append(self.protein.Atoms[an])
        tmp_protein = Molecule("protein")
        tmp_protein.Atoms = Atoms
        tmp_protein.molecule_name = self.protein.molecule_name
        
        PP = ProteinPrepare(tmp_protein,n_terminal=n_terminal,c_terminal=c_terminal)
        PP.run()
        protein = PP.protein
        protein = protein_ring_and_charge_group(protein)
        #protein.mole_name = "_".join([target,mutation])
        return protein

    def get_target_atoms(self,target,mutation):
        if mutation is not None:
            mutation_residue = template_molecule(amino_acid[mutation],mutation)
        else:
            mutation_residue = None
        pre_atoms = []
        target_atoms = []
        post_atoms = []
        flag = False
        for atom in self.protein.Atoms:
            residue_name = f"{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if residue_name == target:
                target_atoms.append(atom)
                flag = True
            else:
                if flag:
                    post_atoms.append(atom)
                else:
                    pre_atoms.append(atom)
        return mutation_residue, target_atoms, pre_atoms,post_atoms

    def old_residue_mutation(self,target,mutation,n_terminal=None,c_terminal=None,):
        mutation_residue,target_atoms,pre_atoms,post_atoms = self.get_target_atoms(target,mutation)
        new_atoms = []
        target_atoms_coord = {atom.atom_name:atom.coordinates for atom in target_atoms}        
        residue_ID = target_atoms[0].residue_ID
        chain_name = target_atoms[0].chain_name
        
        if mutation_residue is not None:
            for atom in mutation_residue.Atoms:
                atom.residue_ID = residue_ID
                atom.chain_name = chain_name
                if atom.atom_name in ["N", "H", "CA", "HA","CB","C", "O","OC1", "OC2","H1", "H2", "H3","HB","HB1","HB2"]:
                    if atom.atom_name in target_atoms_coord:
                        atom.coordinates = target_atoms_coord[atom.atom_name]
                new_atoms.append(atom)
        tmp_protein = Molecule("protein")
        tmp_protein.Atoms = pre_atoms + new_atoms + post_atoms
        tmp_protein.molecule_name = self.protein.molecule_name
        
        PP = ProteinPrepare(tmp_protein,n_terminal=n_terminal,c_terminal=c_terminal)
        PP.run()
        protein = PP.protein
        protein.mole_name = "_".join([target,mutation])
        return protein

    def get_group_template(self,group):
        group_template = modify_groups[group]
        Atoms = [Atom(style="atom") for ii in range(len(group_template))]
        for ii,atom in enumerate(Atoms):
            for kk,vv in group_template[ii].items():
                setattr(atom,kk,vv)
        return Atoms
    
    def residue_modify(self,target,modify,create_3D=True):
        if modify in self._template:
            _label = self._template[modify][0]
            modify_template = self.get_group_template(self._template[modify][1])
            modify_name = modify
        else:
            modify_name = "chg"
            #modify_template = 
            pass
        
        ss = target.split("_")
        if len(ss) == 3:
            if modify not in self._template:
                assert False, f"the modify target set is therr, the modify type must be in pho suf, met, eth, n-met, n-eth, now is {modify}"  
            residue = ss[0]
            residue_id = ss[1]
            chain = ss[2]
            try:
                anchor = f"{_label[residue]}_{target}"
            except:
                anchor = f"{_label['all']}_{target}"
            cutpoint = None
            residue_point = target
        elif len(ss) == 4:
            residue = ss[1]
            residue_id = ss[2]
            chain = ss[3]
            anchor = target
            cutpoint = None
            residue_point = "_".join(ss[1:])
        elif len(ss) == 5:
            residue = ss[1]
            residue_id = ss[2]
            chain = ss[3]
            anchor = "_".join(ss[:4])
            residue_point = "_".join(ss[1:4])
            cutpoint = f"{ss[4]}_{residue_point}"
            
        anchor_atom = None
        cutpoint_atom = None
        for atom in self.protein.Atoms:
            this_atom_name = f"{atom.atom_name}_{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if this_atom_name == anchor:
                anchor_atom = atom
                if cutpoint is None:
                    for an in atom.connectivity:
                        if self.protein.Atoms[an].elem == "H":
                            cutpoint_atom = self.protein.Atoms[an]
                            break
            if this_atom_name == cutpoint:
                cutpoint_atom = atom    
            if anchor_atom is not None and cutpoint_atom is not None:
                break
        an = anchor_atom.connectivity.index(cutpoint_atom.ID)
        anchor_atom.connectivity = [ai for ii,ai in enumerate(anchor_atom.connectivity) if ii != an]
        anchor_atom.bond_type = [ai for ii,ai in enumerate(anchor_atom.bond_type) if ii != an]
        
        remove_ans = self.protein.find_side_componend(cutpoint_atom.ID,anchor_atom.ID)
        remove_ans.append(cutpoint_atom.ID)
        cut_charge = sum([self.protein.Atoms[an].ff_charge for an in remove_ans])
        nn = len(self.protein.Atoms)
        anchor_atom.connectivity.append(nn)
        anchor_atom.bond_type.append("1")
        _ndxs = {atom.atom_name:kk + nn for kk,atom in enumerate(modify_template)}
        for atom in modify_template:
            atom.ID = _ndxs[atom.atom_name]
            atom.connectivity = [_ndxs[an] for an in atom.connectivity]
            atom.residu = residue
            atom.chain_name = chain
            atom.residue_ID = residue_id
            atom.charge_group = anchor_atom.charge_group
        modify_template[0].connectivity.append(anchor_atom.ID)
        modify_template[0].bond_type.append("1")
        if create_3D:
            modify_template[0].ff_charge = modify_template[0].ff_charge + cut_charge
        
        
        Atoms = []
        
        ii = 0
        ids = {}
        for atom in self.protein.Atoms:
            if atom.ID not in remove_ans:
                ids[atom.ID] = ii
                Atoms.append(atom)
                if atom.ID == anchor_atom.ID:
                    for atom_m in modify_template:
                        ii += 1
                        ids[atom_m.ID] = ii
                        Atoms.append(atom_m)
                        
                ii += 1
        
        for atom in Atoms:
            atom.ID = ids[atom.ID]
            atom.connectivity = [ids[an] for an in atom.connectivity]
        
        modify_protein = Molecule(style="protein")
        modify_protein.mole_name = f"{target}_{modify_name}"
        modify_protein.Atoms = Atoms
        modify_protein.create_topols()
        modify_protein.create_intra_nonbond_macromole()
        modify_protein.create_improper()
        
        if create_3D:
            PP = ProteinPrepare(modify_protein)
            PP.create_coor_for_loss_atoms()
            PP.assign_force_field()
            PP.optimize_position_of_loss_atom()
            return PP.protein    
        return modify_protein
    
    def phosphorylation_modify(self,target,create_3D=True):
        _label = {"SER":["OG","HG"],"THR":["OG1","HG1"],"TYR":["OH","HH"]}
        modify_template = self.get_group_template("Pho")
        ss = target.split("_")
        if len(ss) == 3:
            residue = ss[0]
            residue_id = ss[1]
            chain = ss[2]
            anchor = f"{_label[residue][0]}_{target}"
            cutpoint = None
            residue_point = target
        elif len(ss) == 4:
            residue = ss[1]
            residue_id = ss[2]
            chain = ss[3]
            anchor = target
            cutpoint = None
            residue_point = "_".join(ss[1:])
        elif len(ss) == 5:
            residue = ss[1]
            residue_id = ss[2]
            chain = ss[3]
            anchor = "_".join(ss[:4])
            residue_point = "_".join(ss[1:4])
            cutpoint = f"{ss[4]}_{residue_point}"
        if residue not in ["SER","THR","TYR"]:
            assert False, f"residue must be SER, THR or TYR, {residue}"    
        
        anchor_atom = None
        cutpoint_atom = None
        for atom in self.protein.Atoms:
            this_atom_name = f"{atom.atom_name}_{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if this_atom_name == anchor:
                anchor_atom = atom
                if cutpoint is None:
                    for an in atom.connectivity:
                        if self.protein.Atoms[an].elem == "H":
                            cutpoint_atom = self.protein.Atoms[an]
                            break
            if this_atom_name == cutpoint:
                cutpoint_atom = atom    
            if anchor_atom is not None and cutpoint_atom is not None:
                break
        an = anchor_atom.connectivity.index(cutpoint_atom.ID)
        anchor_atom.connectivity = [ai for ii,ai in enumerate(anchor_atom.connectivity) if ii != an]
        anchor_atom.bond_type = [ai for ii,ai in enumerate(anchor_atom.bond_type) if ii != an]
        
        remove_ans = self.protein.find_side_componend(cutpoint_atom.ID,anchor_atom.ID)
        remove_ans.append(cutpoint_atom.ID)
        cut_charge = sum([self.protein.Atoms[an].ff_charge for an in remove_ans])
        nn = len(self.protein.Atoms)
        anchor_atom.connectivity.append(nn)
        anchor_atom.bond_type.append("1")
        _ndxs = {atom.atom_name:kk + nn for kk,atom in enumerate(modify_template)}
        for atom in modify_template:
            atom.ID = _ndxs[atom.atom_name]
            atom.connectivity = [_ndxs[an] for an in atom.connectivity]
            atom.residu = residue
            atom.chain_name = chain
            atom.residue_ID = residue_id
        modify_template[0].connectivity.append(anchor_atom.ID)
        modify_template[0].bond_type.append("1")
        
        
        Atoms = []
        
        ii = 0
        ids = {}
        for atom in self.protein.Atoms:
            if atom.ID not in remove_ans:
                ids[atom.ID] = ii
                Atoms.append(atom)
                if atom.ID == anchor_atom.ID:
                    for atom_m in modify_template:
                        ii += 1
                        ids[atom_m.ID] = ii
                        Atoms.append(atom_m)
                        
                ii += 1
        
        for atom in Atoms:
            atom.ID = ids[atom.ID]
            atom.connectivity = [ids[an] for an in atom.connectivity]
        
        modify_protein = Molecule(style="protein")
        modify_protein.Atoms = Atoms
        modify_protein.create_topols()
        modify_protein.create_intra_nonbond_macromole()
        modify_protein.create_improper()
        
        if create_3D:
            PP = ProteinPrepare(modify_protein)
            PP.create_coor_for_loss_atoms()
            PP.assign_force_field()
            PP.optimize_position_of_loss_atom()
            return PP.protein    
        return modify_protein
    
    def residue_termination(self):
        pass
    
    def residue_ionization(self):
        pass
    
    def acetylation_modify(self,):
        pass
    
    def N_glycosylation_modify(self,):
        pass
    
    def amidation_modify(self):
        pass
    
    def hydroxylation_modify(self):
        pass
    
    def methylation_modify(self):
        pass
    
    def O_glycosylation_modify(self,):
        pass
    
    def ubiquitylation_modify(self):
        pass
    
    def pyrrolidone_carboxylic_acid_modify(self):
        pass
    
    def sulfation_modify(self):
        pass
    
