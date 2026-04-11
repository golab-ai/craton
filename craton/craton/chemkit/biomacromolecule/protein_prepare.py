
import json
from copy import deepcopy
from ...utils.geometry import calc_stru_para
from ...utils import logger
from ...chem.atom import Atom
from ...chem.molecule import Molecule

#from ...force_field import MolForceField as MFF #grasp_force_field

from ..stereo import Stereo
from .protein_utils import optimize_position
from .protein_utils import amino_acid, n_terminal_molecule, c_terminal_molecule, template_molecule
from .protein_utils import TERMINAL_RESIDUE, CUTOFF_RESIDUE, DEFAULT_PROTEIN_FORCEFIELD
from .protein_utils import DEFAULT_COMPEN_TYPING_FILE, DEFAULT_COMPEN_FORCE_FIELD_FILE
      
class ProteinPrepare:
    def __init__(self,protein,n_terminal=None,c_terminal=None) -> None:
        self.protein = protein
        if n_terminal is None:
            self.n_terminal_molecule = n_terminal_molecule
        else:
            self.n_terminal_molecule = template_molecule(amino_acid[n_terminal],n_terminal)
        
        if c_terminal is None:
            self.c_terminal_molecule = c_terminal_molecule
        else:
            self.c_terminal_molecule = template_molecule(amino_acid[c_terminal],c_terminal)
        
        self.protein_name = self.protein.mole_name
        self.residues = []

    def run(self):
        self.get_residue()
        if len(self.proteins_atoms) > 0:
            self.check_residue_connect()
        if len(self.rnadna_atoms) > 0:
            self.check_base_connect()
        self.assign_id_and_check_loss_atoms()
        self.assign_connectivity_and_atom_type()
        self.protein.create_intra_nonbond_macromole()
        self.assign_force_field()
        self.create_coor_for_loss_atoms()
        self.optimize_position_of_loss_atom()
        
    def get_residue(self):
        proteins_atoms = []
        rnadan_atoms = []
        for atom in self.protein.Atoms:
            if atom.residue in ["A","U","G","C","DA","DT","DG","DC"]:
                rnadan_atoms.append(atom)
            else:
                proteins_atoms.append(atom)
        self.proteins_atoms = proteins_atoms
        self.rnadna_atoms = rnadan_atoms

    def check_base_connect(self):
        residues = []
        atom0 = self.rnadna_atoms[0]
        pre_name = f"{atom0.residue}_{atom0.residue_ID}_{atom0.chain_name}"
        tmp = []
        for atom in self.rnadna_atoms:
            res_name = f"{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if res_name == pre_name:
                tmp.append(atom)
            else:
                tmp[0].n_sub_type = ""
                tmp[0].c_sub_type = ""
                tmp[0].sub_residue = tmp[0].residue
                residues.append(tmp)
                
                pre_name = res_name
                tmp = [atom]
        tmp[0].n_sub_type = ""
        tmp[0].c_sub_type = ""
        tmp[0].sub_residue = tmp[0].residue
        residues.append(tmp)
        
        residues = [deepcopy([residues[0][0]])] + residues +[deepcopy([residues[0][0]])]
        residues[0][0].sub_residue  = residues[0][0].residue + "3"
        residues[0][0].c_sub_type = "C"
        
        residues[-1][0].sub_residue  = residues[0][0].residue + "5"
        residues[-1][0].n_sub_type = "N"
        
        for ii in range(1,len(residues)):
            residue = residues[ii]
            pre_c_sub_type = residues[ii-1][0].c_sub_type
            if pre_c_sub_type == "C":
                residue[0].sub_residue = residue[0].residue + "5"
                residue[0].n_sub_type = "N"
            
            H5T = [atom for atom in residue if atom.atom_name == "H5T"]
            H3T = [atom for atom in residue if atom.atom_name == "H3T"]
            if len(H5T) > 0:
                residue[0].sub_residue = residue[0].residue + "5"
                residue[0].n_sub_type = "N"
            if len(H3T) > 0:
                residue[0].sub_residue = residue[0].residue + "3"
                residue[0].c_sub_type = "C"
            
            if residue[0].n_sub_type == "N" and residue[0].c_sub_type == "C":
                residue[0].sub_residue = residue[0].residue + "N"
            
            if residue[0].n_sub_type == "N":
                residues[ii-1][0].c_sub_type == "C"
                residues[ii-1][0].sub_residue = residues[ii-1][0].residue + "3"
                
            else:
                pre_LT = [an for an in residues[ii-1] if an.atom_name == "O3'"][0]
                RT_arr = [an for an in residue if an.atom_name == "P"]
                if len(RT_arr) == 0:
                    residue[0].sub_residue = residue[0].residue + "3"
                    residue[0].c_sub_type = "C"
                else:
                    RT = RT_arr[0]
                    
                    dis = calc_stru_para([pre_LT.coor,RT.coor])
                    if dis > CUTOFF_RESIDUE:
                        logger.warning(f"{pre_LT.residue} {pre_LT.residue_ID} {pre_LT.chain_name} and {RT.residue} {RT.residue_ID} {RT.chain_name} break, the distance is {dis}")
                        logger.warning(f"{pre_LT.residue} {pre_LT.residue_ID} {pre_LT.chain_name} has been changed as C terminal")
                        residues[ii-1][0].c_sub_type == "C"
                        residues[ii-1][0].sub_residue = residues[ii-1][0].residue + "3"    
                        
                        logger.warning(f"{RT.residue} {RT.residue_ID} {RT.chain_name} has been changed as N terminal")
                        residue[0].sub_residue = residue[0].residue + "5"
                        residue[0].n_sub_type = "N"
        self.residues.extend(residues[1:-1])    
                
    def check_residue_connect(self):
        def get_LT_RT(residue_arr):
            template = amino_acid[residue_arr[0].residue]
            LT_name = [an for an,vv in template.items() if "connectivity" in vv and "R*" in vv["connectivity"]]
            RT_name = [an for an,vv in template.items() if "connectivity" in vv and "L*" in vv["connectivity"]]
            if len(LT_name) != 0 and len(RT_name) != 0:
                LT_arr = [atom for atom in residue_arr if atom.atom_name in LT_name]
                RT_arr = [atom for atom in residue_arr if atom.atom_name in RT_name]
                if residue_arr[0].residue not in TERMINAL_RESIDUE:
                    if len(LT_arr) == 0 or len(RT_arr) == 0:
                        return None, None, None, "missing_backbond"
                
                LT = LT_arr[0]
                RT = RT_arr[0]
                RT_C = sum([atom.formal_charge for atom in residue_arr if atom.atom_name in template[RT.atom_name]["connectivity"]]) if RT.atom_name == "C" else None
                RT_H = [atom for atom in residue_arr if atom.atom_name in ["H1","H2","H3","H","1H","2H","3H"]]
                n_label_flag = False
                c_label_flag = False
                if LT.atom_name == "N":
                    if LT.formal_charge == 1:
                        residue_arr[0].n_sub_type = "N"
                        n_label_flag = True  # "n_terminal"
                    else:
                        if len(RT_H) == 2:
                            residue_arr[0].n_sub_type = "H"
                            n_label_flag = True  # "n_terminal"
                            pass
                        elif len(RT_H) == 3:
                            residue_arr[0].n_sub_type = "N"
                            n_label_flag = True  # "n_terminal"
                if RT.atom_name == "C":
                    if RT_C == -1:
                        residue_arr[0].c_sub_type = "C"
                        c_label_flag = True # "c_terminal"
                if n_label_flag and c_label_flag:
                    return LT,RT, RT_C, "double_terminal"
                else:
                    if n_label_flag:
                        return LT, RT, RT_C, "n_terminal"
                    else:
                        if c_label_flag:
                            return LT, RT, RT_C, "c_terminal"
                        else:
                            return LT, RT, RT_C, "normal"
            else:
                if len(LT_name) == 0 and len(RT_name) == 0:
                    return None,None,None,"double_terminal"
                else:
                    if len(LT_name) == 0:
                        RT = [atom for atom in residue_arr if atom.atom_name in RT_name][0]
                        LT = None
                        return LT, RT, None, "n_terminal"
                    else:
                        LT = [atom for atom in residue_arr if atom.atom_name in LT_name][0]
                        RT = None
                        return LT, RT, None, "c_terminal"
       
        def get_terminal_molecule(template,resiude_ID,chain_name):
            atoms = deepcopy(template.Atoms)
            
            atoms[0].c_sub_type = ""
            atoms[0].n_sub_type = ""
            for atom in atoms:
                atom.residue_ID = resiude_ID
                atom.chain_name = chain_name
            return atoms 
        
        ####拆分成残基组合######################################
        labels = []
        residues = []
        atom0 = self.proteins_atoms[0]
        pre_name = f"{atom0.residue}_{atom0.residue_ID}_{atom0.chain_name}"
        tmp = []
        OXT_flag = False
        for atom in self.proteins_atoms:
            #if atom.atom_name[0].isdigit():
            #    atom.atom_name = atom.atom_name[1:] + atom.atom_name[0]
            res_name = f"{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if res_name == pre_name:
                if atom.atom_name == "OXT":
                    atom.atom_name = "OC2"
                    OXT_flag = True
                tmp.append(atom)
                
            else:
                if OXT_flag:
                    try:
                        o_c1 = [ii for ii,aa in enumerate(tmp) if aa.atom_name == "O"][0]
                    except:
                        o_c1 = None
                    if o_c1 is not None:
                        tmp[o_c1].atom_name = "OC1"
                    
                    
                OXT_flag = False
                tmp[0].n_sub_type = ""
                tmp[0].c_sub_type = ""
                LT,RT,RT_C,label = get_LT_RT(tmp)
                if label != "missing_backbond":
                    labels.append([LT,RT,label,tmp[0].residue_ID,tmp[0].chain_name,tmp[0].residue])
                    residues.append(tmp)
                else:
                    logger.warning(f"missing backband atoms: ignore {tmp[0].residue} {tmp[0].residue_ID} {tmp[0].chain_name}")
                pre_name = res_name
                tmp = [atom]
        tmp[0].n_sub_type = ""
        tmp[0].c_sub_type = ""
        LT,RT,RT_C,label = get_LT_RT(tmp)
        if label != "missing_backbond":
            labels.append([LT,RT,label,tmp[0].residue_ID,tmp[0].chain_name,tmp[0].residue])
            residues.append(tmp)
        else:
            logger.warning(f"missing backband atoms: ignore {tmp[0].residue} {tmp[0].residue_ID} {tmp[0].chain_name}")
        
        labels = [[None,None,None,None,None,None]] + labels + [[None,None,None,None,None,None]]
        residues = [None] + residues + [None]
        ##########################################################
        
        ####添加封端###############################################
        new_residues = []
        for ii in range(1,len(labels)):
            residue = residues[ii]
            label = labels[ii]
            pre_label = labels[ii-1]
            if label[2] == "n_terminal":
                if pre_label[2] == "normal":
                    logger.warning(f"C terminal has been added to {pre_label[5]} {pre_label[3]} {pre_label[4]}")
                    new_residues.append(get_terminal_molecule(self.c_terminal_molecule,f"{pre_label[3]}R",pre_label[4]))
                new_residues.append(residue)
            elif label[2] is None:
                if pre_label[2] != "c_terminal":
                    logger.warning(f"C terminal has been added to {pre_label[5]} {pre_label[3]} {pre_label[4]}")
                    new_residues.append(get_terminal_molecule(self.c_terminal_molecule,f"{pre_label[3]}R",pre_label[4]))
            else:
                if pre_label[2] in [None,"c_terminal"]:
                    logger.warning(f"N terminal has been added to {label[5]} {label[3]} {label[4]}")
                    new_residues.append(get_terminal_molecule(self.n_terminal_molecule,f"{label[3]}L",label[4]))
                    new_residues.append(residue)
                    
                else:
                    dis = calc_stru_para([pre_label[1].coor,label[0].coor])
                    if dis > CUTOFF_RESIDUE:
                        logger.warning(f"{pre_label[5]} {pre_label[3]} {pre_label[4]} and {label[5]} {label[3]} {label[4]} break, the distance is {dis}")
                        logger.warning(f"C terminal has been added to {pre_label[5]} {pre_label[3]} {pre_label[4]}")
                        new_residues.append(get_terminal_molecule(self.c_terminal_molecule,f"{pre_label[3]}R",pre_label[4]))
                        
                        logger.warning(f"N terminal has been added to {label[5]} {label[3]} {label[4]}")
                        new_residues.append(get_terminal_molecule(self.n_terminal_molecule,f"{label[3]}L",label[4]))
                        new_residues.append(residue)
                    else:
                        new_residues.append(residue)
        new_residues = self.assign_residue_subtype(new_residues)                
        
        
        self.residues += new_residues

    def CYS_prepare(self,CYS,residues):
        S_Atoms = [atom for sii in CYS for atom in residues[sii] if atom.atom_name == "SG" ]
        SN = len(S_Atoms)
        S_S = {}
        _tmp_1 = set([ii for ii in range(SN)])
        _tmp = []
        if SN > 1:
            for ii in range(SN):
                _tmp.append(ii)
                S1 = S_Atoms[ii]
                for jj in _tmp_1.difference(set(_tmp)):
                    S2 = S_Atoms[jj]
                    dis = calc_stru_para([S1.coor,S2.coor])
                    if dis <= 2.2:
                        S_S[CYS[ii]] = CYS[jj]
                        S_S[CYS[jj]] = CYS[ii]
                        _tmp.append(ii)
                        _tmp.append(jj)
                        break
        return S_S

    def assign_residue_subtype(self,residues):
        CYS = []
        for ii, residue in enumerate(residues):
            sub_residue = residue[0].n_sub_type + residue[0].c_sub_type+residue[0].residue
            if sub_residue == "HIS":
                HD1 = [atom for atom in residue if atom.atom_name == "HD1"]
                HE2 = [atom for atom in residue if atom.atom_name == "HE2"]
                if len(HD1) == 1 and len(HE2) == 1:
                    sub_residue = "HIP"
                else:
                    if len(HD1) == 1:
                        sub_residue = "HID"
                    else:
                        if len(HE2) ==1:
                            sub_residue = "HIS"
                        else:
                            sub_residue = "HIP"
            if sub_residue == "GLU":
                Charge = sum([atom.formal_charge for atom in residue if atom.atom_name in ["OE2","OE1"]])
                OE2 = [atom for atom in residue if atom.atom_name == "HE2"]
                if Charge != 0 or len(OE2) != 1:
                    #sub_residue = "GLH"
                    pass
                else:
                    sub_residue = "GLH"
            if sub_residue == "LYS":
                Charge = sum([atom.formal_charge for atom in residue if atom.atom_name in ["NZ"]])
                HZ3 = [atom for atom in residue if atom.atom_name == "HZ3"]
                if len(HZ3) != 0 or Charge != 0:
                    #sub_residue = "LYN"
                    pass
                else:
                    sub_residue = "LYN"
            if sub_residue == "ASP":
                Charge = sum([atom.formal_charge for atom in residue if atom.atom_name in ["OD2","OD1"]])
                HD2 = [atom for atom in residue if atom.atom_name == "HD2"]
                if len(HD2) != 1 or Charge != 0:
                    #sub_residue = "ASH"
                    pass
                else:
                    sub_residue = "ASH"
            
            #if sub_residue in ["CYS","CCYS","NCYS"]:
            if residue[0].residue == "CYS":
                CYS.append(ii)
            residue[0].sub_residue = sub_residue
                
        S_S = self.CYS_prepare(CYS,residues)
        for sii in CYS:
            if sii not in S_S:
                HG = [atom for atom in residues[sii] if atom.atom_name == "HG"]
                if len(HG) == 0:
                    residues[sii][0].sub_residue = "CYM"
                else:
                    residues[sii][0].sub_residue = "CYS"
            else:
                residues[sii][0].sub_residue = "CYX"
                residues[sii][0].CYX_connectivity = S_S[sii]
        return residues
            
    def assign_id_and_check_loss_atoms(self):
        pre_nn = 0
        _tmp_residues = []
        for ii,residue in enumerate(self.residues):
            residue_name = residue[0].residue
            sub_residue_name = residue[0].sub_residue
            target_atoms = amino_acid[residue_name]["template"][sub_residue_name][0]
            this_atoms = [atom.atom_name for atom in residue]
            _tmp = []
            _dicts = {atom.atom_name:atom for atom in residue}
            for kk, tar in enumerate(target_atoms):
                nn = pre_nn + kk
                if tar in _dicts:
                    atom = _dicts[tar]
                    atom.ID = nn
                    atom.charge_group = ii
                else:
                    _temp = amino_acid[residue_name][tar]
                    atom = Atom("atom")
                    atom.residue = residue_name
                    atom.residue_ID = residue[0].residue_ID
                    atom.chain_name = residue[0].chain_name
                    atom.sub_residue = residue[0].sub_residue
                    atom.c_sub_type = residue[0].c_sub_type
                    atom.n_sub_type = residue[0].n_sub_type
                    atom.element = _temp["element"]
                    atom.atom_name = _temp["atom_name"]
                    atom.ID = nn
                    atom.charge_group = ii
                _tmp.append(atom)
            if not hasattr(_tmp[0],"sub_residue"):
                _tmp[0].sub_residue = residue[0].sub_residue
                _tmp[0].c_sub_type = residue[0].c_sub_type
                _tmp[0].n_sub_type = residue[0].n_sub_type 
            _tmp_residues.append(_tmp)
            pre_nn = nn + 1
        self.residues = _tmp_residues
                    
    def assign_connectivity_and_atom_type(self):
        pre_L = None
        for ii, residue_arr in enumerate(self.residues):
            residue = residue_arr[0].residue
            atom_type_name = amino_acid[residue]["template"][residue_arr[0].sub_residue][1]
            formal_charge = amino_acid[residue]["template"][residue_arr[0].sub_residue][2]
            ff_charge = amino_acid[residue]["template"][residue_arr[0].sub_residue][3]
            connectivity = amino_acid[residue]["template"][residue_arr[0].sub_residue][4]
            bond_type = amino_acid[residue]["template"][residue_arr[0].sub_residue][5]
            
            ID_dict = {atom.atom_name:atom.ID for atom in residue_arr}
            this_L = None
            n_sub_type = residue_arr[0].n_sub_type
            for jj,atom in enumerate(residue_arr):
                atom.atom_type_name = amino_acid[residue][atom.atom_name][atom_type_name]
                atom.formal_charge = amino_acid[residue][atom.atom_name][formal_charge]
                atom.ff_charge = amino_acid[residue][atom.atom_name][ff_charge]
                atom.plate = amino_acid[residue][atom.atom_name]["plate"]
                atom.connectivity = []
                atom.bond_type = []
                if "local" in amino_acid[residue][atom.atom_name]:
                    for ttt in ["has_ring","has_ring_size","has_ring_property","local","partial_formal_charge"]:
                        setattr(atom,ttt,amino_acid[residue][atom.atom_name][ttt])
                    atom.bond_type_aromatic = []
                    atom.connectivity_type = []
                    atom.bond_type_conjugate = []
                
                
                for kk,an in enumerate(amino_acid[residue][atom.atom_name][connectivity]):
                    if an in ID_dict:
                        atom.connectivity.append(ID_dict[an])
                        atom.bond_type.append(amino_acid[residue][atom.atom_name][bond_type][kk])
                        if hasattr(atom,"bond_type_conjugate"):
                            atom.bond_type_aromatic.append(amino_acid[residue][atom.atom_name]["bond_type_aromatic"][kk])
                            atom.connectivity_type.append(amino_acid[residue][atom.atom_name]["connectivity_type"][kk])
                            atom.bond_type_conjugate.append(amino_acid[residue][atom.atom_name]["bond_type_conjugate"][kk])
                            
                    if an == "L*":
                        this_L = [jj,atom.ID]
                    if an == "R*" and n_sub_type == "":
                        if pre_L is not None:
                            atom.connectivity.append(pre_L[1])
                            atom.bond_type.append(amino_acid[residue][atom.atom_name][bond_type][kk])
                            self.residues[ii-1][pre_L[0]].connectivity.append(atom.ID)
                            self.residues[ii-1][pre_L[0]].bond_type.append(amino_acid[residue][atom.atom_name][bond_type][kk])
            pre_L = this_L
        for ii,residue_arr in enumerate(self.residues):
            if residue_arr[0].sub_residue == "CYX":
                this_SG = [atom for atom in residue_arr if atom.atom_name == "SG"][0]
                conn_SG = [atom.ID for atom in self.residues[residue_arr[0].CYX_connectivity] if atom.atom_name == "SG"][0]
                this_SG.connectivity.append(conn_SG)
                this_SG.bond_type.append("1")

        atoms = [atom for residue in self.residues for atom in residue]
        self.protein = Molecule("protein")
        self.protein.molecule_name = self.protein_name 
        self.protein.Atoms = atoms
        self.protein.residu_n = len(self.residues)
        self.protein.create_topols(smiles_flag=False)
        self.protein.create_improper(create_method="atom_type")
        self.protein.create_intra_nonbond_macromole()
    
    def assign_force_field(self):
        from ...force_field import MolForceField as MFF #grasp_force_field
        for atom in self.protein.Atoms:
            atom.nonb_atom_type = atom.atom_type_name
            atom.binc_atom_type = atom.atom_type_name
            atom.atc_atom_type = atom.atom_type_names
        for attr in ["Bonds","Angles","Dihedrals","Impropers","Pair14"]:
            for term in getattr(self.protein,attr,[]):
                for ii in range(1,5):
                    if hasattr(term,f"a{ii}"):
                        setattr(term,f"a{ii}_atom_type",self.protein.Atoms[getattr(term,f"a{ii}")].atom_type_name)
                        setattr(term,f"a{ii}_atom_type_used",self.protein.Atoms[getattr(term,f"a{ii}")].atom_type_name)
            
        #proteins = MFF.grasp_force_field(self.protein,None,DEFAULT_PROTEIN_FORCEFIELD,empi_ff_flag=False,use_scalevdw=False,parallel=False)
        
        proteins = MFF.get_force_field(self.protein,
                                              None,
                                              DEFAULT_PROTEIN_FORCEFIELD,
                                              compensating_at_file=DEFAULT_COMPEN_TYPING_FILE,
                                              compensating_ff_file=DEFAULT_COMPEN_FORCE_FIELD_FILE,
                                              use_scalevdw=False,
                                              empi_ff_flag=True,
                                              parallel=False)
        self.protein = proteins[0]
    
    def create_coor_for_loss_atoms(self):
        loss_atoms = [aa.ID for aa in self.protein.Atoms if not hasattr(aa,"coordinates")]
        self.protein.roots = []
        self.protein.addHs = []
        Stereo.create_coor(self.protein,loss_atoms)
        self.roots,self.addHs = deepcopy(self.protein.roots), deepcopy(self.protein.addHs)
        delattr(self.protein,"roots")
        delattr(self.protein,"addHs")
        
    def optimize_position_of_loss_atom(self):
        optimize_position(self.protein,self.roots,self.addHs)
            


