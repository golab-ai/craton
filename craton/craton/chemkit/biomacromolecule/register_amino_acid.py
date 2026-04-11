import itertools
import json
import os
from copy import deepcopy
from .protein_utils import amino_acid, non_AA_register, non_normal_amino_acid,non_normal_amino_acid_json_f

AA_TEMPT = {
        "N":{"element": "N",
            "atom_name": "N",
            "atom_type_name": "N",
            "connectivity": [
                "H",
                "CA",
                "H1",
                "H2",
                "H3",
                "R*"
            ],
            "bond_type": [
                "1",
                "1",
                "1",
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": None,
            "c_formal_charge": 0,
            "n_formal_charge": 1,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "LT",
            "bond_type_aromatic": ["1","1","1","1","1","1"],
            "connectivity_type": ["S","S","S","S","S","S"],
            "bond_type_conjugate": ["","","","","",""],
            "partial_formal_charge": 0,
            },
        "H":{
            "element": "H",
            "atom_name": "H",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": None,
            "c_formal_charge": 0,
            "n_formal_charge": 0,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "EN",
            "bond_type_aromatic": ["1"],
            "connectivity_type": ["S"],
            "bond_type_conjugate": [""],
            "partial_formal_charge": 0,
        },
        "C":{
            "element": "C",
            "atom_name": "C",
            "atom_type_name": "C",
            "connectivity": [
                "CA",
                "O",
                "OC1",
                "OC2",
                "L*"
            ],
            "bond_type": [
                "1",
                "2",
                "2",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": None,
            "c_formal_charge": 0,
            "n_formal_charge": 0,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "LM",
            "bond_type_aromatic": ["1","2","2","1","1"],
            "connectivity_type": ["S","D","D","S","S"],
            "bond_type_conjugate": ["","","","",""],
            "partial_formal_charge": 0,
        },
        "O": {
            "element": "O",
            "atom_name": "O",
            "atom_type_name": "O",
            "connectivity": [
                "C"
            ],
            "bond_type": [
                "2"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": None,
            "c_formal_charge": 0,
            "n_formal_charge": 0,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "LT",
            "bond_type_aromatic": ["2"],
            "connectivity_type": ["D"],
            "bond_type_conjugate": [""],
            "partial_formal_charge": 0,
        },
        "H1": {
            "element": "H",
            "atom_name": "H1",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": None,
            "plate": "no",
            "ff_charge": None,
            "n_formal_charge": 0,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "c_formal_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "EN",
            "bond_type_aromatic": ["1",],
            "connectivity_type": ["S",],
            "bond_type_conjugate": ["",],
            "partial_formal_charge": 0,
        },
        "H2": {
            "element": "H",
            "atom_name": "H2",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": None,
            "plate": "no",
            "ff_charge": None,
            "n_formal_charge": 0,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "c_formal_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "EN",
            "bond_type_aromatic": ["1",],
            "connectivity_type": ["S",],
            "bond_type_conjugate": ["",],
            "partial_formal_charge": 0,
        },
        "H3": {
            "element": "H",
            "atom_name": "H3",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": None,
            "plate": "no",
            "ff_charge": None,
            "n_formal_charge": 0,
            "n_ff_charge": None,
            "c_ff_charge": None,
            "c_formal_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "EN",
            "bond_type_aromatic": ["1",],
            "connectivity_type": ["S",],
            "bond_type_conjugate": ["",],
            "partial_formal_charge": 0,
        },
        "OC1": {
            "element": "O",
            "atom_name": "OC1",
            "atom_type_name": "O2",
            "connectivity": [
                "C"
            ],
            "bond_type": [
                "2"
            ],
            "formal_charge": None,
            "plate": "no",
            "ff_charge": None,
            "c_formal_charge": 0,
            "c_ff_charge": None,
            "n_ff_charge": None,
            "n_formal_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "LT",
            "bond_type_aromatic": ["2",],
            "connectivity_type": ["D",],
            "bond_type_conjugate": ["",],
            "partial_formal_charge": 0,
        },
        "OC2": {
            "element": "O",
            "atom_name": "OC2",
            "atom_type_name": "O2",
            "connectivity": [
                "C"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": None,
            "plate": "no",
            "ff_charge": None,
            "c_formal_charge": -1,
            "c_ff_charge": None,
            "n_ff_charge": None,
            "n_formal_charge": None,
            "has_ring": [],
            "has_ring_size": [],
            "has_ring_property": [],
            "local": "LT",
            "bond_type_aromatic": ["1",],
            "connectivity_type": ["S",],
            "bond_type_conjugate": ["",],
            "partial_formal_charge": 0,
        }
}

def get_current_position(position_file):
    try:
        with open(position_file, "r") as f:
            return int(f.read().strip())  
    except FileNotFoundError:
        return 000 

def save_position(position_file,position):
    with open(position_file, "w") as f:
        f.write(str(position))  

def _Generate_Char():
    characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    # 使用itertools.product生成所有3个字符的排列组合
    combinations = itertools.product(characters, repeat=3)
    position_file = "data/position.txt"
    current_position = get_current_position(position_file)
    Three_Char = next(itertools.islice(combinations, current_position, current_position+1))
    save_position(position_file,current_position+1)
    # 将元组转化为字符串并打印
    return  ''.join(Three_Char)

def Name_new_AA(Known_Char):
    while True:
        residue_name = Generate_Char()
        if residue_name not in Known_Char:
            return residue_name

def Generate_Char(used):
    characters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    # 使用itertools.product生成所有3个字符的排列组合
    for char in itertools.product(characters, repeat=3):
        char_name = "".join(char)
        if char_name not in used:
            return char_name

def register_amino_acid(molecules):
    new_molecules = []
    reg_aa = [[],[],[]]
    for molecule in molecules:
        if molecule.inchi_key not in non_AA_register[1]:
            register_name_flag = False
            char = molecule.name
            if len(char) == 3:
                if molecule.name not in non_AA_register[0]:
                    register_name_flag = True
            if not register_name_flag:
                char = Generate_Char(non_AA_register[0])
            non_AA_register[0].append(char)
            non_AA_register[1].append(molecule.inchi_key)
            non_AA_register[2].append(molecule.smiles)
            molecule.name = char
            new_molecules.append(assign_atom_name_for_amino_acid(molecule))
    text = create_amino_acid_template(new_molecules)
    dicts = {"registered":non_AA_register}
    dicts.update(non_normal_amino_acid)
    for kk,vv in text.items():
        dicts[kk] = vv
    os.system(f"mv {non_normal_amino_acid_json_f} {non_normal_amino_acid_json_f}.backup")
    
    with open(non_normal_amino_acid_json_f,'w') as outf:
        outf.write(json.dumps(dicts))

def get_ff_charge_of_special_atoms(molecule,tmp):
    if molecule.amino_acid_subtype == "alpha_aa":
        move_charge_atoms = sum([atom.ff_charge for atom in molecule.Atoms if atom.atom_name in ["H2","OC2","H0"]])
        ave_shift_charge = round(move_charge_atoms/4.0,4)
        n_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "N"][0]
        c_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "C"][0]
        o_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "OC1"][0]
        h_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "H1"][0]
        tmp["N"]["ff_charge"] = round(n_charge + ave_shift_charge,4)
        tmp["O"]["ff_charge"] = round(o_charge + ave_shift_charge,4)
        tmp["C"]["ff_charge"] = round(c_charge + ave_shift_charge,4)
        tmp["H"]["ff_charge"] = round(h_charge + ave_shift_charge,4)
        
        #tmp["N"]["ff_charge"] = round(-0.4157 - ave_shift_charge,4)
        #tmp["O"]["ff_charge"] = round(-0.5679 - ave_shift_charge,4)
        #tmp["C"]["ff_charge"] = round(0.5973 - ave_shift_charge,4)
        #tmp["H"]["ff_charge"] = round(0.2719 - ave_shift_charge,4)
        
        n_H_charge_0 = tmp["H"]["ff_charge"] - 0.05
        ave_shift_charge_n = round((0.45 - 2 * n_H_charge_0)/5.0, 4)
        tmp["O"]["n_ff_charge"] = tmp["O"]["ff_charge"]
        tmp["C"]["n_ff_charge"] = tmp["C"]["ff_charge"]
        tmp["N"]["n_formal_charge"] = 1
        tmp["N"]["n_ff_charge"] = round(tmp["N"]["ff_charge"] + 0.53 + ave_shift_charge_n,4)
        tmp["CA"]["n_ff_charge"] = round(tmp["CA"]["ff_charge"] + 0.07 + ave_shift_charge_n,4)
        tmp["H1"]["n_ff_charge"] = round(n_H_charge_0 + ave_shift_charge_n,4)
        tmp["H2"]["n_ff_charge"] = round(n_H_charge_0 + ave_shift_charge_n,4)
        tmp["H3"]["n_ff_charge"] = round(n_H_charge_0 + ave_shift_charge_n,4)
        tmp["H"]["n_ff_charge"] = round(n_H_charge_0 + ave_shift_charge_n,4)
        
        c_O_charge_0 = tmp["O"]["ff_charge"] - 0.24
        ave_shift_charge_c = round((-0.79 - c_O_charge_0)/5.0, 4)
        tmp["H"]["c_ff_charge"] = tmp["H"]["ff_charge"]
        tmp["O"]["c_ff_charge"] = round(c_O_charge_0 + ave_shift_charge_c,4)
        tmp["OC1"]["c_ff_charge"] = round(c_O_charge_0 + ave_shift_charge_c,4)
        tmp["OC2"]["c_ff_charge"] = round(c_O_charge_0 + ave_shift_charge_c,4)
        tmp["OC2"]["c_formal_charge"] = -1
        tmp["C"]["c_ff_charge"] = round(tmp["C"]["ff_charge"] + 0.2 + ave_shift_charge_c,4)
        tmp["N"]["c_ff_charge"] = round(tmp["N"]["ff_charge"] + 0.03 + ave_shift_charge_c,4)
        tmp["CA"]["c_ff_charge"] = round(tmp["CA"]["ff_charge"] - 0.20 + ave_shift_charge_c,4)
        
    elif molecule.amino_acid_subtype == "PRO_like":
        move_charge_atoms = sum([atom.ff_charge for atom in molecule.Atoms if atom.atom_name in ["H","OC2","H0"]]) 
        ave_shift_charge = round(move_charge_atoms/3.0,4)
        n_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "N"][0]
        c_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "C"][0]
        o_charge = [atom.ff_charge for atom in molecule.Atoms if atom.atom_name == "OC1"][0]
        tmp["N"]["ff_charge"] = round(n_charge + ave_shift_charge,4)
        tmp["O"]["ff_charge"] = round(o_charge + ave_shift_charge,4)
        tmp["C"]["ff_charge"] = round(c_charge + ave_shift_charge,4)
        #tmp["N"]["ff_charge"] = round(-0.2548 - ave_shift_charge,4)
        #tmp["O"]["ff_charge"] = round(-0.5748 - ave_shift_charge,4)
        #tmp["C"]["ff_charge"] = round(0.5896 - ave_shift_charge,4)
        
        ave_shift_charge_n = round(0.376/6.0, 4)
        tmp["O"]["n_ff_charge"] = round(tmp["O"]["ff_charge"] + 0.074 + ave_shift_charge_n, 4)
        tmp["C"]["n_ff_charge"] = round(tmp["C"]["ff_charge"] - 0.064 + ave_shift_charge_n, 4)
        tmp["N"]["n_formal_charge"] = 1
        tmp["N"]["n_ff_charge"] = round(tmp["N"]["ff_charge"] + 0.053 + ave_shift_charge_n,4)
        tmp["CA"]["n_ff_charge"] = round(tmp["CA"]["ff_charge"] - 0.063 + ave_shift_charge_n,4)
        tmp["H1"]["n_ff_charge"] = round(0.312 + ave_shift_charge_n,4)
        tmp["H2"]["n_ff_charge"] = round(0.312 + ave_shift_charge_n,4)

        
        c_O_charge_0 = tmp["O"]["ff_charge"] - 0.195
        ave_shift_charge_c = round((-0.7461 - c_O_charge_0)/5.0, 4)
        tmp["O"]["c_ff_charge"] = c_O_charge_0 + ave_shift_charge_c
        tmp["OC1"]["c_ff_charge"] = c_O_charge_0 + ave_shift_charge_c
        tmp["OC2"]["c_ff_charge"] = c_O_charge_0 + ave_shift_charge_c
        tmp["OC2"]["c_formal_charge"] = -1
        tmp["C"]["c_ff_charge"] = tmp["C"]["ff_charge"] + 0.0735 + ave_shift_charge_c
        tmp["N"]["c_ff_charge"] = tmp["N"]["ff_charge"] - 0.0254 + ave_shift_charge_c
        tmp["CA"]["c_ff_charge"] = tmp["CA"]["ff_charge"] - 0.107 + ave_shift_charge_c


def generate_atom_info(molecule,_atom_name,tmp):
    for atom_name in _atom_name:
        atom = _atom_name[atom_name]
        tmp[atom_name] = {
            "element": atom.elem,
            "atom_name": atom_name,
            "atom_type_name": atom.atom_type_name,
            "connectivity": [ molecule.Atoms[an].atom_name for an in atom.connectivity],
            "bond_type": atom.bond_type,
            "formal_charge": atom.formal_charge,
            "plate": atom.plate,
            "c_formal_charge": atom.formal_charge,
            "n_formal_charge": atom.formal_charge,
            "has_ring": atom.has_ring,
            "has_ring_size": atom.has_ring_size,
            "has_ring_property": atom.has_ring_property,
            "local": atom.local,
            "bond_type_aromatic": atom.bond_type_aromatic,
            "connectivity_type": atom.connectivity_type,
            "bond_type_conjugate": atom.bond_type_conjugate,
            "partial_formal_charge": atom.partial_formal_charge,
            }
        if hasattr(atom,"ff_charge"):
            tmp[atom_name]["ff_charge"] = round(atom.ff_charge,4)  #atom.ff_charge,
            tmp[atom_name]["n_ff_charge"] = round(atom.ff_charge,4) #atom.ff_charge,
            tmp[atom_name]["c_ff_charge"] = round(atom.ff_charge,4) #atom.ff_charge,
            tmp[atom_name]["h_ff_charge"] = round(atom.ff_charge,4)
        else:
            tmp[atom_name]["ff_charge"] = None  #atom.ff_charge,
            tmp[atom_name]["n_ff_charge"] = None #atom.ff_charge,
            tmp[atom_name]["c_ff_charge"] = None #atom.ff_charge,
            tmp[atom_name]["h_ff_charge"] = None
    if molecule.amino_acid_subtype == "alpha_aa":
        tmp["N"]["connectivity"] = ["H"] + tmp["N"]["connectivity"] + ["H3","R*"]
        tmp["N"]["bond_type"] = ["1"] + tmp["N"]["bond_type"] + ["1","1"]
        tmp["N"]["bond_type_aromatic"] = ["1"] + tmp["N"]["bond_type_aromatic"] + ["1","1"]
        tmp["N"]["connectivity_type"] = ["S"] + tmp["N"]["connectivity_type"] + ["S","S"]
        tmp["N"]["bond_type_conjugate"] = [""] + tmp["N"]["bond_type_conjugate"] + ["",""]
        tmp["H"] = deepcopy(tmp["H1"])
        tmp["H"]["atom_name"] = "H"
        tmp["H3"] = deepcopy(tmp["H1"])
        tmp["H3"]["atom_name"] = "H3"
        tmp["C"]["connectivity"] = ["O"] + tmp["C"]["connectivity"] + ["L*"]
        tmp["C"]["bond_type"] = ["2"] + tmp["C"]["bond_type"] + ["1"]
        tmp["C"]["bond_type_aromatic"] = ["2"] + tmp["C"]["bond_type_aromatic"] + ["1"]
        tmp["C"]["connectivity_type"] = ["D"] + tmp["C"]["connectivity_type"] + ["S"]
        tmp["C"]["bond_type_conjugate"] = [""] + tmp["C"]["bond_type_conjugate"] + [""]
        tmp["O"] = deepcopy(tmp["OC1"])
        tmp["O"]["atom_name"] = "O"
    elif molecule.amino_acid_subtype == "PRO_like":
        tmp["N"]["connectivity"] = tmp["N"]["connectivity"] + ["H2","R*"]
        ndx = tmp["N"]["connectivity"].index("H")
        tmp["N"]["connectivity"][ndx] = "H1"
        tmp["N"]["bond_type"] = tmp["N"]["bond_type"] + ["1","1"]
        tmp["N"]["bond_type_aromatic"] = tmp["N"]["bond_type_aromatic"] + ["1","1"]
        tmp["N"]["connectivity_type"] = tmp["N"]["connectivity_type"] + ["S","S"]
        tmp["N"]["bond_type_conjugate"] = tmp["N"]["bond_type_conjugate"] + ["",""]
        tmp["H1"] = deepcopy(tmp["H"])
        tmp["H1"]["atom_name"] = "H1"
        tmp["H2"] = deepcopy(tmp["H"])
        tmp["H2"]["atom_name"] = "H2"
        del tmp["H"]
        tmp["C"]["connectivity"] = ["O"] + tmp["C"]["connectivity"] + ["L*"]
        tmp["C"]["bond_type"] = ["2"] + tmp["C"]["bond_type"] + ["1"]
        tmp["C"]["bond_type_aromatic"] = ["2"] + tmp["C"]["bond_type_aromatic"] + ["1"]
        tmp["C"]["connectivity_type"] = ["D"] + tmp["C"]["connectivity_type"] + ["S"]
        tmp["C"]["bond_type_conjugate"] = [""] + tmp["C"]["bond_type_conjugate"] + [""]
        tmp["O"] = deepcopy(tmp["OC1"])
        tmp["O"]["atom_name"] = "O"
    return tmp    

def create_amino_acid_template(molecules):
    template = {}
    
    for molecule in molecules:
        _atom_name = {atom.atom_name:atom for atom in molecule.Atoms}
        if molecule.amino_acid_subtype == "alpha_aa":
            tmp = {"property":[],
                "template":{
                molecule.mole_name:[["N", "H",]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","IM"],
                f"C{molecule.mole_name}":[["N", "H",]+molecule.atom_name_order + [ "C", "OC1","OC2"],"atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type","RT"],
                f"N{molecule.mole_name}":[["N", "H1", "H2","H3"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type","LT"],
                f"H{molecule.mole_name}":[["N", "H1", "H2"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","formal_charge","h_ff_charge","connectivity","bond_type","LT"],
                }}
                
        elif molecule.amino_acid_subtype == "PRO_like":
            tmp = {"property":[],
                "template":{
                molecule.mole_name:[["N"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","IM"],
                f"C{molecule.mole_name}":[["N"]+molecule.atom_name_order + [ "C", "OC1","OC2"],"atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type","RT"],
                f"N{molecule.mole_name}":[["N", "H1", "H2"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type","LT"],
                f"N{molecule.mole_name}":[["N", "H"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","formal_charge","h_ff_charge","connectivity","bond_type","LT"],
                }}
            
        tmp = generate_atom_info(molecule,_atom_name,tmp)
        
        if hasattr(molecule.Atoms[0],"ff_charge"):
            get_ff_charge_of_special_atoms(molecule,tmp)
            zero_charge = sum([vv["ff_charge"] for at,vv in tmp.items() if at in tmp["template"][molecule.mole_name][0]])
            if zero_charge != 0.0:
                tmp["CA"]["ff_charge"] -= zero_charge
            pos_charge = 1 - sum([vv["n_ff_charge"] for at,vv in tmp.items() if at in tmp["template"][f"N{molecule.mole_name}"][0]])
            if pos_charge != 0.0:
                tmp["CA"]["n_ff_charge"] += pos_charge
            neg_charge = -1 - sum([vv["c_ff_charge"] for at,vv in tmp.items() if at in tmp["template"][f"C{molecule.mole_name}"][0]])
            if neg_charge != 0.0:
                tmp["CA"]["c_ff_charge"] += neg_charge

        template[molecule.mole_name] = tmp
    return template


def old_create_amino_acid_template(molecules):
    template = {}
    
    for molecule in molecules:
        _atom_name = {atom.atom_name:atom for atom in molecule.Atoms}
        if molecule.amino_acid_subtype == "alpha_aa":
            tmp = {"property":[],
                "template":{
                molecule.mole_name:[["N", "H",]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","IM"],
                f"C{molecule.mole_name}":[["N", "H",]+molecule.atom_name_order + [ "C", "OC1","OC2"],"atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type","RT"],
                f"N{molecule.mole_name}":[["N", "H1", "H2","H3"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type","LT"],
                }}
            for kk,vv in AA_TEMPT.items():
                tmp[kk] = deepcopy(vv)
                
        elif molecule.amino_acid_subtype == "PRO_like":
            tmp = {"property":[],
                "template":{
                molecule.mole_name:[["N"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","IM"],
                f"C{molecule.mole_name}":[["N"]+molecule.atom_name_order + [ "C", "OC1","OC2"],"atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type","RT"],
                f"N{molecule.mole_name}":[["N", "H1", "H2"]+molecule.atom_name_order + [ "C", "O"],"atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type","LT"],
                }}
            for kk,vv in AA_TEMPT.items():
                if kk not in ["H","H3"]:
                    tmp[kk] = deepcopy(vv)
                tmp["C"]["local"] = "LMC"
                tmp["C"]["connectivity_type"][0] = "eS"
                tmp["N"]["connectivity_type"] = []
            
        for atom_name in molecule.atom_name_order:
            atom = _atom_name[atom_name]
            tmp[atom_name] = {
            "element": atom.elem,
            "atom_name": atom_name,
            "atom_type_name": atom.atom_type_name,
            "connectivity": [ molecule.Atoms[an].atom_name for an in atom.connectivity],
            "bond_type": atom.bond_type,
            "formal_charge": atom.formal_charge,
            "plate": atom.plate,
            "c_formal_charge": atom.formal_charge,
            "n_formal_charge": atom.formal_charge,
            "has_ring": atom.has_ring,
            "has_ring_size": atom.has_ring_size,
            "has_ring_property": atom.has_ring_property,
            "local": atom.local,
            "bond_type_aromatic": atom.bond_type_aromatic,
            "connectivity_type": atom.connectivity_type,
            "bond_type_conjugate": atom.bond_type_conjugate,
            "partial_formal_charge": atom.partial_formal_charge,
            }
            if hasattr(atom,"ff_charge"):
                tmp[atom_name]["ff_charge"] = round(atom.ff_charge,4)  #atom.ff_charge,
                tmp[atom_name]["n_ff_charge"] = round(atom.ff_charge,4) #atom.ff_charge,
                tmp[atom_name]["c_ff_charge"] = round(atom.ff_charge,4) #atom.ff_charge,
            else:
                tmp[atom_name]["ff_charge"] = None  #atom.ff_charge,
                tmp[atom_name]["n_ff_charge"] = None #atom.ff_charge,
                tmp[atom_name]["c_ff_charge"] = None #atom.ff_charge,
        
        if hasattr(molecule.Atoms[0],"ff_charge"):
            get_ff_charge_of_special_atoms(molecule,tmp)
            zero_charge = sum([vv["ff_charge"] for at,vv in tmp.items() if at in tmp["template"][molecule.mole_name][0]])
            if zero_charge != 0.0:
                tmp["CA"]["ff_charge"] -= zero_charge
            pos_charge = 1 - sum([vv["n_ff_charge"] for at,vv in tmp.items() if at in tmp["template"][f"N{molecule.mole_name}"][0]])
            if pos_charge != 0.0:
                tmp["CA"]["n_ff_charge"] += pos_charge
            neg_charge = -1 - sum([vv["c_ff_charge"] for at,vv in tmp.items() if at in tmp["template"][f"C{molecule.mole_name}"][0]])
            if neg_charge != 0.0:
                tmp["CA"]["c_ff_charge"] += neg_charge

        template[molecule.mole_name] = tmp
    return template

def assign_atom_name_for_amino_acid(molecule):
    molecule.make_mole_as_graph()
    _label = ["-","-","A","B","G","D","E","Z","H",
                      "Q","I","K","L","M","N","X",
                      "O","P","R","S","T","U","F",
                      "C","Y","W",
                      ]
    _special_atoms = {
        "c_center":["C","C"],
        "o_center":["O","O"],
        "n_center":["N","N"],
                      }
    tmp_c = []
    for atom in molecule.Atoms:
        if atom.atom_type_name == "c_center":
            tmp = [atom.ID,[an for an in atom.connectivity if molecule.Atoms[an].elem == "O"]]
            for an in atom.connectivity:
                if molecule.Atoms[an].elem == "C":
                    n_center = [ann for ann in molecule.Atoms[an].connectivity if molecule.Atoms[ann].atom_type_name == "n_center"]
                    for ah in molecule.Atoms[n_center[0]].connectivity:
                        if molecule.Atoms[ah].elem == "H":
                            n_center.append(ah)
                    tmp.append(n_center)
            tmp_c.append(tmp)
                    
                    
            
    if len(tmp_c) == 0:
        pass
    
    center_C_arr = tmp_c[0]
    center_C = center_C_arr[0]

    molecule.Atoms[center_C].atom_name = "C"
    molecule.Atoms[center_C].atom_type_name = "C"
    fix_atom_id = [center_C]
    molecule.Atoms[center_C_arr[2][0]].atom_name = "N"
    molecule.Atoms[center_C_arr[2][0]].atom_type_name = "N"
    
    for ii,an in enumerate(center_C_arr[2][1:]):
        molecule.Atoms[an].atom_type_name = "H"
        if len(center_C_arr[2]) == 2:
            molecule.Atoms[an].atom_name = "H"
            molecule.amino_acid_subtype = "PRO_like"
        else:
            molecule.amino_acid_subtype = "alpha_aa"
            molecule.Atoms[an].atom_name = f"H{ii+1}"
    if len(center_C_arr[1]) == 1:
        molecule.Atoms[center_C_arr[1][0]].atom_name = "O"
        molecule.Atoms[center_C_arr[1][0]].atom_type_name = "O"    
    else:
        OH_c = [an for an in center_C_arr[1] if molecule.Atoms[an].atom_type_name == "OH"]
        non_OH_c = [an for an in center_C_arr[1] if molecule.Atoms[an].atom_type_name != "OH"]
        if len(OH_c) == 1:
            molecule.Atoms[non_OH_c[0]].atom_name = "OC1"
            molecule.Atoms[non_OH_c[0]].atom_type_name = "O" 
            molecule.Atoms[OH_c[0]].atom_name = "OC2"
            OH_C_H = [an for an in molecule.Atoms[OH_c[0]].connectivity if molecule.Atoms[an].elem == "H"][0]
            molecule.Atoms[OH_C_H].atom_name = "H0"
        else:
            molecule.Atoms[non_OH_c[0]].atom_name = "OC1"
            molecule.Atoms[non_OH_c[0]].atom_type_name = "O2" 
            molecule.Atoms[non_OH_c[1]].atom_name = "OC2"
            molecule.Atoms[non_OH_c[1]].atom_type_name = "O2"
            
        
    fix_atom_id += center_C_arr[1]
    fix_atom_id += center_C_arr[2]

    _order_dict = {}       
    for atom in molecule.Atoms:
        if atom.ID not in fix_atom_id and atom.elem != "H":
            dist = molecule.calc_bond_distance(center_C,atom.ID)
            if dist not in _order_dict:
                _order_dict[dist] = []
            _order_dict[dist].append(atom.ID)

    atom_name_order = []
    for kk in range(2,2+len(_order_dict)):
        vv = _order_dict[kk]
        _lab_ = _label[kk]
        for ii,an in enumerate(vv):
            if len(vv) > 1:
                lab = f"{_lab_}{str(ii+1)}"
            else:
                lab = _lab_
            molecule.Atoms[an].atom_name = f"{molecule.Atoms[an].elem}{lab}"
            atom_name_order.append(molecule.Atoms[an].atom_name)
            connect_H_atoms = [ann for ann in molecule.Atoms[an].connectivity if molecule.Atoms[ann].elem == "H"]
            if len(connect_H_atoms) > 0:
                for jj,ann in enumerate(connect_H_atoms):
                    lab_H = f"H{lab}"
                    if len(connect_H_atoms) > 1:
                        lab_H += str(jj+1)
                    molecule.Atoms[ann].atom_name = lab_H
                    atom_name_order.append(lab_H)
    molecule.atom_name_order = atom_name_order
    return molecule
