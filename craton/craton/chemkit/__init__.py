from .stereo import Stereo
from .structure import Structure
from .moledit import MolEdit
from .interaction import InteractionModel
from .conformation import MolConformer
from .biomacromolecule import Protein
from ..utils.commons import parallel_run

import csv
import pubchempy as pcp
import re

def save_txt(molecules,fpath,parallel=True):
    if parallel:
        kwds = [{"fpath":fpath} for _ in molecules]
        tmp = parallel_run(write_txt,molecules, kwds=kwds,)
    else:
        for molecule in molecules:
            write_txt(molecule,fpath=fpath)

def atom_txt(atom,molecule):
    text = f"第{atom.ID}个原子是{atom.element}，元素符号{atom.element}，原子序号{atom.number}，原子量{atom.mass}，原子杂化形式{atom.hybrid}，"
    text += f"联接{len(atom.connectivity)}个原子，"
    for ii,an in enumerate(atom.connectivity):
        text += f"与{an}号的{molecule.Atoms[an].elem}原子，通过{atom.bond_type[ii]}联接,"
    if atom.formal_charge != 0:
        text += f"带{atom.formal_charge}电荷，"
        if atom.formal_charge != atom.partial_formal_charge:
            _N = 0
            if abs(atom.partial_formal_charge) == 0.5:
                _N = 2
            elif abs(atom.partial_formal_charge) == 0.33333:
                _N = 3
            if _N != 0:
                text += f"由于共轭，电荷分散到{_N}个原子上，原子带{atom.partial_formal_charge}部分电荷,"
    if hasattr(atom,"_interaction_model"):
        text += f"能形式{','.join(atom._interaction_model)}"
    text += "\n"
        
    return text

def chiral_txt(m):
    chiral_atom = [str(atom.ID) for atom in m.Atoms if atom.chirality_flag]
    text = ""
    if len(chiral_atom) > 0:
        text += f"含有以下手性原子：{','.join(chiral_atom)}\n"

def function_group_txt(m):
    text = "包含以下官能团或环结构：\n"
    for ni,fg in m.function_group_dicts.items():
        if fg[5] == "ring":
            text += f"第{ni}个是环结构：{fg[3]}，包含{'、'.join([str(a) for a in fg[1]])}原子，"
        else:
            text += f"第{ni}个：{fg[3]}，包含{'、'.join([str(a) for a in fg[1]])}原子，"
        for rr in fg[2]:
            text += f"通过{rr[1]}号原子与第{rr[0]}个管能团（{m.function_group_dicts[rr[0]][3]}）中的{rr[2]}号原子相连，"
        text += "\n"
    return text

def torsion_txt(m):
    _label = {"O":["氧","oxygen"],"C":["碳","carbon"],"H":["氢","hydrogen"],"N":["氮","nitrogen"],
              "S":["硫","sulfur"],"P":["磷","phoshporus"],"B":["硼","boron"],"Si":["硅","silicon"],
              "F":["氟","fluorine"],"Cl":["氯","chlorine"],"Br":["溴","bromine"],"I":["碘","iodine"],}
    if len(m.torsions) == 0:
        return ""
    text = "可旋转的键：\n"
    for torsion in m.torsions:
        atom_i = m.Atoms[torsion[1]]
        atom_j = m.Atoms[torsion[2]]
        if atom_i.atom_fg_id == atom_j.atom_fg_id:
            fg_i = atom_i.atom_fg_id
            text += f"第{fg_i}个管能团（{m.function_group_dicts[fg_i][3]}）中，{atom_i.ID}号{_label[atom_i.elem][0]}原子和{atom_j.ID}号{_label[atom_j.elem][0]}原子形成可旋转键\n"
        else:
            fg_i = atom_i.atom_fg_id
            fg_j = atom_j.atom_fg_id
            text += f"第{fg_i}个管能团（{m.function_group_dicts[fg_i][3]}）的{atom_i.ID}号{_label[atom_i.elem][0]}原子，与第{fg_j}个管能团（{m.function_group_dicts[fg_j][3]}）的{atom_j.ID}号{_label[atom_j.elem][0]}原子形成可旋转键\n"
    return text    


def write_txt(m,fpath=".",idx=None):
    text = f"分子smiles: {m.smiles} \n"
    text += f"分子inchi key: {m.inchi_key} \n"
    text += f"分子式: {m.formula} \n"
    text += f"分子量: {m.mass} \n\n"
    text += "原子结构部分：\n"
    for atom in m.Atoms:
        text +=  atom_txt(atom,m)
    text += function_group_txt(m)
    text += torsion_txt(m)
    with open(f"{fpath}/{m.inchi_key}.txt",'w') as outf:
        outf.write(text)
    if idx is not None:
        return "OK", idx

def run_get_moleinfo_from_pubchem(input,input_type="smiles",output_type="name"):
    _name_ = {
        "smiles": "canonical_smiles",
        "nickname": "traditional_iupac_name",
        "iupac_name": "systematic_iupac_name",
        "name": "systematic_iupac_name",
        "weight": "exact_mass",
        "formula": "molecular_formula",
        "mass": "exact_mass",
        "inchikey":"standard_inchikey",
        "inchi_key": "standard_inchikey",
        "inchi": "standard_inchi"
    }
    
    def get_cas_no(smiles):
        cas = []
        ddd = pcp.get_synonyms(smiles, 'smiles')
        for dd in ddd:
            for d in dd["Synonym"]:
                match = re.match('(\d{2,7}-\d\d-\d)',d)
                
                if match:
                    cas.append(match.group(1))
        if len(cas) > 0:
            return ":".join(cas)
        else:
            return None
    
    output_type = output_type.lower()
    
    ccc_ = pcp.get_compounds(input, input_type)
    if len(ccc_) > 0:
        ccc = ccc_[0]
        smiles = ccc.canonical_smiles
        if output_type == "cas_no":
            return get_cas_no(smiles)
        else:
            dicts = {}
            for ur in ccc._record["props"]:
                label_ = ur["urn"]["label"].lower()
                label_ = label_.replace(" ","_")
                _name = ur["urn"]["name"].lower() if "name" in ur["urn"] else ""
                if _name == "":
                    label_name = label_
                else:
                    label_name = _name + "_" + label_
                for aa,bb in ur["value"].items():
                    dicts[label_name] = bb
                #dicts[label_name] = ur["value"]["sval"] if "sval" in ur["value"] else ur["value"]["fval"] if "fval" in ur["value"] else ur["value"]["ival"]
            if output_type in _name_:
                if _name_[output_type] in dicts:
                    return dicts[_name_[output_type]]
                else:
                    return None
            else:
                if output_type in _name_:
                    return dicts[output_type]
                else:
                    return None
    else:
        return None
    
    
def get_moleinfo_from_pubchem(input,input_type="smiles",output_type="name",file_flag=False):
    if file_flag:
        if input[-4:] == ".csv":
            raws = list(csv.reader(open(input)))
            ndx = raws[0].index(input_type)
            input_arr = [rr[ndx] for rr in raws[1:]]
        elif input[-4:] == ".txt":
            with open(input) as inf:
                input_arr = [rr.strip() for rr in inf.readlines()]
    else:
        input_arr = [input]
    
    datas = []
    for rr in input_arr:
        datas.append(run_get_moleinfo_from_pubchem(rr,input_type=input_type,output_type=output_type))
    
    if file_flag:
        results = [["id",input_type,output_type]]
        for ii,ra in enumerate(input_arr):
            results.append([ii,ra,datas[ii]])
        with open("results.csv",'w') as outf:
            writer = csv.writer(outf)
            writer.writerows(results)
        return None
    else:
        return input_arr[0], datas[0] 
