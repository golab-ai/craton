import itertools
from ..chem import FormatMolecule as FM
from ..chemkit.biomacromolecule.protein_utils import assign_force_field
from ..chemkit.structure.structure import protein_ring_and_charge_group
from ..chem.molecule import Molecule
from ..chem.atom import Atom
from ..chem import create_3D
import itertools
from copy import deepcopy
import json
from .. import CRATON_CONFIGURE

amino_acid_json_f = f'{CRATON_CONFIGURE["path"]["template"]}/amino_acid.json'
#from .amino_acid_template import amino_acid
template_amino_acid = json.loads(open(amino_acid_json_f).read())

normal_amino_acid = ["ALA","GLY","SER","THR","LEU","ILE","VAL","ASN","GLN","ARG","HIS","TRP","PHE","TYR","GLU","ASP","LYS","PRO","CYS","MET"]
rna = ["A","U","G","C"]
dna = ["DA","DT","DG","DC"]

def create_template_molecule(template,residue_name):
    atoms = []
    
    infos = template["template"][residue_name]
    idx_dict = {name:ii for ii,name in enumerate(infos[0])}
    
    for ii,name in enumerate(infos[0]):
        atoms.append(Atom("aa"))
        atoms[-1].ID = ii
        atoms[-1].residue = residue_name
        atoms[-1].residue_ID = 0
        
        atoms[-1].atom_name = name
        atoms[-1].element = template[name]["element"]
        atoms[-1].plate = template[name]["plate"]
        
        atoms[-1].atom_type_name = template[name][infos[1]]
        atoms[-1].formal_charge = template[name][infos[2]]
        atoms[-1].ff_charge = template[name][infos[3]]
        if infos[6] == "LT":
            connectivitys = [an for an in  template[name][infos[4]] if an in infos[0] + ["L*"]]
        elif infos[6] == "RT":
            connectivitys = [an for an in  template[name][infos[4]] if an in infos[0] + ["R*"]]
        elif infos[6] == "EN":
            connectivitys = [an for an in  template[name][infos[4]] if an in infos[0]]
        else:
            connectivitys = [an for an in  template[name][infos[4]] if an in infos[0] + ["L*","R*"]]
        atoms[-1].connectivity = [idx_dict[an] if an not in ["L*","R*"] else an for an in  connectivitys]
        atoms[-1].bond_type = [template[name][infos[5]][template[name][infos[4]].index(an)] for an in connectivitys]
        atoms[-1].coordinates = [0.0,0.0,0.0]
        
    rm = Molecule("residue")
    rm.Atoms = atoms
    return rm

def create_template_total_molecule(template,residue_name):
    atoms = []
    
    infos = template["template"][residue_name]
    idx_dict = {name:ii for ii,name in enumerate(infos[0])}
    
    for ii,name in enumerate(infos[0]):
        atoms.append(Atom("aa"))
        atoms[-1].ID = ii
        atoms[-1].residue = residue_name
        atoms[-1].residue_ID = 0
        
        atoms[-1].atom_name = name
        atoms[-1].element = template[name]["element"]
        atoms[-1].plate = template[name]["plate"]
        
        atoms[-1].atom_type_name = template[name][infos[1]]
        atoms[-1].formal_charge = template[name][infos[2]]
        atoms[-1].ff_charge = template[name][infos[3]]
        connectivitys = [an for an in  template[name][infos[4]] if an in infos[0]]
        atoms[-1].connectivity = [idx_dict[an] for an in  connectivitys]
        atoms[-1].bond_type = [template[name][infos[5]][template[name][infos[4]].index(an)] for an in connectivitys]
        atoms[-1].coordinates = [0.0,0.0,0.0]
        
    rm = Molecule("residue")
    rm.Atoms = atoms
    return rm

class MoleculeAssembly:
    def __init__(self) -> None:
        pass
    
    @staticmethod
    def combine_residue(arr,residues):
        mole_name = "-".join(arr)
        Atoms = []
        L_terminal = -1
        for ii,rname in enumerate(arr):
            shift_n = len(Atoms)
            residue = deepcopy(residues[rname])
            for atom in residue.Atoms:
                atom.residue_ID = ii
                #atom["residue"] = rname
                atom.ID += shift_n
                atom.connectivity = [conn+shift_n if conn not in ["L*","R*"] else conn for conn in atom.connectivity]
                if atom.connectivity[-1] == "L*":
                    this_L_terminal = atom.ID
                elif atom.connectivity[-1] == "R*":
                    atom.connectivity[-1] = L_terminal
                    Atoms[L_terminal].connectivity[-1] = atom.ID
            
                Atoms.append(atom)
            L_terminal = this_L_terminal
            
        molecule = Molecule("peptide")
        molecule.Atoms = Atoms
        molecule.mole_name = mole_name
        
        molecule.create_topols()
        molecule.create_improper(create_method="atom_type")
        molecule.create_intra_nonbond_macromole()
        molecule = protein_ring_and_charge_group(molecule)
        molecule = assign_force_field(molecule)
        
        molecule.steps = ["structure","atom type","force field"]
        #molecule = Stru._basic_structure_analyze(molecule)[0]

        create_3D(molecule)
        return molecule
    
    @staticmethod
    def peptide_gen(n,left_cap="ACE",right_cap="NME",terminal_flag=True,templates=None,):
        if templates is None:
            templates = normal_amino_acid
        template_unit = {subaa:create_template_molecule(template_amino_acid[aa],subaa) for aa in template_amino_acid for subaa in template_amino_acid[aa]["template"]}

        _tmp = []
        _tmpn = []
        _tmpc = []
        for templ in templates:
            _tmp.extend([aa for aa in template_amino_acid[templ]["template"] if len(aa)== 3])
            _tmpn.extend([aa for aa in template_amino_acid[templ]["template"] if len(aa) == 4 and aa[0] == "N"])
            _tmpc.extend([aa for aa in template_amino_acid[templ]["template"] if len(aa) == 4 and aa[0] == "C"])

        if n == 1:
            arrs = [_tmp]
        else:
            left_terminal = []
            right_terminal = []
            non_terminal = [[] for __ in range(n-2)]
            
            if left_cap is not None:
                left_terminal.extend(_tmp)
            if right_cap is not None:
                right_terminal.extend(_tmp)
            
            
            if terminal_flag:
                left_terminal.extend(_tmpn)
                right_terminal.extend(_tmpc)
            for ii in range(n-2):
                non_terminal[ii].extend(_tmp)
            arrs = [left_terminal] + non_terminal + [right_terminal]
        assembly = []
        for item in itertools.product(*arrs):
            item = list(item)
            if len(item[0]) == 3:
                item = [left_cap] + item
            if len(item[-1]) == 3:
                item = item + [right_cap]
            assembly.append(item)
        
        molecules = []
        for arr in assembly:
            molecules.append(MoleculeAssembly.combine_residue(arr,template_unit))
        
        FM._convert(molecules,otype="mtx",ofilename=None,opath=None,extra_var="all",parallel=True)
        FM._convert(molecules,otype="mol",ofilename=None,opath=None,extra_var="all",parallel=True)
        return [mole.smiles for mole in molecules]

    @staticmethod
    def dnarna_gen(n,templates=None):
        if templates is None:
            templates = rna
        elif templates == "rna":
            templates = rna
        elif templates == "dna":
            templates = dna
        
        if templates is None:
            templates = normal_amino_acid
        template_unit = {subaa:create_template_molecule(template_amino_acid[aa],subaa) for aa in template_amino_acid for subaa in template_amino_acid[aa]["template"]}

        _tmp = []
        _tmp3 = []
        _tmp5 = []
        _tmpend = []
        for templ in templates:
            _tmp3.extend([aa for aa in template_amino_acid[templ]["template"] if aa[-1] == "3"])
            _tmp5.extend([aa for aa in template_amino_acid[templ]["template"] if aa[-1] == "5"])
            _tmpend.extend([aa for aa in template_amino_acid[templ]["template"] if aa[-1] == "N"])
            _tmp.extend([aa for aa in template_amino_acid[templ]["template"] if aa[-1] not in ["3","5","N"]])
        
        if n == 1:
            molecules = []
            for rname in _tmpend:
                molecule = create_template_total_molecule(template_amino_acid[rname[:-1]],rname)
                molecule.mole_name = rname
                molecule.create_topols()
                molecule.create_improper(create_method="atom_type")
                molecule.create_intra_nonbond_macromole()
                molecule = protein_ring_and_charge_group(molecule)
                molecule = assign_force_field(molecule)
                
                molecule.steps = ["structure","atom type","force field"]
                create_3D(molecule) 
                molecules.append(molecule)
        else:
            arrs = [_tmp5] + [_tmp for __ in range(n-2)] + [_tmp3]

            assembly = []
            for item in itertools.product(*arrs):
                item = list(item)
                assembly.append(item)
        
            molecules = []
            for arr in assembly:
                molecules.append(MoleculeAssembly.combine_residue(arr,template_unit))
        
        
        
        FM._convert(molecules,otype="mtx",ofilename=None,opath=None,extra_var="all",parallel=True)
        FM._convert(molecules,otype="mol",ofilename=None,opath=None,extra_var="all",parallel=True)
        #return [mole.smiles for mole in molecules]
            
    @staticmethod
    def polymer_gen(n,templates=None):
        
        pass        
    
        

class AminoAcid:
    def __init__(self,rc,left_cap=None,right_cap=None,template=None,output_dir="."):
        self.rc = rc
        self.left_cap = [left_cap] if not isinstance(left_cap,list) else left_cap
        self.right_cap = [right_cap] if not isinstance(right_cap,list) else right_cap
        if template is None:
            template = {}
        self.template = template
        self.residues = {}
        self.caps = {}
        self.residues = residues

    def combine_residue(self):
        from ..chemkit import Structure as Stru
        residue_pull = [[rname for rname in self.residues if rname in normal_residues] for ii in range(self.rc)]
        for ii in range(self.rc):
            if ii in self.template:
                residue_pull[ii] = self.template[ii]
        residue_pull = [[cap for cap in self.left_cap]] + residue_pull
        residue_pull += [[cap for cap in self.right_cap]]
        molecules = []
        for _arr_ in itertools.product(*residue_pull):
            arr = list(_arr_)
            if arr[0] is None:
                if arr[1] not in ["LYN","CYM","ASH","GLH"]:
                    del arr[0]
                    arr[0] = f"N{arr[0]}"
                else:
                    continue
                    #arr[0] = "NME"
            if arr[-1] is None:
                if arr[-2] not in ["LYN","CYM","ASH","GLH"]:
                    del arr[-1]
                    arr[-1] = f"C{arr[-1]}"
                else:
                    continue
                    #arr[-1] = "ACE"
            mole_name = "-".join(arr)
            Atoms = []
            L_terminal = -1
            for ii,rname in enumerate(arr):
                shift_n = len(Atoms)
                residue = deepcopy(self.residues[rname])
                for atom in residue:
                    atom["residue_ID"] = ii
                    #atom["residue"] = rname
                    atom["ID"] += shift_n
                    atom["connectivity"] = [conn+shift_n if conn not in ["L*","R*"] else conn for conn in atom["connectivity"]]

                    if atom["connectivity"][-1] == "L*":
                        this_L_terminal = atom["ID"]
                    elif atom["connectivity"][-1] == "R*":
                        atom["connectivity"][-1] = L_terminal
                        Atoms[L_terminal]["connectivity"][-1] = atom["ID"]
                
                    Atoms.append(atom)
                L_terminal = this_L_terminal
            
            keys = list(Atoms[0].keys())
            
            
            datas = {kk:[] for kk in keys}
            for an in Atoms:
                for kk in keys:
                    datas[kk].append(an[kk])
            datas["elements"] = datas["element"]
            del datas["element"]
            molecule = FM._create_molecule(datas)[0]
            molecule.mole_name = mole_name
            molecule.create_improper(create_method="atom_type")
            molecule = Stru._basic_structure_analyze(molecule)[0]
            
            create_3D(molecule)
            molecules.append(molecule)
            
            
        FM._convert(molecules,otype="mtx",ofilename=None,opath=None,extra_var="all",parallel=True)
        FM._convert(molecules,otype="mol",ofilename=None,opath=None,extra_var="all",parallel=True)
        return [mole.smiles for mole in molecules]

def _run():
    __residue_label = {
        "normal":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"formal_charge",
                  "ff_charge":"ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "c_terminal":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"c_formal_charge",
                      "ff_charge":"c_ff_charge","plate":"plate","connectivity_extra":["R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "n_terminal":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"n_formal_charge",
                      "ff_charge":"n_ff_charge","plate":"plate","connectivity_extra":["L*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},

        "HID":{"element":"element","connectivity":"connectivity","bond_type":"HID_bond_type","formal_charge":"HID_formal_charge",
               "ff_charge":"HID_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "HIP":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"HIP_formal_charge",
               "ff_charge":"HIP_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "LYN":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"LYN_formal_charge",
               "ff_charge":"LYN_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "CYM":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"CYM_formal_charge",
               "ff_charge":"CYM_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "CYX":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"CYX_formal_charge",
               "ff_charge":"CYX_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "ASH":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"ASH_formal_charge",
               "ff_charge":"ASH_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "GLH":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"GLH_formal_charge",
               "ff_charge":"GLH_ff_charge","plate":"plate","connectivity_extra":["L*","R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "CHID":{"element":"element","connectivity":"connectivity","bond_type":"HID_bond_type","formal_charge":"HID_c_formal_charge",
               "ff_charge":"HID_c_ff_charge","plate":"plate","connectivity_extra":["R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "NHID":{"element":"element","connectivity":"connectivity","bond_type":"HID_bond_type","formal_charge":"HID_n_formal_charge",
               "ff_charge":"HID_n_ff_charge","plate":"plate","connectivity_extra":["L*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "CHIP":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"HIP_c_formal_charge",
               "ff_charge":"HIP_c_ff_charge","plate":"plate","connectivity_extra":["R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "NHIP":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"HIP_n_formal_charge",
               "ff_charge":"HIP_n_ff_charge","plate":"plate","connectivity_extra":["L*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "CCYX":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"CYX_c_formal_charge",
               "ff_charge":"CYX_c_ff_charge","plate":"plate","connectivity_extra":["R*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        "NCYX":{"element":"element","connectivity":"connectivity","bond_type":"bond_type","formal_charge":"CYX_n_formal_charge",
               "ff_charge":"CYX_n_ff_charge","plate":"plate","connectivity_extra":["L*"],"atom_name":"atom_name","atom_type_name":"atom_type_name"},
        }
    residues = {}

    for kname,vv in amino_acid.items():
        tmp = []
        rname = vv[1]
        r_atoms = template_amino_acid[rname]
        _ids = {kk:ii for ii,kk in enumerate(vv[2])}
        if kname in __residue_label:
            attr_label = __residue_label[kname]
        else:
            if len(kname) == 4:
                if kname[0] == "C":
                    attr_label = __residue_label["c_terminal"]
                elif kname[0] == "N":
                    attr_label = __residue_label["n_terminal"]
            else:
                attr_label = __residue_label["normal"]
        for an in vv[2]:

            _cc_ = cc = [[cn,ii] for ii,cn in enumerate(r_atoms[an][attr_label["connectivity"]]) if cn in vv[2]+attr_label["connectivity_extra"]]
            #atom = deepcopy(r_atoms[an])
            atom = {}
            atom["ID"] = _ids[an]
            atom["residue"] = rname
            atom["element"] = r_atoms[an][attr_label["element"]]
            atom["atom_name"] = r_atoms[an][attr_label["atom_name"]]
            atom["atom_type_name"] = r_atoms[an][attr_label["atom_type_name"]]
            atom["connectivity"] = [_ids[rr[0]] if rr[0] in _ids else rr[0] for rr in _cc_]
            atom["bond_type"] = [r_atoms[an][attr_label["bond_type"]][rr[1]] for rr in _cc_]
            atom["formal_charge"] = r_atoms[an][attr_label["formal_charge"]]
            atom["ff_charge"] = r_atoms[an][attr_label["ff_charge"]]
            atom["plate"] = r_atoms[an][attr_label["plate"]]
            atom["coordinates"] = [0.0,0.0,0.0]
            tmp.append(atom)
        residues[kname] = tmp
    normal_residues = ["ALA","GLY","SER","THR","LEU","ILE","VAL","ASN","GLN","ARG",
                    "HIS","TRP","PHE","TYR","GLU","ASP","LYS","PRO","CYS","MET",
                    "HID","HIP","LYN","CYM","ASH","GLH"]
    ##normal_residues = ["SER","THR","PHE","TYR","ASN","GLN","ARG",]
