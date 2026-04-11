#from ..chemistry.constants import DEFAULT_TOPOLS_TO_MTX_ATTRIBUTES, DEFAULT_MOLECULE_TO_MTX_ATTRIBUTES
from copy import deepcopy
from .db_parse import covert_frag_to_dict
from ...utils import logger


transfer_attr= {
    "ID":"No","element":"elem","formal_charge":"formal_charge","connectivity":"connect","bond_type":"bond_type_old","coordinates":"coor",
    "atom_number":"number","atom_name":"atom_nae","local":"position","plate":"plate","score":"score",
    "bond_type_aromatic":"bond_type","connectivity_type":"bond_type_detail",
    "has_ring":"inring","has_ring_size":"ring_size","has_ring_property":"ring_prop",
    "residue":"residu","residue_id":"residu_number","chain_name":"chain_name","charge_group":"charge_group",
    "partial_formal_charge":"partial_formal_charge","point_charge":"ff_charge",
    "atom_type_name":"atom_type_name","nonb_atom_type":"nonb_atom_type","binc_atom_type":"binc_atom_type","atc_atom_type":"atc_atom_type","atom_type_ID":"atom_type_ID",
    "binc_style":"binc_style","binc_parameters":"binc_para","binc_score":"binc_score","binc_str":"binc_tag",
    "pstyle":"pstyle","parameters":"parameter","pscore":"pscore","ptag":"tag",
    "a1":"a1","a2":"a2","a3":"a3","a4":"a4","value":"value","style":"style","para":"para",
    "smiles":"smiles","inchi_key":"inchi_key","inchi":"inchi","molecule_name":"mole_name",
    "iupac_name":"iupac_name","nick_name":"nick_name","drug_name":"drug_name","internal_name":"internal_name",
    "moecule_mass":"mass","formula":"formula","heavy_atoms":"heavy_atoms","net_charge":"net_charge","multiple":"multi","element_set":"element_set",
    "function_group":"function_group","function_group_label":"function_group_label",
    "torsions":"scan_term","torsion_number":"torsion_number","constrain_term":"constrain_term","constrain_value":"constrain_value",
    "rings":"ring_dict","ring_number":"ring_number","ring_size":"ring_size","ring_property":"ring_property","ring_blocks":"ring_stru",
    "energy":"energy","dipole":"dipole","inertia":"inertia","density":"density",
    "force":"force","hessian":"hessian","frequency":"freq",
}

def transfer_data_type(dd):
    def _single(d):
        if d == "None" or d is None:
            return None
        try:
            v = float(d)
            try:
                v = int(d)
                return v
            except:
                return float(d)
        except:
            return d

    if isinstance(dd,list):
        _tmp = []
        for d in dd:
            _tmp.append(_single(d))
        return _tmp
    else:
        return _single(dd)


def reorder_topology_atoms(topol,improper_flag = False):
    
    topol = sorted(topol,key=lambda pp:pp.a1)
    return topol

def reorder_list(molecule):
    for atom in molecule.Atoms:
        if hasattr(atom,"connectivity"):
            conn_order = sorted(atom.connectivity)
            #atom.connectivity = conn_order
            for attr in ["bond_type","bond_type_aromatic"]:#,"connectivity_type","bond_type_conjugate"]:
                if hasattr(atom,attr):
                    
                    new_order = [getattr(atom,attr)[atom.connectivity.index(an)] for an in conn_order]
                    setattr(atom, attr, new_order)
            atom.connectivity = conn_order
        if hasattr(atom,"has_ring"): 
            if len(getattr(atom,"has_ring")) > 0:
                has_ring_order = sorted(atom.has_ring)
                for attr in ["has_ring_size","has_ring_property"]:
                    new_order = [getattr(atom,attr)[atom.has_ring.index(an)] for an in has_ring_order]
    #for term in ["Bonds","Angles","Dihedrals","Impropers",
    #                   "Pair12","Pair13","Pair14","Pair1n","AlteredPair",]:
    #    if hasattr(molecule,term):
    #        setattr(molecule,term,reorder_topology_atoms(getattr(molecule,term),improper_flag=True if term == "Improers" else False))

    for attr in ["torsions","ring_blocks"]:
        if hasattr(molecule,attr):
            setattr(molecule,attr,sorted(getattr(molecule,attr)))
        

LIST_TOPOL_ATTRIBUTES = ["connectivity","bond_type","coordinates","bond_type_aromatic","connectivity_type",
                        "has_ring","has_ring_size","has_ring_property",
                        "binc_style","binc_parameter","binc_score","binc_str","binc_count",
                        "fix_parameter","parameter",
                        ]

LIST_MOLECULE_ATTRIBUTES = ["function_group","constrain_term","ring_size","ring_property",
                            "dipole","inertia","force","hessian","frequency","source","steps"]
                            #"element_set"]
LIST_2D_MOLECULE_ATTRIBUTES = ["torsions","ring_blocks"]
DICT_MOLECULE_ATTRIBUTES = ["rings","count","associated_data","charge_group"]

DEFAULT_TOPOLS_TO_MTX_ATTRIBUTES = {
                                  "Atoms":
                                  [
                                  "ID","element","mass","formal_charge","connectivity","bond_type","coordinates",
                                  "atom_name","residue","residue_ID","chain_name","local","plate","score","bond_type_aromatic","connectivity_type",
                                  "has_ring","has_ring_size","has_ring_property",
                                  "atom_fg","atom_fg_id","atom_fg_tag",
                                  "atom_cluster","atom_cluster_id","atom_cluster_tag",
                                  "point_charge","esp_charge","mulliken_charge","ff_charge","partial_formal_charge",
                                  "atom_type_name","nonb_atom_type","binc_atom_type","atc_atom_type","atom_type_ID",
                                  "binc_style","binc_parameter","binc_score","binc_str","binc_count",
                                  "pstyle","parameter","fix_parameter","pscore","ptag","pcount",
                                 ],
                                  "Bonds":["a1","a2","value","style",
                                            "a1_atom_type","a2_atom_type",
                                            "a1_atom_type_used","a2_atom_type_used",
                                            "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Angles":["a1","a2","a3","value","style","value_a",
                                            "a1_atom_type","a2_atom_type","a3_atom_type",
                                            "a1_atom_type_used","a2_atom_type_used","a3_atom_type_used",
                                            "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Dihedrals":["a1","a2","a3","a4","value","value_a","style",
                                                "a1_atom_type","a2_atom_type","a3_atom_type","a4_atom_type",
                                                "a1_atom_type_used","a2_atom_type_used","a3_atom_type_used","a4_atom_type_used",
                                                "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Impropers":["a1","a2","a3","a4","value","value_a","style",
                                                "a1_atom_type","a2_atom_type","a3_atom_type","a4_atom_type",
                                                "a1_atom_type_used","a2_atom_type_used","a3_atom_type_used","a4_atom_type_used",
                                                "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Pair12":["a1","a2","value","style",
                                            "a1_atom_type","a2_atom_type",
                                            "a1_atom_type_used","a2_atom_type_used",
                                            "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Pair13":["a1","a2","value","style",
                                            "a1_atom_type","a2_atom_type",
                                            "a1_atom_type_used","a2_atom_type_used",
                                            "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Pair14":["a1","a2","value","style",
                                            "a1_atom_type","a2_atom_type",
                                            "a1_atom_type_used","a2_atom_type_used",
                                            "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "Pair1n":["a1","a2","value","style",
                                            "a1_atom_type","a2_atom_type",
                                            "a1_atom_type_used","a2_atom_type_used",
                                            "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "AlteredPair":["a1","a2","value","style",
                                                "a1_atom_type","a2_atom_type",
                                                "a1_atom_type_used","a2_atom_type_used",
                                                "atom_type_used_name","pstyle","parameter","fix_parameter","pscore","pcount","ptag"],
                                  "constrain":["a1","a2","a3","a4","fix_value","value","value_a"]
                                 }
DEFAULT_MOLECULE_ATTR_TO_CSV = [
    "smiles","inchi_key","inchi","molecule_name",
    "iupac_name","nick_name","drug_name","internal_name",
    "molecule_mass","formula","heavy_atoms","net_charge","multiple",
    "function_group","function_group_label",
    "torsions","torsion_number",
    "rings","ring_number","ring_size","ring_property","ring_blocks",
    "charge_group",
    "associated_data"
    #"energy","dipole","inertia","density",
    #"constrain_term","constrain_value",
]

DEFAULT_MOLECULE_TO_MTX_ATTRIBUTES = DEFAULT_MOLECULE_ATTR_TO_CSV + [
    "energy","force","hessian","frequency","density","steps","conform_type",
    #"elem_frag","sssr_frag","csf_frag","rsf_frag",
    #"chain_frag","scaffold_frag","tf_frag",
    #"seco_frag","sketch_frag","fragments"
    ]
TOPOLS = ["Atoms","Bonds","Angles","Dihedrals","Impropers","Pair12","Pair13","Pair14","Pair1n","AlteredPair","constrain"]
_ALL_ATTRS = {"molecule": DEFAULT_MOLECULE_TO_MTX_ATTRIBUTES,}
for term in DEFAULT_TOPOLS_TO_MTX_ATTRIBUTES:
    _ALL_ATTRS[term] = DEFAULT_TOPOLS_TO_MTX_ATTRIBUTES[term]

_DEFAULT_ATTRS = {
            "molecule": ["molecule_name","inchi_key","smiles",],
            "Atoms":["ID","element","connectivity","bond_type","formal_charge"],
            "Bonds":["a1","a2",],
            "Angles":["a1","a2","a3"],
            "Dihedrals": ["a1","a2","a3","a4"],
            "Impropers": ["a1","a2","a3","a4"],
            "Pair12":["a1","a2",],
            "Pair13":["a1","a2",],
            "Pair14":["a1","a2",],
            "Pair1n":["a1","a2",],
            "AlteredPair":["a1","a2",],
            "constrain":["a1","a2","a3","a4","fix_value"]
        }

FRAGMENT_TYPE = ["elem_frag","csf_frag", "rsf_frag", "sssr_frag", "scaffold_frag", "tf_frag"]

class MtxData:
    def __init__(self) -> None:
        pass

    def _parse(self,input_script,extra_var=None):
        script = input_script.splitlines()
        infos = self._split_script(script)
        return self._parse_script(infos)

    def _split_script(self,script):
        infos = []
        tmp = []
        for line in script:
            if line == "":
                infos.append(tmp)
                tmp = []
            else:
                tmp.append(line.strip())
        return infos

    def _parse_script(self,script):
        datas = {"Molecule":{}}
        for info in script:
            if "@" in info[0]:
                title,typ = info[0].strip().split("@")
            else:
                title = info[0].strip()
                typ = ""
            if title in DEFAULT_TOPOLS_TO_MTX_ATTRIBUTES:
                datas[title] = self._parse_topolgy(info)
            else:
                if len(info) == 1:
                    continue
                if title not in FRAGMENT_TYPE:
                    datas["Molecule"][title] = self._parse_molecule_attr(info,typ)
                else:
                    datas["Molecule"][title] = self._parse_frags(info)   
        return datas
    
    def _parse_frags(self,info):
        arr = []
        _tmp_ = info[1].split()
        attrs = []
        list_attr = []
        list_2d_attr = []
        for rr in _tmp_:
            if "@" in rr:
                s0,s1 = rr.split("@")
            else:
                s0 = rr.strip()
                s1 = ""
            attrs.append(s0)
            if s1 == "list":
                list_attr.append(s0)
            elif s1 == "list_2d":
                list_2d_attr.append(s0)
        for rr in info[2:]:
            _tmp = {}
            ss = rr.split()
            for jj,attr in enumerate(attrs):
                if ss[jj] == "None" or ss[jj] is None:
                    continue
                if attr in list_2d_attr:
                    sss = ss[jj].split("::")
                    _tmp[attr] = []
                    for s in sss[1:-1]:
                        _tmp[attr].append(transfer_data_type(s.split(":")))
                elif attr in list_attr:
                    _tmp[attr] = transfer_data_type(ss[jj].split(":")[1:-1])
                else:
                    _tmp[attr] = transfer_data_type(ss[jj])
            arr.append(_tmp)
        return arr

    def _parse_molecule_attr(self,info,typ):
        if typ == "list":
            _tmp = []
            if info[1] == ":":
                return _tmp
            for rr in info[1:]:
                _tmp.extend(transfer_data_type(rr.split(":")[1:-1]))
        elif typ == "list_2d":
            _tmp = []
            if info[1] == "::":
                return _tmp
            for rr in info[1:]:
                _tmp.append(transfer_data_type(rr.split(":")[1:-1]))
        elif typ == "dict":
            _tmp = {}
            if info[1] == "::":
                return _tmp
            for rr in info[1:]:
                tt,ss = rr.split()
                if ":" not in ss:
                    _tmp[tt] = transfer_data_type(ss.strip())
                else:
                    _tmp[tt] = transfer_data_type(ss.split(":")[1:-1])
        elif typ == "script":
            _tmp = "\n".join(info[1:])
        else:
            _tmp = transfer_data_type(info[1].strip())
        return _tmp

    def _parse_topolgy(self,info):
        arr = []
        _tmp_ = info[1].split()
        attrs = []
        list_attr = []
        for rr in _tmp_:
            if "@" in rr:
                s0,s1 = rr.split("@")
            else:
                s0 = rr.strip()
                s1 = ""
            attrs.append(s0)
            if s1 == "list":
                list_attr.append(s0)
        for rr in info[2:]:
            _tmp = {}
            ss = rr.split()
            for jj,attr in enumerate(attrs):
                if ss[jj] == "None" or ss[jj] is None:
                    continue
                if attr in ["bond_type","bond_type_aromatic"]:
                    _tmp[attr] = ss[jj].split(":")[1:-1]
                else:
                    if attr in list_attr:
                        _tmp[attr] = transfer_data_type(ss[jj].split(":")[1:-1])
                    else:
                        _tmp[attr] = transfer_data_type(ss[jj])
            arr.append(_tmp)
        return arr

    def _convert(self,molecule,extra_var=None):
        molecule_attrs = "default" if extra_var is None else extra_var
        if molecule_attrs == "default":
            molecule_attrs = _DEFAULT_ATTRS
        elif molecule_attrs == "all":
            molecule_attrs = _ALL_ATTRS
        elif molecule_attrs == "all-noncoor":
            molecule_attrs = deepcopy(_ALL_ATTRS)
            index = molecule_attrs["Atoms"].index("coordinates")
            del molecule_attrs["Atoms"][index]
            for term in TOPOLS[1:]:
                index = molecule_attrs[term].index("value")
                del molecule_attrs[term][index]

        #texts = []
        #for molecule in molecules:
        reorder_list(molecule)
        return self._convert_molecule(molecule,molecule_attrs,TOPOLS)
        #texts.append(self._convert_molecule(molecule,molecule_attrs,TOPOLS))

    def _read_file(self):
        pass

    def _convert_molecule(self,molecule,molecule_attrs,TOPOLS):
        text = "molecule_name\n%s\n\n" %(molecule.molecule_name)
        
        for tt in TOPOLS:
            if tt in molecule_attrs and molecule_attrs[tt] is not None and hasattr(molecule,tt):
                text += self._convert_molecule_topol(molecule,tt,molecule_attrs[tt])
        for attr in molecule_attrs["molecule"]:
                text += self._convert_molecule_attr(molecule,attr)
        return text

    def _convert_molecule_topol(self,molecule,topol,topol_attrs):
        text = "%s\n"%topol
        tmp = [f"{attr}@list" if attr in LIST_TOPOL_ATTRIBUTES else f"{attr}" for attr in topol_attrs]
        text += f'{" ".join(tmp)}\n'
        for term in getattr(molecule,topol):
            for attr in topol_attrs:
                value = getattr(term,attr,None)
                if value is None:
                    text += f"{value} "
                else:
                    if attr in LIST_TOPOL_ATTRIBUTES:
                        if not isinstance(value,list):
                            logger.error(f"the attribution of {attr} is not a list {value}: {attr} is setted as list")
                            text += f":{value}: "
                        else:
                            if len(value) > 0:
                                text += f':{":".join([str(vv) for vv in value])}: '
                            else:
                                text += ": "
                    else:
                        text += f'{value} '
                #text += (f'{":".join([str(vv) for vv in value])} ' if len(value) > 0 else "None ") if isinstance(value,list) else f'{value} '
            text += "\n"
        text += "\n"
        return text

    def _convert_list_2d(self,attr,value):
        text = "%s@list_2d\n"%attr
        if value is None:
            text += "\n"
            return text 
        
        if len(value) == 0:
            text += "::\n\n"
            return text
        
        for rr in value:
            text += f":{rr}:\n" if not isinstance(rr,list) \
                else f":{':'.join([str(r) for r in rr])}:\n"
        text += '\n'
        return text

    def _convert_list(self,attr,value):
        text = "%s@list\n"%attr
        if value is None:
            text += "\n"
            return text

        if len(value) == 0:
            text += ":\n\n"
            return text

        nn = int(len(value)/10)
        if nn == 0:
            text += f':{":".join([str(vv) for vv in value])}:\n\n'
            return text
        for ii in range(nn + 1):
            text += f':{":".join([str(vv) for vv in value[ii*10:(ii+1)*10]])}:\n'
        if len(value) > nn * 10:
            text += f':{":".join([str(vv) for vv in value[(nn+1)*10:]])}:\n'
            text += "\n"
            return text
        else:
            text += "\n"
            return text
        
    def _convert_dict(self,attr,value):
        text = "%s@dict\n"%attr
        if value is None:
            text += "\n"
            return text
        
        if len(value) == 0:
            text += "::\n\n"
            return text

        keys = sorted([kk for kk in value])
        for kk in keys:
            vv = value[kk]
            if isinstance(vv,list):
                text += f"{kk} :{':'.join([str(r) for r in vv])}:\n"
            else:
                text += f"{kk} {vv}\n"
        text += "\n"
        return text

    def _convert_frag(self,molecule):
        fragments = covert_frag_to_dict(molecule)
        __frag_attr = ["frag_type","index","components","offsets","connects",
                        "smiles","inchi","inchi_key","label","ef_type"]
        text = ""
        for ftype in FRAGMENT_TYPE:
            text += "%s@dict\n"%ftype
            this_frags = fragments[ftype]
            if len(this_frags) == 0:
                text += "::\n\n"
                continue
            this_attrs = [attr for attr in __frag_attr if attr in this_frags[0]]
            for attr in this_attrs:
                if attr not in ["offsets","connects"]:
                    text += f"{attr}@list " if attr == "components" else f"{attr} "
                else:
                    text += f"{attr}@list_2d "
            text += "\n"
            for frag in this_frags:
                for attr in this_attrs:
                    if attr in ["offsets","connects"]:
                        text_tmp = "::"
                        for rr in frag[attr]:
                            text_tmp += f"{':'.join([str(r) for r in rr])}::"
                        text += "%s "%text_tmp
                    elif attr in ["components"]:
                        text += f":{':'.join([str(r) for r in frag[attr]])}: "
                    else:
                        text += f"{frag[attr] } "
                text += "\n"
            text += "\n"
        return text    

    def _convert_molecule_attr(self,molecule,attr):
        if attr == "fragments":
            return self._convert_frag(molecule)
        else:
            value = getattr(molecule,attr,None)
            if attr in LIST_2D_MOLECULE_ATTRIBUTES:
                text = self._convert_list_2d(attr,value)
            elif attr in LIST_MOLECULE_ATTRIBUTES:
                text = self._convert_list(attr,value)
            elif attr in DICT_MOLECULE_ATTRIBUTES:
                text = self._convert_dict(attr,value)
            else:
                if value is None:
                    text = "%s\n\n"%attr
                else:
                    text = "%s\n%s\n\n"%(attr,value)
            return text

    
        

        
