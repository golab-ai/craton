import functools
import itertools
import json
import os
import re
import networkx as nx
from typing import Iterable, List, Any
from copy import deepcopy
from .. import CRATON_CONFIGURE
DEFAULT_TYPING_RULE = CRATON_CONFIGURE["ForceFieldSetting"]["DEFAULT_TYPING_FILE"]
from ..utils import logger

def _convert_equtable_(inputf,outputf):
    with open(inputf) as inf:
        tmp = [rr.strip().split() for rr in inf.readlines()]
        convert_dict  = {rr[0]:rr[1] for rr in tmp}
    text = "#######################################################\n"
    text += "###### The following section for EquivalTable #########\n"
    text += "EquivalTable:\n"
    text += "#Type      ATC      NONB       BINC      BOND      A_C      A_S      T_C      T_S      O_C      O_S \n"

    for aa,bb in convert_dict.items():
        text += f"{aa}      {bb}      {bb}      {bb}      {bb}      {bb}      {bb}      {bb}      {bb}      {bb}      {bb}\n"

    with open(outputf,'w') as outf:
        outf.write(text)

def check_kuohao(es, kuohao=["[", "]"]):
    n = len(es)
    n1 = es.count(kuohao[0])
    n2 = es.count(kuohao[1])
    if n1 != n2:
        raise Exception("Unmatched brackets")
    else:
        kuohao_n = [i for i in range(0, n) if es[i] in kuohao]
        node_arr = []
        json_str = ""
        for i in range(0, len(kuohao_n) - 1):
            if kuohao_n[i + 1] - kuohao_n[i] != 1:
                node_arr.append(es[kuohao_n[i] + 1 : kuohao_n[i + 1]])
                if es[kuohao_n[i]] == kuohao[0]:
                    json_str = json_str + ",[" + '"%s"' % es[kuohao_n[i] + 1 : kuohao_n[i + 1]]
                elif es[kuohao_n[i]] == kuohao[1]:
                    json_str += "]"
            else:
                json_str += "]"
        json_str = json_str[1:] + "]"
        tree = json.loads(json_str)
        return tree

class AtomTypeTXT:
    _bondType = {"-": "1", "~": "ar", "=": "2", "#": "3", "*": "-1"}  # aromatic bond type  # any bond type

    def __init__(self):

        self.At_Name = ""
        self.Expression = ""
        self.Hybrid = ""
        self.Arom = ""
        self.Charge = ""
        self.At_Nodes: {str: [str, [int], [float]]} = {}  # {index: [element, [partner_index], [partner_bond_type]]}
        self.At_Attrs: {str: {str: [str]}} = {}  # {index: {'coor': ['coor', ...], 'elem':['elem', ...]}
        self.At_Connect: [int, int, float] = []  # [[atom1, atom2, bond_type]]

        self._at_n = 0
        self._m = 0

    def read_script(self, script):
        ss = script[0].strip().split()
        self.At_Name = ss[1]
        self.Hybrid = ss[2]
        if "a" in ss[3]:
            self.Arom = "a"
            if len(ss[3]) == 2:
                self.Charge = ss[3][1]
        else:
            if ss[3] != "0":
                self.Charge = ss[3]
        self.Expression = script[1].strip()
        express = "[" + script[1].strip() + "]"

        tree = check_kuohao(express)

        self.interpret_tree(tree, 0)

        for line in script[2:]:
            self.interpret_attrs(line)

        self.create_connect()

        for aa, bb in self.At_Attrs.items():
            if "elem" not in bb.keys():
                self.At_Attrs[aa]["elem"] = [self.At_Nodes[aa][0]]
        
        extend_depth = self.get_extend_depth([int(an) for an in self.At_Attrs.keys()],self.At_Connect)
        return {"At_Name":self.At_Name,"Expression":self.Expression,"Hybrid":self.Hybrid, "Arom":self.Arom,"Charge":self.Charge,
                "At_Attrs":self.At_Attrs,"At_Connect":self.At_Connect,"Extend_depth":extend_depth}

    def get_extend_depth(self,atoms,connect):
        G = nx.Graph()
        G.add_nodes_from(atoms)
        G.add_edges_from([(bond[0], bond[1]) for bond in connect])
        extend_depth = 0
        for an in atoms:
            if an != 0:
                nn = len(nx.shortest_path(G, source=0, target=an))
                if nn > extend_depth:
                    extend_depth = nn 
        return extend_depth

    def interpret_attrs(self, scr):
        def __is_float(str):
            s = str.split(".")
            if len(s) > 2:
                return False
            else:
                for si in s:
                    if not si.isdigit():
                        return False
                return True

        string = scr.strip().split(":")
        nnn = str(int(string[0].split()[1]) - 1)
        self.At_Attrs[nnn] = {}
        attrs = string[1].split(",")
        for attr in attrs:
            kkk = attr.split("=")
            if kkk[0].strip() in ["coor","ring"]:
                self.At_Attrs[nnn][kkk[0].strip()] = []
                for aaa in kkk[1].strip().split():
                    if aaa.isdigit():
                        self.At_Attrs[nnn][kkk[0].strip()].append(int(aaa))
                    else:
                        self.At_Attrs[nnn][kkk[0].strip()].append(aaa)
            elif kkk[0].strip() in ["charge"]:
                self.At_Attrs[nnn][kkk[0].strip()] = []
                for aaa in kkk[1].strip().split():
                    try:
                        self.At_Attrs[nnn][kkk[0].strip()].append(float(aaa))
                    except:
                        self.At_Attrs[nnn][kkk[0].strip()].append(aaa)
            else:
                self.At_Attrs[nnn][kkk[0].strip()] = kkk[1].strip().split()

    def interpret_tree(self, r, m):
        if self._at_n == 0:
            self.At_Nodes["0"] = [r[0], [], []]
        else:
            self.At_Nodes[str(self._at_n)] = [r[0][1:], [], []]
            self.At_Nodes[str(m)][1].append(self._at_n)
            self.At_Nodes[str(m)][2].append(self._bondType[r[0][0]])
            m = self._at_n

        self._at_n += 1
        for r0 in r[1:]:
            self.interpret_tree(r0, m)

    def create_connect(self):
        for a, b in self.At_Nodes.items():
            if len(b[1]) != 0:
                for i in range(len(b[1])):
                    self.At_Connect.append([int(a), b[1][i], b[2][i]])


@functools.lru_cache(maxsize=None)
def _parse_decorated_atom_type(string):
    regex = re.compile(r"^([^~@^]+)(.*)$")
    match = re.search(regex, string)
    if match is None:
        raise Exception("Invalid atom type name: " + string)
    at_base = match.group(1)
    tag = match.group(2)

    if tag.endswith("^"):
        tag_conj = "^"
        tag = tag.rstrip("^")
    else:
        tag_conj = ""

    if tag.endswith("@"):
        tag_endo = "@"
        tag = tag.rstrip("@")
    else:
        tag_endo = ""

    if tag:
        regex = re.compile(r"^(~[0-9]+=*)$")
        match = re.search(regex, tag)
        if match is None:
            raise Exception("Invalid atom type name: " + string)

    return at_base, tag, tag_endo, tag_conj

class DecoratedAtomType:
    def __init__(self, string):
        self.at_base, self.tag_ring, self.tag_endo, self.tag_conj = _parse_decorated_atom_type(string)

    def is_subset_of(self, type2):
        if self.at_base != type2.at_base:
            return False
        if self.tag_conj != type2.tag_conj:
            return False
        if type2.tag_ring != "" and self.tag_ring != type2.tag_ring:
            return False
        if type2.tag_endo != "" and self.tag_endo != type2.tag_endo:
            return False

        return True

    @property
    def tag(self):
        return f"{self.tag_ring}{self.tag_endo}{self.tag_conj}"

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.at_base} ({self.tag_ring}) ({self.tag_endo}) ({self.tag_conj})>"

    def __str__(self):
        return f"{self.at_base}{self.tag_ring}{self.tag_endo}{self.tag_conj}"

    @staticmethod
    def is_pair_subset_of(pair1, pair2):
        deco1, deco2 = pair1
        ref1, ref2 = pair2
        return (
            deco1.is_subset_of(ref1)
            and deco2.is_subset_of(ref2)
            or deco2.is_subset_of(ref1)
            and deco1.is_subset_of(ref2)
        )

class TypingDefine:
    def __init__(self,atversion=None,style=None,ctx=None,check_at_types=True,convert_flag=False):
        self.at_types = {}

        self.procedure_tree = []
        self._at_arrs = []
        self._procedure_tree_script = ""


        if atversion is None:
            self.atversion = DEFAULT_TYPING_RULE
        else:
            self.atversion = atversion
        if style is None:
            self.style = "JSON" if atversion.find(".json") != -1 else "TXT"
        else:
            self.style = style

        self._at_file = os.path.basename(atversion)

        if self.style == "JSON":
            self.read_atom_type_json(self.atversion,convert_flag=convert_flag)
        else:
            self.read_atom_type_txt(self.atversion,convert_flag=convert_flag)


        if check_at_types:
            self._check_at_types()

    def _create_equ_table(self,at,arr=None):  
        if arr is None:
            arr = [at] * 11
        data = {"at":at,
            "atc":arr[1],
            "nonb":arr[2],
            "binc":arr[3],
            "bond":arr[4],
            "a_c":arr[5],
            "a_s":arr[6],
            "d_c":arr[7],
            "d_s":arr[8],
            "i_c":arr[9],
            "i_s":arr[10]
            }
        return data

    def _check_at_types(self):
        # if there is no equivalent table, create a default one
        #if len(self.equ_table) == 0:
        #    self.equ_table = {at: TypeEquivalence(at) for at in self.at_types}
        ats_tree = []
        for vv in self.at_types.values():
            ats_tree.extend(vv["children"])
        #ats_tree = [at0 for vv in self.at_types.values() for at0 in vv["children"] ]
        ats_define = [at for at,define in self.at_types.items()]

        _at_tree_undefined = set(ats_tree).difference(set(ats_define))
        if _at_tree_undefined:
            logger.error(f"Undefined atom type appearing in hierarchical tree: {_at_tree_undefined}")

        # warn on charged atom types
        _charge_missing = [at for at in self.at_types if at[-1] in ("+", "-") and at not in self.at_types["root"]["primitive_charge"]]
        if self.at_types["root"]["primitive_charge"] and _charge_missing:
            logger.error(f"PrimitiveCharge missing for atom types: {_charge_missing}")

        # validate the integrity of typing rule
        ats_equ = []
        for at ,equ in self.at_types["root"]["equ_table"].items():
            for __,aa in equ.items():
                ats_equ.append(aa)

        _at_equ_undefined = set(ats_equ) - set(ats_define)
        if _at_equ_undefined:
            logger.error(f"Undefined atom type appearing in equivalent table: {_at_equ_undefined}")

        def __check_attrs_value(attr):
            err = ""
            for aa,bb in attr.items():
                if aa in ["coor","ring"]:
                    for bbb in bb:
                        if not isinstance(bbb,int):
                            if bbb != "?":
                                err += "%s is not int, " %aa
                elif aa in ["charge"]:
                    for bbb in bb:
                        if not isinstance(bbb,float):
                            if bbb != "?":
                                err += "%s is not float, " %aa
            return err


        attrs_setting_error = []
        for at,define in self.at_types.items():
            for __,attr in define["At_Attrs"].items():
                err = __check_attrs_value(attr)
                if err != "":
                    logger.error(f"{at} {err}")
                    attrs_setting_error.append(at)

        if _at_equ_undefined or _at_tree_undefined or _charge_missing or attrs_setting_error:
            raise Exception("Invalid typing rule")

    def _special_section_generate(self):
        self.origin_at_types = deepcopy(self.at_types)


        def _is_at_in_list(at, lst):
            for i in lst:
                if at == i:
                    return True
                if i.endswith("*") and at.startswith(i[:-1]):
                    return True
            return False

        # expand wildcards
        #self.ats_has_improper.extend([at for at in self.at_types if _is_at_in_list(at, self.ats_has_improper)])
        #self.ats_linear_center.extend([at for at in self.at_types if _is_at_in_list(at, self.ats_linear_center)])
        #self.ats_degenerated_improper.extend(
        #    [at for at in self.at_types if _is_at_in_list(at, self.ats_degenerated_improper)]
        #)
        self.at_types["root"]["has_improper"].extend([at for at in self.at_types if _is_at_in_list(at, self.at_types["root"]["has_improper"])])
        self.at_types["root"]["linear_center"].extend([at for at in self.at_types if _is_at_in_list(at, self.at_types["root"]["linear_center"])])
        self.at_types["root"]["degenerated_improper"].extend(
            [at for at in self.at_types if _is_at_in_list(at, self.at_types["root"]["degenerated_improper"])]
        )

        # planer dihedral will be degenerated
        self.at_types["root"]["degenerated_dihedral"] = list(dict.fromkeys(self.at_types["root"]["degenerated_dihedral"]
                                                                            + self.at_types["root"]["planer_dihedral"]))
        # parse decorations
        def _parse_decorations(lst):
            new = []
            for i in lst[:]:
                if "$" in i:
                    at1, at2 = i.split("$")
                    deco1 = DecoratedAtomType(at1)
                    deco2 = DecoratedAtomType(at2)
                    new.append((deco1, deco2))
                else:
                    deco = DecoratedAtomType(i)
                    new.append(deco)
            lst.clear()
            lst.extend(new)

        _parse_decorations(self.at_types["root"]["conjugation"])
        _parse_decorations(self.at_types["root"]["linear_center"])
        _parse_decorations(self.at_types["root"]["degenerated_improper"])
        _parse_decorations(self.at_types["root"]["degenerated_dihedral"])
        _parse_decorations(self.at_types["root"]["planer_dihedral"])

    def read_atom_type_json(self,atf,convert_flag=False):
        """
        root node是一个特殊的Node, 它是根节点，很多特殊定义的性质都内含在root node中，例如：
            equ_table: {}, 等价表
            rings_special: [3,4,5,6,ar]，表示环的原子类型中要添加的环前缀
            primitive_charge: {"o_1-":0.5}, 带电原子电荷分配的数值
            has_improper: [], 能形成Improper的原子类型
            conjugation: [], 可形成共轭的原子类型
            linear_center: [], 线性的原子类型
            degenerated_improper:[], improper参数的i_s变成通配类型，如“Y”
            degenerated_dihedral: [], 二面角参数的d_s原子变成通配类型，如“Y”
            planer_dihedral: [], 平面的二面角
            AtomTypeDecoChar:{"RING":"~","ENDO":"@","CONJ":"^","DEG":"Y"}, 特殊修饰符号
        以上可以在root node里统一设置，也可在具体的原子类型中设置，
        还可以在具体原子类型中通过decorated的方法设置修饰后的原子类型,如decorated_has_improper
        """
        at_types = json.loads(open(atf).read())
        _special_section = {"primitive_charge","has_improper","conjugation","linear_center",
                            "degenerated_improper","degenerated_dihedral","planer_dihedral","same"}
        if "equ_table" not in at_types["root"]:
            at_types["root"]["equ_table"] = {}
        for attr in _special_section:
            if attr not in at_types["root"]:
                if attr in ["primitive_charge"]:
                    at_types["root"][attr] = {}
                else:
                    at_types["root"][attr] = []
        for at,define in at_types.items():
            if at not in ["root"]:
                for attr in _special_section:
                    if attr in define:
                        if define[attr]:
                            if attr in ["primitive_charge","same"]:
                                at_types["root"][attr][at] = define[attr]
                            else:
                                at_types["root"][attr].append(at)
                    if f"decorated_{attr}" in define:
                        if define[f"decorated_{attr}"]:
                            if attr in ["primitive_charge","same"]:
                                at_types["root"][attr] = dict(at_types["root"][attr],**define[attr])
                                #at_types["root"][attr][at] = define[attr]
                            else:
                                at_types["root"][attr].extend(define[f"decorated_{attr}"])
                if "equ_table" not in define or define["equ_table"] is None:
                    define["equ_table"] = self._create_equ_table(at)
                at_types["root"]["equ_table"][at] = define["equ_table"]
        self.at_types = at_types
        if convert_flag:
            return
        self._special_section_generate()

    def _find_at_in_tree(self, at, ats, parent_ats):
        flag = False
        if ats[0] == at:
            self.parent = parent_ats[0]
            if len(ats) > 1:
                self.children = []
                for rr in ats[1:]:
                    self.children.append(rr[0])
            if len(parent_ats) > 2:
                self.sibling = []
                for rr in parent_ats[1:]:
                    if rr[0] != at:
                        self.sibling.append(rr[0])
            flag = True
        for ats0 in ats[1:]:
            if flag:
                break
            self._find_at_in_tree(at, ats0, ats)

    def _find_parent_children(self, at, tree):
        self.parent = "root"
        self.children = []
        self.sibling = []
        self._find_at_in_tree(at, tree, tree)

    def _create_parent_children(self):
        for at,define in self.at_types.items():
            self._find_parent_children(at, self.procedure_tree)
            define["parent"] = self.parent
            define["children"] = self.children
            define["sibling"] = self.sibling

    def _create_at_types_txt(self):
        for i in range(len(self._at_arrs)):
            name = self._at_arrs[i][0].strip().split()[1]
            ATT = AtomTypeTXT()
            self.at_types[name] = ATT.read_script(self._at_arrs[i])
            #self.at_types[name] = atom_type_txt(self._at_arrs[i])

    def _create_produce_tree_txt(self):
        procedure_tree_script = "(root" + self._procedure_tree_script + ")"
        self.procedure_tree = check_kuohao(procedure_tree_script, kuohao=["(", ")"])

    def read_atom_type_txt(self, atf, convert_flag=False,ctx=None):
        """
        读取txt文件，转换成json。后续可能不再维护！！！！
        """
        local = ctx is None or not ctx.remote
        if local:
            with open(atf) as f:
                lines = f.read().splitlines()
        else:
            pass

        self.at_types = {"root":{"At_Name":"root","Expression":"*","Hybrid":"*","elem":["?"],
                                 "At_Attrs":{0:{"coor":["?"],"elem":["?"]}},
                                 "At_Connect":[],
                                 "has_improper":[],"equ_table":{},"primitive_charge":{},
                                 "rings_special":[],"conjugation":[],"linear_center":[],
                                 "degenerated_improper":[],"degenerated_dihedral":[],
                                 "planer_dihedral":[],"same":{},
                                 "AtomTypeDecoChar":{"RING":"~","ENDO":"@","CONJ":"^","DEG":"Y"},
                                 "_at_file":self._at_file,
                                 }
                        }
        sections = {
            "Define":None,
            "ProcedureTree":None,
            "EquivalTable":None,
            "PrimitiveCharge":None,
            "Improper": self.at_types["root"]["has_improper"],
            "Conjugation": self.at_types["root"]["conjugation"],
            "JoinedRing": self.at_types["root"]["rings_special"],
            "LinearCenter": self.at_types["root"]["linear_center"],
            "DegeneratedImproper": self.at_types["root"]["degenerated_improper"],
            "DegeneratedDihedral": self.at_types["root"]["degenerated_dihedral"],
            "PlanerDihedral": self.at_types["root"]["planer_dihedral"],
            "Same":self.at_types["root"]["same"],
        }

        lines = [line.split("!")[0].strip() for line in lines if not line.startswith("#")]
        lines = [line for line in lines if line != ""]

        flag = "no"
        for i in range(0, len(lines)):
            if flag == "no":
                s = lines[i].split()[0].rstrip(":")
                if s in sections:
                    flag = s
                elif s[0].isupper():
                    raise Exception("Invalid section in type define: %s" % s)
                if s == "Define":
                    arr = [lines[i]]
            elif flag == "Define":
                if lines[i].startswith("End"):
                    flag = "no"
                    self._at_arrs.append(arr)
                else:
                    arr.append(lines[i])
            elif flag == "ProcedureTree":
                if lines[i].startswith("End"):
                    flag = "no"
                else:
                    self._procedure_tree_script += lines[i].strip()
            elif flag == "EquivalTable":
                if lines[i].startswith("End"):
                    flag = "no"
                else:
                    arr = lines[i].strip().split()
                    self.at_types["root"]["equ_table"][arr[0]] = self._create_equ_table(arr[0],arr)
            elif flag == "PrimitiveCharge":
                if lines[i].startswith("End"):
                    flag = "no"
                else:
                    at, charge = lines[i].strip().split()
                    self.at_types["root"]["primitive_charge"][at] = float(charge)
            elif flag in sections:
                if lines[i].startswith("End"):
                    flag = "no"
                else:
                    for a in lines[i].strip().split():
                        sections[flag].append(a)

        self._create_at_types_txt()
        self._create_produce_tree_txt()
        self._create_parent_children()
        for at,eqt in self.at_types["root"]["equ_table"].items():
            self.at_types[at]["equ_table"] = eqt
        if convert_flag:
            self.at_types["root"]["equ_table"] = {}
            return
        self._special_section_generate()

class TypingEngine:
    def __init__(self, atversion=None, extend_depth=8, ctx=None,check_at_types=True):

        if atversion is None:
            atversion = DEFAULT_TYPING_RULE

        self._at_file = os.path.basename(atversion)

        ats = TypingDefine(atversion, ctx=ctx, check_at_types=check_at_types)
        self.at_types = ats.at_types  # Read from local file or remote path

        self.origin_at_types = ats.origin_at_types

    def assign_mole_at(
        self,
        molecule,
        atoms_arr=None,
        ignore_existing=False,
        ignore_ff_existing=False,
    ):
        if "force field" in molecule.steps and not ignore_existing:
            return molecule
        if not ignore_existing and "atom type" in molecule.steps:
            return molecule
        if atoms_arr is None:
            atoms_arr = [i for i in range(len(molecule.Atoms))]
        for a in atoms_arr:
            self.search_at_tree(a, molecule.Atoms)
            molecule.Atoms[a].ats_tree = self.ats_tree
            if self.ats_tree[-1] in self.at_types["root"]["same"]:
                this_at = self.at_types["root"]["same"][self.ats_tree[-1]]
                #molecule.Atoms[a].atom_type_name = self.at_types["root"]["same"][self.ats_tree[-1]]
            else:
                this_at = self.ats_tree[-1]
            molecule.Atoms[a].atom_type_name = this_at
            molecule.Atoms[a].atom_type_name_arom = self.at_types[this_at]["Arom"]
            molecule.Atoms[a].atom_type_name_charge = self.at_types[this_at]["Charge"]
            molecule.Atoms[a].atom_type_name_hybrid = self.at_types[this_at]["Hybrid"]
            molecule.Atoms[a].plate = "yes" if this_at in self.at_types["root"]["has_improper"] else "no"
        return molecule

    def search_at_tree(self, a, Atoms):
        self.ats_tree, self.find_flag = [], False
        for at in self.at_types["root"]["children"]:
            self.find_at(a, at, Atoms)
            if len(self.ats_tree) != 0:
                break

    def find_at(self, a, at, Atoms):
        if not self.find_flag and self.match_at_a(at, a, Atoms):
            self.ats_tree.append(at)
            for at0 in self.at_types[at]["children"]:
                self.find_at(a, at0, Atoms)
            self.find_flag = True

    def match_at_a(self, at, a, Atoms):
        this_connect = self.at_types[at]["At_Connect"]
        this_attrs = self.at_types[at]["At_Attrs"]
        extend_depth = self.at_types[at]["Extend_depth"]
        if not self.match_at_attrs(Atoms[a], this_attrs["0"]):
            return False

        match_atom = [[] for i in range(len(this_attrs))]
        match_atom[0].append(a)

        def get_candidate_atom(a, n, Atoms):
            arr = [a]
            tmp = set()
            for i in range(0, n):
                res = set(arr) - tmp
                tmp = set(arr)
                for aa in res:
                    arr += getattr(Atoms[aa], "connect")
            return set(arr)

        candidate_atom = get_candidate_atom(a, extend_depth + 1, Atoms)
        
        for i in range(1, len(this_attrs)):
            for aa in candidate_atom:
                if self.match_at_attrs(Atoms[aa], this_attrs[str(i)]):
                    match_atom[i].append(aa)
            if len(match_atom[i]) == 0:
                return False
        
        for arr in itertools.product(*match_atom):
            if len(set(arr)) != len(this_attrs):
                continue
            for conn in this_connect:
                aa, bb, cc = arr[conn[0]], arr[conn[1]], conn[2]
                #if bb not in Atoms[aa].connect or (
                #    cc != "-1" and (Atoms[aa].bond_type[Atoms[aa].connect.index(bb)] != cc
                #                    and Atoms[aa].bond_type_aromatic[Atoms[aa].connect.index(bb)] != cc )
                #):
                #    break
                if bb not in Atoms[aa].connect or (
                    cc != "-1" and Atoms[aa].bond_type_aromatic[Atoms[aa].connect.index(bb)] != cc ):
                    break
            else:
                return True
        return False

    def match_at_attrs(self, atom, attrs):
        for attr, value in attrs.items():
            if "?" in value:
                continue
            if attr == "coor":
                if len(atom.connect) not in value:
                    return False
            elif attr == "formal_charge":
                if atom.formal_charge not in value:
                    return False
            elif attr == "charge":
                # cfl@2023-08-18 match multi charge sets
                __flag = False
                for cc in value:
                    if abs(atom.primitive_formal_charge - cc) <= 0.01:
                        __flag = True
                        break
                if not __flag:
                    return False
            elif attr == "elem":
                if atom.elem not in value:
                    return False
            elif attr == "ring":
                if len(atom.ring_size) == 0 or min(atom.ring_size) not in value:
                    return False
            elif attr == "arom":
                # consider the atom as aromatic even if the smallest ring is not
                if all(prop not in value for prop in atom.ring_prop):
                    return False
            elif attr == "conjugation":
                if not set(atom.bond_type_conjugate).intersection({"J", "eJ", "fj", "bJ", "j"}):
                    return False
            elif attr == "local":
                if atom.local not in value:
                    return False
            elif attr == "res":
                if atom.residu not in value:
                    return False
            elif attr == "hybrid":
                if atom.hybrid not in value:
                    return False
            elif attr == "btype":
                if len(set(atom.connectivity_type) & set(value)) == 0:
                    return False
        return True

    def get_electron_pairs(self, at_string):
        """
        Returns 4 for sp3, 3 for sp2, 2 for sp1
        """
        deco = DecoratedAtomType(at_string)
        
        #for d in self.ats_degenerated_dihedral:

        for d in self.at_types["root"]["degenerated_dihedral"]:
            if isinstance(d, DecoratedAtomType) and deco.is_subset_of(d):
                return 3

        _elem_ep = {
            "C": 0,
            "N": 1,
            "P": 1,
            "O": 2,
            "S": 2,
        }
        at = deco.at_base
        try:
            conn = int(at.split("_")[1][0])
        except:
            conn = 8
        elem = at[:2].rstrip("_").capitalize()
        ep = _elem_ep.get(elem, 0) + conn
        #ep = _elem_ep.get(elem, 0) + int(at[2])
        if at.endswith("-"):
            ep += 1
        elif at.endswith("+"):
            ep -= 1
        return ep


class AtomTypeDecoChar:
    RING = "~"
    ENDO = "@"
    CONJ = "^"
    DEG = "Y"

#default_typer = TypingEngine()

def assign_atom_type_to_term(
    molecule,
    this_terms=None,
    atom_arr=None,
    atom_type_rule=None, 
    equtable_flag=True,
    create_improper=True,
    raise_invalid_charge=False,
    raise_unsupported_coordination=False,
    ignore_existing = False,
):
    """
    将原子类型转换，并分配到各个term中，
    输入：
        molecule: 分子对象
        equ_table: 原子类型等价表
        rings_special: 环结构
        ats_conjugation: 共轭结构
        ats_linear_center: 线性结构
        ats_degenerated_dihedral: 降解的二面角
        ats_degenerated_improper: 降解的improper
    输出：
        为每个作用项中的每个原子分配到可以用来提取力场参数的原子类型
    """
    def _check(raise_invalid_charge,raise_unsupported_coordination,primitive_charge_flag=False):
    
        
        if primitive_charge_flag :
            try:
                ff_charge_total = sum(atom.ff_charge_base for atom in molecule.Atoms)
            except:
                ff_charge_total = sum(atom.ff_charge for atom in molecule.Atoms)
            if abs(molecule.net_charge - ff_charge_total) > 1e-4:
            
                msg = (
                    f"Invalid DEF `{atom_type_rule['_at_file']}` for molecule {molecule.name}: "
                    f"Primitive charge ({ff_charge_total}) not equal to net charge ({molecule.net_charge})"
                )
                if raise_invalid_charge:
                    raise Exception(msg)
                else:
                    (msg)

        unsupported = [atom.atom_type_name for atom in molecule.Atoms if atom.atom_type_name.endswith("_")]
        if unsupported:
            msg = f"Incomplete DEF `{atom_type_rule['_at_file']}` for molecule {molecule.name}: Unsupported coordination {unsupported}"
            if raise_unsupported_coordination:
                raise Exception(msg)
            else:
                logger.error(msg)

    def _in_same_smallest_ring(m, *ids):
        try:
            smallest_ring_sizes = [min(m.Atoms[i].ring_size) for i in ids]
        except ValueError:
            return False
        if len(set(smallest_ring_sizes)) > 1:
            return False
        for k, v in m.ring_dict.items():
            ring = v[:-1]
            if len(ring) == smallest_ring_sizes[0] and all(i in ring for i in ids):
                return True
        return False

    if atom_arr is None:
        atom_arr = [atom for atom in molecule.Atoms]
    else:
        atom_arr = [molecule.Atoms[ii] for ii in atom_arr]
    if "force field" in molecule.steps:
        return molecule

    if not ignore_existing and "atom type" in molecule.steps:
        return molecule

    if atom_type_rule is None:
        atom_type_rule = {
                        "rings_special":None,
                        "equ_table":None,
                        "conjugation":None,
                        "linear_center":None,
                        "degenerated_dihedral":None,
                        "degenerated_improper":None,
                        "primitive_charge":None
                        }
    
    if not equtable_flag:
        atom_type_rule["equ_table"] = None

    if create_improper:
        molecule.create_improper(create_method="mix")

    if atom_type_rule["primitive_charge"]:
        for atom in atom_arr:
            atom.ff_charge_base = atom_type_rule["primitive_charge"].get(atom.atom_type_name, 0)
    else:
        for atom in atom_arr:
            atom.ff_charge_base = 0.0

    _check(raise_invalid_charge,raise_unsupported_coordination,primitive_charge_flag=True if atom_type_rule["primitive_charge"] else False)


    # decorate atom types with ring information
    if atom_type_rule["rings_special"]:
        for atom in atom_arr:
            if not atom.ring_size:
                continue
            min_size = min(atom.ring_size)
            if str(min_size) not in atom_type_rule["rings_special"]:
                for size, prop in zip(atom.ring_size, atom.ring_prop):
                    if size == min_size and prop in atom_type_rule["rings_special"]:
                        break
                else:
                    continue
            tag = atom_type_rule["AtomTypeDecoChar"]["RING"] + str(min_size)
            if {"bD", "eD"}.intersection(atom.connectivity_type):
                tag += "="
            elif atom.atom_type_name == "c_3n":
                # special cases: n_3im+, n_3gu+
                neighbors = [molecule.Atoms[i] for i in atom.connect]
                ats_set = {DecoratedAtomType(neigh.atom_type_name).at_base for neigh in neighbors}
                if {"n_3im+", "n_3gu+"}.intersection(ats_set):
                    tag += "="
            atom.atom_type_name = atom.atom_type_name + tag
    if atom_type_rule["equ_table"]:
        for atom in atom_arr:
            deco = DecoratedAtomType(atom.atom_type_name)
            atom.nonb_atom_type = atom_type_rule["equ_table"][deco.at_base]["nonb"] + deco.tag
            atom.atc_atom_type = atom_type_rule["equ_table"][deco.at_base]["atc"] + deco.tag
            atom.binc_atom_type = atom_type_rule["equ_table"][deco.at_base]["binc"] + deco.tag

    else:
        for atom in atom_arr:
            atom.nonb_atom_type = atom.atom_type_name
            atom.atc_atom_type = atom.atom_type_name
            atom.binc_atom_type = atom.atom_type_name
    _label = {
        "Bonds": ["bond", "bond"],
        "Angles": ["a_s", "a_c", "a_s"],
        "Dihedrals": ["d_s", "d_c", "d_c", "d_s"],
        "Impropers": ["i_c", "i_s", "i_s", "i_s"],
        "Pair1n": ["nonb", "nonb"],
        "Pair12": ["nonb", "nonb"],
        "Pair13": ["nonb", "nonb"],
        "Pair14": ["nonb", "nonb"],
    }
    if this_terms is None:
        total_terms = {typ:getattr(molecule,typ,[]) for typ in _label}
    else:
        total_terms = {typ:[getattr(molecule,typ)[ii] for ii in this_terms[typ]] if typ in this_terms else [] for typ in _label}
        
        
    for typ,terms in total_terms.items():
        for term in terms:
            for i in range(1, 5):
                key = "a" + str(i)
                if not hasattr(term, key):
                    continue
                idx = getattr(term, key)
                attr_at = key + "_atom_type"
                attr_at_used = key + "_atom_type_used"
                at = molecule.Atoms[idx].atom_type_name
                setattr(term, attr_at, at)
                if atom_type_rule["equ_table"] is None:
                    at_used = at
                else:
                    deco = DecoratedAtomType(at)
                    attr_equ = _label[typ][i - 1]
                    at_used = atom_type_rule["equ_table"][deco.at_base][attr_equ]
                    #at_used = getattr(equ_table[deco.at_base], attr_equ)
                    if attr_equ not in ("i_c" "i_s"):
                        at_used += deco.tag
                setattr(term, attr_at_used, at_used)

    # ring
    if atom_type_rule["rings_special"]:
        for term in total_terms["Bonds"]:
            if (
                DecoratedAtomType(term.a1_atom_type_used).tag_ring
                and DecoratedAtomType(term.a2_atom_type_used).tag_ring
            ):
                if _in_same_smallest_ring(molecule, term.a1, term.a2):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
        for term in total_terms["Angles"]:
            if (
                DecoratedAtomType(term.a1_atom_type_used).tag_ring
                and DecoratedAtomType(term.a2_atom_type_used).tag_ring
                and DecoratedAtomType(term.a3_atom_type_used).tag_ring
            ):
                if _in_same_smallest_ring(molecule, term.a1, term.a2, term.a3):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
        for term in total_terms["Dihedrals"]:
            if (
                DecoratedAtomType(term.a2_atom_type_used).tag_ring
                and DecoratedAtomType(term.a3_atom_type_used).tag_ring
            ):
                if _in_same_smallest_ring(molecule, term.a2, term.a3):
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
    # conjugation
    if atom_type_rule["conjugation"]:
        for term in total_terms["Bonds"]:
            a1, a2 = term.a1, term.a2
            deco1 = DecoratedAtomType(term.a1_atom_type_used)
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            if molecule.Atoms[a1].bond_type[molecule.connectivity[a1].index(a2)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco1, deco2), pair) for pair in atom_type_rule["conjugation"]):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
        for term in total_terms["Angles"]:
            a1, a2, a3 = term.a1, term.a2, term.a3
            deco1 = DecoratedAtomType(term.a1_atom_type_used)
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            if molecule.Atoms[a1].bond_type[molecule.connectivity[a1].index(a2)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco1, deco2), pair) for pair in atom_type_rule["conjugation"]):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
            if molecule.Atoms[a2].bond_type[molecule.connectivity[a2].index(a3)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in atom_type_rule["conjugation"]):
                    term.a2_atom_type_used = (
                        term.a2_atom_type_used.rstrip(atom_type_rule["AtomTypeDecoChar"]["CONJ"]) + atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    )
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
        for term in total_terms["Dihedrals"]:
            a2, a3 = term.a2, term.a3
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            if molecule.Atoms[a2].bond_type[molecule.connectivity[a2].index(a3)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in atom_type_rule["conjugation"]):
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
    # linear dihedrals
    if atom_type_rule["linear_center"]:
        for term in total_terms["Dihedrals"]:
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            if any(deco2.is_subset_of(deco) for deco in atom_type_rule["linear_center"]) or any(
                deco3.is_subset_of(deco) for deco in atom_type_rule["linear_center"]
            ):
                term.a1_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a2_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a3_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
    # dihedral wildcard
    if atom_type_rule["degenerated_dihedral"]:
        deco_types = [i for i in atom_type_rule["degenerated_dihedral"] if isinstance(i, DecoratedAtomType)]
        deco_pairs = [i for i in atom_type_rule["degenerated_dihedral"] if isinstance(i, Iterable)]
        for term in total_terms["Dihedrals"]:
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            # degenerate endocyclic torsions
            if deco2.tag_endo and deco3.tag_endo:
                term.a1_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                continue
            try:
                if any(deco2.is_subset_of(deco) for deco in deco_types) and any(
                    deco3.is_subset_of(deco) for deco in deco_types
                ):
                    raise StopIteration
                if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in deco_pairs):
                    raise StopIteration
            except StopIteration:
                term.a1_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                continue
    # improper wildcard
    if atom_type_rule["degenerated_improper"]:
        for term in total_terms["Impropers"]:
            deco1 = DecoratedAtomType(term.a1_atom_type_used)
            if any(deco1.is_subset_of(deco) for deco in atom_type_rule["degenerated_improper"]):
                term.a2_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a3_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
    return molecule


def old_assign_atom_type_to_term(
    molecule,
    atom_type_rule=None, 
    equtable_flag=True,
    create_improper=True,
    raise_invalid_charge=False,
    raise_unsupported_coordination=False,
    ignore_existing = False,
):
    """
    将原子类型转换，并分配到各个term中，
    输入：
        molecule: 分子对象
        equ_table: 原子类型等价表
        rings_special: 环结构
        ats_conjugation: 共轭结构
        ats_linear_center: 线性结构
        ats_degenerated_dihedral: 降解的二面角
        ats_degenerated_improper: 降解的improper
    输出：
        为每个作用项中的每个原子分配到可以用来提取力场参数的原子类型
    """
    def _check(raise_invalid_charge,raise_unsupported_coordination,primitive_charge_flag=False):
    
        
        if primitive_charge_flag :
            ff_charge_total = sum(atom.ff_charge_base for atom in molecule.Atoms)
            if abs(molecule.net_charge - ff_charge_total) > 1e-4:
            
                msg = (
                    f"Invalid DEF `{atom_type_rule['_at_file']}` for molecule {molecule.name} {molecule.smiles}: "
                    f"Primitive charge ({ff_charge_total}) not equal to net charge ({molecule.net_charge})"
                )
                if raise_invalid_charge:
                    raise Exception(msg)
                else:
                    (msg)

        unsupported = [atom.atom_type_name for atom in molecule.Atoms if atom.atom_type_name.endswith("_")]
        if unsupported:
            msg = f"Incomplete DEF `{atom_type_rule['_at_file']}` for molecule {molecule.name}: Unsupported coordination {unsupported}"
            if raise_unsupported_coordination:
                raise Exception(msg)
            else:
                logger.error(msg)

    def _in_same_smallest_ring(m, *ids):
        try:
            smallest_ring_sizes = [min(m.Atoms[i].ring_size) for i in ids]
        except ValueError:
            return False
        if len(set(smallest_ring_sizes)) > 1:
            return False
        for k, v in m.ring_dict.items():
            ring = v[:-1]
            if len(ring) == smallest_ring_sizes[0] and all(i in ring for i in ids):
                return True
        return False

    if "force field" in molecule.steps:
        return molecule

    if not ignore_existing and "atom type" in molecule.steps:
        return molecule

    if atom_type_rule is None:
        atom_type_rule = {
                        "rings_special":None,
                        "equ_table":None,
                        "conjugation":None,
                        "linear_center":None,
                        "degenerated_dihedral":None,
                        "degenerated_improper":None,
                        "primitive_charge":None
                        }
    
    if not equtable_flag:
        atom_type_rule["equ_table"] = None

    if create_improper:
        molecule.create_improper(create_method="mix")

    if atom_type_rule["primitive_charge"]:
        for atom in molecule.Atoms:
            atom.ff_charge_base = atom_type_rule["primitive_charge"].get(atom.atom_type_name, 0)
    else:
        for atom in molecule.Atoms:
            atom.ff_charge_base = 0.0

    _check(raise_invalid_charge,raise_unsupported_coordination,primitive_charge_flag=True if atom_type_rule["primitive_charge"] else False)


    # decorate atom types with ring information
    if atom_type_rule["rings_special"]:
        for atom in molecule.Atoms:
            if not atom.ring_size:
                continue
            min_size = min(atom.ring_size)
            if str(min_size) not in atom_type_rule["rings_special"]:
                for size, prop in zip(atom.ring_size, atom.ring_prop):
                    if size == min_size and prop in atom_type_rule["rings_special"]:
                        break
                else:
                    continue
            tag = atom_type_rule["AtomTypeDecoChar"]["RING"] + str(min_size)
            if {"bD", "eD"}.intersection(atom.connectivity_type):
                tag += "="
            elif atom.atom_type_name == "c_3n":
                # special cases: n_3im+, n_3gu+
                neighbors = [molecule.Atoms[i] for i in atom.connect]
                ats_set = {DecoratedAtomType(neigh.atom_type_name).at_base for neigh in neighbors}
                if {"n_3im+", "n_3gu+"}.intersection(ats_set):
                    tag += "="
            atom.atom_type_name = atom.atom_type_name + tag
    if atom_type_rule["equ_table"]:
        for atom in molecule.Atoms:
            deco = DecoratedAtomType(atom.atom_type_name)
            atom.nonb_atom_type = atom_type_rule["equ_table"][deco.at_base]["nonb"] + deco.tag
            atom.atc_atom_type = atom_type_rule["equ_table"][deco.at_base]["atc"] + deco.tag
            atom.binc_atom_type = atom_type_rule["equ_table"][deco.at_base]["binc"] + deco.tag

    else:
        for atom in molecule.Atoms:
            atom.nonb_atom_type = atom.atom_type_name
            atom.atc_atom_type = atom.atom_type_name
            atom.binc_atom_type = atom.atom_type_name
    _label = {
        "Bonds": ["bond", "bond"],
        "Angles": ["a_s", "a_c", "a_s"],
        "Dihedrals": ["d_s", "d_c", "d_c", "d_s"],
        "Impropers": ["i_c", "i_s", "i_s", "i_s"],
        "Pair1n": ["nonb", "nonb"],
        "Pair12": ["nonb", "nonb"],
        "Pair13": ["nonb", "nonb"],
        "Pair14": ["nonb", "nonb"],
    }
    for typ in _label:
        for term in getattr(molecule, typ, []):
            for i in range(1, 5):
                key = "a" + str(i)
                if not hasattr(term, key):
                    continue
                idx = getattr(term, key)
                attr_at = key + "_atom_type"
                attr_at_used = key + "_atom_type_used"
                at = molecule.Atoms[idx].atom_type_name
                setattr(term, attr_at, at)
                if atom_type_rule["equ_table"] is None:
                    at_used = at
                else:
                    deco = DecoratedAtomType(at)
                    attr_equ = _label[typ][i - 1]
                    at_used = atom_type_rule["equ_table"][deco.at_base][attr_equ]
                    #at_used = getattr(equ_table[deco.at_base], attr_equ)
                    if attr_equ not in ("i_c" "i_s"):
                        at_used += deco.tag
                setattr(term, attr_at_used, at_used)

    # ring
    if atom_type_rule["rings_special"]:
        for term in getattr(molecule, "Bonds", []):
            if (
                DecoratedAtomType(term.a1_atom_type_used).tag_ring
                and DecoratedAtomType(term.a2_atom_type_used).tag_ring
            ):
                if _in_same_smallest_ring(molecule, term.a1, term.a2):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
        for term in getattr(molecule, "Angles", []):
            if (
                DecoratedAtomType(term.a1_atom_type_used).tag_ring
                and DecoratedAtomType(term.a2_atom_type_used).tag_ring
                and DecoratedAtomType(term.a3_atom_type_used).tag_ring
            ):
                if _in_same_smallest_ring(molecule, term.a1, term.a2, term.a3):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
        for term in getattr(molecule, "Dihedrals", []):
            if (
                DecoratedAtomType(term.a2_atom_type_used).tag_ring
                and DecoratedAtomType(term.a3_atom_type_used).tag_ring
            ):
                if _in_same_smallest_ring(molecule, term.a2, term.a3):
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["ENDO"]
    # conjugation
    if atom_type_rule["conjugation"]:
        for term in getattr(molecule, "Bonds", []):
            a1, a2 = term.a1, term.a2
            deco1 = DecoratedAtomType(term.a1_atom_type_used)
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            if molecule.Atoms[a1].bond_type[molecule.connectivity[a1].index(a2)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco1, deco2), pair) for pair in atom_type_rule["conjugation"]):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
        for term in getattr(molecule, "Angles", []):
            a1, a2, a3 = term.a1, term.a2, term.a3
            deco1 = DecoratedAtomType(term.a1_atom_type_used)
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            if molecule.Atoms[a1].bond_type[molecule.connectivity[a1].index(a2)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco1, deco2), pair) for pair in atom_type_rule["conjugation"]):
                    term.a1_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
            if molecule.Atoms[a2].bond_type[molecule.connectivity[a2].index(a3)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in atom_type_rule["conjugation"]):
                    term.a2_atom_type_used = (
                        term.a2_atom_type_used.rstrip(atom_type_rule["AtomTypeDecoChar"]["CONJ"]) + atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    )
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
        for term in getattr(molecule, "Dihedrals", []):
            a2, a3 = term.a2, term.a3
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            if molecule.Atoms[a2].bond_type[molecule.connectivity[a2].index(a3)] == "1":
                if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in atom_type_rule["conjugation"]):
                    term.a2_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
                    term.a3_atom_type_used += atom_type_rule["AtomTypeDecoChar"]["CONJ"]
    # linear dihedrals
    if atom_type_rule["linear_center"]:
        for term in getattr(molecule, "Dihedrals", []):
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            if any(deco2.is_subset_of(deco) for deco in atom_type_rule["linear_center"]) or any(
                deco3.is_subset_of(deco) for deco in atom_type_rule["linear_center"]
            ):
                term.a1_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a2_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a3_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
    # dihedral wildcard
    if atom_type_rule["degenerated_dihedral"]:
        deco_types = [i for i in atom_type_rule["degenerated_dihedral"] if isinstance(i, DecoratedAtomType)]
        deco_pairs = [i for i in atom_type_rule["degenerated_dihedral"] if isinstance(i, Iterable)]
        for term in getattr(molecule, "Dihedrals", []):
            deco2 = DecoratedAtomType(term.a2_atom_type_used)
            deco3 = DecoratedAtomType(term.a3_atom_type_used)
            # degenerate endocyclic torsions
            if deco2.tag_endo and deco3.tag_endo:
                term.a1_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                continue
            try:
                if any(deco2.is_subset_of(deco) for deco in deco_types) and any(
                    deco3.is_subset_of(deco) for deco in deco_types
                ):
                    raise StopIteration
                if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in deco_pairs):
                    raise StopIteration
            except StopIteration:
                term.a1_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                continue
    # improper wildcard
    if atom_type_rule["degenerated_improper"]:
        for term in getattr(molecule, "Impropers", []):
            deco1 = DecoratedAtomType(term.a1_atom_type_used)
            if any(deco1.is_subset_of(deco) for deco in atom_type_rule["degenerated_improper"]):
                term.a2_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a3_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
                term.a4_atom_type_used = atom_type_rule["AtomTypeDecoChar"]["DEG"]
    return molecule