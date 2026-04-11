"""
create by CFL 2023 oct. 9


力场参数有两种记录方式：
一是以字典的形式，程序运行中通过该方式,该字典为三级：
    第一级为相互作用方法名称，如atomtype, bondterm, angleterm, dihedralterm, improperterm, binc等
        还有一些特殊的key，如special_bond(12，13，14pair的作用方式), qmmodel(电荷的模型),
        combination_rule（混合规则）, equ_table(是否使用等价表)
    第二级为参数的名称：如c_3$c_3, c_3$c_3$c_3等
    第三级为具体的参数内容：主要包含以下key: name(名称), style(函数形式), para(参数值), tag(参数的标记符), sorce(参数的评价)
    例如：
    {"atomtype":{"c_3":{"name":"c_3","style":"LJ12_6","para":[3.5000,0.2000],"tag":"V","sorce":"nan"},
                 "c_4":{"name":"c_4","style":"LJ12_6","para":[3.6000,0.0500],"tag":"V","sorce":"nan"}
                },
     "bondterm":{"c_3$c_3":{"name":"c_3$c_3","style":"harmonic","para":[1.500,300.0],"tag":"V","sorce":"nan"},
                 "c_3$c_4":{"name":"c_3$c_4","style":"harmonic","para":[1.600,200.0],"tag":"V","sorce":"nan"},
                },
     "angleterm":{"c_3$c_3$c_3":{"name":"c_3$c_3$c_3","style":"harmonic","para":[120.0,30.0],"tag":"V","sorce":"nan"},
                 },
     "dihedralterm":{"c_3$c_3$c_3$c_3":{"c_3$c_3$c_3":{"name":"c_3$c_3$c_3$c_3","style":"amber","para":[0.0,3.0,0.0,0.0],"tag":"V","sorce":"nan"},
                },
    }

二是以文本文件的形式保存，力场参数的存储形式，
    文本文件中包含一些通用的设定special_bond，qmmodel，combination_rule, equ_table
    剩下的每行记录一条参数。每行分为以下几类字段，一是相互作用名称，二是组成该参数的原子类型（该部分与相互作用类型有关）
    三是函数形式，四是具体的参数（参数的数目与函数形式有关），剩下的包括tag, sorce,和count(有多少个数据拟合得到的该参数)
    范例可参考cadd/models/force_field/0.ff
"""
import json
from copy import deepcopy

from .. import CRATON_CONFIGURE
from ..utils import logger

def get_hybrid_froce_field(atf,fff1,fff2):
    """
    元素、杂化、净电荷、环大小(~)、芳香性(a)、环内外(@)、共轭(^)
    元素_杂化 (elem_hy):
        h_s, h_sp(氢键给体), 
         c_s, c_sp, c_sp2, c_sp3, c_su,
         o_sp2, o_sp3
         n_sp, n_sp2, n_sp3, n_spd2
         s_sp2, s_sp3, s_spd, s_spd2
         p_s, p_sp3, p_sp3d2
         b_sp3
         si_sp3
         f_s, cl_s, br_s, i_s
    芳香性  (arom): a (c_sp2, o_sp3, o_sp3-, n_sp2, n_sp3, n_sp2+, n_sp3+, n_sp3-, s_sp3有)
    净电荷 (charge): -, + (o_sp3-, n_sp2+, n_sp3+, n_sp3-,s_sp3-,b_sp3-,f_s-, cl_s-, br_s-, i_s-,na_s+, k_s+)
    环大小  (rz): ~3,4,5(除h_s, h_sp, halogen,metal ion, c_s, o_sp2, n_sp, s_sp2, p_s外)
    环内外 (ex): @ (除h_s, h_sp, halogen,metal ion, c_s, o_sp2, n_sp, s_sp2, p_s外)
    共轭 (ju): ^(c_sp,c_sp2, n_sp2)

    order:  elem_hy, arom, charge, rz, ex, ju

        atomtype: 元素-杂化 + 芳香性 + 净电荷, c_sp2
        bondtype: 元素-杂化 + 芳香性 + 环内外 + 共轭, c_sp2 c_sp2
       angletype: 元素-杂化 + 环大小 + 环内外 + 共轭, X c_sp2 X
    dihedraltype: 无素-杂化 + 芳香性 + 环内外 + 共轭， X c_sp2 c_sp3 X
    impropertype: 通用值, X X X X
            binc: 元素-杂化 + 芳香性 + 净电荷, c_sp2 c_sp3
    """
    
    #_term_attr = {
    #                "atomtype":["elem_hy","arom","charge"],
    #                "bondterm":["elem_hy","arom","ex","ju"],
    #                "angleterm":["elem_hy","rz","ex","ju"],
    #                "dihedralterm":["elem_hy","arom","ex","ju"],
    #                "binc":["elem_hy","arom","charge"],
    #              }
    _term_attr = {
                    "atomtype":["elem_hy","arom","charge"],
                    #"bondtype":["elem_hy","arom","ex","ju"],
                    "bondtype":["elem_hy","ju"],
                    "angletype":["elem_hy","rz","ex","ju"],
                    #"dihedraltype":["elem_hy","arom","ex","ju"],
                    "dihedraltype":["elem_hy","ju"],
                    "binc":["elem_hy","arom","charge"],
                  }
    
    from .atom_type import TypingDefine
    ATT = TypingDefine(atversion=atf)
    at_dicts = {}
    for term,vv in ATT.at_types.items():
        if term != "root":
            at_dicts[term] = {"hy":vv["Hybrid"],"arom":vv["Arom"],"charge":vv["Charge"]}
    ff_dicts = ForceField.read_files(fff1)

    def _get_atom_type_attr(at):
        rat = at.split("@")[0]
        rat = rat.split("^")[0]
        rat = rat.split("~")[0]
        elem = rat.split("_")[0]
        if rat not in at_dicts:
            return None
        attrs = at_dicts[rat]
        attrs["elem"] = elem
        attrs["elem_hy"] = f"{elem}_{attrs['hy']}"
        attrs["rz"] = "~" + at.split("~")[1][0] if "~" in at else ""
        attrs["ex"] = "@" if "@" in at else ""
        attrs["ju"] = "^" if "^" in at else ""
        return attrs

    _template = {"name": "X$X$X$X","fix_parameter": [],"ptag": "empi","pscore": "nan","pcount": "nan"}
    _term_ = {
                    "atomtype":"LJ12_6",
                    "bondterm":"harmonic",
                    "angleterm":"harmonic",
                    "dihedralterm":"amber",
                    "binc":"binc",
                  }

    empi_dicts = {
                  "atomtype":{},"bondterm":{},"angleterm":{},"dihedralterm":{},"binc":{}
                  }

    for term in ["atomtype","angleterm"]:
        for rr,vv in ff_dicts[term].items():
            at = rr if term == "atomtype" else rr.split("$")[1]
            attrs = _get_atom_type_attr(at)
            if attrs is None:
                continue
            eat = "".join([attrs[ss] for ss in _term_attr[term]])
            name = eat if term == "atomtype" else f"X${eat}$X"
            if name not in empi_dicts[term]:
                empi_dicts[term][name] = deepcopy(_template)
                if term == "atomtype":
                    empi_dicts[term][name]["mass"] = vv["mass"]
                empi_dicts[term][name]["pstyle"] = _term_[term]
                empi_dicts[term][name]["parameter"] = [[],[]]
                empi_dicts[term][name]["name"] = name
            empi_dicts[term][name]["parameter"][0].append(vv["parameter"][0])
            empi_dicts[term][name]["parameter"][1].append(vv["parameter"][1])

    for term in ["bondterm","binc","dihedralterm"]:
        for rr,vv in ff_dicts[term].items():
            ats = rr.split("$")
            if ats[1] == "X":
                continue
            attrs_0 = _get_atom_type_attr(ats[1]) if term == "dihedralterm" else _get_atom_type_attr(ats[0])
            attrs_1 = _get_atom_type_attr(ats[2]) if term == "dihedralterm" else _get_atom_type_attr(ats[1])
            if attrs_0 is None or attrs_1 is None:
                continue
            eat_0 = "".join([attrs_0[ss] for ss in _term_attr[term]])
            eat_1 = "".join([attrs_1[ss] for ss in _term_attr[term]])
            _tmp = [f"X${eat_0}${eat_1}$X",f"X${eat_1}${eat_0}$X"] if term == "dihedralterm" else [f"{eat_0}${eat_1}",f"{eat_1}${eat_0}"]

            names = list(set(_tmp).intersection(empi_dicts[term].keys()))
            x = 1
            if len(names) == 0:
                name = _tmp[0]
                empi_dicts[term][name] = deepcopy(_template)
                empi_dicts[term][name]["pstyle"] = _term_[term]
                empi_dicts[term][name]["parameter"] = [[]] if term == "binc" else [[],[]] if term == "bondterm" else [[],[],[],[],[],[],[],[]]
                empi_dicts[term][name]["name"] = name
            else:
                name = names[0]
                if name == _tmp[1] and term == "binc":
                    x = -1
            
            for ii in range(len(empi_dicts[term][name]["parameter"])):
                empi_dicts[term][name]["parameter"][ii].append(vv["parameter"][ii] * x)

    for term, items in empi_dicts.items():
        for name, item in items.items():
            if term == "binc":
                ss = name.split("$")
                if ss[0] == ss[1]:
                    para = [0.0000]
                else:
                    para = [round(sum(rr)/len(rr),4) for rr in item["parameter"]]
                empi_dicts[term][name]["parameter"] = para
            else:
                para = [sum(rr)/len(rr) for rr in item["parameter"]]
                empi_dicts[term][name]["parameter"] = para

    empi_dicts["general"] = {"use_nnff": False,"combination_rule": "LB","special_bond": ["None","None",0.8333,"None","None",0.5],"qmodel": "binc","equ_table": "yes"}
    empi_dicts["improperterm"] = {"X$X$X$X": {"name": "X$X$X$X","pstyle": "amber","fix_parameter": [],"parameter": [6.0],"ptag": "G","pscore": "nan","pcount": "nan"},}
    for term,vv in empi_dicts.items():
        print(term,len(vv))
    for rr in empi_dicts["dihedralterm"].keys():
        print(rr)
    with open(fff2,'w') as outf:
        outf.write(json.dumps(empi_dicts))

def convert_scalce_vdw(f1,f2):
    scalce_at = []
    with open(f1) as f:
        lines = f.readlines()
    with open(f2,'w') as outf:
        for line in lines:
            if line[:8] != "atomtype":
                outf.write(line)
            else:
                ss = line.split()
                factor = float(ss[6])
                if factor != 1.0:
                    scalce_at.append([ss[1],factor])
                del ss[6]
                outf.write("atomtype ")
                for s in ss[1:]:
                    outf.write("%15s " %s)
                outf.write("\n")
        outf.write("\n")
        for rr in scalce_at:
            outf.write("scalevdw@o_2w ")
            outf.write("%15s " % rr[0])
            outf.write("LJ12_6 ")
            outf.write("%15.4f " % rr[1])
            outf.write("fit ")
            outf.write("nan ")
            outf.write("nan \n")


_topol_to_ff_term = {
            "binc":"binc",
            "Atoms":"atomtype",
            "Bonds": "bondterm",
            "Angles": "angleterm",
            "Dihedrals": "dihedralterm",
            "Impropers": "improperterm",
            "Pairs": "pairwise",
            "Pair1n": "pairwise",
            "Pair14": "pairwise14",
            "Pair13": "pairwise13",
            "Pair12": "pairwise12",
        }

_Pair_terms = ["Pair1n", "Pair14", "Pair13", "Pair12"]
_Special_index = {"Pair12": 0, "Pair13": 1, "Pair14": 2}

def get_term_name(arr, ptyp):
    """
    考虑了名称的对称性，如o_1$c_3$n_3$c_4与c_4$n_3$c_3$o_1等价
    返回排序后的固定值
    """

    if ptyp == "improperterm":
        arr1 = sorted(arr[1:])
        return arr[0] + "$" + "$".join(arr1)
    else:
        arr = min(arr,list(reversed(arr)))
        return "$".join(arr)


def convert_dihe_4param_to_8param(ff):
    with open(ff) as inf:
        lines = inf.readlines()
    with open("new_empirical.ff",'w') as outf:
        for line in lines:
            rr = line.strip().split()
            if len(rr) > 1:
                if rr[0] != "dihedralterm":
                    outf.write(line)
                else:
                    outf.write("%s %s %s %s %s %s %s 0.0000 %s 180.0 %s 0.0000 %s 180.0 %s %s %s\n"%(rr[0],rr[1],rr[2],rr[3],rr[4],rr[5],rr[6],rr[7],rr[8],rr[9],rr[10],rr[11],rr[12]))
            else:
                outf.write(line)

class ForceField:
    def __init__(self,qmodel=None,use_nnff=False):
        self.qmodel = qmodel
        self.use_nnff = use_nnff

    ########read force field parameters file#############
    @staticmethod
    def read_files(ff,use_scalevdw,style=None):
        if style is None:
            style = "JSON" if ff.find(".json") != -1 else "TXT"
        if style == "JSON":
            jsondata = json.loads(open(ff).read())
        else:
            jsondata = ForceField.read_txt_file(ff)
        jsondata["general"]["use_scalevdw"] = use_scalevdw
        return jsondata

    @staticmethod
    def read_txt_file(ff):
        """
        读入txt格式的力场参数文件
        """
        with open(ff) as f:
            lines = f.read().splitlines()

        general_setting = ["special_bond", "combination_rule", "qmodel", "equ_table"]
        term_labels = {"atomtype": 1, "scalevdw":1,"binc": 2,
                 "bondterm": 2, "angleterm": 3, "dihedralterm": 4, "improperterm": 4,
                 "pairwise12":2,"pairwise13":2,"pairwise14":2,"pairwise":2,
                 "bondbondterm":3,"bondangleterm":3,
                 "bonddihedralterm":4,"angledihedralterm":4,"angledihedralangleterm":4,
                 "angleangleterm":4,
                 }
        
        ffjson = {"general":{"use_nnff": False}}

        for rr in lines:
            if rr.strip() != "" and rr[0].strip() != "#":
                try:
                    ss = [s.strip() for s in rr.strip().split()]
                    ss[0] = ss[0].lower()
                    if ss[0].find("scalevdw") != -1:
                        vdw_at = ss[0].split("@")[1]
                        ss = ["scalevdw"] + ss[1:] 
                        

                    if ss[0] not in general_setting and ss[0] not in term_labels:
                        raise Exception("Invalid line in file: %s" % rr)

                    if ss[0] in general_setting:
                        if ss[0] == "special_bond":
                            ffjson["general"]["special_bond"] = [float(ss[i]) if ss[i] != "None" else ss[i] for i in range(1,len(ss))]
                        else:
                            ffjson["general"][ss[0]] = ss[1]
                        continue
                    

                    if ss[0] not in ffjson.keys():
                        ffjson[ss[0]] = {}
                    
                    if ss[0] == "scalevdw":
                        if vdw_at not in ffjson[ss[0]]:
                            ffjson[ss[0]][vdw_at] = {}

                    term = ss[0]
                    start = term_labels[ss[0]]
                    atoms = ss[1 : start + 1]
                    name = get_term_name(atoms,term)

                    if name in ffjson[term]:
                        logger.warning("Duplicated ff line in %s, %s used" % (ff, rr))

                    _tmp = {"name":name,"pstyle":ss[start + 1],"fix_parameter":[],"parameter":[],"fit_parameter":[]}

                    if ss[0] == "atomtype":
                        _tmp["mass"] = float(ss[start + 2])
                        start = start + 3
                    else:
                        start = start + 2

                    for k in range(start, len(ss) - 3):
                        if ss[k][-1] == "*":
                            _tmp["fix_parameter"].append(k - start)
                            ss[k] = ss[k].strip("*")
                        elif ss[k][-1] == "~":
                            _tmp["fit_parameter"].append(k - start)
                            ss[k] = ss[k].strip("~")
                        _tmp["parameter"].append(float(ss[k]))

                    _tmp["ptag"] = ss[-3]
                    _tmp["pscore"] = float(ss[-2]) if ss[-2] != "nan" else "nan"
                    _tmp["pcount"] = int(ss[-1]) if ss[-1] != "nan" else "nan"
                    if term != "scalevdw":
                        ffjson[term][name] = _tmp
                    else:
                        ffjson[term][vdw_at][name] = _tmp
                except:
                    logger.error("Parse FF error: %s" % rr)
                    raise


        if ffjson["general"]["qmodel"] not in ("None", "binc", "atc"):
            raise Exception("Invalid qmodel. Should be None, binc, atc")

        return ffjson  

    #####分配力场参数########################
    @staticmethod
    def assign_para(molecule,this_ff,empi_ff = None,this_terms=None,return_ff=True):
        """
        力场参数分配到分子每个term项,总入口：
        input:
            molecule: Molecule object
            this_ff: force field parameter datas
        output:
            None
        """
        special_terms = ["Pair12", "Pair13", "Pair14"]

        for i in range(0, 3):
            if this_ff["general"]["special_bond"][0][i] == "None" and this_ff["general"]["special_bond"][1][i] == "None":
                if hasattr(molecule, special_terms[i]):
                    delattr(molecule, special_terms[i])

        
        if this_ff["general"]["qmodel"] not in [None, "None"]:
            ForceField.assign_charge_para(molecule,this_ff,empi_ff=empi_ff,this_terms=this_terms)
        
        ForceField.assign_ff_para(molecule,this_ff,empi_ff=empi_ff,this_terms=this_terms)
        if hasattr(molecule,"binc_loss_para_items"):
            molecule.loss_para_items["binc"] = molecule.binc_loss_para_items
            delattr(molecule,"binc_loss_para_items")
        if hasattr(molecule,"loss_charge_pair"):
            molecule.loss_para_items["pair_charge"] = molecule.loss_charge_pair
            delattr(molecule,"loss_charge_pair")


    ####分配非电荷参数########################
    @staticmethod
    def get_Y_term_name(ats,typ):
        if typ == "Angles":
            return [f"Y${ats[1]}$Y"]
        elif typ == "Dihedrals":
            return [f"Y${ats[1]}${ats[2]}$Y",f"Y${ats[2]}${ats[1]}$Y"]
        elif typ == "Impropers":
            return ["Y$Y$Y$Y"]
        else:
            return [ ]

    @staticmethod
    def get_empi_term_name(item,atoms,typ):
        _term_attr = {
                    "atomtype":["elem_hy","arom","charge"],
                    #"bondtype":["elem_hy","arom","ex","ju"],
                    "bondtype":["elem_hy","ju"],
                    "angletype":["elem_hy","rz","ex","ju"],
                    #"dihedraltype":["elem_hy","arom","ex","ju"],
                    "dihedraltype":["elem_hy","ju"],
                    "binc":["elem_hy","arom","charge"],
                  }
        def _get_empi_attrs(at,atom):
            rz = ""
            if len(atom.ring_size) > 0:
                if min(atom.ring_size) <= 5:
                    rz=f"~{min(atom.ring_size)}" 
            return {"elem_hy":f"{atom.elem.lower()}_{atom.atom_type_name_hybrid}","arom":atom.atom_type_name_arom,"charge":atom.atom_type_name_charge,
                    "rz":rz,"ex": "@" if "@" in at else "", "ju": "^" if "ju" in at else ""}


        if typ == "Atoms":
            return [f"{item.elem.lower()}_{item.atom_type_name_hybrid}{item.atom_type_name_arom}{item.atom_type_name_charge}"]
        if typ == "Bonds":
            _tmp0 = _get_empi_attrs(item.a1_atom_type_used,atoms[item.a1])
            _tmp1 = _get_empi_attrs(item.a2_atom_type_used,atoms[item.a2])
            emat0 = "".join([_tmp0[rr] for rr in _term_attr["bondtype"]])
            emat1 = "".join([_tmp1[rr] for rr in _term_attr["bondtype"]])
            return [f"{emat0}${emat1}",f"{emat1}${emat0}"]
        if typ == "Angles":
            _tmp = _get_empi_attrs(item.a2_atom_type_used,atoms[item.a2])
            emat = "".join([_tmp[rr] for rr in _term_attr["angletype"]])
            return [f"X${emat}$X"]
        if typ == "Dihedrals":
            _tmp0 = _get_empi_attrs(item.a2_atom_type_used,atoms[item.a2])
            _tmp1 = _get_empi_attrs(item.a3_atom_type_used,atoms[item.a3])
            emat0 = "".join([_tmp0[rr] for rr in _term_attr["dihedraltype"]])
            emat1 = "".join([_tmp1[rr] for rr in _term_attr["dihedraltype"]])
            return [f"X${emat0}${emat1}$X",f"X${emat1}${emat0}$X"]
        if typ == "Impropers":
            return ["X$X$X$X"]
        if typ == "binc":
            _tmp0 = _get_empi_attrs(item[0][1],atoms[item[0][0]])
            _tmp1 = _get_empi_attrs(item[1][1],atoms[item[1][0]])
            emat0 = "".join([_tmp0[rr] for rr in _term_attr["binc"]])
            emat1 = "".join([_tmp1[rr] for rr in _term_attr["binc"]])
            return [f"{emat0}${emat1}",f"{emat1}${emat0}"]
        return [ ]

    @staticmethod
    def assign_single_ff_para(term,item,parameter):
        if term == "Atoms":
            item.pstyle = parameter["pstyle"]
            for kk,vv in parameter.items():
                if kk not in ["pstyle","name"]:
                    setattr(item,kk,vv)
        else:
            for kk,vv in parameter.items():
                if kk not in ["name"]:
                    setattr(item,kk,vv)

    @staticmethod
    def assign_ff_para(molecule, this_ff, empi_ff=None, this_terms=None):
        """
        分配非charge参数
        注意：vdw参数优先考虑给定的pair参数，如果缺失则通过combine rule生成
        output:
            datas：使用到的参数和缺失的参数
        """
        if this_terms is None:
            this_terms = {term:[i for i in range(len(getattr(molecule,term)))] 
                          for term in _topol_to_ff_term.keys() if hasattr(molecule,term)}
        _null_parameter = {
                    "Atoms":{"pstyle":"LJ12_6","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Bonds":{"pstyle":"harmonic","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Angles":{"pstyle":"harmonic","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Dihedrals":{"pstyle":"amber","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Impropers":{"pstyle":"amber","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Pair12":{"pstyle":"LJ12_6","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Pair13":{"pstyle":"LJ12_6","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Pair14":{"pstyle":"LJ12_6","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Pair1n":{"pstyle":"LJ12_6","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"},
                    "Pairs":{"pstyle":"LJ12_6","fix_parameter":[],"parameter":None,"ptag":"null","pscore":"nan","pcount":"nan"}
                   }
        #from_empi = {"Atoms":[],"Bonds":[],"Angles":[],"Dihedrals":[],"Impropers":[]}
        #loss_para = {"Atoms":[],"Bonds":[],"Angles":[],"Dihedrals":[],"Impropers":[]}
        #loss_para_items = {"Atoms":[],"Bonds":[],"Angles":[],"Dihedrals":[],"Impropers":[]}
        from_empi = {}
        loss_para = {}
        loss_para_items = {}
        for term in ["Atoms","Bonds","Angles","Dihedrals","Impropers","Pair12","Pair13","Pair14","Pair1n","Pairs"]:
            if hasattr(molecule,term) and term in this_terms:
                terms = getattr(molecule,term)
                item = _topol_to_ff_term[term]
                pp = this_ff[item] if item in this_ff else {}
                pp_items = set(pp.keys())
                empi_pp = None
                if empi_ff is not None:
                    if item in empi_ff:
                        empi_pp = empi_ff[item]
                        empi_pp_items = set(empi_pp.keys())
                    else:
                        empi_pp = None
                for ii in this_terms[term]:
                    tt = terms[ii]
                    _name = list(set(tt.atom_type_used_names).intersection(pp_items))
                    parameter = deepcopy(pp[_name[0]]) if len(_name) > 0 else None
                    if len(_name) > 1:
                        logger.warning(f"there multi parameter match this term, the {_name[0]} has been used")
                    if parameter is None:
                        _name = list(set(ForceField.get_Y_term_name(tt.atom_type_used,term)).intersection(pp_items))
                        parameter = deepcopy(pp[_name[0]]) if len(_name) > 0 else None
                    name = _name[0] if len(_name) > 0 else min(tt.atom_type_used_names)
                    if parameter is None:
                        if empi_pp is not None and term not in ["Pair12","Pair13","Pair14","Pair1n","Pairs"]:
                            _empi_name = list(set(ForceField.get_empi_term_name(tt,molecule.Atoms,term)).intersection(empi_pp_items))
                            parameter = deepcopy(empi_pp[_empi_name[0]]) if len(_empi_name) > 0 else None
                            if parameter is not None:
                                if term not in from_empi:
                                    from_empi[term] = []
                                from_empi[term].append(name)
                            #else:
                            #    if name not in loss_para[term]:
                            #        loss_para[term].append(name)
                        
                    if parameter is None:
                        if term not in loss_para:
                            loss_para[term] = []
                        loss_para[term].append(name)
                        if term not in loss_para_items:
                            loss_para_items[term] = []
                        loss_para_items[term].append(ii)
                        parameter = deepcopy(_null_parameter[term])
                    ForceField.assign_single_ff_para(term,tt,parameter)

                    tt.atom_type_used_name = name
                    parameter["name"] = name
                    tt._ff_parameter = parameter
        for tt,vv in from_empi.items():
            if tt not in _Pair_terms:
                logger.warning(f"{tt} parameter from empi froce field: {set(vv)}")
        molecule.loss_para_items = {}
        for tt,vv in loss_para.items():
            if tt not in _Pair_terms:
                logger.warning(f"{tt} parameter loss: {set(vv)}")
                molecule.loss_para_items[tt] = loss_para_items[tt]
        #####通过combine_rule生成缺失的pair vdw参数
        if this_ff["general"]["combination_rule"] is not None:
            ForceField.assign_vdw_to_pair(molecule,{term:pair for term,pair in this_terms.items() if term in _Pair_terms},
                                          combination_rule=this_ff["general"]["combination_rule"],special_bond=this_ff["general"]["special_bond"][3:6])

        if hasattr(molecule,"loss_pair_para"):
            for kk,vv in molecule.loss_pair_para.items():
                molecule.loss_para_items[kk] = vv
            delattr(molecule,"loss_pair_para")

    @staticmethod
    def assign_vdw_to_pair(molecule, pair_terms, combination_rule="LB", special_bond=[None, None, 0.5]):
        """
        根据conbine rule计算出pair的参数
        """
        #pair_terms = [term for term in _Pair_terms if hasattr(molecule, term)]
        loss_pair_para = {}
        for term,indxs in pair_terms.items():
            if term != "Pair1n":
                scale_factor = special_bond[_Special_index[term]]
            else:
                scale_factor = 1
            if scale_factor not in ["None",None]:
                items = getattr(molecule, term)
                for ii in indxs:
                    pp = items[ii]
                    _tmp_pp = getattr(pp,"parameter",None)
                    if _tmp_pp is None:
                        para1 = molecule.Atoms[pp.a1].parameter if hasattr(molecule.Atoms[pp.a1],"parameter")  else None
                        para2 = molecule.Atoms[pp.a2].parameter if hasattr(molecule.Atoms[pp.a2],"parameter")  else None
                        if para1 is not None and para2 is not None:
                            if combination_rule in ["LB", "lb"]:
                                sigma = para1[0] / 2.0 + para2[0] / 2.0
                                espi = scale_factor * (para1[1] * para2[1]) ** 0.5
                                pp.parameter = [sigma, espi]
                                pp.pstyle = molecule.Atoms[pp.a1].pstyle
                                pp.ptag = molecule.Atoms[pp.a1].ptag
                                pp.fix_parameter = []
                                pp.pscore = molecule.Atoms[pp.a1].vdw_score
                                pp.pcount = molecule.Atoms[pp.a1].pcount
                                pp.combination_rule = combination_rule
                                pp.scale_factor = scale_factor
                        else:
                            if term not in loss_pair_para:
                                loss_pair_para[term] = []
                            loss_pair_para[term].append(ii)
        molecule.loss_pair_para = loss_pair_para
                
    #####分配电荷参数########################
    @staticmethod
    def assign_charge_para(molecule, this_ff,empi_ff=None,this_terms=None):

        """
        按一定的模式，把电荷分配到每个原子中去
        """
        qmodel = this_ff["general"]["qmodel"]
        if this_terms is None:
            #this_terms = {atom.ID:[i for i in range(len(atom.connectivity))] for atom in molecule.Atoms}
            this_terms={qmodel:None}
            this_terms["pair_charge"] = {term:[ii for ii in range(len(getattr(molecule,term,[]) ))]for term in _Pair_terms}
        else:
            if qmodel not in this_terms:
                return
        
        if qmodel == "binc":
            if empi_ff is not None and "binc" in empi_ff:
                empi_binc_ff = empi_ff["binc"]
            else:
                empi_binc_ff = None
            ForceField.get_binc_charge(molecule,this_ff["binc"],empi_binc_ff=empi_binc_ff,this_terms=this_terms["binc"])
        elif qmodel == "atc":
           ForceField.get_atc_charge(molecule,this_ff["atc"],this_terms=this_terms["atc"])
        elif qmodel in ("esp", "am1bcc", "nn", "nnv2", "nnbcc","manual"):
            ForceField.get_X_charge(molecule,qmodel,this_ff[qmodel])
        else:
            raise Exception("Invalid qmodel. Should be binc, atc, esp, am1bcc, nn, nnv2, nnbcc")
        ForceField.assign_charge_to_pair(molecule,special_bond=this_ff["general"]["special_bond"][0:3],this_terms=this_terms["pair_charge"])

    @staticmethod
    def get_binc_charge(molecule, binc_ff,empi_binc_ff=None,this_terms=None):
        """
        计算binc charge
        """
        binc_from_empi = []
        binc_loss_para = []
        binc_loss_para_items = {}
        if this_terms is None:
            this_terms = {atom.ID:[i for i in range(len(atom.connectivity))] for atom in molecule.Atoms}
        
        for idx,rr in this_terms.items():
            atom = molecule.Atoms[idx]
            if not hasattr(atom,"binc_parameter"):
                nn = len(atom.connectivity)
                atom.binc_parameter = [None for i in range(nn)]
                atom.binc_score = [None for i in range(nn)]
                atom.binc_tag = [None for i in range(nn)]
                atom.binc_style = [None for i in range(nn)]
                atom.binc_count = [None for i in range(nn)]
                atom._ff_binc_parameter = [None for i in range(nn)]

            for ii in rr:
                atom2 = molecule.Atoms[atom.connect[ii]]
                name1 = atom.binc_atom_type + "$" + atom2.binc_atom_type
                name2 = atom2.binc_atom_type + "$" + atom.binc_atom_type
                if name1 in binc_ff:
                    name = name1
                    parameter = deepcopy(binc_ff[name1])
                    factor = 1.0
                elif name2 in binc_ff:
                    name = name2
                    parameter = deepcopy(binc_ff[name2])
                    factor = -1.0
                else:
                    name = min([name1,name2])
                    parameter = None

                if parameter is None and empi_binc_ff is not None:
                    _empi_name1,_empi_name2 = ForceField.get_empi_term_name([[atom.ID,atom.binc_atom_type],[atom2.ID,atom2.binc_atom_type]],molecule.Atoms,"binc")
                    if _empi_name1 in empi_binc_ff:
                        parameter = deepcopy(empi_binc_ff[_empi_name1])
                        factor = 1.0
                    elif _empi_name2 in empi_binc_ff:
                        parameter = deepcopy(empi_binc_ff[_empi_name2])
                        factor = -1.0
                    else:
                        parameter = None
                    if parameter is not None:
                        if name not in binc_from_empi:
                            binc_from_empi.append(name)
                    #else:
                    #    if name not in binc_loss_para:
                    #        binc_loss_para.append(name)
                    
                if parameter is None: 
                    if idx not in binc_loss_para_items:
                        binc_loss_para_items[idx] = []
                    binc_loss_para_items[idx].append(ii)

                    binc_loss_para.append(name)
                    parameter = {"parameter":None,"ptag":"null","pscore":"nan","pstyle":"binc","pcount":"nan"}

                parameter["name"] = name
                atom.binc_parameter[ii] = parameter["parameter"][0] * factor if parameter["parameter"] is not None else None
                atom.binc_tag[ii] = parameter["ptag"]
                atom.binc_score[ii] = parameter["pscore"]
                atom.binc_style[ii] = parameter["pstyle"]
                atom.binc_count[ii] = parameter["pcount"]
                atom._ff_binc_parameter[ii] = parameter
        #if len(binc_from_empi) > 0:
        #    logger.warning(f"binc parameter from empi: {binc_from_empi}")
        if len(binc_loss_para) > 0:
            logger.warning(f"binc loss parameter: {set(binc_loss_para)}")

        for idx,rr in this_terms.items():
            atom = molecule.Atoms[idx]
            #for atom in molecule.Atoms:
            if "None" not in atom.binc_parameter and None not in atom.binc_parameter:
                
                atom.point_charge = sum(atom.binc_parameter) + atom.point_charge_base
            else:
                atom.point_charge = None 
        molecule.binc_loss_para_items = binc_loss_para_items

    @staticmethod
    def get_atc_charge(molecule,atc_ff,this_terms=None):
        if this_terms is None:
            this_terms = [atom.ID for atom in molecule.Atoms]
        for ii in this_terms:
            atom = molecule.Atoms[ii]
            atom.point_charge = atc_ff[atom.atc_atom_type]["parameters"][0] if atom.atc_atom_type in atc_ff else None

    @staticmethod
    def get_X_charge(molecule,ctype,charge_ff):
        if molecule.mole_name in charge_ff:
            for atom, q in zip(molecule.Atoms, charge_ff[molecule.mole_name]):
                atom.point_charge = q
        else:
            for atom in molecule.Atoms:
                atom.point_charge = None        
    
    @staticmethod
    def assign_charge_to_pair(molecule, special_bond=[None, None, 0.8333],this_terms=None):
        """
        分配charge到pair中
        """
        if this_terms is None:
            this_terms = {term:[ii for ii in range(len(getattr(molecule,term,[]) ))]for term in _Pair_terms}
        
        loss_charge_pair = {}
        #pair_terms = [term for term in _Pair_terms if hasattr(molecule, term)]
        for term, arrs in this_terms.items():
            if term != "Pair1n":
                scale_factor = special_bond[_Special_index[term]]
            else:
                scale_factor = 1
            if scale_factor not in [None, "None"]:
                items = getattr(molecule, term,[])
                for ii in arrs:
                    aa = items[ii]
                    aa.charge_parameter = []
                    aa.charge_scale_factor = scale_factor
                    for atom in [aa.a1, aa.a2]:
                        if molecule.Atoms[atom].point_charge is not None and molecule.Atoms[atom].point_charge != "None":
                            aa.charge_parameter.append(scale_factor**0.5 * molecule.Atoms[atom].point_charge)
                        else:
                            if term not in loss_charge_pair:
                                loss_charge_pair[term] = []
                            loss_charge_pair[term].append(ii)
                            aa.charge_parameter.append(molecule.Atoms[atom].point_charge)
        molecule.loss_charge_pair = loss_charge_pair
    ##########################################
