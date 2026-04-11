import itertools

from copy import deepcopy
from typing import Iterable, List, Any

import networkx as nx
import numpy as np
import string


from .atom import VS, Atom
from .topology import Angle, Bond, Constrain, Dihedral, Improper, Pair
from .elements import ELEMENT_ORDER
from . import vs_gen

#from .mole_analy import calc_dipole, calc_inertia, find_center
#from changban.chemistry.molecule_mechanics import force_field
#from changban.chemistry.molecule_mechanics.atom_type import AtomTypeDecoChar, DecoratedAtomType
#from changban.chemistry.molecule_mechanics.ff_calc import calculator
#from changban.chemistry.structure import vs_gen


#from changban.chemistry.structure.information import get_formula, get_inchikey
from .constants import ATOM_DEFAULT_ATTRIBUTES
from .format.smiles_parse import SmilesData
        
class _Molecule:

    attrs_topol = [
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        "Pair12",
        "Pair13",
        "Pair14",
        "Pair1n",
        "AlteredPair",
    ]

    def __init__(self, style=""):
        self.style = style
        self.Atoms: List[Atom] = []
        self._rings = {}
        self._rings_defined = False
        self._ob_mol = None
        self._renew_inchikey_flag = False
        self._inchikey_3d_structure_flag = False

    #####Graph network
    def make_mole_as_graph(self, ignore_elems=None):
        self.G = nx.Graph()
        if ignore_elems is None:
            ignore_elems = []
        atoms = [i for i in range(len(self.Atoms)) if self.Atoms[i].elem not in ignore_elems]
        self.G.add_nodes_from(atoms)
        self.G.add_edges_from([(bond.a1, bond.a2) for bond in self.Bonds if bond.a1 in atoms and bond.a2 in atoms])

    def calc_bond_distance(self, p0, p1):
        if not hasattr(self,"G"):
            self.make_mole_as_graph()
        return len(nx.shortest_path(self.G, source=p0, target=p1))

    ####编辑分子      
    def add_atom(self, style="aa"):
        atom_n = len(self.Atoms)
        self.Atoms.append(Atom(style, atom_n)) 

    def create_atoms(self, n: int) -> List[Atom]:
        """
        create Atoms:
            n: int, atom count of molecule
        """
        self.Atoms = [Atom("aa",i) for i in range(n)]


    #def create_atoms(molecule, n):
    #    """
    #    生成self.Atoms,包含n个Atom对象
    #    输入：
    #        n: int, 原子的个数
    #    输出：
    #        产生self.Atoms
    #    """
    #    for i in range(n):
    #        molecule.add_atom("aa")
        #molecule.Atoms = [Atom("aa") for i in range(n)]

    def give_atoms_info(self, arr):
        """
        读取原子的属性，输出一个字典。是get_atoms_info的逆操作
        ###TODO: 该方法需要UT###
        输入：
            arr: List[str]，需要提取的原子的属性，（必须在self.__defaultatomatts中)
        输出：
            包括原子属性的字典
        """
        datas = {}
        for aa in arr:
            datas[aa] = []
            for i in range(len(self.Atoms)):
                datas[aa].append(getattr(self.Atoms[i], aa))
        return datas

    def get_atoms_info(self, informations):
        """
        从一个字典中读取信息，分配到每个原子。
        通常情况下，从SMILES, mol, pdb, gjf等方式生成一个Molecule对象，需要调用该方法得到分子最基本的属性
        如每个原子的元素，连接性，键类型和formal charge
        ###TODO: 该方法需要UT###
        输入：
            datas: dict, 记录每个原子的属性，key为原子属性（key 在self.__defaultatomatts中的才被使用到）
                     value 为list，size为原子的数目
        输出：
            为每个原子分配相应的属性值。通常情况下是基本的元素，连接性，键类型和formal charge
        """
        _special_term = [  # noqa
            "bonds",
            "charge",
        ]
        if "atom_count" not in informations:
            informations["atom_count"] = len(informations["elements"])
        self.create_atoms(informations["atom_count"])
        for kk, vv in informations.items():
            if kk in ATOM_DEFAULT_ATTRIBUTES:
                if kk == "charge":
                    for term,vvv in vv.items():
                        for i in range(len(self.Atoms)):

                        #dd = {}
                            setattr(self.Atoms[i],f"{term}_charge",vvv[i])
                            #dd[term] = vv[term][i]
                        #setattr(self.Atoms[i], "charge", dd)
                elif kk == "formal_charge":
                    if len(vv) > 0 and isinstance(vv[0],list):
                        for tt in vv:
                            setattr(self.Atoms[tt[0]], "formal_charge", tt[1])
                    else:
                        for ii,v in enumerate(vv):
                            setattr(self.Atoms[ii],"formal_charge", v)
                elif kk in ["ID","mass","ring_size","cis_trans","chirality"]:
                    pass
                else:
                    if len(vv) != len(self.Atoms):
                        vv = [None for ii in range(len(self.Atoms))]
                    for i in range(len(self.Atoms)):
                        setattr(self.Atoms[i], kk, vv[i])
        if "connect" not in informations.keys():
            if hasattr(self, "Bonds"):
                self.connectivity_from_bonds()

    def get_mole_info(self, informations):
        for key, value in informations.items():
            if key not in ATOM_DEFAULT_ATTRIBUTES and key.islower() \
                and key not in ["script","formula","heavy_atoms","mass","molecule_mass",
                            "net_charge","constrain","constrain_term","constrain_value"]: #删除了torsion_number
                setattr(self, key, value)
            elif key == "confID":
                setattr(self, key, value)
            if key == "constrain" and value:
                self.create_constrain([value[0][0]+[value[0][1]]])
            if key == "constrain_term" and value not in [None,"None"]:
                self.create_constrain([value + [informations["constrain_value"]]])

    def remove_atom(self):
        pass

    def add_bond(self):
        pass

    def remove_bond(self):
        pass

    #####topology生成#################
    #生成topol时，顺序按item.atoms排序#
    # 以保持同一分子的顺序固定         #
    ##################################
    def create_topols(self,smiles_flag=True):
        """
        create the topols of molecule, including Bonds, Angles, Dihedras, and check the smiles, inchi_key, inchi
        每次运用该，Bonds, Angles, Dihedrals会重新生成
        NOTE: the Impropers is created by the function of create_improper
        """
        self.Bonds = []
        self.Angles =[]
        self.Dihedrals =[] 
        self.Impropers = []
        
        bond_name, angle_name, dihedral_name = [], [], []
        for ii in range(0, len(self.Atoms)):
            for jj in self.Atoms[ii].connect:
                if (f"{ii}-{jj}" not in bond_name) and (f"{jj}-{ii}" not in bond_name):
                    self.Bonds.append(Bond("bond", ii,jj))
                    bond_name.append(self.Bonds[-1].name)
                for kk in self.Atoms[jj].connect:
                    if kk not in [ii, jj]:
                        if (f"{ii}-{jj}-{kk}" not in angle_name) and (f"{kk}-{jj}-{ii}" not in angle_name):
                            self.Angles.append(Angle("angle", ii,jj,kk))
                            angle_name.append(self.Angles[-1].name)
                        for ll in self.Atoms[kk].connect:  # noqa
                            if ll not in [ii, jj, kk]:
                                if (f"{ii}-{jj}-{kk}-{ll}" not in dihedral_name) and (f"{ll}-{kk}-{jj}-{ii}" not in dihedral_name):
                                    self.Dihedrals.append(Dihedral("dihedral", ii,jj,kk,ll))
                                    dihedral_name.append(self.Dihedrals[-1].name) 
                                    
        for term in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
            if len(getattr(self, term)) == 0:
                delattr(self, term)
            else:
                setattr(self,term,sorted(getattr(self,term), key=lambda item:item.atoms))

        flag = False
        if not hasattr(self,"_molecule_name"):
            flag = True
            self._molecule_name = self.get_formula()
        if smiles_flag:
            if not hasattr(self,"_smiles") and len(self.Atoms) < 500:
                self.get_inchikey()
        if flag:
            self._molecule_name = self._inchi_key

    def create_improper(self, create_method="classic"):
        """
        注意：通常情况下，create_topols不会生成Impropers
        输入：
            每个原子有正确的connect属性
            create_method: str, 三种不同的生成Impropers的方法
                           "classic" 传统的方法，即含有双键的三连接的原子有improper项
                           "atom_type"，根据atom_type的plate属性确定是否生成improper项
                           "mixture",混合上面两种方法
        """
        centers = []
        if create_method == "classic":
            for atom in self.Atoms:
                if "2" in atom.bond_type and len(atom.connectivity) == 3:
                    centers.append(atom)
        elif create_method == "atom_type":
            for atom in self.Atoms:
                if atom.plate == "yes" and len(atom.connectivity) == 3:
                    centers.append(atom)
        elif create_method == "mixture" or create_method == "mix":
            for atom in self.Atoms:
                if len(atom.connectivity) == 3 and ("2" in atom.bond_type or atom.plate == "yes"):
                    centers.append(atom)

        if len(centers) != 0:
            self.Impropers = sorted([Improper("improper",atom.ID,*atom.connectivity,
                                              atom_numbers=[self.Atoms[an].atom_number for an in atom.connectivity]
                                              ) for atom in centers], 
                                    key = lambda item:item.atoms)


    def create_intra_nonbond(self):
        """
        生成非键拓扑结构
        包括
           self.Pair12, 1 bond pairs
           self.Pair13, 2 bond pairs
           self.Pair14, 3 bond pairs
           self.Pair1n, > 3 bond pairs
        """
        pairs12, pairs13, pairs14, pairs1n = self.get_pairs_set()
        self.Pair12 = sorted([Pair("pair", *pair) for pair in pairs12], key = lambda item:item.atoms)
        self.Pair13 = sorted([Pair("pair", *pair) for pair in pairs13], key = lambda item:item.atoms)
        self.Pair14 = sorted([Pair("pair", *pair) for pair in pairs14], key = lambda item:item.atoms)
        self.Pair1n = sorted([Pair("pair", *pair) for pair in pairs1n], key = lambda item:item.atoms)
        for term in ["Pair12", "Pair13", "Pair14", "Pair1n"]:
            if len(getattr(self, term)) == 0:
                delattr(self, term)
            
    def create_intra_nonbond_macromole(self):
        """
        只生成Pair14
        """
        pairs14 = set()
        for dihe in getattr(self,"Dihedrals",[]):
            pairs14.add(tuple(sorted([dihe.a1,dihe.a4])))
        #pairs12, pairs13, pairs14, pairs1n = self.get_pairs_set()
        self.Pair14 = sorted([Pair("pair", *pair) for pair in pairs14], key = lambda item:item.atoms)

    def get_pairs_set(self):
        pairs12 = set()
        pairs13 = set()
        pairs14 = set()
        pairs1n = set(itertools.combinations(range(len(self.Atoms)), 2))
        for bond in getattr(self, "Bonds", []):
            a1, a2 = bond.a1, bond.a2
            pairs12.add(tuple(sorted([a1, a2])))
            for neigh1 in self.connectivity[a1]:
                if neigh1 != a2:
                    pairs13.add(tuple(sorted([neigh1, a2])))
            for neigh2 in self.connectivity[a2]:
                if neigh2 != a1:
                    pairs13.add(tuple(sorted([a1, neigh2])))
                for neigh1 in self.connectivity[a1]:
                    if neigh1 != a2 and neigh2 != a1 and neigh1 != neigh2:
                        pairs14.add(tuple(sorted([neigh1, neigh2])))
        pairs13 = pairs13.difference(pairs12)
        pairs14 = pairs14.difference(pairs13).difference(pairs12)
        pairs1n = pairs1n.difference(pairs14).difference(pairs13).difference(pairs12)
        return list(sorted(pairs12)), list(sorted(pairs13)), list(sorted(pairs14)), list(sorted(pairs1n))

    def get_dummy_1234n_pairs(self):
        dummy = []
        for aa in self.Atoms:
            if aa.atom_type_name == "_D":
                dummy.append(aa.No)
        for an in dummy:
            for an0 in self.Atoms[an].connect:
                if an not in self.Atoms[an0].connect:
                    self.Atoms[an0].connect.append(an)
                    self.Atoms[an0].bond_type.append("1")
                    self.Atoms[an0].bond_type_aromatic.append("1")
        pairs12 = set()
        pairs13 = set()
        pairs14 = set()
        pairs1n = set([tuple(sorted([an, nn])) for an in dummy for nn in range(len(self.Atoms)) if nn != an])
        for an in dummy:
            for neigh0 in self.Atoms[an].connect:
                pairs12.add(tuple(sorted([an, neigh0])))
                for neigh1 in self.Atoms[neigh0].connect:
                    if neigh1 != an:
                        pairs13.add(tuple(sorted([neigh1, an])))
                        for neigh2 in self.Atoms[neigh1].connect:
                            if neigh2 != an and neigh2 != neigh0:
                                pairs14.add(tuple(sorted([neigh2, an])))
        pairs13 = pairs13.difference(pairs12)
        pairs14 = pairs14.difference(pairs13).difference(pairs12)
        pairs1n = pairs1n.difference(pairs14).difference(pairs13).difference(pairs12)
        return list(sorted(pairs12)), list(sorted(pairs13)), list(sorted(pairs14)), list(sorted(pairs1n))

    def create_constrain(self, arrs):
        """
        生成constrain
        输入：
            arrs: List[[int and float]], 受限制的原子及数值，如[[1,2,1.4],[1,2,3,4, 120.0]]，
                  表示有两个受限制的项，第一个限制1和2号原子的距离为1.4 A，第二个限制1，2，3和4号原子组成的二面角为120.0度
        """
        if not hasattr(self, "constrain"):
            self.constrain = []
        for arr in arrs:
            self.constrain.append(Constrain(arr[:-1], arr[-1]))
            if hasattr(self,"Dihedral"):
                if len(self.Dihedrals) > 0 and hasattr(self.Dihedrals[0], "para"):
                    for dihe in self.Dihedrals:
                        if self.constrain[-1].name in dihe.names:
                            self.constrain[-1].para = dihe.para
                            self.constrain[-1].para_style = dihe.style
                            break
        if len(self.constrain) == 0:
            delattr(self, "constrain")

    def connectivity_from_bonds(self):
        """
        从self.Bonds中产生原子的connect属性
        ###TODO: 该方法需要UT###
        输入：
            预先知道self.Bonds
        输出：
            为每个原子分配或更新connect属性,
        """
        for aa in self.Atoms:
            aa.connect = []
        for bb in self.Bonds:
            if bb.a2 not in self.Atoms[bb.a1].connect:
                self.Atoms[bb.a1].connect.append(bb.a2)
            if bb.a1 not in self.Atoms[bb.a2].connect:
                self.Atoms[bb.a2].connect.append(bb.a1)
            if hasattr(bb, "bond_type"):
                for a in [[bb.a1, bb.a2], [bb.a2, bb.a1]]:
                    n = self.Atoms[a[1]].connect.index(a[2])
                    if hasattr(self.Atoms[a[1]], "bond_type"):
                        nn = n - len(self.Atoms[a[1]].bond_type)
                        if nn < 0:
                            self.Atoms[a[1]].bond_type[n] = bb.bond_type
                        else:
                            for i in range(0, nn):
                                self.Atoms[a[1]].bond_type.append(None)
                            self.Atoms[a[1]].bond_type.append(bb.bond_type)
                    else:
                        self.Atoms[a[1]].bond_type = []
                        for i in range(0, n):
                            self.Atoms[a[1]].bond_type.append(None)
                        self.Atoms[a[1]].bond_type.append(bb.bond_type)

    ########################################################

    def update_topol_value(self):
        """
        更新每个拓扑项的值，即调用Bond, Angle, Dihedral等对象中的calc_vaules方法计算
        通常情况下，当一个分子的3D坐标发生变化时，需要用该方法进行更新
        ###TODO: 该方法需要UT###
        输入：
        输出：
        """
        this_terms = [term for term in self.__dict__.keys() if term[0].isupper()]
        del this_terms[this_terms.index("Atoms")]
        if hasattr(self, "constrain"):
            this_terms.append("constrain")
        
        for term in this_terms:
            for aa in getattr(self, term):
                this_coor = []
                if hasattr(aa, "a1"):
                    this_coor.append(self.Atoms[aa.a1].coor)
                if hasattr(aa, "a2"):
                    this_coor.append(self.Atoms[aa.a2].coor)
                if hasattr(aa, "a3"):
                    this_coor.append(self.Atoms[aa.a3].coor)
                if hasattr(aa, "a4"):
                    this_coor.append(self.Atoms[aa.a4].coor)
                aa.calc_value(this_coor)

    def get_formula(self):
        element_number = {element:len([1 for atom in self.Atoms if atom.elem == element]) for element in set(self.elements)}
        formula = ""
        for element in ELEMENT_ORDER:
            formula += (f"{element}{str(element_number[element])}" if element_number[element] > 1 else f"{element}") if element in element_number else ""
        return formula

    def get_inchikey(self):
        rdk = SmilesData()
        rdk._convert(self,extra_var = {"normalization":True,"structure_3d": self._inchikey_3d_structure_flag})
        self._inchi_key = rdk.inchi_key
        self._smiles = rdk.smiles
        self._inchi = rdk.inchi

    def get_null_atom_type(self):
        alphabet = list(string.ascii_lowercase)
        elems = {}
        for aa in self.Atoms:
            if aa.elem not in elems.keys():
                elems[aa.elem] = 1
                aa.atom_name = "%s1" % aa.elem
            else:
                elems[aa.elem] += 1
                aa.atom_name = "%s%d" % (aa.elem, elems[aa.elem])
            n = int(aa.No / 26)
            m = aa.No % 26
            aa.atom_type_null = "%s%s" % (alphabet[n], alphabet[m])

    def reorder_atom_number(self,smiles_flag=True):
        match_dict = {}
        atoms = []
        for ii, atom in enumerate(self.Atoms):
            match_dict[atom.No] = ii
            atom.No = ii
            atoms.append(atom)
        for atom in atoms:
            bond_type = []
            connectivity = []
            bond_type_aromatic = []
            connectivity_type = []
            for ii, an in enumerate(atom.connectivity):
                if an in match_dict:
                    connectivity.append(match_dict[an])
                    bond_type.append(atom.bond_type[ii])
                    if hasattr(atom, "bond_type_aromatic"):
                        bond_type_aromatic.append(atom.bond_type_aromatic[ii])
                    if hasattr(atom, "connectivity_type"):
                        connectivity_type.append(atom.connectivity_type[ii])
            atom.connectivity = connectivity
            atom.bond_type = bond_type
            if len(bond_type_aromatic) > 0:
                atom.bond_type_aromatic = bond_type_aromatic
            if len(connectivity_type) > 0:
                atom.connectivity_type = connectivity_type
        self.Atoms = atoms
        self.atom_count = len(self.Atoms)
        self.create_topols(smiles_flag=smiles_flag)
        self.bond_count = len(self.Bonds)
        self.create_intra_nonbond()
        return match_dict

    def reorder_atoms(self, indexes_map: {int: int}):
        """
        For FEP. Merge two ligands into one.
        indexes_map = {old_id: new_id}
        """
        indexes_new = [v for k, v in sorted(indexes_map.items())]
        if len(indexes_new) != len(self.Atoms) or set(indexes_new) != set(range(len(self.Atoms))):
            raise Exception("Number of atoms and order not match")
        for i, atom in enumerate(self.Atoms):
            atom.No = indexes_new[i]
            atom.connect = [indexes_new[i] for i in atom.connect]
        self.Atoms = list(sorted(self.Atoms, key=lambda x: x.No))
        for attr in self.attrs_topol:
            for term in getattr(self, attr, []):
                try:
                    term.a1 = indexes_new[term.a1]
                    term.a2 = indexes_new[term.a2]
                    term.a3 = indexes_new[term.a3]
                    term.a4 = indexes_new[term.a4]
                except AttributeError:
                    pass

    def check_vs(self):
        """
        检查分子里面是否有虚原子
        ###TODO: 该方法需要UT###
        输入：
        输出：
        """
        vs_arr = []
        self.Vss = []
        for aa in self.Atoms:
            if aa.elem == "EP":
                patoms = [aa.connect[0]]
                for ii in self.Atoms[patoms[0]].connect:
                    if ii != aa.No:
                        patoms.append(ii)
                if self.Atoms[patoms[0]].elem == "N":
                    if len(patoms) - 1 == 2:
                        v = VS("style1", aa.No, patoms, [0.5, -0.03])
                        v.coor = aa.coor
                if self.Atoms[patoms[0]].elem == "O":
                    if len(patoms) - 1 == 2:
                        v = VS("style2", aa.No, patoms, [0.5, 0.02])
                        v.coor = aa.coor
                self.Vss.append(v)
                vs_arr.append(aa.No)
        if len(vs_arr) == 0:
            delattr(self, "Vss")
        else:
            tmp_dict = {"Bonds": [], "Angles": [], "Dihedrals": [], "Impropers": []}
            for term in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
                if hasattr(self, term):
                    rrs = getattr(self, term)
                    for rr in rrs:
                        atom_arr = []
                        for i in range(1, 5):
                            if hasattr(rr, "a%d" % i):
                                atom_arr.append(getattr(rr, "a%d" % i))
                        if len(set(atom_arr).intersection(vs_arr)) == 0:
                            tmp_dict[term].append(rr)
                    setattr(self, term, tmp_dict[term])
        for vs in vs_arr:
            self.Atoms[vs].connect = []
            self.Atoms[vs].bond_type = []
            for aa in self.Atoms:
                flag = False
                for ii in range(len(aa.connect)):
                    if aa.connect[ii] == vs:
                        flag = True
                        break
                if flag:
                    del aa.connect[ii]
                    del aa.bond_type[ii]

    def create_vs_coor(self, aa, setting_type="vs_setting"):
        acoor = [aa.coor]
        for ca in aa.connect:
            acoor.append(self.Atoms[ca].coor)
        vparas = getattr(aa, setting_type)["paras"]
        vstyle = getattr(aa, setting_type)["style"]
        vbinc = getattr(aa, setting_type)["binc"]
        vcoor = vs_gen.vs_generate_coor(vstyle, acoor, vparas)
        vcoor = np.array(vcoor)
        sp = vcoor.shape
        vcoor = list(vcoor)
        if len(sp) == 1:
            vcoor = [list(vcoor)]
        else:
            vcoor = list(vcoor)
        return vcoor, vstyle, vparas, vbinc

    def create_vs_atom(self, aa, vs_info, m1_vs, m2_vs):
        for rr in vs_info[0]:
            self.add_atom(style="vs")
            self.Atoms[-1].coor = rr
            self.Atoms[-1].elem = "Bq"
            self.Atoms[-1].atom_type_name = m1_vs[0]
            self.Atoms[-1].nonb_atom_type = m1_vs[0]
            self.Atoms[-1].binc_atom_type = m1_vs[0]
            self.Atoms[-1].para = [0.0, 0.0]
            if m2_vs[0] is not None:
                self.Atoms[-1].atom_type_name_m2 = m2_vs[0]
                self.Atoms[-1].nonb_atom_type_m2 = m2_vs[0]
                self.Atoms[-1].binc_atom_type_m2 = m2_vs[0]
                self.Atoms[-1].mass_m2 = 0.0
                self.Atoms[-1].para_m2 = [0.0, 0.0]
            self.Atoms[-1].ff_charge = m1_vs[1]
            aa.ff_charge = aa.ff_charge - m1_vs[1]
            if m2_vs[1] is not None:
                self.Atoms[-1].ff_charge_m2 = m2_vs[1]
                aa.ff_charge_m2 = aa.ff_charge_m2 - m2_vs[1]
            patoms = [aa.No]
            for an in aa.connect:
                if self.Atoms[an].elem != "D":
                    if hasattr(self.Atoms[an], "atom_type_name_m2"):
                        if not self.Atoms[an].atom_type_name_m2 == "_D":
                            patoms.append(an)
                    else:
                        patoms.append(an)
            # patoms = [aa.No] + [an for an in aa.connect if self.Atoms[an].elem != "D"]
            v = VS(vs_info[1], self.Atoms[-1].No, patoms, vs_info[2])
            v.coor = aa.coor
            self.Vss.append(v)
            self.rfe_exclusions.append([aa.No, self.Atoms[-1].No])
            tmp_an = [aa.No, self.Atoms[-1].No]
            for an0 in aa.connect:
                if an0 not in tmp_an and self.Atoms[an0].elem != "D":
                    self.rfe_exclusions.append([self.Atoms[-1].No, an0])
                    tmp_an.append(an0)
                    for an1 in self.Atoms[an0].connect:
                        if an1 not in tmp_an and self.Atoms[an1].elem != "D":
                            self.rfe_exclusions.append([self.Atoms[-1].No, an1])
                            tmp_an.append(an1)
                            for an2 in self.Atoms[an1].connect:
                                if an2 not in tmp_an and self.Atoms[an2].elem != "D":
                                    self.rfe_exclusions.append([self.Atoms[-1].No, an2])
                                    tmp_an.append(an2)
                                    pp = self.create_Pairs([self.Atoms[-1].No, an2])
                                    self.Pair14.append(pp)

    def create_Pairs(self, term):
        PP = Pair("pair", *term)
        PP.style = "LJ12_6"
        PP.a1_atom_type = self.Atoms[term[0]].atom_type_name
        PP.a1_atom_type_used = self.Atoms[term[0]].nonb_atom_type
        PP.a2_atom_type = self.Atoms[term[1]].atom_type_name
        PP.a2_atom_type_used = self.Atoms[term[1]].nonb_atom_type
        PP.name = f"{PP.a1_atom_type_used}${PP.a2_atom_type_used}"
        PP.scale_factor = 0.5
        PP.charge_para = [self.Atoms[term[0]].ff_charge * 0.8333**0.5, self.Atoms[term[1]].ff_charge * 0.8333**0.5]
        PP.vdw_para = [
            (self.Atoms[term[0]].para[0] + self.Atoms[term[1]].para[1]) * 0.5,
            0.5 * (self.Atoms[term[0]].para[1] * self.Atoms[term[1]].para[1]) ** 0.5,
        ]  # noqa
        if (
            hasattr(self.Atoms[term[0]], "ff_charge_m2")
            or hasattr(self.Atoms[term[0]], "ff_charge_m2")
            or hasattr(self.Atoms[term[1]], "para_m2")
            or hasattr(self.Atoms[term[1]], "para_m2")
        ):  # noqa
            c1 = (
                self.Atoms[term[0]].ff_charge_m2
                if hasattr(self.Atoms[term[0]], "ff_charge_m2")
                else self.Atoms[term[0]].ff_charge
            )
            c2 = (
                self.Atoms[term[1]].ff_charge_m2
                if hasattr(self.Atoms[term[1]], "ff_charge_m2")
                else self.Atoms[term[1]].ff_charge
            )
            vdw1 = self.Atoms[term[0]].para_m2 if hasattr(self.Atoms[term[0]], "para_m2") else self.Atoms[term[0]].para
            vdw2 = self.Atoms[term[1]].para_m2 if hasattr(self.Atoms[term[1]], "para_m2") else self.Atoms[term[1]].para
            PP.charge_para_m2 = [c1 * 0.8333**0.5, c2 * 0.8333**0.5]
            PP.para_m2 = [(vdw1[0] + vdw2[0]) * 0.5, 0.5 * (vdw1[1] * vdw2[1]) ** 0.5]
        return PP

    def create_altered_pair(self, term):
        PP = Pair("pair", *term)
        PP.a1_atom_type = self.Atoms[term[0]].atom_type_name
        PP.a1_atom_type_used = self.Atoms[term[0]].nonb_atom_type
        PP.a2_atom_type = self.Atoms[term[1]].atom_type_name
        PP.a2_atom_type_used = self.Atoms[term[1]].nonb_atom_type
        a1_ff_charge = self.Atoms[term[0]].ff_charge
        a1_ff_charge_m2 = self.Atoms[term[0]].ff_charge_m2 \
            if hasattr(self.Atoms[term[0]], "ff_charge_m2") else self.Atoms[term[0]].ff_charge
        a2_ff_charge = self.Atoms[term[1]].ff_charge
        a2_ff_charge_m2 = self.Atoms[term[1]].ff_charge_m2 \
            if hasattr(self.Atoms[term[1]], "ff_charge_m2") else self.Atoms[term[1]].ff_charge
        PP.ff_charge = [0.8333, a1_ff_charge, a2_ff_charge]
        PP.ff_charge_m2 = [0.8333, a1_ff_charge_m2, a2_ff_charge_m2]
        a1_vdw = self.Atoms[term[0]].para
        a2_vdw = self.Atoms[term[1]].para
        a1_vdw_m2 = self.Atoms[term[0]].para_m2 if hasattr(self.Atoms[term[0]], "para_m2") else self.Atoms[term[0]].para
        a2_vdw_m2 = self.Atoms[term[1]].para_m2 if hasattr(self.Atoms[term[1]], "para_m2") else self.Atoms[term[1]].para
        PP.para = [0.5, (a1_vdw[0] + a2_vdw[0]) * 0.5, (a1_vdw[1] * a2_vdw[1]) ** 0.5]
        PP.para_m2 = [0.5, (a1_vdw_m2[0] + a2_vdw_m2[0]) * 0.5, (a1_vdw_m2[1] * a2_vdw_m2[1]) ** 0.5]
        if hasattr(self, "AlteredPairs") and self.AlteredPairs:
            self.AlteredPairs.append(PP)
        else:
            self.AlteredPairs = [PP]
        return PP

    def create_vs(self):
        """
        生成虚原子
        ###TODO: 该方法需要UT###
        输入：
        输出：
        """
        self.Vss = []
        if not hasattr(self, "rfe_exclusions"):
            self.rfe_exclusions = [[], []]
        for aa in self.Atoms:
            vcoor = None
            vcoor_m2 = None
            if hasattr(aa, "vs_setting"):
                # self.add_atom(style='vs')
                vcoor, vstyle, vparas, vbinc = self.create_vs_coor(aa)
            if hasattr(aa, "vs_setting_m2"):
                vcoor_m2, vstyle_m2, vparas_m2, vbinc_m2 = self.create_vs_coor(aa, setting_type="vs_setting_m2")
            if vcoor is not None:
                if hasattr(aa, "atom_type_name_m2"):
                    if aa.atom_type_name_m2 == aa.atom_type_name:
                        if hasattr(aa, "ff_charge_m2"):
                            self.create_vs_atom(aa, [vcoor, vstyle, vparas], ["Bq", vbinc], ["Bq", vbinc])
                        else:
                            self.create_vs_atom(aa, [vcoor, vstyle, vparas], ["Bq", vbinc], [None, None])
                    else:
                        if vcoor_m2 is not None:
                            self.create_vs_atom(aa, [vcoor, vstyle, vparas], ["Bq", vbinc], ["Bq_dummy", 0.0])
                            self.rfe_exclusions[1].append(self.Atoms[-1].No)
                            self.create_vs_atom(
                                aa, [vcoor_m2, vstyle_m2, vparas_m2], ["Bq_dummy", 0.0], ["Bq", vbinc_m2]
                            )
                            self.rfe_exclusions[0].append(self.Atoms[-1].No)
                            self.rfe_exclusions.append([self.Atoms[-2].No, self.Atoms[-1].No])
                        else:
                            self.create_vs_atom(aa, [vcoor, vstyle, vparas], ["Bq", vbinc], ["Bq_dummy", 0.0])
                            self.rfe_exclusions[1].append(self.Atoms[-1].No)
                else:
                    self.create_vs_atom(aa, [vcoor, vstyle, vparas], ["Bq", vbinc], [None, None])
            elif vcoor_m2 is not None:
                self.create_vs_atom(aa, [vcoor_m2, vstyle_m2, vparas_m2], ["Bq_dummy", 0.0], ["Bq", vbinc_m2])
                self.rfe_exclusions[0].append(self.Atoms[-1].No)
        if len(self.Vss) == 0:
            delattr(self, "Vss")
                    
    def find_side_componend(self, n0, n1, ne=None):
        """
        查找某个flexible bond一边的所有原子
        通过调用self.run_find_side_componend递归完成
        输入：
            n0,n1: int, flexible bond的两个原子
        输出：
            self.tmp_arr: list
        """
        tmp_set = {n1, n0}
        if ne is not None:
            tmp_set.add(ne)
        self.run_find_side_componend(n0, tmp_set)
        tmp_set.remove(n0)
        tmp_set.remove(n1)
        if ne is not None:
            tmp_set.remove(ne)
        return list(tmp_set)
    
    def run_find_side_componend(self, n, tmp_set):
        """
        self.find_side_componend调用的方法
        """
        for aa in self.Atoms[n].connect:
            if aa not in tmp_set:
                tmp_set.add(aa)
                self.run_find_side_componend(aa, tmp_set)

     
