from ..utils.geometry import calc_stru_para

######################################################################################
# 除Constrain, Improper外，所有topol中：                                              #
#    a1,a2,a3,a4 按编号顺序排列                                                       #
#    atom type 顺序与a1 a2 a3 a4一致                                                  #
#    atom type name 按字母顺序排序，即可能与a1 a2 a3 a4的atom type 顺序不一致（或能反序) #
#    Improper,如果给定了atom number, 则按原子量大小排序，但H原子在最后面                 #
#             如果否则按a2 a3 a4 编号排序                                              #
######################################################################################


class Bond:
    """
    定义一根键
    a1: int,第一个原子的编号
    a2: int, 第二个原子的编号
    value: float, 键的距离
    style: 键的力场函数形式，如harmonic, class2等
    para: [float], 键的力场参数
    isfitting: str, 该参数是否被拟合
    a1_atom_type: str, 第一个原子的原子类型 （原始的原子类型）
    a2_atom_type: str, 第二个原子的原子类型 （原始的原子类型）
    a1_atom_type_used: 第一个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a2_atom_type_used: 第二个原子的实际用到的原子类型 （等价表映射后的原子类型）
    """

    def __init__(self, style, atom1, atom2):
        self.s = style
        self.style = style
        if atom1 <= atom2:
            self.a1 = atom1
            self.a2 = atom2
        else:
            self.a1 = atom2
            self.a2 = atom1

    def calc_value(self, this_coor):
        """
        计算键的键长
        输入：
            this_coor: [[float],[float]], 两个原子的坐标
        输出：
            self.value: float, 键的长度，单位A
        """
        self.value = calc_stru_para(this_coor)

    def __repr__(self):
        return f"<Bond: {self.a1} {self.a2}>"

    def __eq__(self, other):
        return {self.a1, self.a2} == {other.a1, other.a2}

    @property
    def str(self):
        return f"<Bond: {self.a1} {self.a2}>"

    @property
    def atoms(self):
        return [self.a1, self.a2]
    
    @property
    def name(self):
        return "-".join([str(an) for an in self.atoms])
    
    @property
    def names(self):
        return [self.name,"-".join([str(an) for an in reversed(self.atoms)])]

    @property
    def atom_type(self):
        return [self.a1_atom_type,self.a2_atom_type]

    @property
    def atom_type_used(self):
        return [self.a1_atom_type_used,self.a2_atom_type_used]

    @property
    def atom_type_names(self):
        return ["$".join(self.atom_type),"$".join(reversed(self.atom_type))]
    
    @property
    def atom_type_name(self):
        return min(self.atom_type_names)

    @property
    def atom_type_used_names(self):
        return ["$".join(self.atom_type_used),"$".join(reversed(self.atom_type_used))]

    @property
    def atom_type_used_name(self):
        if not hasattr(self,"_atom_type_used_name"):
            return min(self.atom_type_used_names)
        else:
            if self._atom_type_used_name is None:
                return min(self.atom_type_used_names)
            else:
                return self._atom_type_used_name
            
    @atom_type_used_name.setter
    def atom_type_used_name(self,ss):
        self._atom_type_used_name = ss

    def reset(self,atoms,atom_type=None,atom_type_used=None):
        self.a1, self.a2 = atoms
        if atom_type is not None:
            self.a1_atom_type, self.a2_atom_type = atom_type
        if atom_type_used is not None:
            self.a1_atom_type_used, self.a2_atom_type_used = atom_type_used
        #    self.name = "$".join(atom_type_used)
        #else:
        #    self.name = "$".join(atom_type)

    def get_type(self, mol):
        #try:
        #    return mol.Atoms[self.a1].bond_type_aromatic[mol.Atoms[self.a1].connect.index(self.a2)]
        #except:
        return mol.Atoms[self.a1].bond_type[mol.Atoms[self.a1].connect.index(self.a2)]

    def get_type_aromatic(self, mol):
        #try:
        return mol.Atoms[self.a1].bond_type_aromatic[mol.Atoms[self.a1].connect.index(self.a2)]
        #except:
        #   return mol.Atoms[self.a1].bond_type[mol.Atoms[self.a1].connect.index(self.a2)]

    def get_type_detail(self, mol):
        return mol.Atoms[self.a1].connectivity_type[mol.Atoms[self.a1].connect.index(self.a2)]

class Angle:
    """
    定义一个键角
    a1: int,第一个原子的编号
    a2: int, 第二个原子的编号，角中心的原子
    a3: int, 第二个原子的编号
    value: float, 键角的角度，用弧度表示，特殊情况下还会包括其他的数值，详细见calc_value方法中的说明
    value_a: float, 键角的角度，用度表示
    style: 键角的力场函数形式，如harmonic, class2等
    para: [float], 键角的力场参数
    isfitting: str, 该参数是否被拟合
    a1_atom_type: str, 第一个原子的原子类型 （原始的原子类型）
    a2_atom_type: str, 第二个原子的原子类型 （原始的原子类型）
    a3_atom_type: str, 第三个原子的原子类型 （原始的原子类型）
    a1_atom_type_used: 第一个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a2_atom_type_used: 第二个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a3_atom_type_used: 第三个原子的实际用到的原子类型 （等价表映射后的原子类型）
    """

    def __init__(self, style, atom1, atom2, atom3):
        self.s = style
        self.style = style
        self.a2 = atom2
        if atom1 <= atom3:
            self.a1 = atom1
            self.a3 = atom3
        else:
            self.a1 = atom3
            self.a3 = atom1

    def calc_value(self, this_coor):
        """
        计算键角的结构性质，体数值与力场函数形式有关
        输入：
            this_coor: [[float],[float],[float]], 三个原子的坐标
        输出：
            self.value_a: float, 键角的角度，用度表示
            self.value:
                        float, 键角，单位弧度
                        float(r12), 1号和2号原子间的距离（style=class2）
                        float(r12), 2号和3号原子间的距离 （style=class2）
                        float(r13), 1号和3号原子间的距离 （style=urey_bradley, charmm, or sdk）
        """
        self.value = calc_stru_para(this_coor)
        if hasattr(self,"pstyle") and self.pstyle == "class2":
            r12 = calc_stru_para([this_coor[0], this_coor[1]])
            r23 = calc_stru_para([this_coor[1], this_coor[2]])
            self.value = [self.value, r12, r23]
        elif self.style in ["urey_bradley", "charmm", "sdk"]:
            r13 = calc_stru_para([this_coor[0], this_coor[2]])
            self.value = [self.value, r13]
        else:
            self.value_r = self.value * 3.1415926 / 180.0

    def __repr__(self):
        return f"<Angle: {self.a1} {self.a2} {self.a3}>"

    def __eq__(self, other):
        return self.a2 == other.a2 and {self.a1, self.a3} == {other.a1, other.a3}

    @property
    def str(self):
        return f"<Angle: {self.a1} {self.a2} {self.a3}>"

    @property
    def atoms(self):
        return [self.a1, self.a2, self.a3]
    
    @property
    def name(self):
        return "-".join([str(an) for an in self.atoms])
    
    @property
    def names(self):
        return [self.name,"-".join([str(an) for an in reversed(self.atoms)])]

    @property
    def atom_type(self):
        return [self.a1_atom_type,self.a2_atom_type,self.a3_atom_type]
    
    @property
    def atom_type_used(self):
        return [self.a1_atom_type_used,self.a2_atom_type_used,self.a3_atom_type_used]
    
    @property
    def atom_type_names(self):
        return ["$".join(self.atom_type),"$".join(reversed(self.atom_type))]

    @property
    def atom_type_name(self):
        return min(self.atom_type_names)

    @property
    def atom_type_used_names(self):
        return ["$".join(self.atom_type_used),"$".join(reversed(self.atom_type_used))]

    @property
    def atom_type_used_name(self):
        if not hasattr(self,"_atom_type_used_name"):
            return min(self.atom_type_used_names)
        else:
            if self._atom_type_used_name is None:
                return min(self.atom_type_used_names)
            else:
                return self._atom_type_used_name
            
    @atom_type_used_name.setter
    def atom_type_used_name(self,ss):
        self._atom_type_used_name = ss

    def reset(self,atoms,atom_type=None,atom_type_used=None):
        self.a1, self.a2, self.a3 = atoms
        if atom_type is not None:
            self.a1_atom_type, self.a2_atom_type, self.a3_atom_type = atom_type
        if atom_type_used is not None:
            self.a1_atom_type_used, self.a2_atom_type_used, self.a3_atom_type_used = atom_type_used
        #    self.name = "$".join(atom_type_used)
        #else:
        #    self.name = "$".join(atom_type)

class Dihedral:
    """
    定义一个二面角
    a1: int,第一个原子的编号，与2号原子相连
    a2: int, 第二个原子的编号，中间的原子
    a3: int, 第二个原子的编号，中间的原子
    a4: int,第四个原子的编号，与3号原子相连
    value: float, 二面角的角度，用弧度表示，特殊情况下还会包括其他的数值，详细见calc_value方法中的说明
    value_a: float, 二面角的角度，用度表示
    style: 键角的力场函数形式，如fourier，amber，opls, class2等
    para: [float], 键角的力场参数
    isfitting: str, 该参数是否被拟合
    a1_atom_type: str, 第一个原子的原子类型 （原始的原子类型）
    a2_atom_type: str, 第二个原子的原子类型 （原始的原子类型）
    a3_atom_type: str, 第三个原子的原子类型 （原始的原子类型）
    a4_atom_type: str, 第四个原子的原子类型 （原始的原子类型）
    a1_atom_type_used: 第一个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a2_atom_type_used: 第二个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a3_atom_type_used: 第三个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a4_atom_type_used: 第四个原子的实际用到的原子类型 （等价表映射后的原子类型）
    """

    def __init__(self, style, atom1, atom2, atom3, atom4):
        self.s = style
        self.style = style
        if atom1 <= atom2:
            self.a1 = atom1
            self.a2 = atom2
            self.a3 = atom3
            self.a4 = atom4
        else:
            self.a1 = atom4
            self.a2 = atom3
            self.a3 = atom2
            self.a4 = atom1

    def calc_value(self, this_coor):
        """
        计算二面角的结构性质，具体数值与力场函数形式有关
        输入：
            this_coor: [[float],[float],[float],[float]], 四个原子的坐标
        输出：
            self.value_a: float, 二面角的角度，用度表示
            self.value:
                        float, 二面角的度数，单位弧度
                        float(r12), 1号和2号原子间的距离（style=class2）
                        float(r12), 2号和3号原子间的距离 （style=class2）
                        float(r34), 3号和4号原子间的距离 （style=class2）
                        float(a123), 1号,2号和3号原子形成角的度数（弧度） （style=class2）
                        float(a234), 2号,3号和4号原子形成角的度数（弧度） （style=class2）
                        float(a123_a), 1号,2号和3号原子形成角的度数 （style=class2）
                        float(a234_a), 2号,3号和4号原子形成角的度数 （style=class2）
        """
        self.value = calc_stru_para(this_coor)

        # linear dihedral
        a123 = calc_stru_para([this_coor[0], this_coor[1], this_coor[2]])
        a234 = calc_stru_para([this_coor[1], this_coor[2], this_coor[3]])
        self.is_linear = a123 > 179 or a234 > 179

        if hasattr(self,"pstyle") and self.pstyle == "class2":
            a123_r = a123 * 3.1415926 / 180.0
            a234_r = a234 * 3.1415926 / 180.0
            r12 = calc_stru_para([this_coor[0], this_coor[1]])
            r23 = calc_stru_para([this_coor[1], this_coor[2]])
            r34 = calc_stru_para([this_coor[2], this_coor[3]])
            self.value = [self.value, r12, r23, r34, a123, a234]
        else:
            self.value_r = self.value * 3.1415926 / 180.0

    def __repr__(self):
        return f"<Dihedral: {self.a1} {self.a2} {self.a3} {self.a4}>"

    def __eq__(self, other):
        return self.atoms == other.atoms or list(reversed(self.atoms)) == other.atoms

    @property
    def str(self):
        return f"<Dihedral: {self.a1} {self.a2} {self.a3} {self.a4}>"

    @property
    def atoms(self):
        return [self.a1, self.a2, self.a3, self.a4]
    
    @property
    def atom_type(self):
        return [self.a1_atom_type,self.a2_atom_type,self.a3_atom_type,self.a4_atom_type]
    
    @property
    def atom_type_used(self):
        return [self.a1_atom_type_used,self.a2_atom_type_used,self.a3_atom_type_used,self.a4_atom_type_used]

    @property
    def atom_type_names(self):
        return ["$".join(self.atom_type),"$".join(reversed(self.atom_type))]
    
    @property
    def atom_type_name(self):
        return min(self.atom_type_names)
    
    @property
    def atom_type_used_names(self):
        return ["$".join(self.atom_type_used),"$".join(reversed(self.atom_type_used))]

    @property
    def atom_type_used_name(self):
        if not hasattr(self,"_atom_type_used_name"):
            return min(self.atom_type_used_names)
        else:
            if self._atom_type_used_name is None:
                return min(self.atom_type_used_names)
            else:
                return self._atom_type_used_name
            
    @atom_type_used_name.setter
    def atom_type_used_name(self,ss):
        self._atom_type_used_name = ss

    @property
    def name(self):
        return "-".join([str(an) for an in self.atoms])
    
    @property
    def names(self):
        return [self.name,"-".join([str(an) for an in reversed(self.atoms)])]

    def reset(self,atoms,atom_type=None,atom_type_used=None):
        self.a1, self.a2, self.a3, self.a4 = atoms
        if atom_type is not None:
            self.a1_atom_type, self.a2_atom_type, self.a3_atom_type, self.a4_atom_type = atom_type
        if atom_type_used is not None:
            self.a1_atom_type_used, self.a2_atom_type_used, self.a3_atom_type_used, self.a4_atom_typer_used = atom_type_used
        #    self.name = "$".join(atom_type_used)
        #else:
        #    self.name = "$".join(atom_type)

class Improper:
    """
    定义一个improper二面角
    a1: int,第一个原子的编号, 中心原子
    a2: int, 第二个原子的编号
    a3: int, 第二个原子的编号
    a4: int,第四个原子的编号
    value: float, 二面角的角度，用弧度表示，特殊情况下还会包括其他的数值，详细见calc_value方法中的说明
    value_a: float, 二面角的角度，用度表示
    style: 键角的力场函数形式，如fourier，amber，opls, harmonic, cvff, charmm,class2等
    para: [float], 键角的力场参数
    isfitting: str, 该参数是否被拟合
    a1_atom_type: str, 第一个原子的原子类型 （原始的原子类型）
    a2_atom_type: str, 第二个原子的原子类型 （原始的原子类型）
    a3_atom_type: str, 第三个原子的原子类型 （原始的原子类型）
    a4_atom_type: str, 第四个原子的原子类型 （原始的原子类型）
    a1_atom_type_used: 第一个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a2_atom_type_used: 第二个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a3_atom_type_used: 第三个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a4_atom_type_used: 第四个原子的实际用到的原子类型 （等价表映射后的原子类型）
    """

    def __init__(self, style, atom1, atom2, atom3, atom4,atom_numbers=None,ignore_order=False):
        self.s = style
        self.style = style
        self.a1 = atom1
        if ignore_order:
            self.a2,self.a3,self.a4 = atom2,atom3,atom4
        else:
            if atom_numbers is None:
                self.a2,self.a3,self.a4 = sorted([atom2,atom3,atom4])
            else:
                self.a2,self.a3,self.a4 = self.__order_by_atom_number([[atom_numbers[0],atom2],[atom_numbers[1],atom3],[atom_numbers[2],atom4]])

    def __order_by_atom_number(self,_tmp_):
            _tmp_ = sorted(_tmp_)
            return [rr[1] for rr in _tmp_ if rr[0] != 1] + [rr[1] for rr in _tmp_ if rr[0] == 1]

    def calc_value(self, this_coor):
        """
        计算二面角的结构性质，具体数值与力场函数形式有关
        输入：
            this_coor: [[float],[float],[float],[float]], 四个原子的坐标
        输出：
            self.value_a: float, 键角的角度，用度表示
            self.value:
                        float, improper的二面角，单位弧度
                        float(a213), 2号,1号和3号原子形成角的度数（弧度） （style=class2）
                        float(a314), 3号,1号和4号原子形成角的度数（弧度） （style=class2）
                        float(a214), 2号,1号和4号原子形成角的度数（弧度） （style=class2）
                        float(a213_a), 2号,1号和3号原子形成角的度数 （style=class2）
                        float(a314_a), 3号,1号和4号原子形成角的度数 （style=class2）
                        float(a214_a), 2号,1号和4号原子形成角的度数 （style=class2）
        """
        # for fourier_2n, the function is K(1+cos(2*aijkl-pi)), i is center atom
        # if you want use the function of K(1+cos(2*aijkl)), k is center atom
        self.value = calc_stru_para(this_coor)
        self.value_r = self.value * 3.1415926 / 180.0
        if hasattr(self,"pstyle") and self.pstyle == "class2":
            # there are some error, may be
            value_ = calc_stru_para([this_coor[1], this_coor[0], this_coor[2], this_coor[3]])
            value_r_ = value_ * 3.1415926 / 180.0

            a213 = calc_stru_para([this_coor[1], this_coor[0], this_coor[2]])
            a314 = calc_stru_para([this_coor[2], this_coor[0], this_coor[3]])
            a214 = calc_stru_para([this_coor[1], this_coor[0], this_coor[3]])
            a213_r = a213 * 3.1415926 / 180.0
            a314_r = a314 * 3.1415926 / 180.0
            a214_r = a214 * 3.1415926 / 180.0
            self.value = [value_, a213, a314, a214]

    def __repr__(self):
        return f"<Improper: {self.a1} {self.a2} {self.a3} {self.a4}>"

    def __eq__(self, other):
        return self.a1 == other.a1 and {self.a2, self.a3, self.a4} == {other.a2, other.a3, other.a4}

    @property
    def str(self):
        return f"<Improper: {self.a1} {self.a2} {self.a3} {self.a4}>"

    @property
    def atoms(self):
        return [self.a1, self.a2, self.a3, self.a4]

    @property
    def name(self):
        return "-".join([str(an) for an in reversed(self.atoms)])

    @property
    def names(self):
        return [self.name,
                "-".join([str(an) for an in [self.a1,self.a2,self.a4,self.a3]]),
                "-".join([str(an) for an in [self.a1,self.a3,self.a2,self.a4]]),
                "-".join([str(an) for an in [self.a1,self.a3,self.a4,self.a2]]),
                "-".join([str(an) for an in [self.a1,self.a4,self.a2,self.a3]]),
                "-".join([str(an) for an in [self.a1,self.a4,self.a3,self.a2]]),
                ]

    @property
    def atom_type(self):
        return [self.a1_atom_type,self.a2_atom_type,self.a3_atom_type,self.a4_atom_type]
    
    @property
    def atom_type_used(self):
        return [self.a1_atom_type_used,self.a2_atom_type_used,self.a3_atom_type_used,self.a4_atom_type_used]
    
    @property
    def atom_type_names(self):
        return [
                "$".join(self.atom_type),
                "$".join([self.a1_atom_type,self.a2_atom_type,self.a4_atom_type,self.a3_atom_type]),
                "$".join([self.a1_atom_type,self.a3_atom_type,self.a2_atom_type,self.a4_atom_type]),
                "$".join([self.a1_atom_type,self.a3_atom_type,self.a4_atom_type,self.a2_atom_type]),
                "$".join([self.a1_atom_type,self.a4_atom_type,self.a2_atom_type,self.a3_atom_type]),
                "$".join([self.a1_atom_type,self.a4_atom_type,self.a3_atom_type,self.a2_atom_type]),
                ]
    
    @property
    def atom_type_name(self):
        return min(self.atom_type_names)

    @property
    def atom_type_used_names(self):
        return [
                "$".join(self.atom_type_used),
                "$".join([self.a1_atom_type_used,self.a2_atom_type_used,self.a4_atom_type_used,self.a3_atom_type_used]),
                "$".join([self.a1_atom_type_used,self.a3_atom_type_used,self.a2_atom_type_used,self.a4_atom_type_used]),
                "$".join([self.a1_atom_type_used,self.a3_atom_type_used,self.a4_atom_type_used,self.a2_atom_type_used]),
                "$".join([self.a1_atom_type_used,self.a4_atom_type_used,self.a2_atom_type_used,self.a3_atom_type_used]),
                "$".join([self.a1_atom_type_used,self.a4_atom_type_used,self.a3_atom_type_used,self.a2_atom_type_used]),
                ]

    @property
    def atom_type_used_name(self):
        if not hasattr(self,"_atom_type_used_name"):
            return min(self.atom_type_used_names)
        else:
            if self._atom_type_used_name is None:
                return min(self.atom_type_used_names)
            else:
                return self._atom_type_used_name
            
    @atom_type_used_name.setter
    def atom_type_used_name(self,ss):
        self._atom_type_used_name = ss


    def reset(self,atoms,atom_type=None,atom_type_used=None):
        self.a1, self.a2, self.a3, self.a4 = atoms
        if atom_type is not None:
            self.a1_atom_type, self.a2_atom_type, self.a3_atom_type, self.a4_atom_type = atom_type
        if atom_type_used is not None:
            self.a1_atom_type_used, self.a2_atom_type_used, self.a3_atom_type_used, self.a4_atom_typer_used = atom_type_used
            #self.name = "$".join(atom_type_used)
        #else:
        #    self.name = "$".join(atom_type)

class Pair:
    """
    定义一个pair相互作用对
    a1: int,第一个原子的编号
    a2: int, 第二个原子的编号
    value: float, 两个原子间的距离
    style: 键的力场函数形式，如LJ12_6等
    para: [float], pair的vdw参数
    charge_para: [float], 两个原子的点电荷
    isfitting: str, 该参数是否被拟合
    a1_atom_type: str, 第一个原子的原子类型 （原始的原子类型）
    a2_atom_type: str, 第二个原子的原子类型 （原始的原子类型）
    a1_atom_type_used: 第一个原子的实际用到的原子类型 （等价表映射后的原子类型）
    a2_atom_type_used: 第二个原子的实际用到的原子类型 （等价表映射后的原子类型）
    """

    def __init__(self, style, atom1, atom2):
        self.s = style
        if atom1 <= atom2:
            self.a1 = atom1
            self.a2 = atom2
        else:
            self.a1 = atom2
            self.a2 = atom1

    def calc_value(self, this_coor):
        """
        计算两个原子间的距离
        输入：
            this_coor: [[float]], 两个原子的坐标
        输出：
            self.value: float, 键的长度，单位A
        """
        self.value = calc_stru_para(this_coor)

    def __repr__(self):
        return f"<Pair: {self.a1} {self.a2}>"

    def __eq__(self, other):
        return {self.a1, self.a2} == {other.a1, other.a2}

    @property
    def str(self):
        return f"<Pair: {self.a1} {self.a2}>"

    @property
    def atoms(self):
        return [self.a1, self.a2]
    
    @property
    def name(self):
        return "-".join([str(an) for an in self.atoms])
    
    @property
    def names(self):
        return [self.name,"-".join([str(an) for an in reversed(self.atoms)])]

    @property
    def atom_type(self):
        return [self.a1_atom_type,self.a2_atom_type]
    
    @property
    def atom_type_used(self):
        return [self.a1_atom_type_used,self.a2_atom_type_used]
    
    @property
    def atom_type_names(self):
        return ["$".join(self.atom_type),"$".join(reversed(self.atom_type))]
    
    @property
    def atom_type_name(self):
        return min(self.atom_type_names)

    @property
    def atom_type_used_names(self):
        return ["$".join(self.atom_type_used),"$".join(reversed(self.atom_type_used))]

    @property
    def atom_type_used_name(self):
        if not hasattr(self,"_atom_type_used_name"):
            return min(self.atom_type_used_names)
        else:
            if self._atom_type_used_name is None:
                return min(self.atom_type_used_names)
            else:
                return self._atom_type_used_name
            
    @atom_type_used_name.setter
    def atom_type_used_name(self,ss):
        self._atom_type_used_name = ss

    def reset(self,atoms,atom_type=None,atom_type_used=None):
        self.a1, self.a2 = atoms
        if atom_type is not None:
            self.a1_atom_type, self.a2_atom_type = atom_type
        if atom_type_used is not None:
            self.a1_atom_type_used, self.a2_atom_type_used = atom_type_used
            #self.name = "$".join(atom_type_used)
        #else:
        #    self.name = "$".join(atom_type)


class Constrain:
    def __init__(self, atoms, fix_value):
        #self.atoms = atoms
        self.fix_value = fix_value
        if len(atoms) == 2:
            self.a1 = atoms[0]
            self.a2 = atoms[1]
            self.style = "bond"
        elif len(atoms) == 3:
            self.a1 = atoms[0]
            self.a2 = atoms[1]
            self.a3 = atoms[2]
            self.style = "angle"
        elif len(atoms) == 4:
            self.a1 = atoms[0]
            self.a2 = atoms[1]
            self.a3 = atoms[2]
            self.a4 = atoms[3]
            self.style = "dihedral"

    def calc_value(self, this_coor):
        self.value = calc_stru_para(this_coor)
        if len(this_coor) > 2:
            self.value_r = self.value * 3.1415926 / 180.0

    @property
    def atoms(self):
        return [getattr(self,f"a{ii}") for ii in range(1,5) if hasattr(self,f"a{ii}")]

    @property
    def name(self):
        return "-".join([str(an) for an in self.atoms])
    
    @property
    def names(self):
        return [self.name,"-".join([str(an) for an in reversed(self.atoms)])]
