
# Physical constants
HF_TO_KCAL_PER_MOL = 627.509
BOHR_TO_ANGSTROM = 0.529177249
ONE_4PI_EPS0 = 138.93545522028575  # kJ/mol * nm / e^2

AVOGADRO_CONSTANT = 6.02214076e23  # /mol
BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
MOLAR_GAS_CONSTANT = 8.3144598  # AVOGADRO_CONSTANT * BOLTZMANN_CONSTANT  # 8.314462 J/mol/K
DALTON_CONSTANT = 1.6605390666e-27  # kg
CAL_TO_J = 4.184

#molecule strucutre constants
HBOND_DIST_MAX = 4.1
HBOND_ANGLE_MIN = 100.0

# default file setting


# position: L: link; T: terminal; M: medium, C: cycle; F: fused; B: bonded,
#           S: spiro
# 
# TA (terminal atom), CE(end atom of chain), CM(medium atom of chain), CER(end atom chain wiht ring)
# CMR(midium atom of chain with ring), R(ring atom), RC(ring atom with chain), RH(head atom of ring bridge),
# RF(medium atom of fused ring),RFM(medium atom of fused ring), RS(spiro atom)
#position_TA = ["TA"]
#position_chain = ["TA", "CE", "CM", "CER", "CMR"]
#position_chaint = position_TA + position_chain
#position_ring_fused = ["RFM", "RF", "RFMH", "RFH", "RFMC", "RFC"]
#position_ring_bridge = ["RH", "RHC", "RFMH", "RFH", "RSH"]
#position_ring = ["R","RC","RH","RHC","RFM","RFMH","RFMC","RF",
#                "RFH","RFC","RS","RSH","RSC",]
#position = position_chaint + position_ring

LOCAL_EN = ["EN","LT","LTC",]
LOCAL_LINK = ["EN", "LT", "LM", "LTC", "LMC","LBC"]
LOCAL_LINK_T = ["LMC","LM","LT","LTC","EN","LBC"]
LOCAL_CYCLE_F = ["CFML","CFM", "CF", "CFMB", "CFB", "CFL"]
LOCAL_CYCLE_B = ["CFB","CB", "CBL", "CFMB", "CSB"]
LOCAL_CYCLE = ["C","CB","CL","CS","CF","CBL","CFM",
                "CFB","CFL","CSB","CSL","CFMB","CFML",]
LOCAL = LOCAL_LINK + LOCAL_CYCLE

# bond type: S (single bond in chain), 
# D (double bond in chain), T(triple bond in chain), M(multiple bond in chain)
# J(conjugation bond), O (coordinated bond), 
# s (single bond in ring), d (double bond in ring), 
# t (triple bond in ring), # noqa
# r(aromatic bond), j (conjugation bond in ring), 
# fs (fused single bond),fd(fused double bond),fr (fused aromatic bond), # noqa
# fj(fused conjugation bond), ft (fused triple bond),
# bs (single bridge bond),bd(double bridge bond),es(exocyclic single bond), # noqa
# ed(exocyclic double bond)
#single_bond = ["S", "eS", "fs", "bS", "s"]
#double_bond = ["D", "eD", "fd", "bD", "d"]
#triple_bond = ["T", "ft", "t"]
#arom_bond = ["r", "fr"]
#conju_bond = ["J", "eJ", "fj", "bJ", "j"]
#multi_bond = double_bond + triple_bond  # +arom_bond
#break_bond = ["S", "eS", "bS", "br", "J"]

#connectivity: s: sigma; c: cycle; o: out; f: fused; h: head;
#              p: pi; d: detal; a: aromatic; u: conju; 
CONNECT_TYPE = {
    "SINGLE_CONNECT": ["s","cs","os","fcs","hs"],
    "DOUBLE_CONNECT": ["p","op","fcp","hp","cp"],
    "TRIPLE_CONNECT": ["d","fcd","cd"],
    #"MULTI_CONNECT":[], 
    "AROMATIC_CONNECT": ["a","fa","ha"],
    "CONJU_CONNECT": ["u","ou","fcu","hu","cu"],
    "CUT_CONNECT": ["s","os","hs","ha","u"]
}

#single_bond = ["S", "eS", "fs", "bS", "s"]
#double_bond = ["D", "eD", "fd", "bD", "d"]
#triple_bond = ["T", "ft", "t"]
#arom_bond = ["r", "fr"]
#conju_bond = ["J", "eJ", "fj", "bJ", "j"]
#multi_bond = double_bond + triple_bond  # +arom_bond
#break_bond = ["S", "eS", "bS", "br", "J"]

SINGLE_BOND = ["S", "eS", "fs", "bS", "s"]
DOUBLE_BOND = ["D", "eD", "fd", "bD", "d"]
TRIPLE_BOND = ["T", "ft", "t"]
AROM_BOND = ["r", "fr"]
CONJU_BOND = ["J", "eJ", "fj", "bJ", "j"]
MULTI_BOND = DOUBLE_BOND + TRIPLE_BOND  # +arom_bond
BREAK_BOND = ["S", "eS", "bS", "br", "J"]
RING_BOND = ["s","fs","fd","d","r","fr","j","fj","t","ft"]
CHAIN_BOND = ["S","eS","bS","D","eD","bD","T","J","eJ","bJ"]
# break_score=[1,2,3,4,5,6,7,8,9,10,11,12,13]
# old setting
# score_type=[[1,4,7,8,10],[2,3,5,6,9,11,12,13],[30,31,32],[50,51,52,53,54,55],[60,61,62,63,64,65,66,67]]
score_type = [
    [1, 4, 7, 8, 10],
    [2, 3, 5, 6, 9, 11, 12, 13],
    [30, 31, 32],
    [50, 51, 52, 53, 54, 55],
    [60, 61, 62, 63, 64, 65, 66, 67],
]
small_ring = [3, 4, 5, 6, 7, 8, 9]


#single_connect_atom_or_group = ["H", "F", "Cl", "Br"]
SINGLE_CONNECT_ATOM_OR_GROUP = ["H", "F", "Cl", "Br"]
# single_side_chain_type=[2,9,11,20,26,27]
# double_side_chain_type=[3,12,13,21,28,29]
# one_terminal_chain_type=[4,8,15,22,25,31]
# two_terminal_chain_type=[6,18,19,24,34,35]
# one_terminal_one_side_chain_type=[5,14,16,17,23,30,32,33]
# terminal_chain_type = one_terminal_chain_type+two_terminal_chain_type+one_terminal_one_side_chain_type
# side_chain_type = single_side_chain_type+double_side_chain_type+one_terminal_one_side_chain_type

DEFAULT_TOPOLS_TO_MTX_ATTRIBUTES = {"Atoms":["ID","element","atom_number","formal_charge","atom_name","local","plate",
                                  "connectivity","bond_type","coordinates","bond_type_aromatic","connecticity_type",
                                  "has_ring","has_ring_size","has_ring_property",
                                  "partial_formal_charge","point_charge","charge_group",
                                  "residue","residue_id",
                                  "atom_type_name","nonb_atom_type","binc_atom_type","atc_atom_type","atom_type_ID",
                                  "binc_str","binc_parameters","binc_style","binc_score",
                                  "vdw_str","parameters","pstyle"
                                 ],
                                  "Bonds":["a1","a2","value","style","para"],
                                  "Angles":["a1","a2","a3","value","style","para"],
                                  "Dihedrals":["a1","a2","a3","a4","value","style","para"],
                                  "Impropers":["a1","a2","a3","a4","value","style","para"],
                                  "Pair12":["a1","a2","value","style","para"],
                                  "Pair13":["a1","a2","value","style","para"],
                                  "Pair14":["a1","a2","value","style","para"],
                                  "Pair1n":["a1","a2","value","style","para"],
                                  "AlteredPair":["a1","a2","value","style","para"],
                                 }

DEFAULT_MOLECULE_ATTR_TO_CSV = [
    "smiles","inchi_key","inchi","name","mole_name",
    "iupac_name","nick_name","drug_name","internal_name",
    "mass","formula","heavy_atoms","net_charge","multiple","element_set",
    "function_group","function_group_label",
    "torsions","torsion_number","constrain_term","constrain_value",
    "rings","ring_number","ring_size","ring_property","ring_blocks",
    "energy","dipole","inertia","density",
    "source","count","topol_label","frag_type",
]

DEFAULT_MOLECULE_TO_MTX_ATTRIBUTES = DEFAULT_MOLECULE_ATTR_TO_CSV + [
    "force","hessain","frequency",
    "elem_frag","sssr_frag","csf_frag","rsf_frag",
    "chain_frag","scaffold_frag","tf_frag",
    "seco_frag","sketch_frag",
    ]


OLD_ATOM_DEFAULT_ATTRIBUTES = [
        "No",  # ID int, 原子编号，如 0 1 2 3等
        "elem",  # element,str, 元素符号, 如 C,H,O，N等
        "atom_number",  # 同number
        "number", #atom_number, 原子序数，如碳6，氧8
        "mass",  # atom_mass, float, 原子量
        "name",  # atom_name, str， 命称，如C1，C2, H1,O1等
        "formal_charge",  # 净电荷，如1，2，-1，-2等
        "primitive_formal_charge",  #partial_formal_charge float, 原子的净电荷在共轭体系上的平均值
        "atom_type_name",  # str，原子类型的名称，如c_4h3, o_1等
        "nonb_atom_type",  # str, 提取vdw参数时用到的原子类型
        "binc_atom_type",  # str, 提取binc参数时用到的原子类型
        "atc_atom_type",  # str, 提取atc参数时用到的原子类型
        "atom_type_ID",  # int, 原子类型的编号，lammps软件用到
        "mole_ID",  # molecule_ID int, 原子所在分子的分子编号，lammps软件用到
        "mole_type",  #molecule_type str, 原子所在分子的分子类型，gromacs软件会用到
        "residu",  #residue str, 原子所在残基的类型
        "residu_number",  # residue_id int, 原子所在残基的编号
        "charge_group",  # int, 原子所在的chrage group，gromacs软件专用
       
        "charge",  # dict, 原子的各种partial charge, 如{“esp": 0.123, "mullien": 0.231}
        # 力场参数相关
        "ff_charge",  # point_charge float, 原子的力场电荷
        "ff_charge_base",  #point_charge_base float, Assign原子力场电荷时的基准值
        "binc_style",  # List[str], binc参数的类型
        "binc_tag",  # binc_str List[str], binc的参数的标记。因一个原子会有多个相关联的binc参数，所以是个list
        "binc_score",  # List[float] binc参数的评价
        "binc_para",  # binc_parameters List[float], binc的参数
        #"style",  # pstyle str, vdw的函数形式
        "para",  # vdw_parameters List[float] vdw的参数
        "tag",  # vdw_str str, vdw参数的标记
        "pscore",  #vdw_score vdw参数的评价分数
        # #####
        "local",  # str, 原子的位置描述，详细说明见48-50行,之前是position
        "inring",  #has_ring [str], 原子所在环的标记符
        "ring_size",  #has_ring_size [int],原子所在环的大小，如[6]表示原子在一个6元环中，
        #              [6,5]表示原子在一个6元环和5元环组成的稠环中，如果原子不在任何环中为[]
        "ring_prop",  #has_ring_property [str],原子所在环的环性质, 数组size与ring_size一样。环的性质为：nonar, ar1,ar2,ar3,ar4,ar5
        "coor",  # coordinates [float, float, float], 原子的坐标，单位A,
        "connect",  # connectivity [int], 和该原子相连的其他原子，如[2,3,4]，表示该原子与2，3，4号原子相连
        "bond_type",  # [str], 该原子的键类型，列表的size与connect一致。
        #              如果connect为[2,3,4],bond_type为["ar","ar","1"],表示与2号，3号，4号原子分别形成芳香，芳香和单键
        "bond_type_detail",  # connectivity_type [str] 更详细的键类型表示方式,详细说明见51-55行
        "bond_type_old",  # bond_type_aromatic [str] 以纯粹的单，双，三键的方式记录键类型。即键类型中不包括ar键类型
        "plate",  # yes or no， 该原子是否在一个平面中，决定是否能形成improper项
        "vlocity",  # 原子的速度，（目前没有用到）
        "force",  # 原子受到的力，（目前没有用到)
        "chirality",  # 原子的手性，（目前没有用到）
        "cis_trans",  # 原子的顺反异构，（目前没有用到）
        "lone_pair",  # 原子孤对电子对的数目，（目前没有用到）
        "hydrogen_bond",  # 原子在形成氢键时是受体还是供体，（目前没有用到）
        "polaribility",  # 原子的极性，（目前没有用到）
        "hydrophibic_philic",  #hydrophibic and hydrophilic 原子的亲水或憎水性，（目前没有用到）
                ]

#note: 目前bond_type不再因为aromatic而改变，等于以前的bond_type_old。新增加的bond_type_aromatic等于以前aromatic后的bond_type
#      之前的style变成pstyle
ATOM_DEFAULT_ATTRIBUTES = ["No", "elem",              "number",     "mass",   "primitive_formal_charge",
                           "ID", "element","elements","atom_number","atom_mass","atom_name","partial_formal_charge",
    "mole_ID",   "mole_type",    "residu", "residu_number","ff_charge",   "ff_charge_base","binc_tag","binc_para",
    "molecle_ID","molecule_type","residue","residue_id",   "point_charge","point_charge_base","binc_str","binc_parameters",
    "para",          "tag",    "pscore",   "position", "inring",                  "ring_prop",        "coor",       "connect",     "bond_type_detail",
    "vdw_parameters","vdw_str","vdw_score","local",    "has_ring","has_ring_size","has_ring_property","coordinates","connectivity","connectivity_type",
    
    "formal_charge", "atom_type_name","nonb_atom_type","binc_atom_type","atc_atom_type","atom_type_ID","charge_group","charge","binc_style","binc_score",
    "pstyle","bond_type","bond_type_aromatic","plate","vlocity","atom_force","chirality","cis_trans","lone_pair",
    "hydrogen_bond","polaribility","hydrophibic_philic",
    "esp_charge","mulliken_charge","am1bcc_charge","residue_ID","chain_name"
    ]
#"ring_size" = has_ring_size

MOLECULE_DEFAULT_ATTRIBUTES = {
        # 分子的基本组成
        "Atoms": "arr",  # List[Atom]， 分子中所有原子的列表
        "Bonds": "arr",  # List[Bond],  分子中所有的键
        "Angles": "arr",  # List[Angle]， 分子中所有的键角
        "Dihedrals": "arr",  # List[Dihedras], 分子中所有的二面角
        "Impropers": "arr",  # List[Improper]， 分子中所有的improper项
        "constrain": "arr",  # List[Constrain], 分子中受限制的自由度
        "Pair12": "arr",  # list[Pair], 1-2 pair
        "Pair13": "arr",  # list[Pair], 1-3 pair
        "Pair14": "arr",  # list[Pair], 1-4 pair
        "Pair1n": "arr",  # list[Pair], 1-n pair
        #"AlteredPairs": "arr",
        "element_set": "set", # 包含的元素种类

        # 原子基本属性的组合
        "elements": "arr",  # List[str], 分子中所有原子的元素符号(atom.elem)形成的一个列表
        "coordinates": "arr",  # List[[float]], 分子中所有原子的坐标(atom.coor)形成的一个列表
        "connectivity": "arr",  # List[[int]], 分子中所有原子的连接关系(atom.connect)形成的一个列表
        "bond_type": "arr",  # List[[str]], 分子中所有原子的键类型(atom.bond_type)形成的一个列表
        "formal_charge": "arr",  # List[int], 分子中所有原子的formal_charge形成的一个列表
        "charges": "arr", # {str: [float]}, 分子中所有原子的点电荷(aa.charge)形成的一个字典。charges["esp"], charges["am1bcc"],
        "esp_charge": "arr",  # List[float], 分子中所有原子的esp电荷(aa.charge["esp"])形成的一个列表
        "mulliken_charge": "arr",  # List[float], 分子中所有原子的mulliken电荷(aa.charge["mulliken"])形成的一个列表
        # 与flexible torsion相关的属性
        "scan_term": "arr",  # torsions List[],分子中所有的flexible torsion。由conformation_analysis.py中的方法得到。以下三个属性由该属性得到
        "torsion_number": "int",  # int, 分子中flexible torsion的数目
        "constrain_term": "arr",  # List[[int]], 该分子对象目前被限制的自由度，如[1,2,3,4]表示对1，2，3，4四个原子形成的二面角进行限制
        "constrain_value": "arr",  # float, 自由度被限制的数值
        # 与环相关的属性
        "ring_dict": "dict",  # rings dict, 分子中环相关的信息，从topol_structure.py中的方法得到。以下三个属性以及原子对象的环属性由该属性得到
        #                      如{"R6-1-3-6":[1,2,3,4,5,6,"ar1"],"R5-10-12-14":[10,11,12,13,14,15,"nonar"]}
        #                      其中key为某一个标记符，value[:-1]为该环所包含的原子,value[-1]表示该环的性质
        "ring_stru": "arr",  #ring_block List[List[int]]
        "ring_number": "int",  # int，分子中sssr环的数目
        "ring_size": "arr",  # List[int]， 分子中所有sssr环的大小
        "ring_property": "arr",  # List[str], 分子中所有sssr环的性质
        ""
        # 管能团相关的属性
        "function_group": "arr",  # List[str]，分子中所包含的管能团种类。由function_class.py中的方法得到
        "function_group_label": "str",  # str，分子中所包含的管能团生成的一个标记符
        # 分子的整体属性
        "net_charge": "int",  # int, 分子的净电荷
        "multi": "int",  # multiple int, 多重度
        "mass": "float",  # molecule_mass float, 分子量
        "formula": "str",  # str, 分子式
        "id": "int",  # ID int, 分子类型的编号，lammps软件会用到
        "smiles": "str",  # str, smiles字符串
        "inchi": "str",  # str, inchi字符串
        # 命名系统
        "inchi_key": "str",  # str, inchi_key
        "name": "str",  # str, 分子的名称，通常为 inchi_key
        "mole_name": "str",  # molecule_name str, alias of name
        "iupac_name": "str",  # str, iupac_name
        "nick_name": "str",  # str, nick_name
        "drug_name": "str",  # str, drug_name
        "internal_name": "str",  # str, 内部名称
        # 分子的热力学性质
        "energy": "float",  # float, 能量
        "force": "arr",  # List[float],所有原子受力的值。按原子编号和x,y,z的顺序排列
        "volicity": "arr",  # List[float], 所有原子的速度。按原子编号和x,y,z的顺序排列
        "hessian": "arr",  # List[float], hessian矩阵。按原子编号组合和xx,xy,xz,yx,yy,yz,zx,zy,zz的顺序形成的一个下三角矩阵的数组
        "freq": "arr",  # List[float]，频率的值。从小到大的排序的数组
        "dipole": "arr",  # List[float], 偶极矩，按x, y, z顺序排列
        "inertia": "arr",  # List[float], 转运惯量，按x, y, z顺序排列
        "density": "float",  # float, 密度，单位g/ml
        # QM 计算相关
        "qm_method": "str",  # str, QM计算的方法
        "qm_basis_set": "str",  # str, QM 计算的基组
    }


md_para = {
    "gromacs": {
        "integrator":"integrator",
        "timestep": "dt",
        "nsteps": "nsteps",
        "emtol": "emtol",
        "temperature": "ref_t",
        "temperature_coupl": "Tcoupl",
        "temperature_group": "tc_grps",
        "temperature_tau": "tau_t",
        "pressure": "ref_p",
        "pressure_coupl": "Pcoupl",
        "pressure_coupl_type": "pcoupltype",
        "pressure_tau": "tau_p",
        "pressure_compressibility": "compressibility",
        "rlist": "rlist",
        "cutoff_scheme": "cutoff_scheme",
        "coulomb_type": "coulombtype",
        "coul_cut": "rcoulomb",
        "vdw_cut": "rvdw",
        "vdw_long": "DispCorr",
        "velocity_generate": "gen_vel",
        "velocity_generate_temp": "gen_temp",
        "constraints": "constraints",
        "constraint-algorithm": "constraint-algorithm",
        "nstxout": "nstxout",
        "nstvout": "nstvout",
        "nstfout": "nstfout",
        "nstxout-compressed": "nstxout-compressed",
        "nstenergy": "nstenergy",
        "nstlog": "nstlog",
        "free_energy": "free_energy",
        "init_lambda_state": "init_lambda_state",
        "delta_lambda": "delta_lambda",
        "calc_lambda_neighbors": "calc_lambda_neighbors",
        "couple-moltype": "couple-moltype",
        "couple-lambda0": "couple-lambda0",
        "couple-lambda1": "couple-lambda1",
        "couple-intramol": "couple-intramol",
        "restraint-lambdas": "restraint_lambdas",
        "fep_lambdas": "fep_lambdas",
        "vdw_lambdas": "vdw_lambdas",
        "vdw_lambdas_2": "vdw_lambdas",
        "coul_lambdas": "coul_lambdas",
        "bonded_lambdas": "bonded_lambdas",
        "restraint_lambdas": "restraint_lambdas",
        "restraint_lambdas": "restraint_lambdas",
        "mass_lambdas": "mass_lambdas",
        "temperature_lambdas": "temperature_lambdas",
        "sc-alpha": "sc-alpha",
        "sc-coul": "sc-coul",
        "sc-power": "sc-power",
        "sc-sigma": "sc-sigma",
        "nstdhdl": "nstdhdl",
        "disre": "disre",
        "disre-fc": "disre-fc"
    },
}
gromacs_ftype_trans = {
    "bond": {
        "harmonic": [1, [0, 0.1], [1, 836.8]],
        "morse": [3, 10, [0,0.1],[1,4.184],[2,10.0]],
        "tabulate_2":[9,],
    },
    "angle": {
        "harmonic": [1, [0, 1.0], [1, 8.368]],
    },
    "improper": {
        "amber": [4, [0, 4.184]],
    },
    "dihedral": {
        "amber": [9, [0, 4.184], [1, 4.184], [2, 4.184], [3, 4.184]],
    },
    "pair": {
        "LJ12_6": [1, [0, 4.184], [1, 0.1]],
    },
}


#elem_table = elements.get_legacy_elements_table()