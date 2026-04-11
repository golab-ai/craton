import os,sys
import itertools
def identify_same_mole(moles):
    model_mole = {}
    count = []
    mole_type_n = 1
    former_type_name = ""
    for m in range(0, len(moles)):
        flag = 0
        for type_name in model_mole.keys():
            if (
                moles[str(m)].Elem_dict == model_mole[type_name].Elem_dict
                and moles[str(m)].Connect_dict == model_mole[type_name].Connect_dict
            ):
                if type_name == former_type_name:
                    count[-1][1] += 1
                else:
                    count.append([type_name, 1])
                    former_type_name = type_name
                flag += 1
        if flag == 0:
            model_mole["mole_type_" + str(mole_type_n)] = moles[str(m)]
            count.append(["mole_type_" + str(mole_type_n), 1])
            mole_type_n += 1
        elif flag > 1:
            sys.exit()
    return (model_mole, count)


def match_atom_number(
    m1,
    m2,
    pointp=[
        "elem",
        "formal_charge",
        "connectivity",
        "local",
    ],  # "atom_type_name"],
    bondp=["bond_type"],
    remove_H=True,
    more_connect=True,
):

    g1 = create_graphy(m1, pointp, bondp, remove_H, more_connect)
    g2 = create_graphy(m2, pointp, bondp, remove_H, more_connect)

    match_table = match_graphy(g1, g2)
    if remove_H:
        match_table = match_Hatom(m1, m2, match_table)
    return match_table


def match_Hatom(m1, m2, match_table):
    used_no = []
    for aa in m2.Atoms:
        if aa.elem in ["H", "F", "Cl", "Br"]:
            ac2 = aa.connect[0]
            ac1 = match_table[ac2]
            for ii in m1.Atoms[ac1].connect:
                if m1.Atoms[ii].elem == aa.elem:
                    if ii not in used_no:
                        used_no.append(ii)
                        match_table[aa.No] = ii
                        break
    return match_table


def get_formula(elems, position=None):
    __default_order = [
        "C",
        "H",
        "O",
        "N",
        "S",
        "P",
        "F",
        "Cl",
        "Br",
        "I",
        "B",
        "Si",
        "As",
        "Te",
        "Na",
        "K",
    ]
    __default_position_order = [
        "RFM",
        "RFMH",
        "RFMC",
        "RF",
        "RFH",
        "RFC",
        "RS",
        "RSH",
        "RSC",
        "R",
        "RC",
        "RH",
        "RHC",
        "CMR",
        "CER",
        "CM",
        "CE",
        "TA",
    ]
    formula = ""
    if position is None:
        for ee in __default_order:
            if ee in elems.keys():
                formula += ee
                if elems[ee] > 1:
                    formula += str(elems[ee])
        return formula
    else:
        for ee in __default_order:
            for pp in __default_position_order:
                this_label = ee + "$" + pp
                if this_label in position.keys():
                    formula += this_label
                    if position[this_label] > 1:
                        formula += str(position[this_label])
        return formula

def remove_hydrogen_atoms(element_dict,aa_connect_dict):
    nonh_connect_dict = {}
    terminal_atom_arr = []
    H_atom_arr = []
    for an,conn in aa_connect_dict.items():
        if element_dict[an] == "H":
            H_atom_arr.append(an)
        else:
            _tmp = []
            for aa in conn:
                if element_dict[aa] != "H":
                    _tmp.append(aa)
            if len(_tmp) < 2:
                terminal_atom_arr.append(an)
            nonh_connect_dict[an] = _tmp
    return nonh_connect_dict, terminal_atom_arr, H_atom_arr
            


def create_graphy(m, pointp, bondp, remove_H, more_connect):
    g = {"point": {}, "bond": {}}

    element_dict = {}
    aa_connect_dict = {}
    for aa in m.Atoms:
        element_dict[aa.No] = aa.elem
        aa_connect_dict[aa.No] = aa.connect
    if remove_H:
        nonh_connect_dict, terminal_atom_arr, H_atom_arr = remove_hydrogen_atoms(element_dict, aa_connect_dict)
        connect_dict = {}
        for aa, bb in nonh_connect_dict.items():
            tmp = []
            for iii in bb:
                if iii not in H_atom_arr:
                    tmp.append(iii)
            connect_dict[aa] = tmp
    else:
        H_atom_arr = []
        connect_dict = aa_connect_dict

    bond_dict = {}
    for aa, bb in connect_dict.items():
        bond_dict[aa] = [
            bb,
        ]
        for pp in bondp[1:]:
            tmp = []
            for ii in bb:
                tmp.append(getattr(m.Atoms[aa], pp)[m.Atoms[aa].connect.index(ii)])
            bond_dict[aa].append(tmp)
    g["bond"] = bond_dict

    for aa in m.Atoms:
        if aa.No not in H_atom_arr:
            label = ""
            for pp in pointp:
                if pp == "connectivity":
                    this_elems = {}
                    if more_connect:
                        this_position = {}
                    for ii in aa.connect:
                        if m.Atoms[ii].elem not in this_elems.keys():
                            this_elems[m.Atoms[ii].elem] = 0
                        this_elems[m.Atoms[ii].elem] += 1
                        if more_connect:
                            this_label = m.Atoms[ii].elem + "$" + m.Atoms[ii].local
                            if this_label not in this_position.keys():
                                this_position[this_label] = 0
                            this_position[this_label] += 1
                    if more_connect:
                        sss = get_formula(this_elems, position=this_position)
                    else:
                        sss = get_formula(this_elems)
                    label += f"{sss}:"
                else:
                    label += "%s:" % getattr(aa, pp)
            label = label[:-1]
            if label not in g["point"].keys():
                g["point"][label] = []
            g["point"][label].append(aa.No)
    return g


def match_graphy(g1, g2):
    import time
    print(time.time())
    point1, bond1 = g1["point"], g1["bond"]
    point2, bond2 = g2["point"], g2["bond"]
    point1_combine,point2_combine = {},{}
    point1_arr, point2_arr = [],[]
    point1_order,point2_order = [],[]
    point2_combine, point2_arr, point1_order = {}, [], []
    for aa,b1 in point1.items():
        if aa in point2:
            b2 = point2[aa]
            nn = min([len(b1),len(b2)])
            point1_combine[aa] = list(itertools.combinations(b1,nn))
            point2_combine[aa] = list(itertools.permutations(b2,nn))

 
    for aa, bb in point1_combine.items():
        point1_arr.append(point1_combine[aa])
        point2_arr.append(point2_combine[aa])

    #print(point1_arr)
    #print("##############")
    #print(point2_arr)
    for aa, bb in bond1.items():
        bnn = len(bb)
        break

    larger_match_table = []
    nn = 0
    for arr1 in itertools.product(*point1_arr):
        print("####arr1:",arr1)
        tmp1 = []
        for rr in arr1:
            tmp1.extend(rr)
        print("@@@tmp1:",tmp1)


        for arr2 in itertools.product(*point2_arr):
            print("   ###################")
            print("   @@@@arr2:",arr2)
            tmp2 = []
            for rr in arr2:
                tmp2.extend(rr)
            print("   @@@tmp2:",tmp2)
            
            tmp_match_table, tmp_bond, match_table = {}, {}, {}
            
            for ii,an in enumerate(tmp2):
                tmp_match_table[an] = tmp1[ii]
            print("   tmp_match_table:",tmp_match_table)
            tmp_match_table1 = {vv:kk for kk,vv in tmp_match_table.items()}

            for kk, vv in tmp_match_table.items():
                connect1,connect2 = [], []
                for jj in bond2[kk][0]:
                    _an = tmp_match_table[jj] if jj in tmp_match_table else None
                    if _an is not None:
                        connect2.append(_an)
                for jj in bond1[vv][0]:
                    if jj in tmp_match_table1:
                        connect1.append(jj)
                print("      kk,vv,connect1,connect2:",kk,vv,connect1,connect2)
                #if set(connect1) == set(connect2):
                if set(connect1).issubset(set(connect2)) or set(connect2).issubset(set(connect1)):
                    match_table[kk] = vv
            print("   match_table:",match_table)
            if len(match_table) == nn:
                larger_match_table.append(match_table)
            elif len(match_table) > nn:
                larger_match_table = []
                larger_match_table.append(match_table)
                nn = len(match_table)
    _tmp_ = [[kk,vv] for kk ,vv in larger_match_table[0].items()]
    _tmp_ = sorted(_tmp_,key=lambda x:x[0])
    return larger_match_table

def match_graphy_all(g1, g2):
    point1, bond1 = g1["point"], g1["bond"]
    point2, bond2 = g2["point"], g2["bond"]
    
    point2_combine, point2_arr, point1_order = {}, [], []
    for aa, bb in point2.items():
        point2_combine[aa] = list(itertools.permutations(bb))

    # set1 = set(point1.keys()) - set(point2.keys())
    # set2 = set(point2.keys()) - set(point1.keys())

    for aa, bb in point1.items():
        point1_order.extend(bb)
        point2_arr.append(point2_combine[aa])

    for arr in itertools.product(*point2_arr):
        tmp = []
        for rr in arr:
            tmp.extend(rr)
        flag = True
        match_table, tmp_bond = {}, {}
        for aa, bb in bond1.items():
            bnn = len(bb)
            break

        for i in range(len(tmp)):
            match_table[tmp[i]] = point1_order[i]
        for kk, vv in match_table.items():
            tmp_bond[vv] = [
                [],
            ]
            for jj in bond2[kk][0]:
                tmp_bond[vv][0].append(match_table[jj])
            for iii in range(1, bnn):
                tmp_bond[vv].append(bond2[kk][iii])
        for aaa, bbb in tmp_bond.items():
            if set(bond1[aaa][0]) != set(bbb[0]):
                flag = False
                break
        if flag:
            break
    if flag:
        return match_table
    return None