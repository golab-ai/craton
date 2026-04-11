from copy import deepcopy
from collections import OrderedDict
import itertools
from ...utils.geometry import calculate_distance,calculate_angle,calculate_dihedral,calc_stru_para
from .ring import find_cyclo
from .ring import cyclo_blocks
from .chain import search_chain
from .structure import find_mole
from ...chem.elements import get_bonded_distance, Element,get_bonded_type_distance

def determine_bonded(e1, e2, a, b):  # 两个原子间是否形成键
    d = calculate_distance(a, b)
    if d > 3.0:
        return False
    else:
        d_r = get_bonded_distance(e1, e2)
        #文献上给定的是d_r + 0.4, rkdit生成的C1=CC=C1,1-3原子间的距离会稍微小于这一值
        #if d < d_r + 0.4 and d > 0.8:
        if d < d_r + 0.35 and d > 0.8:
            return True
        else:
            return False

def connected_atoms(elements,coordinates):
    an = len(elements)
    connectivitys = {ii:[] for ii in elements}
    for ii in range(an):
        for jj in range(ii + 1, an):
            if determine_bonded(elements[ii],elements[jj],coordinates[ii],coordinates[jj]):
                connectivitys[ii].append(jj)
                connectivitys[jj].append(ii)
    return connectivitys

def over_connected_check(elements,coordinates,connectivitys):
    for ii,elem in elements.items():
        if elem in ["C","N","S","P"] and len(connectivitys[ii]) > 4:
            tmp = [[jj, calculate_distance(coordinates[ii],coordinates[jj])] for jj in connectivitys[ii]]
            tmp = sorted(tmp,key=lambda x:x[1])
            connectivitys[ii] = [tmp[nn][0] for nn in range(4)]
            for rr in tmp[4:]:
                index = connectivitys[rr[0]].index(ii)
                del connectivitys[rr[0]][index]
    return connectivitys

def single_bond_case(elem,conn,three_N_flag = False):
    flag = False
    if elem in ["H","F","Cl","Br","I"] and conn == 1:
        flag = True
    elif elem in ["O","S"] and conn == 2:
        flag = True
    elif elem in ["C","N"] and conn == 4:
        flag = True
    elif elem in ["N","P"] and conn == 3:
        flag = three_N_flag
    return flag

def create_initial_bond_types(elements,connectivitys):
    used = []
    bond_types = {}
    for ii,conn in connectivitys.items():
        for jj in conn:
            if [ii,jj] not in used and [jj,ii] not in used:
                used.append([ii,jj])
                bond_types[f"{ii}-{jj}"] = "N"
    return bond_types
    
def classify_based_on_element_connects(elements,connectivitys,ring_atoms):
    _types = {}
    for ii,elem in elements.items():
        ss = f"{elem}{len(connectivitys[ii])}"
        if ss in ["N3","P3"]:
            #if ii not in ring_atoms:
            for rr in connectivitys[ii]:
                if elements[rr] in ["O","N"] and len(connectivitys[rr]) == 1:
                    ss += "A"
                    break
        if ss not in _types:
            _types[ss] = []
        _types[ss].append(ii)
    return _types

def get_bond_name(ii,jj,all_bonds):
    return list(set([f"{ii}-{jj}",f"{jj}-{ii}"]).intersection(set(all_bonds)))[0]

def scale_length_bond(ii,connects,elements,coordinates,bt="2"):
    tmp = {}
    for jj in connects:
        dist = calculate_distance(coordinates[ii],coordinates[jj])
        double_bond_length = get_bonded_type_distance(elements[ii], elements[jj], bt)
        if double_bond_length is None:
            tmp[jj] = 0
        else:
            tmp[jj] = dist / double_bond_length
    return tmp

def N3A(ii,elements,coordinates,connectivitys,bond_types,all_type_bonds,conju_length_bond_types):
    def _carboxy_():
        for jj in conn:
            if elements[jj] == "C":
                for kk in connectivitys[jj]:
                    if elements[kk] == "O" and conju_length_bond_types[get_bond_name(ii,jj,all_type_bonds)] >= 1.5:
                        return True
        return False

    conn = connectivitys[ii]
    length_bond_type = [conju_length_bond_types[get_bond_name(ii,jj,all_type_bonds)] for jj in conn]
    elems = [elements[jj] for jj in conn]
    if set(length_bond_type) == set([1]):
        for jj in conn:
            bond_types[get_bond_name(ii,jj,all_type_bonds)] = 1
        return bond_types
    
    if elems.count("C") >=2:
        if _carboxy_():
            for jj in conn:
                bond_types[get_bond_name(ii,jj,all_type_bonds)] = 1
            return bond_types

    sa = -1
    
    for jj in conn:
        if single_bond_case(elements[jj],len(connectivitys[jj])):
            bond_types[get_bond_name(ii,jj,all_type_bonds)] = 1
            sa = jj
            break
    if sa == -1:
        for jj in conn:
            if elements[jj] in ["C","H","P","S"] and conju_length_bond_types[get_bond_name(ii,jj,all_type_bonds)] == 1:
                bond_types[get_bond_name(ii,jj,all_type_bonds)] = 1
                sa = jj
                break
    if sa == -1:
        tmp0 = scale_length_bond(ii,connectivitys[ii],elements,coordinates)
        tmp = sorted([[aa,bb] for aa,bb in tmp0.items()],key=lambda x:x[1])
        bond_types[get_bond_name(ii,tmp[-1][0],all_type_bonds)] = 1
        sa = tmp[-1][0]
    for jj in conn:
        if jj != sa:
            bond_types[get_bond_name(ii,jj,all_type_bonds)] = 2
    return bond_types    

def SPClBr34(ii,elements,coordinates,connectivitys,bond_types,all_bonds):
    nn = len(connectivitys[ii])
    if elements[ii] == "S" and len(connectivitys[ii]) == 3:
        double_bond_number = 1
    elif elements[ii] in ["S","N","Se"]:
        double_bond_number = 2
    elif elements[ii] in ["P","Cl","Br","As"]:
        double_bond_number = 1
    
    tmp0 = scale_length_bond(ii,connectivitys[ii],elements,coordinates)
    tmp = {}
    for aa,bb in tmp0.items():
        if bb == 0:
            tmp[aa] = [1,0]
        else:
            if single_bond_case(elements[aa],len(connectivitys[aa]),three_N_flag=False):
                tmp[aa] = [1,bb]
            else:
                tmp[aa] = [0,bb]

    tn = nn - sum([tt[0] for jj,tt in tmp.items()])

    if tn > double_bond_number:
        tmp_tmp = [[jj,tt[1]] for jj,tt in tmp.items() if tt[0] == 0]
        tmp_tmp = sorted(tmp_tmp,key=lambda x:x[1])
        for kk,jj in enumerate([rr[0] for rr in tmp_tmp]):
            if kk < double_bond_number:
                tmp[jj][0] = 2
            else:
                tmp[jj][0] = 1
    else:
        for jj,tt in tmp.items():
            if tt[0] == 0:
                tmp[jj][0] = 2
        if tn < double_bond_number:
            tmp_tmp = [[jj,tt[1]] for jj,tt in tmp.items() if tt[0] == 1]
            tmp_tmp = sorted(tmp_tmp,key=lambda x:x[1],reverse=True)
            for kk,jj in enumerate([rr[0] for rr in tmp_tmp]):
                if kk >= double_bond_number:
                    tmp[jj][0] = 2

    for jj,tt in tmp.items():
        bond_types[get_bond_name(ii,jj,all_bonds)] = tt[0]
    return bond_types

def get_linear_plane_in_ring(ii,coordinates,connectivitys,ring_atoms):
    coors = []
    connects = connectivitys[ii]
    for jj in connects:
        if jj in ring_atoms:
            pass

def get_linear_plane_single(ii,
                            elements,
                            coordinates,
                            connectivitys,
                            ring_atoms,
                            small_ring_atoms,
                            length_bond_types,
                            conju_length_bond_types,
                            all_bond_types):
    
    def _judge_SP(this_tha,ii,jj,elements,conju_length_bond_types,all_bond_types):
        if abs(this_tha) >= 165 or abs(this_tha) <= 15: #原始值 175，5,曾改为170，10
            if elements[ii] in ["C"] or elements[jj] in ["C"]:
                if conju_length_bond_types[get_bond_name(ii,jj,all_bond_types)] != 1:
                    return "SP"
                else:
                    return "NP"
            else:
                return "SP"
        return "NP"
            
    connects = connectivitys[ii]
    coors = [coordinates[connects[0]],coordinates[ii]]
    for rr in connects[1:]:
        coors.append(coordinates[rr])

    if len(coors) > 4:
        return "M"
    tha = calc_stru_para(coors)
    if len(coors) <= 2:
        return "UN"
    else:
        if len(coors) == 4:
            if abs(tha) >= 160: ### 原始值175,后改为168
                return "P"
            else:
                if ii in small_ring_atoms:
                    return "UP"
                if ii in ring_atoms:
                    conns = [jj for jj in connectivitys[ii]]
                    #if len(set(conns).difference(set(ring_atoms))) > 0:
                    #    return "NP"
                    #else:
                    sp_arrs = []
                    for jj in conns:
                        for kk in connectivitys[jj]:
                            if kk != ii:
                                for ll in conns:
                                    if ll != jj and ll != kk:
                                        if jj in ring_atoms and kk in ring_atoms and ll in ring_atoms:
                                            this_tha = calc_stru_para([coordinates[ll],coordinates[ii],coordinates[jj],coordinates[kk]])
                                            if _judge_SP(this_tha,ii,jj,elements,conju_length_bond_types,all_bond_types) == "SP":
                                                return "SP"
                    return "NP"
                return "NP"
        else:
            if abs(tha) >= 175: ### 原始值175
                return  "L"
            else:
                if len(connectivitys[connects[-1]]) == 1 and len(connectivitys[connects[0]]) == 1:
                    return "JT"
                if ii in small_ring_atoms:
                    return "UP"
                if ii in ring_atoms:
                    for jj in connectivitys[connects[-1]]:
                        if jj != ii and jj in ring_atoms and jj != connects[0]:
                            this_coors = deepcopy(coors)
                            this_coors.append(coordinates[jj])
                            this_tha = calc_stru_para(this_coors)
                            tmp =  _judge_SP(this_tha,ii,connects[-1],elements,conju_length_bond_types,all_bond_types)
                            if tmp == "SP":
                                return "SP"
                    for jj in connectivitys[connects[0]]:
                        if jj != ii and jj in ring_atoms and jj != connects[-1]:
                            this_coors = deepcopy(coors)
                            this_coors = [coordinates[jj]] + this_coors
                            this_tha = calc_stru_para(this_coors)
                            tmp =  _judge_SP(this_tha,ii,connects[0],elements,conju_length_bond_types,all_bond_types)
                            if tmp == "SP":
                                return "SP"
                    return "NP"
                else:
                    return "UP"

def assign_single_bond(atoms,connectivitys,bond_types,all_bond_types):
    for ii in atoms:
        for jj in connectivitys[ii]:
            bn = get_bond_name(ii,jj,all_bond_types)
            if bond_types[bn] == "N":
                bond_types[bn] = 1
    return bond_types

def assign_acid_bond_type(atoms,elements,coordinates,connectivitys,bond_types,all_bonds,ss,conju_length_bond_types):
    if ss == "N3A":
        for ii in atoms:
            bond_types = N3A(ii,elements,coordinates,connectivitys,bond_types,all_bonds,conju_length_bond_types)
        return bond_types
    
    for ii in atoms:
        bond_types = SPClBr34(ii,elements,coordinates,connectivitys,bond_types,all_bonds)
    return bond_types

def get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types):
    unassign = []
    for ii,conn in connectivitys.items():
        tmp = [bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in conn]
        if "N" in tmp:
            unassign.append(ii)
    return unassign

def assign_plane_line(style_atoms,connectivitys,bond_types,all_bond_types):
    for bn,bt in bond_types.items():
        if bt == "N":
            
            ii,jj = [int(an) for an in bn.split("-")]
            
            if style_atoms[ii] == "P" and style_atoms[jj] == "L":
                bond_types[bn] = 2
                for kk in connectivitys[ii]:
                    if kk != jj:
                        bond_types[get_bond_name(ii,kk,all_bond_types)] = 1
                #bond_types[get_bond_name(jj,connectivitys[jj][0],all_bond_types)] = 2
            elif style_atoms[jj] == "P" and style_atoms[ii] == "L":
                bond_types[bn] = 2
                for kk in connectivitys[jj]:
                    if kk != ii:
                        bond_types[get_bond_name(jj,kk,all_bond_types)] = 1
                #bond_types[get_bond_name(ii,connectivitys[ii][0],all_bond_types)] = 2
    return bond_types
        
def assign_linear_bond(elements,connectivitys,bond_types,style_atoms,all_bond_types,formal_charges):
    linear_atoms = [ii for ii,st in style_atoms.items() if st == "L"]
    for ii in linear_atoms:
        conn = connectivitys[ii]
        elems = [elements[jj] for jj in conn]
        sts = [style_atoms[jj] for jj in conn]
        if elements[ii] == "N":
            if "UN" in sts:
                if elems[sts.index("UN")] in ["N","C"]:
                    bond_types[get_bond_name(ii,conn[sts.index("UN")],all_bond_types)] = 3
                    formal_charges[conn[sts.index("UN")]] = 3 - Element.get(elements[conn[sts.index("UN")]]).valents[0]
                formal_charges[ii] = 1
        else:
            conn = connectivitys[ii]
            elems = [elements[jj] for jj in conn]
            #if len(set(["O","S","Se"]).intersection(set(elems))) > 0:
            #    for jj in conn:
            #        bond_types[get_bond_name(ii,jj,all_bond_types)] = 2
            #else:
            if "UN" in sts:
                if elems[sts.index("UN")] == "N":
                    bond_types[get_bond_name(ii,conn[sts.index("UN")],all_bond_types)] = 3
    return bond_types,formal_charges

def assign_three_n(elements,connectivitys,bond_types,style_atoms,all_bond_types,formal_charges):
    for ii,elem in elements.items():
        if len(connectivitys[ii]) == 2:
            if elem == "N":
                conn = connectivitys[ii]
                elems = [elements[jj] for jj in conn]
                sts = [style_atoms[jj] for jj in conn]
                if elems == ["N","N"] and "UN" in sts:
                    for jj in conn:
                        bond_types[get_bond_name(ii,jj,all_bond_types)] = 2

    return bond_types,formal_charges

def assign_bond_for_plane_linear(elements,connectivitys,bond_types,style_atoms,all_bond_types):
    unassign_atoms = get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types)

    while 1:
        for ii in unassign_atoms:
            flag = False
            conn = connectivitys[ii]
            bts = [bond_types[get_bond_name(ii,jj,all_bond_types)] if bond_types[get_bond_name(ii,jj,all_bond_types)] != "N" else 0 for jj in conn]
            if style_atoms[ii] == "P":
                if bts.count(1) == 2:
                    flag = True
                    nn = 2
                elif bts.count(2) > 0:
                    flag = True
                    nn = 1
            if style_atoms[ii] == "L":
                if bts.count(0) == 1:
                    nn = 4 - sum(bts)
                    flag = True
            if flag:
                for kk,bt in enumerate(bts):
                    if bt == 0:
                        bond_types[get_bond_name(ii,conn[kk],all_bond_types)] = nn
        tmp = unassign_atoms = get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types)
        if len(set(tmp).difference(unassign_atoms)) == 0:
            break
        else:
            unassign_atoms = tmp
    return bond_types

def judge_plane_C_without_double_bond(atoms,connectivitys,bond_types_tmp,style_atoms,all_bond_types):
    nn = 0
    for ii in atoms:
        if style_atoms[ii] in ["P","SP","L"]:
            this_bond_types = [bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]]
            if "N" in this_bond_types:
                nn += 1
            else:
                kk = this_bond_types.count(2)
                if kk != 1:
                    nn += 1
    return nn

def get_unassign_bond_types_and_atoms(ring,connectivitys,bond_types,all_bond_types):
    unassign_atoms = list(set(ring).intersection(set(get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types))))
    ex_ring = [ii for rr in ring for ii in connectivitys[rr] if ii not in ring]
    total_ring = ring + ex_ring
    bond_types_ring = {}
    connectivitys_ring = {}
    for ii,conn in connectivitys:
        connectivitys_ring[ii] = []
        for jj in conn:
            if jj in total_ring:
                connectivitys_ring[ii].append(jj)
    for bn,bt in bond_types.items():
        ii,jj = [int(aa) for aa in bn.split("-")]
        if ii in total_ring and jj in total_ring:
            bond_types_ring[bn] = bt
    unassign_bond_types = [bn for bn,bt in bond_types_ring.items() if bt == "N"]
    return unassign_atoms, unassign_bond_types,connectivitys_ring,bond_types_ring

def sing_ring_bond_type(accu_ring_atoms,connectivitys,bond_types,style_atoms,all_bond_types):
    #total_ring_atoms = accu_ring_atoms + ring
    unassign_atoms = list(set(accu_ring_atoms).intersection(set(get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types))))
    
    #unassign_atoms, unassign_bond_types,connectivitys_ring,bond_types_ring = get_unassign_bond_types_and_atoms(ring,connectivitys,bond_types,all_bond_types)

    if len(unassign_atoms) == 0:
        nn = judge_plane_C_without_double_bond(accu_ring_atoms,connectivitys,bond_types,style_atoms,all_bond_types)
        if nn == 0:
            return [bond_types]
        else:
            return []
        
    unassign_bond_types = []
    for ii in unassign_atoms:
        for jj in connectivitys[ii]:
            bn = get_bond_name(ii,jj,all_bond_types)
            if bn not in unassign_bond_types:
                if bond_types[bn] == "N":
                    unassign_bond_types.append(bn)

    tmp = [[[bn,1],[bn,2]] for bn in unassign_bond_types]
    correct_combine = []
    for combine in itertools.product(*tmp):
        this_bond_types_ring = deepcopy(bond_types)
        for rr in combine:
            this_bond_types_ring[rr[0]] = rr[1]
        nn = judge_plane_C_without_double_bond(accu_ring_atoms,connectivitys,this_bond_types_ring,style_atoms,all_bond_types)
        if nn == 0:
            correct_combine.append(this_bond_types_ring)
    return correct_combine

def get_combinations(arrs,reverse=False):
    if reverse:
        tmp_tmp = []
        for nn in range(len(arrs),0,-1):
            for rr in itertools.combinations(arrs,nn):
                tmp_tmp.append(rr)
        tmp_tmp.append([])
    else:
        tmp_tmp = [[]]
        for nn in range(len(arrs)):
            for rr in itertools.combinations(arrs,nn+1):
                tmp_tmp.append(rr)
    return tmp_tmp

def get_ring_special_atoms(unassign_atoms,
                           ring_atoms,
                           ring_connect_atoms,
                           style_atoms,
                           length_bond_types,
                           conju_length_bond_types,
                           elements,
                           coordinates,
                           connectivitys,
                           all_bond_types):
        ##ring的排列
    #block_ring_arrs = []
    #for rr in itertools.permutations(block_ring,len(block_ring)):
    #    block_ring_arrs.append(list(rr))
    ##环内平面的C原子，但不存在双键
    tmp1 = []
    tmp2 = []
    for ii in unassign_atoms:
        if elements[ii] == "C":
            if style_atoms[ii] == "SP":
                elems = [elements[jj] for jj in connectivitys[ii] if jj in ring_atoms]
                if set(elems) == set(["C"]):
                    if 2 not in [conju_length_bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]]:
                        tmp1.append(ii)
                elif "N" in elems:
                    if 2 not in [length_bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]]:
                        tmp1.append(ii)
    ps_plane_C = get_combinations(tmp1,reverse=True)
    ##环内连接性为3的N原子，变成N正离子，可以形成双键
    tmp = [ii for ii in ring_atoms if elements[ii] == "N" and len(connectivitys[ii]) == 3]
    ion_N =get_combinations(tmp)
    
    ##环外可能形成双键的O，S，N，C原子：
    tmp = []
    for ii in ring_connect_atoms:
        if elements[ii] in ["O","S","N","C"] and len(connectivitys[ii]) < Element.get(elements[ii]).valents[0]:
            for jj in connectivitys[ii]:
                if jj in unassign_atoms:
                    nn = length_bond_types[get_bond_name(ii,jj,all_bond_types)]
                    if nn == 2:
                        tmp.append(jj)
    ex_ring_double = get_combinations(tmp,reverse=True)

    #得到不存在双键的连接性为2的N原子：
    tmp = []
    for ii in ring_atoms:
        if ii in unassign_atoms:
            if elements[ii] == "N" and len(connectivitys[ii]) == 2:
                tmp_tmp = [conju_length_bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]]
                if tmp_tmp == [1,1]:
                    tmp.append(ii)
    confused_N = get_combinations(tmp,reverse=True)

    #得到可能破坏共轭结构的同一平面的N原子的组合
    tmp = [ii for ii in unassign_atoms if elements[ii] == "N" and len(connectivitys[ii]) == 2]
    for rr in get_combinations(tmp):
        if rr not in confused_N:
            confused_N.append(rr)
    p_atoms = [ii for ii in ring_atoms if style_atoms[ii] == "P"]
    return ps_plane_C,ion_N, ex_ring_double,confused_N,p_atoms

def assign_ring_special_atoms(nc,
                              ion,
                              ex_ring,
                              con,
                              ring_atoms,
                              unassign_atoms_tmp,
                              connectivitys,
                              bond_types_tmp,
                              p_atoms,
                              all_bond_types):
    #处理环中平面的CH2基团
    for ii in nc:
        for jj in connectivitys[ii]:
            bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] = 1

        index_n = unassign_atoms_tmp.index(ii)
        del unassign_atoms_tmp[index_n]
    #处理N正离子
    for ii in ion:
        for jj in connectivitys[ii]:
            #if jj in ring_atoms:
            bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] = "N"
            if jj in p_atoms:
                for kk in connectivitys[jj]:
                    if kk != ii and kk in ring_atoms:
                        bond_types_tmp[get_bond_name(jj,kk,all_bond_types)] = "N"
        unassign_atoms_tmp.append(ii)
    #处理环外双键
    if ex_ring is None:
        pass
    else:
        for ii in ex_ring:
            for jj in connectivitys[ii]:
                if jj in ring_atoms:
                    bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] = 1
                else:
                    bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] = 2
            index_n = unassign_atoms_tmp.index(ii)
            del unassign_atoms_tmp[index_n]
        for ii in unassign_atoms_tmp:
            if ii not in ion:
                for jj in connectivitys[ii]:
                    if jj not in ring_atoms:
                        bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] = 1

    ## 将这种N原子的键设成单键
    for ii in con:
        for jj in connectivitys[ii]:
            bond_types_tmp[get_bond_name(ii,jj,all_bond_types)] = 1
        index_n = unassign_atoms_tmp.index(ii)
        del unassign_atoms_tmp[index_n]
    return bond_types_tmp,unassign_atoms_tmp

def assign_ring_block_double_bond(block_ring,
                                  elements,
                                  coordinates,
                                  connectivitys,
                                  bond_types,
                                  style_atoms,
                                  all_bond_types,
                                  length_bond_types,
                                  conju_length_bond_types,
                                  find_all=False):
    def _back(ri,ring_n,bond_types_dict,correct_count):
        flag = False
        for aa in range(ring_n + 1):

            if ri == 0:
                flag = True
                break
            if correct_count[ri - 1] == len(bond_types_dict[ri - 1]) - 1:
                del bond_types_dict[ri-1]
                del correct_count[ri-1]
            else:
                correct_count[ri-1] += 1
                break
            ri -= 1
        return flag,ri,bond_types_dict,correct_count

    plane_label = ["P","SP","L"]
    ring_atoms = [ii for ring in block_ring for ii in ring]
    ring_connect_atoms = []
    for ii in ring_atoms:
        for jj in connectivitys[ii]:
            if jj not in ring_atoms:
                ring_connect_atoms.append(jj)
    unassign_atoms = list(set(ring_atoms).intersection(set(get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types))))
    unassign_atoms = [ii for ii in unassign_atoms if style_atoms[ii] not in ["UP"]]
    ps_plane_C,ion_N,ex_ring_double,confused_N,p_atoms = get_ring_special_atoms(unassign_atoms,
                                                             ring_atoms,
                                                             ring_connect_atoms,
                                                             style_atoms,
                                                             length_bond_types,
                                                             conju_length_bond_types,
                                                             elements,
                                                             coordinates,
                                                             connectivitys,
                                                             all_bond_types)
    
    ring_n = len(block_ring)
    total_correct_bond_types = []
    for nc_ii,nc in enumerate(ps_plane_C):
        for ex_ring_ii,ex_ring in enumerate(ex_ring_double):
            for ion_ii,ion in enumerate(ion_N):
                for con_ii,con in enumerate(confused_N):
                    #循环迭代每种assign_N_single_bond，找出正确或最接近正确的结果
                    bond_types_tmp = deepcopy(bond_types)
                    unassign_atoms_tmp = deepcopy(unassign_atoms)
                    bond_types_tmp,unassign_atoms_tmp = assign_ring_special_atoms(nc,
                                                                                  ion,
                                                                                  ex_ring,
                                                                                  con,
                                                                                  ring_atoms,
                                                                                  unassign_atoms_tmp,
                                                                                  connectivitys,
                                                                                  bond_types_tmp,
                                                                                  p_atoms,
                                                                                  all_bond_types)

                    block_ring_unassign_atoms = [list(set(ring).intersection(set(unassign_atoms_tmp))) for ring in block_ring]

                    bond_types_dict = {0:[]}
                    correct_count = {0:0}
                    accu_ring_atoms = [itertools.chain(*block_ring_unassign_atoms[:ii+1]) for ii,rr in enumerate(block_ring_unassign_atoms)]
                    accu_ring_atoms = [list(set(arr)) for arr in accu_ring_atoms]
                    
                    this_bond_types_tmp = deepcopy(bond_types_tmp)
                    tmp = sing_ring_bond_type(accu_ring_atoms[0],connectivitys,this_bond_types_tmp,style_atoms,all_bond_types)
                    if len(tmp) > 0:
                        for rr in tmp:
                            bond_types_dict[0].append(rr)
                    else:
                        continue
                    ri = 1
                    flag = False
                    while 1:
                        if ri == ring_n:
                            if not find_all:
                                return [bond_types_dict[ri-1][0]]
                            else:
                                total_correct_bond_types.extend(bond_types_dict[ri-1])
                                ri -= 1
                                flag,ri,bond_types_dict,correct_count = _back(ri,ring_n,bond_types_dict,correct_count)
                                if flag:
                                    break
                        this_bond_types_tmp = deepcopy(bond_types_dict[ri-1][correct_count[ri-1]])
                        tmp = sing_ring_bond_type(accu_ring_atoms[ri],connectivitys,this_bond_types_tmp,style_atoms,all_bond_types)
                        if len(tmp) > 0:
                            correct_count[ri] = 0
                            bond_types_dict[ri] = []
                            for rr in tmp:
                                bond_types_dict[ri].append(rr)
                            ri += 1
                        else:
                            flag,ri,bond_types_dict,correct_count = _back(ri,ring_n,bond_types_dict,correct_count)
                        if flag:
                            break
                    
                #if len(total_correct_bond_types) > 0:
                #    return total_correct_bond_types
            #if len(total_correct_bond_types) > 0:
            #    return total_correct_bond_types
        #if len(total_correct_bond_types) > 0:
        #    return total_correct_bond_types
    if len(total_correct_bond_types) > 0:
        return total_correct_bond_types
    return [bond_types]

def assign_ring_double_bond(rings,
                            ring_block_components,
                            elements,coordinates,
                            connectivitys,
                            bond_types,
                            style_atoms,
                            all_bond_types,
                            length_bond_types,
                            conju_length_bond_types,
                            find_all=False
                            ):
    block_rings = []
    for rr in ring_block_components:
        block_rings.append([])
        for r0 in rr:
            block_rings[-1].append(rings[r0][:-1])
        #block_rings[-1] = sorted(block_rings[-1], key=lambda x:len(x))
    if len(block_rings) == 0:
        return [bond_types]
    ring_bond_types_dict = {
                            -1: [deepcopy(bond_types)],
                            }
    for ii,block_ring in enumerate(block_rings):
        ring_bond_types_dict[ii] = []
        for bond_types in ring_bond_types_dict[ii-1]:
            bond_types_arr = assign_ring_block_double_bond(block_ring,
                                                            elements,
                                                            coordinates,
                                                            connectivitys,
                                                            bond_types,
                                                            style_atoms,
                                                            all_bond_types,
                                                            length_bond_types,
                                                            conju_length_bond_types,
                                                            find_all=find_all
                                                            )
            ring_bond_types_dict[ii].extend(bond_types_arr)
    return ring_bond_types_dict[len(block_rings) - 1]

def assign_single_ring_block_double_bond(atoms,elements,coordinates,connectivitys,bond_types_tmp,style_atoms):
    plane_label = ["P","SP","L"]
    all_bond_types = list(bond_types_tmp.keys())
    this_unassign_atoms = list(set(atoms).intersection(set(get_unassign_bond_type_atoms(connectivitys,bond_types_tmp,all_bond_types))))
    if len(this_unassign_atoms) == len(atoms):
        block_arrs = [atoms]
    else:
        block_arrs = get_block_chain(this_unassign_atoms,connectivitys)

    for arr in block_arrs:
        tmp_arrs = []
        unbond_atoms = [jj for jj in connectivitys[arr[0]] if bond_types_tmp[(get_bond_name(arr[0],jj,all_bond_types))] == "N"]
        for jj in unbond_atoms:
            this_bond_types_tmp = deepcopy(bond_types_tmp)
            this_bond_types_tmp[get_bond_name(arr[0],jj,all_bond_types)] = 2
            this_bond_types_tmp = get_ring_conjugate_bond_type(block_arrs,elements,coordinates,connectivitys,this_bond_types_tmp,style_atoms,plane_label,all_bond_types,this_unassign_atoms) 

    return bond_types_tmp

def get_bond_type_based_on_length(e1,e2,coor1,coor2,conju_flag=False):
    dist = calculate_distance(coor1,coor2)
    s = get_bonded_type_distance(e1, e2, "1") # single_bond_length
    
    d = get_bonded_type_distance(e1,e2,"2") # double_bond_length
    t = get_bonded_type_distance(e1,e2,"3") # triple_bond_length
    if conju_flag:
        conju = get_bonded_type_distance(e1, e2, "1.5")
        if conju is not None:
            if dist > (s+conju) * 0.5:
                return 1
            else:
                if d is not None:
                    #if dist > (conju+d)*0.5:
                    if dist > conju - 0.01:
                        return 1.5
                    else:
                        return 2
                return 1.5
        else:
            if d is not None:
                if dist > (s+d)*0.5:
                    return 1
                else:
                    return 2
    
    if d is not None:
        if dist > (s+d)*0.5:
            return 1
        else:
            if t is not None:
                if dist > (d+t)*0.5:
                    return 2
                else:
                    return 3
            return 2
    else:
        return 1
        
def judge_style_of_unassign_atoms(elements,
                                  coordinates,
                                  connectivitys,
                                  bond_types,
                                  style_atoms,
                                  all_bond_types,
                                  length_bond_types,
                                  conju_length_bond_types):
    unassign_atoms = get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types)
    for ii in unassign_atoms:
        if style_atoms[ii] in ["UP","JT","L"]:
            conn = [jj for jj in connectivitys[ii] if jj in unassign_atoms]
            elems = [elements[jj] for jj in conn]
            elems += [elements[ii]]
            if elems == ["C","C","C"]:
                tmp = [conju_length_bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in conn]
            else:
                tmp = [length_bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in conn]
            #tmp = [get_bond_type_based_on_length(elements[ii],elements[jj],coordinates[ii],coordinates[jj]) for jj in conn]
            nn = sum([1 if n >= 2 else 0 for n in tmp])
            if nn >= 2:
                tmp_tmp = [conju_length_bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in conn]
                #tmp_tmp = [get_bond_type_based_on_length(elements[ii],elements[jj],coordinates[ii],coordinates[jj],conju_flag=True) for jj in conn]
                nnn = sum([1 if n == 2 else 0 for n in tmp_tmp])
                if nnn >= 2:
                    style_atoms[ii] = "CU"
                else:
                    style_atoms[ii] = "P" if style_atoms[ii] != "L" else "L"
            elif nn == 1:
                style_atoms[ii] = "P" if style_atoms[ii] != "L" else "L"
            else:
                style_atoms[ii] = "NP" if style_atoms[ii] != "L" else "L"

    return style_atoms

def assgin_cum(elements,connectivitys,bond_types,style_atoms,all_bond_types,formal_charges):
    unassign_atoms = get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types)
    for ii in unassign_atoms:
        if style_atoms[ii] == "CU":
            for jj in connectivitys[ii]:
                bond_types[get_bond_name(ii,jj,all_bond_types)] = 2
                if style_atoms[jj] == "UN":
                    formal_charges[jj] = 2 - Element.get(elements[ii]).valents[0]
            oxide = sum([bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]])
            formal_charges[ii] = oxide - Element.get(elements[ii]).valents[0]

    return bond_types, formal_charges

def get_block_chain(atoms,connectivitys):
    tmp_connects = {}
    for ii in atoms:
        tmp_connects[ii] = []
        for jj in connectivitys[ii]:
            if jj in atoms:
                tmp_connects[ii].append(jj)
    arrs = find_mole(tmp_connects)
    block_arrs = []
    for arr in arrs:
        terminal_atoms = [ii for ii in arr if len(tmp_connects[ii]) == 1]

        block_arrs.extend(search_chain(terminal_atoms,{ii:rr for ii,rr in tmp_connects.items() if ii in arr}))
    return block_arrs

def get_conjugate_bond_type(block_arrs,
                            elements,
                            coordinates,
                            connectivitys,
                            bond_types,
                            style_atoms,
                            plane_label,
                            all_bond_types,
                            unassign_atoms,
                            length_bond_types,
                            conju_length_bond_types
                            ):
    for arr in block_arrs:
        for ii in arr:
            
            if style_atoms[ii] in plane_label:#["P","SP","L"]:
                
                bts_ii = [bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]]
                unbond_atoms = [jj for jj in connectivitys[ii] if bond_types[(get_bond_name(ii,jj,all_bond_types))] == "N"]
                if 2 in bts_ii or 3 in bts_ii:
                    for jj in connectivitys[ii]:
                        if bond_types[get_bond_name(ii,jj,all_bond_types)] != "N":
                            bond_types[get_bond_name(ii,jj,all_bond_types)] == 1
                else:
                    double_atoms = [jj for jj in unbond_atoms if length_bond_types[get_bond_name(ii,jj,all_bond_types)] >= 2]
                    single_atoms = [jj for jj in unbond_atoms if length_bond_types[get_bond_name(ii,jj,all_bond_types)] == 1]
                    tmp = scale_length_bond(ii,double_atoms,elements,coordinates)
                    double_bond_level = sorted([[jj,tt] for jj,tt in tmp.items()],key=lambda x:x[1])

                    flag = True
                    for rr in double_bond_level:
                        bts_jj = [bond_types[get_bond_name(rr[0],jj,all_bond_types)] for jj in connectivitys[rr[0]]]
                        if 2 in bts_jj or 3 in bts_jj:
                            bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 1
                        else:
                            if flag:
                                if style_atoms[ii] in ["L"]:
                                    bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 3
                                else:
                                    bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 2
                                flag = False
                            else:
                                bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 1
                    for jj in single_atoms:
                        bond_types[get_bond_name(ii,jj,all_bond_types)] = 1
    return bond_types

def get_ring_conjugate_bond_type(arr,elements,coordinates,connectivitys,bond_types,style_atoms,plane_label,all_bond_types):
    for ii in arr:
        if style_atoms[ii] in plane_label:#["P","SP","L"]:
            bts_ii = [bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in connectivitys[ii]]
            unbond_atoms = [jj for jj in connectivitys[ii] if bond_types[(get_bond_name(ii,jj,all_bond_types))] == "N"]
            if 2 in bts_ii or 3 in bts_ii:
                for jj in connectivitys[ii]:
                    #if bond_types[get_bond_name(ii,jj,all_bond_types)] != "N":
                    if jj in unbond_atoms:
                        bond_types[get_bond_name(ii,jj,all_bond_types)] = 1
            else:
                tmp = scale_length_bond(ii,unbond_atoms,elements,coordinates)
                double_bond_level = sorted([[jj,tt] for jj,tt in tmp.items()],key=lambda x:x[1])
                flag = True
                for nn,rr in enumerate(double_bond_level):
                    bts_jj = [bond_types[get_bond_name(rr[0],jj,all_bond_types)] for jj in connectivitys[rr[0]]]
                    if 2 in bts_jj or 3 in bts_jj:
                        bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 1
                    else:
                        if flag:
                            if style_atoms[ii] in ["L"]:
                                bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 3
                                flag = False
                            else:
                                if rr[0] not in arr and elements[rr[0]] in ["C","P"]:
                                    if nn == len(double_bond_level) - 1:
                                        bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 2
                                        flag = False
                                    else:
                                        bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 1
                                else:
                                    bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 2
                                    flag = False
                        else:
                            bond_types[get_bond_name(ii,rr[0],all_bond_types)] = 1
    return bond_types

def assign_conjugate_chain(elements,
                           coordinates,
                           connectivitys,
                           bond_types,
                           style_atoms,
                           all_bond_types,
                           ring_atoms,
                           length_bond_types,
                           conju_length_bond_types
                           ):
    plane_label = ["P","L"]
    unassign_atoms = get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types)
    ring_unassign_atoms = []
    chain_unassign_atoms = []
    for ii in unassign_atoms:
        if ii not in ring_atoms:
            chain_unassign_atoms.append(ii)
        else:
            ring_unassign_atoms.append(ii)
    bond_types = get_ring_conjugate_bond_type(ring_unassign_atoms,elements,coordinates,connectivitys,bond_types,style_atoms,plane_label,all_bond_types)
    
    for bn,bt in bond_types.items():
        if bt == "N":
            ii,jj = [int(aa) for aa in bn.split("-")]
            dist = calculate_distance(coordinates[ii],coordinates[jj])
            double_bond_length = get_bonded_type_distance(elements[ii], elements[jj], "2")
            triple_bond_length = get_bonded_type_distance(elements[ii], elements[jj], "3")
            if style_atoms[ii] in plane_label or style_atoms[jj] in plane_label:
                if triple_bond_length is not None and dist <= triple_bond_length * 1.02:
                    bond_types[bn] = 3
                elif dist <= double_bond_length * 1.02:
                    bond_types[bn] = 2
                else:
                    if (elements[ii] == "S" and style_atoms[ii] == "UN") or (elements[jj] == "S" and style_atoms[jj] == "UN"):
                        if dist <= double_bond_length * 1.04:
                            bond_types[bn] = 2
    
    #block_arrs = get_block_chain(chain_unassign_atoms,connectivitys)
    block_arrs = [chain_unassign_atoms]
    bond_types = get_conjugate_bond_type(block_arrs,
                                         elements,
                                         coordinates,
                                         connectivitys,
                                         bond_types,
                                         style_atoms,
                                         plane_label,
                                         all_bond_types,
                                         chain_unassign_atoms,
                                         length_bond_types,
                                         conju_length_bond_types
                                         )
    return bond_types

def assign_charge_before_add_H(elements,connectivitys,bond_types,all_bond_types):
    for bn,bt in bond_types.items():
        if bt == "N":
            bond_types[bn] = 1
    formal_charge = {}
    for ii,elem in elements.items():
        if elem in ["N"]:
            conn = connectivitys[ii]
            oxide = sum([bond_types[get_bond_name(ii,jj,all_bond_types)] for jj in conn])
            if oxide > 3 and oxide != 5:
                formal_charge[ii] = oxide - 3
                for jj in conn:
                    if elements[jj] == "O":
                        formal_charge[jj] = -1
                    #elif elements[jj] == "N" and bond_types[get_bond_name(ii,jj,all_bond_types)] == 2:
                    #    formal_charge[jj] = -1
                    #elif elements[jj] == "C" and bond_types[get_bond_name(ii,jj,all_bond_types)] == 3:
                    #    formal_charge[jj] = -1
    
    return formal_charge    

def get_all_bond_type_based_on_length(elements,coordinates,connectivitys,all_bond_types):
    length_bond_types = {}
    conju_length_bond_types = {}
    for bn in all_bond_types:
        ii,jj = [int(aa) for aa in bn.split("-")]
        length_bond_types[bn] =  get_bond_type_based_on_length(elements[ii],elements[jj],coordinates[ii],coordinates[jj],conju_flag=False)
        conju_length_bond_types[bn] = get_bond_type_based_on_length(elements[ii],elements[jj],coordinates[ii],coordinates[jj],conju_flag=True)
    return length_bond_types,conju_length_bond_types

def check_double_NP_bond(elements,coordinates,connectivitys,bond_types,style_atoms,length_bond_types,all_bond_types):
    for bt,bn in bond_types.items():
        ii,jj = [int(aa) for aa in bt.split("-")]
        if bn == 1:
            if style_atoms[ii] == "NP" or style_atoms[jj] == "NP": 
                bns_ii  = [bond_types[get_bond_name(ii,kk,all_bond_types)] for kk in connectivitys[ii]]
                bns_jj  = [bond_types[get_bond_name(jj,kk,all_bond_types)] for kk in connectivitys[jj]]
                if set(bns_ii) == set([1]) and set(bns_jj) == set([1]):
                    #e1 = elements[ii]
                    #e2 = elements[jj]
                    #coor1 = coordinates[ii]
                    #coor2 = coordinates[jj]
                    #dist = calculate_distance(coor1,coor2)
                    #d = get_bonded_type_distance(e1, e2, "d") # single_bond_length
                    if length_bond_types[get_bond_name(ii,jj,all_bond_types)] not in [1,None]:
                        bond_types[bt] = length_bond_types[get_bond_name(ii,jj,all_bond_types)]
                else:
                    e1 = elements[ii]
                    e2 = elements[jj]
                    if e1 == "C" and e2 == "C":
                        bn_2_ii = connectivitys[ii][bns_ii.index(2)] if 2 in bns_ii else -1
                        bn_2_jj = connectivitys[jj][bns_jj.index(2)] if 2 in bns_jj else -1
                        bn_2_ii_elem = elements[bn_2_ii] if bn_2_ii != -1 else "NONE"
                        bn_2_jj_elem = elements[bn_2_jj] if bn_2_jj != -1 else "NONE"
                        if bn_2_ii_elem == "N" or bn_2_jj_elem == "N":
                            coor1 = coordinates[ii]
                            coor2 = coordinates[jj]
                            if calculate_distance(coor1,coor2) <= get_bonded_type_distance(e1, e2, "2") * 1.02:
                                bond_types[bt] = 2
                                if bn_2_ii != -1:
                                    bond_types[get_bond_name(ii,bn_2_ii,all_bond_types)] = 1
                                if bn_2_jj != -1:
                                    bond_types[get_bond_name(jj,bn_2_jj,all_bond_types)] = 1
            elif style_atoms[ii] == "SP" and style_atoms[jj] == "SP":
                bns_ii  = [bond_types[get_bond_name(ii,kk,all_bond_types)] for kk in connectivitys[ii]]
                bns_jj  = [bond_types[get_bond_name(jj,kk,all_bond_types)] for kk in connectivitys[jj]]    
                if 2 not in bns_ii and 2 not in bns_jj:
                    if length_bond_types[get_bond_name(ii,jj,all_bond_types)] not in [1,None]:
                        bond_types[bt] = length_bond_types[get_bond_name(ii,jj,all_bond_types)]
        elif (style_atoms[ii] == "L" and style_atoms[jj] == "UN") or (style_atoms[jj] == "L" and style_atoms[ii] == "UN"):
                if bn == 2:
                    if  length_bond_types[get_bond_name(ii,jj,all_bond_types)] == 3:
                        bond_types[bt] = 3
    return bond_types

def map_connectivity(elements,coordinates):
    connectivitys = connected_atoms(elements,coordinates)
    connectivitys = over_connected_check(elements,coordinates,connectivitys)
    return connectivitys

def map_bond_type(elements,coordinates,connectivitys,find_all=False):
    cyclo = find_cyclo(connectivitys)
    ring_atoms = list(set([ii for ring in cyclo for ii in ring[:-1]]))
    rings,ring_blocks,ring_blcok_components = cyclo_blocks(cyclo,connectivity=connectivitys,ring_atoms=ring_atoms)

    small_ring_atoms = []
    for rn,ring in rings.items():
        if len(ring) in [4,5]:
            small_ring_atoms.extend(ring)
    small_ring_atoms = set(small_ring_atoms)

    bond_types = create_initial_bond_types(elements,connectivitys)
    all_bond_types = list(bond_types.keys())
    length_bond_types,conju_length_bond_types = get_all_bond_type_based_on_length(elements,coordinates,connectivitys,all_bond_types)
    elem_conn = classify_based_on_element_connects(elements,connectivitys,ring_atoms)
    
    style_atoms = {ii:get_linear_plane_single(ii,elements,coordinates,connectivitys,ring_atoms,small_ring_atoms,length_bond_types,conju_length_bond_types,all_bond_types) for ii in elements}
    
    formal_charges = {ii:"N" for ii in range(len(elements))}
    ##饱和的情况 
    for ss in ["H1","F1","Cl1","Br1","I1","O2","S2","C4","N4"]:
        if ss in elem_conn:
            bond_types = assign_single_bond(elem_conn[ss],connectivitys,bond_types,all_bond_types)
    
    ##酸结构
    for ss in ["Cl4","Br4","S4","P4","Cl3","Br3","S3","N3A","P3A"]:
        if ss in elem_conn:
            bond_types = assign_acid_bond_type(elem_conn[ss],elements,coordinates,connectivitys,bond_types,all_bond_types,ss,conju_length_bond_types)
    #除去酸后，剩下的饱和的N,P
    for ss in ["N3","P3"]:
        if ss in elem_conn:
            bond_types = assign_single_bond(elem_conn[ss],connectivitys,bond_types,all_bond_types)
    #剩下没有指定键类型的原子
    #unassign_atoms = get_unassign_bond_type_atoms(connectivitys,bond_types,all_bond_types)
    #style_atoms = {ii:get_linear_plane_single(ii,coordinates,connectivitys) for ii in unassign_atoms}
    #线性的原子组合，暂时被关闭
    #bond_types, formal_charges = assign_linear_bond(elements,connectivitys,bond_types,style_atoms,all_bond_types,formal_charges)
    #处理N=N=N结构
    #bond_types, formal_charges = assign_three_n(elements,connectivitys,bond_types,style_atoms,all_bond_types,formal_charges)

    ###下面暂时被关闭了
    #bond_types = assign_plane_line(style_atoms,connectivitys,bond_types,all_bond_types)

    #非平面的原子的单键
    np_atoms = [ii for ii,st in style_atoms.items() if st == "NP"]
    bond_types = assign_single_bond(np_atoms,connectivitys,bond_types,all_bond_types)
    #迭代平面和线性原子,对于平面原子，如果有两个单键，则剩下的为双键，如果有一个双键，则剩下的为单键，
    #                 对于线性原子，如果一个为三键，剩下的为单键，反之亦然
    bond_types = assign_bond_for_plane_linear(elements,connectivitys,bond_types,style_atoms,all_bond_types)
    #("after plane linear",change_bond_types_by_ii_add_1(bond_types))
    #rings
    
    bond_types_arr = assign_ring_double_bond(rings,
                                         ring_blcok_components,
                                         elements,
                                         coordinates,
                                         connectivitys,
                                         bond_types,
                                         style_atoms,
                                         all_bond_types,
                                         length_bond_types,
                                         conju_length_bond_types,
                                         find_all=find_all)

    #再一次迭代平面和线性原子
    complete_bond_types_arr = []
    complete_formal_charges_arr = []

    for bond_types in bond_types_arr:
        bond_types = assign_bond_for_plane_linear(elements,connectivitys,bond_types,style_atoms,all_bond_types)
        style_atoms = judge_style_of_unassign_atoms(elements,
                                                coordinates,
                                                connectivitys,
                                                bond_types,
                                                style_atoms,
                                                all_bond_types,
                                                length_bond_types,
                                                conju_length_bond_types)
        # cumene类型的键（两个双键直接相连）
        bond_types,formal_charges = assgin_cum(elements,connectivitys,bond_types,style_atoms,all_bond_types,formal_charges)
        np_atoms = [ii for ii,st in style_atoms.items() if st == "NP"]
        bond_types = assign_single_bond(np_atoms,connectivitys,bond_types,all_bond_types)
        bond_types = assign_bond_for_plane_linear(elements,connectivitys,bond_types,style_atoms,all_bond_types)
        #conjugate chain atoms
        bond_types = assign_conjugate_chain(elements,
                                            coordinates,
                                            connectivitys,
                                            bond_types,
                                            style_atoms,
                                            all_bond_types,
                                            ring_atoms,
                                            length_bond_types,
                                            conju_length_bond_types
                                            )
        
        bond_types = check_double_NP_bond(elements,coordinates,connectivitys,bond_types,style_atoms,length_bond_types,all_bond_types)
        
        #
        formal_charges = assign_charge_before_add_H(elements,connectivitys,bond_types,all_bond_types)
        complete_bond_types_arr.append(bond_types)
        complete_formal_charges_arr.append(formal_charges)

    return complete_bond_types_arr,complete_formal_charges_arr

def coordinate_to_bond_type(elements,coordinates,connectivity=None,method="all",find_all=False):
    elements = {ii:elem for ii,elem in enumerate(elements)}
    coordinates = {ii:coor for ii,coor in enumerate(coordinates)}
    if connectivity is None:
        connectivitys = map_connectivity(elements,coordinates)
    else:
        connectivitys = {ii:conn for ii,conn in enumerate(connectivity)}
    molecules = find_mole(connectivitys)
    if len(molecules) == 1:
        if method == "all":
            bond_type, formal_charge = map_bond_type(elements,coordinates,connectivitys,find_all=find_all)
    
        return connectivitys, bond_type, formal_charge
    else:
        return {},{},{}
    
def change_bond_types_by_ii_add_1(bond_types):
    tmp = {}
    for aa,bb in bond_types.items():
        s = aa.split("-")
        ss = f"{int(s[0])+1}-{int(s[1])+1}"
        tmp[ss] = bb
    return tmp
# todo for caofl
def arom_to_singledouble(info):
    pass

def molecule_coordinate_to_bond_type(molecule):
    connectivity,bond_type,formal_charge = coordinate_to_bond_type(molecule.elements,molecule.coordinates,find_all=False)


if __name__ == "__main__":
    from Chem import Chem
    import sys
    molecule = Chem.molecule_create(sys.argv[1])[0]
    if len(sys.argv) > 2:
        if sys.argv[2] == "y":
            find_all = True
        else:
            find_all = False
    else:
        find_all = False
    connectivity,bond_type,forma_charge = coordinate_to_bond_type(molecule.elements,molecule.coordinates,find_all=find_all)


