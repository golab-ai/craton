#!/usr/bin/env python
"""
"""
from collections import OrderedDict
from copy import deepcopy
from ...chem.constants import SINGLE_CONNECT_ATOM_OR_GROUP as SCAOG
from ...chem.constants import LOCAL_LINK,LOCAL_CYCLE,LOCAL_CYCLE_B,LOCAL_CYCLE_F
from ...chem.constants import MULTI_BOND, SINGLE_BOND, CONJU_BOND
from ...chem.constants import HBOND_ANGLE_MIN, HBOND_DIST_MAX
from ...chem.constants import RING_BOND, CHAIN_BOND
from ...utils.geometry import calc_stru_para
from dataclasses import dataclass, field
from ...chem.group import Group

def find_mole_old(
    connect_dict,
):
    molecules = []
    for ii,conn in connect_dict.items():
        tmp = [ii] + connect_dict[ii]
        in_flag = False
        for arr in molecules:
            if len(set(tmp).intersection(set(arr))) > 0:
                arr.extend(tmp)
                in_flag = True
                break
        if not in_flag:
            molecules.append(tmp)
    
    while 1:
        n = len(molecules)
        m_tmp = []
        for rr in molecules:
            in_flag = False
            for rr1 in m_tmp:
                if len(set(rr).intersection(set(rr1))) > 0:
                    rr1.extend(rr)
                    in_flag = True
            if not in_flag:
                m_tmp.append(rr)
        molecules = deepcopy(m_tmp)
        if len(molecules) == n:
            break

    molecules = [list(set(arr)) for arr in molecules]
    return molecules

def find_mole(
    connect_dict,
):
    def rearrange_dict(dict, del_arr):
        arr = [i for i in range(0, len(dict)) if i not in del_arr]
        new_dict = {}
        for i in range(0, len(arr)):
            new_dict[i] = dict[arr[i - 1]]
        return new_dict
    this_an = list(connect_dict.keys())
    for an,conn in connect_dict.items():
        connect_dict[an] = [aj for aj in conn if aj in this_an]
    mole = {0: []}
    mole_n = 0
    for i in this_an:
        arr = [i] + connect_dict[i]
        samefrag = []
        for j in range(0, len(mole)):
            if len(set(mole[j]).intersection(set(arr))) > 0:
                mole[j] = list(set(mole[j]).union(set(arr)))
                samefrag.append(j)
        if len(samefrag) == 0:
            mole[mole_n] = arr
            mole_n += 1
            mole[mole_n] = []
        elif len(samefrag) > 1:
            samefrag.sort()
            for k in range(1, len(samefrag)):
                mole[samefrag[0]] = list(set(mole[samefrag[0]]).union(set(mole[samefrag[k]])))
            del samefrag[0]
            mole = rearrange_dict(mole, samefrag)
            mole_n = mole_n - len(samefrag)
            mole[mole_n] = []
    for i in range(0, len(mole)):
        mole[i] = sorted(mole[i])
    del mole[len(mole)-1]
    return [mm for __,mm in mole.items() if len(mm) > 0]

def get_2d_connectivity(molecule):
    element_number = {ii:elem for ii,elem in enumerate(molecule.elements)}
    all_connect = {ii:connect for ii,connect in enumerate(molecule.connectivity)}
    return element_number, all_connect

def remove_hydrogen_atoms(element_number,all_connect):
    reduce_connect = deepcopy(all_connect)
    H_atom_arr = []
    for item in all_connect.keys():
        if element_number[item] in SCAOG:
            del reduce_connect[item]
            H_atom_arr.append(item)
    for item in reduce_connect.keys():
        tmp = reduce_connect[item][:]
        n = len(reduce_connect[item])
        for i in range(0, n):
            if reduce_connect[item][i] in H_atom_arr:
                tmp.remove(reduce_connect[item][i])
        reduce_connect[item] = tmp[:]
    terminal_atom_arr = []
    for item in reduce_connect.keys():
        if len(reduce_connect[item]) == 1:
            terminal_atom_arr.append(int(item))
    return reduce_connect, terminal_atom_arr, H_atom_arr

def assign_structure_info(molecule, rings,ring_blocks,ring_block_components):
    #rings = {}
    #ring_stru = []
    #for arr in cyclos:
    #    n = len(arr) - 1
    #    nn = int(n / 2)
    #    ss = f"R{n}-{arr[0]}_{arr[nn]}_{arr[-1]}"
    #    if ss in rings.keys():
    #        ss = f"R{n}-{arr[0]}_{arr[nn]}--{randint(0, 100)}_{arr[-1]}"
    #    rings[ss] = arr
#
    #tmp = []
    #for aa, bb in rings.items():
    #    tmp.append(bb[:-1])
    #    ring_stru = combine_arr(tmp)

    for i in range(0, len(molecule.Atoms)):
        inring = []
        ring_size = []
        ring_prop = []
        for aa, bb in rings.items():
            if i in bb:
                inring.append(aa)
                ring_size.append(len(bb) - 1)
                ring_prop.append(bb[-1])
        molecule.Atoms[i].has_ring = inring
        molecule.Atoms[i].has_ring_size = ring_size
        molecule.Atoms[i].has_ring_property = ring_prop

    molecule.rings = rings
    molecule.ring_blocks = ring_blocks
    molecule.ring_block_components = ring_block_components

def assign_local_info(molecule,reduce_connect):
    cycloset = []
    for __, cyclo in molecule.ring_dict.items():
        cycloset.extend(cyclo[:-1])

    cycloset = set(cycloset)

    for atom_n,atom in enumerate(molecule.Atoms):
        if atom.elem in SCAOG:
            atom.local = "EN"
            continue

        connect_atom_number_in_ring = len(set(reduce_connect[atom_n]) & cycloset)
        connect_atom_number_out_ring = len(reduce_connect[atom_n]) - connect_atom_number_in_ring
        
        if atom_n not in cycloset:
            if len(reduce_connect[atom_n]) == 1:
                if connect_atom_number_in_ring == 0:
                    atom.local = "LT"
                else:
                    atom.local = "LTC"
            else:
                if connect_atom_number_in_ring == 0:
                    atom.local = "LM"
                else:
                    atom.local = "LBC" if connect_atom_number_in_ring > 1 else "LMC"
            continue

        if len(atom.inring) == 1:
            if connect_atom_number_in_ring == 2:
                atom.local = "CL" if connect_atom_number_out_ring != 0 else "C"
            else:
                atom.local = "CBL" if connect_atom_number_out_ring != 0 else "CB"
            continue

        bridge = 0
        for connect_n in atom.connect:
            if connect_n in cycloset and len(set(molecule.Atoms[connect_n].inring) & set(atom.inring)) == 0:
                bridge += 1
        if connect_atom_number_in_ring - bridge == 2:
            atom.local = "CFM" if bridge == 0 and connect_atom_number_out_ring == 0 else ("CFMB" if bridge != 0 else "CFML")
        elif connect_atom_number_in_ring - bridge == 3:
            atom.local = "CF" if bridge == 0 and connect_atom_number_out_ring == 0 else ("CFB" if bridge != 0 else "CFL")
        else:
            atom.local = "CS" if bridge == 0 and connect_atom_number_out_ring == 0 else ("CSB" if bridge != 0 else "CSL")

def get_atom_flag(x):
    if x in LOCAL_LINK:
        return "c"
    else:
        if x in LOCAL_CYCLE_F:
            return "f"
        else:
            if x in LOCAL_CYCLE_B:
                return "b"
            else:
                return "x"

def get_connect(bn,str,c1,c2):
    if str == "bb":
        str = "nr" if len(set(c1).intersection(c2)) == 0 else "rr"
    index_label = [["cc"],["cb","cf","cx","bc","fc","xc"],["ff"],["nr"],
             ["xx","bf","bx","fb","fx","xf","xb","rr"]]
    for i,arr in enumerate(index_label):
        if str in arr:
            index = i + 1
            break
    types = {
            1: ["S", "D", "T", "J", "M"], #1
            2: ["eS", "eD", "eM", "es", "eM"],#2
            3: ["fs", "fd", "ft", "fr", "fm"],#3
            4: ["bS", "bD", "bm", "br", "bm"],#4
            5: ["s", "d", "t", "r", "m"],#5
            }
    btype = ["1", "2", "3", "ar"]
    if bn not in btype:
        return types[index][-1]
    else:
        return types[index][btype.index(bn)]
    
def get_aromatic_connect(molecule):
    arom_atom_dict = {}
    for __, ring in molecule.rings.items():
        if ring[-1] in ["ar1", "ar2"]:
            for ai in ring[:-1]:
                if ai not in arom_atom_dict:
                    arom_atom_dict[ai] = []
                arom_atom_dict[ai] = list(set(arom_atom_dict[ai]).union(set(ring[:-1])))

    for atom in molecule.Atoms:
        #if not hasattr(atom, "bond_type_old"):
        #    atom.bond_type_old = deepcopy(atom.bond_type)
        if atom.ID in arom_atom_dict:
            atom.bond_type_aromatic = []
            for i in range(len(atom.connectivity)):
                if atom.connectivity[i] in arom_atom_dict[atom.No]:
                    atom.bond_type_aromatic.append("ar")
                else:
                    atom.bond_type_aromatic.append(atom.bond_type[i])
        else:
            atom.bond_type_aromatic = deepcopy(atom.bond_type)

def assign_connect_type(molecule):
    get_aromatic_connect(molecule)
    for atom_n, atom in enumerate(molecule.Atoms):
        atom_cyclo = atom.inring
        atom_str = get_atom_flag(atom.local)
        atom.connectivity_type = []
        for atom_m, conn_an in enumerate(atom.connect):
            atom_conn = molecule.Atoms[conn_an]
            atom_conn_cyclo = atom_conn.inring
            atom_conn_str = get_atom_flag(atom_conn.local)
            bn = atom.bond_type_aromatic[atom_m]
            atom.connectivity_type.append(get_connect(bn,f"{atom_str}{atom_conn_str}",atom_cyclo,atom_conn_cyclo))

def assign_conjugate_info(molecule):
    for atom in molecule.Atoms:
        atom.bond_type_conjugate = ["" for _ in atom.connect]
    conju_atoms = []
    for ii,atom in enumerate(molecule.Atoms):
        conju_flag = False
        if set(atom.connectivity_type) & set(MULTI_BOND) == set():
            continue
        for jj,cbn in enumerate(atom.connectivity_type):
            atom1 = molecule.Atoms[atom.connectivity[jj]]
            if cbn in SINGLE_BOND and len(set(atom1.connectivity_type) & set(MULTI_BOND)) > 0:
                conju_flag = True
            if conju_flag:
                break
        
        if conju_flag:
            conju_atoms.append(ii)

    for ii,atom in enumerate(molecule.Atoms):
        for jj,an in enumerate(atom.connectivity):
            if an in conju_atoms:
                if atom.connectivity_type[jj] in CHAIN_BOND:
                    atom.bond_type_conjugate[jj] = "J"
                elif atom.connectivity_type[jj] in RING_BOND:
                    atom.bond_type_conjugate[jj] == "j"
                
def assign_conjugate_info_old(molecule):
    for atom in molecule.Atoms:
        atom.bond_type_conjugate = ["" for _ in atom.connect]
    conju_atoms = []
    for i in range(len(molecule.Atoms)):
        a = molecule.Atoms[i]
        conju_flag = False
        if set(a.connectivity_type) & set(MULTI_BOND) == set():
            continue
        for j in range(len(a.connectivity_type)):
            if (
                a.connectivity_type[j] in SINGLE_BOND
                and len(
                    set(molecule.Atoms[a.connect[j]].connectivity_type)
                    & set(MULTI_BOND + SINGLE_BOND)
                )
                > 0
            ):
                if a.connectivity_type[j] in ["S", "eS", "bS"]:
                    # if a.bond_type_detail[j] in ["S"]:
                    a.bond_type_conjugate[j] = "J"
                else:
                    a.bond_type_conjugate[j] = "j"
                conju_flag = True
                conju_atoms.append(i)
        if conju_flag:
            for j in range(len(a.connectivity_type)):
                # if a.bond_type_detail in chem_const.single_bond + chem_const.double_bond:
                if a.connectivity_type[j] in ["D", "eD", "bD", "T"]:
                    # if a.bond_type_detail[j] in ["D", "T"]:
                    a.bond_type_conjugate[j] = "J"
                elif a.connectivity_type[j] in ["fd", "d", "ft", "t"]:
                    a.bond_type_conjugate[j] = "j"
    for conju_atom in conju_atoms:
        a = molecule.Atoms[conju_atom]
        for jj in range(len(a.connectivity)):
            if a.connectivity_type[jj] in CONJU_BOND:
                aa = molecule.Atoms[a.connectivity[jj]]
                if aa.connectivity_type[aa.connectivity.index(conju_atom)] not in CONJU_BOND:
                    aa.bond_type_conjugate[aa.connectivity.index(conju_atom)] = a.bond_type_conjugate[jj]

def divide_charge(molecule):
    """
    Sum out formal charges between bonded opposite charges
    Average formal charges between several uni-coordinated O/S on C/S/P
    Average formal charges between several tri-coordinated N on C
    Spread formal charges to all heteroatoms in a conjugated system (TODO Only N-containing 5-member aromatic rings are considered)
    """
    charge_group = []
    for atom in molecule.Atoms:
        atom.primitive_formal_charge = atom.formal_charge
    for bond in molecule.Bonds:
        atom1, atom2 = molecule.Atoms[bond.a1], molecule.Atoms[bond.a2]
        if {atom1.formal_charge, atom2.formal_charge} == {1, -1}:
            atom1.primitive_formal_charge = 0
            atom2.primitive_formal_charge = 0
            flag1 = atom1.elem == "N" and atom2.elem == "O" \
                    and len(atom1.connectivity) == 3 and "2" in atom1.bond_type 
            flag2 = atom2.elem == "N" and atom1.elem == "O" \
                    and len(atom2.connectivity) == 3 and "2" in atom2.bond_type
            if flag1 or flag2:
                atom1.bond_type[atom1.connectivity.index(atom2.ID)] = "2"
                atom2.bond_type[atom2.connectivity.index(atom1.ID)] = "2"
                ######Changed by CFL on 2025.2.27#######
                atom1.formal_charge = 0
                atom2.formal_charge = 0
    _processed = []
    for atom in molecule.Atoms:
        if atom in _processed:
            continue
        if atom.elem in ["O", "S"] and len(atom.connect) == 1 and atom.primitive_formal_charge == -1:
            root = molecule.Atoms[atom.connect[0]]
            if root.elem not in ["C", "S", "P"]:
                continue
            neighbors = [
                molecule.Atoms[i]
                for i in root.connect
                if molecule.Atoms[i].elem in ["O", "S"] and len(molecule.connectivity[i]) == 1
            ]
            if len(neighbors) > 1:
                charge = sum(neigh.primitive_formal_charge for neigh in neighbors)
                for neigh in neighbors:
                    neigh.primitive_formal_charge = charge / len(neighbors)
                _processed.extend(neighbors)
                if charge > 0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group.append([atom.ID for atom in neighbors]+[root.ID,label])

    #####改该函数必测试：[O-]CC(N1N=N[NH+]=C(N)1)C([O-])=O,result:[10, 11, 9], [5, 4, 3, 2, 0, 1], [8]
    rn_used = []
    _processed = []
    for atom in molecule.Atoms:
        if atom in _processed:
            continue
        if atom.elem in ["N"] and len(atom.connect) == 3 and atom.primitive_formal_charge == 1:
            for idx_neigh in atom.connect:
                root = molecule.Atoms[idx_neigh]
                if root.elem not in ["C"]:
                    if atom.formal_charge> 0:
                        label = "positive"
                    else:
                        label = "negative"
                    charge_group.append([atom.ID,label])
                    continue
                neighbors = [
                    molecule.Atoms[i]
                    for i in root.connect
                    if molecule.Atoms[i].elem in ["N"] and len(molecule.connectivity[i]) == 3
                ]
                if len(neighbors) > 1:
                    for ring in atom.has_ring:
                        ss = ring.split("_")
                        if ring[-3:] in ["ar1","ar2"] and ring[1] in ["5","6"]:
                            members = [molecule.Atoms[an] for an in molecule.rings[ring][:-1]]
                            neighbors += [atom for atom in members if atom.elem == "N"]
                            rn_used.append(ring)
                    neighbors = list(set(neighbors))
                    charge = sum(neigh.primitive_formal_charge for neigh in neighbors)
                    for neigh in neighbors:
                        neigh.primitive_formal_charge = charge / len(neighbors)
                    _processed.extend(neighbors)
                    if charge > 0:
                        label = "positive"
                    else:
                        label = "negative"
                    charge_group.append([atom.ID for atom in neighbors]+[root.ID,label])
                    break
    
    for rn,ring in molecule.ring_dict.items():
        if len(ring) != 6 or ring[-1] not in ["ar1", "ar2"]:
            continue
        if rn not in rn_used:
            members = [molecule.Atoms[i] for i in ring[:-1]]
            nitrogens = [atom for atom in members if atom.elem == "N"]
            charge = sum(atom.primitive_formal_charge for atom in nitrogens)
            if charge != 0:
                if charge > 0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group.append([atom.ID for atom in nitrogens])
            for atom in nitrogens:
                atom.primitive_formal_charge = charge / len(nitrogens)

    _tmp = []
    for rr in charge_group:
        _tmp.extend(rr)
    for atom in molecule.Atoms:
        if atom.formal_charge != 0:
            if atom.ID not in _tmp:
                if  atom.formal_charge> 0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group.append([atom.ID,label])

    _charge_group = {}
    for ii,vv in enumerate(charge_group):
        _charge_group[f"charge_{ii}"] = vv
    molecule.charge_group = _charge_group

def flexible_torsion(molecule,idx=None):
    ignore_atom_style = ["dummy"]
    if len(molecule.Atoms) < 2:
        return None
    for bb in molecule.Bonds:
        a1 = bb.a1
        a2 = bb.a2
        if (
            molecule.Atoms[a1].bond_type[molecule.Atoms[a1].connect.index(a2)] == "1"
            and molecule.Atoms[a1].s not in ignore_atom_style
            and molecule.Atoms[a2].s not in ignore_atom_style
        ):
            if molecule.Atoms[a1].connectivity_type[molecule.Atoms[a1].connect.index(a2)] in [
                "S",
                "eS",
                "bS",
                "br",
                "J",
                "eJ",
                "bJ",
            ]:
                if molecule.Atoms[a1].local != "EN" and molecule.Atoms[a2].local != "EN":
                    if molecule.Atoms[a1].local not in ["LT", "LTC"] and molecule.Atoms[a2].local not in ["LT", "LTC"]:
                        bb.flexible = "yes"
                    else:
                        if len(molecule.Atoms[a1].connect) > 1 and len(molecule.Atoms[a2].connect) > 1:
                            cea = None
                            if molecule.Atoms[a1].local in ["LT", "LTC"]:
                                cea = a1
                            else:
                                if molecule.Atoms[a2].local in ["LT", "LTC"]:
                                    cea = a2
                            if cea is not None:
                                if molecule.Atoms[cea].elem != "C":
                                    bb.flexible = "yes"
                                else:
                                    elems = []
                                    for ca in molecule.Atoms[cea].connect:
                                        if ca != a1 and ca != a2:
                                            elems.append(molecule.Atoms[ca].elem)
                                    if len(set(elems)) > 1:
                                        bb.flexible = "yes"   
    if idx is not None:     
        return molecule,idx

def scan_torsion(molecule, ignore_alkane=False,idx=None):
    """
    确定scan_torsion。
    输入：
        m: Molecule
        ignore_alkane: True or False。是否忽略烷基链部分的二面角
    输出：
        scan_term: List[[int]]可以被旋转的二面角
    """
    ignore_atom_style = ["dummy"]
    scan_term = []
    if len(molecule.Atoms) < 4:
        molecule.torsions = []
    for bb in molecule.Bonds:
        if not hasattr(bb, "flexible") or bb.flexible != "yes":
            continue
        a1 = bb.a1
        a2 = bb.a2
        if "3" in molecule.Atoms[a1].bond_type or "3" in molecule.Atoms[a2].bond_type:
            continue
        #ats = [molecule.Atoms[a1].atom_type_name, molecule.Atoms[a2].atom_type_name]
        neighbors1 = [molecule.Atoms[i] for i in molecule.Atoms[a1].connect if i != a2 and molecule.Atoms[i].s not in ignore_atom_style]
        atom0 = neighbors1[0]
        for a in neighbors1[1:]:
            if a.mass > atom0.mass:
                atom0 = a
        #ats.append(atom0.atom_type_name)
        neighbors2 = [molecule.Atoms[i] for i in molecule.Atoms[a2].connect if i != a1 and molecule.Atoms[i].s not in ignore_atom_style]
        atom3 = neighbors2[0]
        for a in neighbors2[1:]:
            if a.mass > atom3.mass:
                atom3 = a
        #ats.append(atom3.atom_type_name)
        #if ignore_alkane and set(ats).issubset({"c_4", "c_4h", "c_4h2", "c_4h3", "c_4h4", "h_1"}):
        #    continue
        scan_term.append([atom0.No, a1, a2, atom3.No])
    molecule.torsions = scan_term
    if idx is not None:
        return molecule,idx
        #return scan_term

def determine_hydrogen_bond(a1, a2, a3):
    distance = calc_stru_para([a1, a2])
    angle = calc_stru_para([a1, a2, a3])
    if distance <= HBOND_DIST_MAX:
        if angle >= HBOND_ANGLE_MIN:
            return "yes"
        else:
            return "no"
    else:
        return "no"

def assign_hybrid(molecule,idx=None):
    _label = {
            "C":{0:"s",1:"sp",2:"sp",3:"sp2",4:"sp3"},
            "H":{0:"s",1:"s",},
            "O":{0:"s",1:"sp2",2:"sp3"},
            "N":{0:"s",1:"sp",2:"sp2",3:"sp3",4:"sp3+",5:"sp3d2"},
            "S":{0:"s",1:"sp2",2:"sp3",3:"sp2d",4:"sp3d2"},
            "P":{0:"s",1:"sp",2:"sp2",3:"sp3",4:"sp3+",5:"sp3d2"},
            "F":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            "Cl":{0:"s",1:"s",2:"sp3d",3:"sp3t",4:"sp3q"},
            "Br":{0:"s",1:"s",2:"sp3d",3:"sp3t",4:"sp3q"},
            "Si":{0:"s",1:"s",2:"sp",3:"sp2",4:"sp3"},
            "B":{0:"s",1:"s",2:"sp2",3:"sp3",4:"sp3-"},
            "I":{0:"s",1:"s",2:"sp3d",3:"sp3t",4:"sp3q",5:"sp3f"},
            "Na":{-1:"s",0:"s",1:"s",2:"s",3:"s",4:"s"},
            "K":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            "Ca":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            "Mg":{0:"s",1:"s",2:"s",3:"s",4:"s"},
            }
    for atom in molecule.Atoms:
        bond_number = len(atom.connectivity)
        bond_number -= atom.formal_charge

        if atom.element in ["N","P"]:
            if len(atom.connectivity) == 4:
                bond_number = 4
            else:
                o_num = sum([int(bn) - molecule.Atoms[atom.connectivity[jj]].formal_charge 
                             for jj,bn in enumerate(atom.bond_type) 
                             if molecule.Atoms[atom.connectivity[jj]].element == "O"])
                if o_num >= 2:
                    bond_number = 5
                    #bond_number += o_num 
                    #if bond_number >= 4:
                    #    bond_number = 5
        atom.hybrid = _label[atom.element][bond_number]
    if idx is not None:
        return molecule,idx
    else:
        return molecule
    
def _update_molecule_topol_value(molecule,idx=None):
    molecule.update_topol_value()
    if idx is not None:
        return molecule, idx
    else:
        return molecule

def old_assign_hybrid(molecule,idx=None):
    for atom in molecule.Atoms:
        if atom.elem in ["H","F","Cl","Br","I"]:
            if len(atom.connectivity) == 1:
                atom.hybrid = "s"
            else:
                atom.hybrid = "sp3"
        elif atom.elem in ["O","C"]:
            if "3" in atom.bond_type:
                atom.hybrid = "sp"
            elif "2" in atom.bond_type:
                atom.hybrid = "sp2"
            else:
                atom.hybrid = "sp3"
        elif atom.elem in ["N"]:
            if "3" in atom.bond_type:
                atom.hybrid = "sp"
            elif "2" in atom.bond_type:
                atom.hybrid = "sp2"
            else:
                atom.hybrid = "sp3"
        elif atom.elem in ["S"]:
            if len(atom.connectivity) - atom.formal_charge <= 2:
                if "2" in atom.bond_type:
                    atom.hybrid = "sp2"
                else:
                    atom.hybrid = "sp3"
            else:
                atom.hybrid = "sp3"
        elif atom.elem in ["P"]:
            if len(atom.connectivity) - atom.formal_charge <= 3:
                if "2" in atom.bond_type:
                    atom.hybrid = "sp2"
                else:
                    atom.hybrid = "sp3"
            else:
                atom.hybrid = "sp3"
    if idx is not None:
        return molecule,idx

def protein_ring_and_charge_group(protein):
    residues = {}
    for atom in protein.Atoms:
        res_name = f"{atom.residue}_{atom.residue_ID}_{getattr(atom,'chain_name','A')}"
        if res_name not in residues:
            residues[res_name] = []
        residues[res_name].append(atom.ID)

    Groups = [Group(group_str=kk,
                    group_name=kk.split("_")[0],
                    #group_idx=int(kk.split("_")[1]),
                    group_idx=kk.split("_")[1],
                    group_chain_name=kk.split("_")[2],
                    atoms=[an for an in vv],
                    net_charge=sum([protein.Atoms[an].formal_charge for an in vv])
                    ) 
            for kk,vv in residues.items()]

    ring_label = {
        "TYR":["CG","CD1","CD2","CE1","CE2","CZ"],
        "TRP":["CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"],
        "HIS":["CG","ND1","CD2","CE1","NE2",],
        "HID":["CG","ND1","CD2","CE1","NE2",],
        "HIE":["CG","ND1","CD2","CE1","NE2",],
        "PHE":["CG","CD1","CD2","CE1","CE2","CZ"],
        "DNA_nonar":["C1'","C2'","C3'","C4'","O4'"],
        "DNA_ar1":["N1","C2","N3","C4","C5","C6"],
        "DNA_ar2":["C4","C5","N7","C8","N9"],
        }
    charge_label={
        "ARG":["NE","CZ","NH1","NH2","HE","HH11","HH12","HH21","HH22"],
        "HIS":["ND1","NE2","HD1","HE2",],
        "HIP":["ND1","NE2","HD1","HE2",],
        "LYS":["NZ","HZ1","HZ2","HZ3"],
        "GLU":["CD","OE1","OE2"],
        "ASP":["CG","OD1","OD2"],
    }

    rings = {}
    for res in Groups:
        if res.group_name =="TRP":
            rings[f"{res.group_str}_5"] = [an for an in res.atoms if protein.Atoms[an].name in ["CG","CD1","CD2","NE1","CE2"]] + ["ar2"]
            rings[f"{res.group_str}_6"] = [an for an in res.atoms if protein.Atoms[an].name in ["CD2","CE2","CE3","CZ2","CZ3","CH2"]] + ["ar1"]
        elif res.group_name in ["TYR", "PHE",]:
            rings[res.group_str] = [an for an in res.atoms if protein.Atoms[an].name in ring_label[res.group_name]] + ["ar1"]
        elif res.group_name in ["HIS","HID","HIE"]:
            rings[res.group_str] = [an for an in res.atoms if protein.Atoms[an].name in ring_label[res.group_name]] + ["ar2"]
        elif res.group_name in ["A","G","DA","DG","A3","G3","DA3","DG3","A5","G5","DA5","DG5","AN","GN","DAN","DGN"]:
            rings[f"{res.group_str}_nonar"] = [an for an in res.atoms if protein.Atoms[an].name in ring_label["DNA_nonar"]] + ["nonar"]
            rings[f"{res.group_str}_5"] = [an for an in res.atoms if protein.Atoms[an].name in ring_label["DNA_ar2"]] + ["ar2"]
            rings[f"{res.group_str}_6"] = [an for an in res.atoms if protein.Atoms[an].name in ring_label["DNA_ar1"]] + ["ar1"]
        elif res.group_name in ["U","C","DT","DC","U3","C3","DT3","DC3","U5","C5","DT5","DC5","UN","CN","DTN","DCN"]:
            rings[f"{res.group_str}_nonar"] = [an for an in res.atoms if protein.Atoms[an].name in ring_label["DNA_nonar"]] + ["nonar"]
            rings[f"{res.group_str}_6"] = [an for an in res.atoms if protein.Atoms[an].name in ring_label["DNA_ar1"]] + ["ar1"]    
            
    charge_group = {}
    for res in Groups:
        if res.group_name in ["ARG", "HIS", "LYS", "GLU", "ASP","HIP"]:
            net_charge = sum([protein.Atoms[an].point_charge for an in res.atoms])
            if abs(net_charge) > 0.001:
                if net_charge > 0.0:
                    label = "positive"
                else:
                    label = "negative"
                charge_group[res.group_str] = [an for an in res.atoms if protein.Atoms[an].name in charge_label[res.group_name]] + [label]
        if res.group_name in ["A3","U3","G3","C3","DA3","DT3","DG3","DC3","A","U","G","C","DA","DT","DG","DC"]:
            charge_group[res.group_str] = [an for an in res.atoms if protein.Atoms[an].name in ["P","O1P","O2P"]] + ["negative"]
        
            
    protein.rings = rings
    protein.charge_group = charge_group
    protein.Groups = Groups
    get_aromatic_connect(protein)

    return protein

def get_fragment_inchi_key(molecule,ans,ignore_double_bond=False):
    from ...chem.molecule import Molecule
    atom_id = {}
    for ii in range(len(ans)):
        atom_id[ans[ii]] = ii
    nn = len(atom_id)

    fm = Molecule("frag")
    fm.create_atoms(nn)
    # fm.create_atoms(len(frag["components"]))
    if not ignore_double_bond:
        for a, b in atom_id.items():
            fm.Atoms[b].elem = molecule.Atoms[a].elem
            fm.Atoms[b].formal_charge = molecule.Atoms[a].formal_charge
            fm.Atoms[b].coor = molecule.Atoms[a].coor
            fm.Atoms[b].connect = []
            fm.Atoms[b].bond_type = []
            for jj,an in enumerate(molecule.Atoms[a].connect):
                if an in atom_id:
                    fm.Atoms[b].connect.append(atom_id[an])
                    fm.Atoms[b].bond_type.append(molecule.Atoms[a].bond_type[jj])
                else:
                    fm.add_atom()
                    fm.Atoms[nn].elem = "H"
                    fm.Atoms[nn].formal_charge = 0
                    fm.Atoms[nn].coor = molecule.Atoms[an].coor
                    fm.Atoms[nn].connect = [b]
                    fm.Atoms[nn].bond_type = ["1"]
                    fm.Atoms[b].connect.append(nn)
                    fm.Atoms[b].bond_type.append("1")
                    nn += 1
        fm.create_topols()
        return fm.inchi_key,fm
    else:
        for a,b in atom_id.items():
            fm.Atoms[b].elem = molecule.Atoms[a].elem
            fm.Atoms[b].formal_charge = molecule.Atoms[a].formal_charge
            fm.Atoms[b].coor = molecule.Atoms[a].coor
            fm.Atoms[b].connect = []
            fm.Atoms[b].bond_type = []
            for jj,an in enumerate(molecule.Atoms[a].connect):
                if an in atom_id:
                    fm.Atoms[b].connect.append(atom_id[an])
                    if molecule.Atoms[a].bond_type[jj] in ["2","3"]:
                        fm.Atoms[b].bond_type.append("1")
                        fm.add_atom()
                        fm.Atoms[nn].elem = "H"
                        fm.Atoms[nn].formal_charge = 0
                        fm.Atoms[nn].coor = molecule.Atoms[an].coor
                        fm.Atoms[nn].connect = [b]
                        fm.Atoms[nn].bond_type = ["1"]
                        fm.Atoms[b].connect.append(nn)
                        fm.Atoms[b].bond_type.append("1")
                        nn += 1
                    else:
                        fm.Atoms[b].bond_type.append(molecule.Atoms[a].bond_type[jj])
                else:
                    fm.add_atom()
                    fm.Atoms[nn].elem = "H"
                    fm.Atoms[nn].formal_charge = 0
                    fm.Atoms[nn].coor = molecule.Atoms[an].coor
                    fm.Atoms[nn].connect = [b]
                    fm.Atoms[nn].bond_type = ["1"]
                    fm.Atoms[b].connect.append(nn)
                    fm.Atoms[b].bond_type.append("1")
                    nn += 1
        fm.create_topols()
        return fm.inchi_key,fm

def ring_atom_order(molecule,ring_block):
    for rname in ring_block:
        pass

def check_chiral_atom(molecule,idx=None):
    for atom in molecule.Atoms:
        chirality_flag = False
        if atom.elem == "C" and len(atom.connect) == 4:
            if len(atom.has_ring) == 0:
                side_fragments = [get_fragment_inchi_key(molecule,[ani for ani in molecule.find_side_componend(an,atom.ID)+[an]])[1] 
                                  for an in atom.connect]
                if len(set([frag.inchi_key for frag in side_fragments])) == 4:
                    chirality_flag = True
            else:
                side_fragments = [get_fragment_inchi_key(molecule,[ani for ani in molecule.find_side_componend(an,atom.ID)+[an]])[1] 
                                  for an in atom.connect if len(molecule.Atoms[an].has_ring) == 0]
                if len(set([frag.inchi_key for frag in side_fragments])) == 2:
                    chirality_flag = True
        
        atom.chirality_flag = chirality_flag

    if idx is not None:
        return molecule, idx
    else:
        return molecule


