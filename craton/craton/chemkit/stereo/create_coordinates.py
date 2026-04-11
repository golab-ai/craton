from copy import  deepcopy
import numpy as np
from ...chem.elements import get_bonded_type_distance
from ...utils.geometry import *
#tetrahedron, triangle, line, change_bond, calculate_dihedral, rotation_dihedral,change_angle,line_to_line_angle

#molecule = deepcopy(molecule_target)
#connects = [atom.connect for atom in molecule.Atoms if atom.elem not in ["H"]]

hybrid = {"C4": tetrahedron,"N4": tetrahedron,"N3": triangle,"O2": tetrahedron,
          "P4": tetrahedron,"P3": tetrahedron,"P2": tetrahedron,
          "S4": tetrahedron,"S2": tetrahedron,
          "C3":triangle,"N2":triangle,"O1":triangle,
          "C2":line,"N1":line,
            "H1": "end",
            "Cl1":line,
            "F1":line,
            "Br1":line,
            "I1":line,
         }

plate_type = ["C3","N2"]
line_type = ["C2","N1"]

class Coordinate:
    def __init__(self,molecule,unatoms=None):
        self.molecule = molecule
        if unatoms is None:
            self.molecule.Atoms[0].coordinates = [0.0000,0.0000,0.0000]
            self.unatoms = [n for n in range(1,len(self.molecule.Atoms))]
        else:
            self.unatoms = unatoms
        

    def create_coordinates(self):
        pass

    def run_create_coordinates(self,center):
        for an in self.molecule.Atoms[center].connectivity:
            pass

def create_coordinates(center,exist_vertexs,predict_atoms):
    distance = get_bonded_type_distance(center.element,"H","1")
    if not hasattr(center,"hybrid"):
        center.hybrid = f"{center.elem}{len(center.connect)}"
    #this_hybrid = f"{center.elem}{len(center.connect)}"
    run_type = hybrid[center.hybrid]
    #run_type = hybrid[f"{center.elem}{len(center.connect)}"]
    if exist_vertexs is None:
        #vertexs = center_(center.coor,distance)
        new_coors = []
    else:
        #new_coors = [change_bond(center.coor,atom.coor,atom.coor,
        #             np.sqrt(sum([(center.coor[ii] - atom.coor[ii])**2 for ii in range(3)]))-distance) 
        #             for atom in exist_vertexs]
        new_coors = [change_bond(center.coor,atom.coor,atom.coor,
                     distance - np.sqrt(sum([(center.coor[ii] - atom.coor[ii])**2 for ii in range(3)]))) 
                     for atom in exist_vertexs]
    vertexs = run_type(center.coor,distance,new_coors)
    
    for jj,conn in enumerate(predict_atoms):
        
        dis = get_bonded_type_distance(center.elem,conn.elem,center.bond_type[center.connectivity.index(conn.ID)])
        #conn.coordinates = change_bond(center.coor,vertexs[jj],vertexs[jj],distance  - dis )
        conn.coordinates = change_bond(center.coor,vertexs[jj],vertexs[jj],dis - distance )
    
def old_chain_conformation(molecule,center,predict_atoms,exist_vertexs):
    if exist_vertexs is not None:
        arr = get_group_size_order(molecule,center.ID, delete_special=[atom.ID for atom in [an[0] for an in exist_vertexs]])
        tmp1 = [atom for atom in predict_atoms if atom.ID == arr[-1][0]][0]
        for rr in exist_vertexs:
            bond_type1 = f"{center.element}{len(center.connect)}"
            bond_type2 = f"{rr[0].element}{len(rr[0].connect)}"
            if bond_type1 in plate_type and bond_type2 in plate_type:
                tha = calculate_dihedral(molecule.Atoms[rr[1][0]].coor,rr[0].coor,center.coor,tmp1.coor)
                deltha = 180 - tha
                for conn in predict_atoms:
                    conn.coordinates = rotation_dihedral(rr[0].coor,center.coor,conn.coordinates,deltha)
                return
        
        tmp = exist_vertexs[0]
        
        tha = calculate_dihedral(molecule.Atoms[tmp[1][0]].coor,tmp[0].coor,center.coor,tmp1.coor)
        deltha = 180 - tha
        for conn in predict_atoms:
            conn.coordinates = rotation_dihedral(tmp[0].coor,center.coor,conn.coordinates,deltha)            
        
        
def get_group_size_order(molecule,kk,delete_special=None):
    arr = []
    if delete_special is None:
        delete_special = [-1]
    for an in molecule.Atoms[kk].connect:
        if an not in delete_special:
            tmp = 0 if molecule.Atoms[an].element in ["H"] else 1
            tmp += len([1 for aa in molecule.Atoms[an].connect if molecule.Atoms[aa].element not in ["H"]])
            arr.append([an,tmp])
    arr = sorted(arr,key=lambda x:x[1])
    return arr

def refine_conformation(molecule,center_ii,predict_vertexs,exist_vertexs):
    _connect_ = []
    for kk in exist_vertexs:
        arr = get_group_size_order(molecule,kk,delete_special = [center_ii])
        _connect_.append(get_group_size_order(molecule,kk,delete_special = [center_ii])[-1])

    center = molecule.Atoms[center_ii]
    predict_atoms = [molecule.Atoms[kk] for kk in predict_vertexs]
    exist_atoms = [[molecule.Atoms[kk],_connect_[ii]] for ii,kk in enumerate(exist_vertexs)]

    if exist_atoms is not None:
        arr = get_group_size_order(molecule,center.ID, delete_special=[atom.ID for atom in [an[0] for an in exist_atoms]])
        tmp1 = [atom for atom in predict_atoms if atom.ID == arr[-1][0]][0]
        for rr in exist_atoms:
            bond_type1 = f"{center.element}{len(center.connect)}"
            bond_type2 = f"{rr[0].element}{len(rr[0].connect)}"
            if bond_type1 in plate_type and bond_type2 in plate_type:
                tha = calculate_dihedral(molecule.Atoms[rr[1][0]].coor,rr[0].coor,center.coor,tmp1.coor)
                deltha = 180 - tha
                for conn in predict_atoms:
                    conn.coordinates = rotation_dihedral(rr[0].coor,center.coor,conn.coordinates,deltha)
                return
        
        tmp = exist_atoms[0]
        
        tha = calculate_dihedral(molecule.Atoms[tmp[1][0]].coor,tmp[0].coor,center.coor,tmp1.coor)
        deltha = 180 - tha
        for conn in predict_atoms:
            conn.coordinates = rotation_dihedral(tmp[0].coor,center.coor,conn.coordinates,deltha)

def old_ring_conformation(molecule,center_ii,atoms_list,used,predict_vertexs,ring_angle):
    predict_ring_atoms =[molecule.Atoms[kk] for kk in list(set(atoms_list).intersection(set(predict_vertexs)))]
    exist_ring_atoms = [molecule.Atoms[kk] for kk in list(set(atoms_list).intersection(set(used)))]
    center = molecule.Atoms[center_ii]

    atoms = [molecule.Atoms[kk] for kk in atoms_list]

    if len(predict_ring_atoms) > 0:
        this_atom = predict_ring_atoms[0]
        this_ii = atoms.index(this_atom)

        angle_atom = atoms[this_ii - 2]
        delangle = line_to_line_angle(this_atom.coor,center.coor,angle_atom.coor,center.coor) - ring_angle
        move_atoms = list(set(molecule.find_side_componend(this_atom.ID, center.ID)).intersection(set(used)))
        for an in [this_atom.ID] + move_atoms:
            molecule.Atoms[an].coor = change_angle(angle_atom.coor,center.coor,this_atom.coor,molecule.Atoms[an].coor,-delangle)
        
        if len(exist_ring_atoms) >= 3:
            tha = calculate_dihedral(atoms[this_ii - 3].coor,atoms[this_ii - 2].coor,atoms[this_ii - 1].coor,this_atom.coor)
            deltha = 180 - tha
            for an in [this_atom.ID] + move_atoms:
                molecule.Atoms[an].coor = rotation_dihedral(atoms[this_ii - 2].coor,atoms[this_ii - 1].coor,molecule.Atoms[an].coor,deltha)

def ring_conformation(molecule,atoms_list,ring_exist_atoms,pre_used):
    if len(atoms_list) <= 8:
        ring_angle = 180.0 - 360.0/len(atoms_list)
    #if len(ring_exist_atoms) > 0:


    #for ii in range(1,len(atoms_list)-1):
    for ii in range(1,len(atoms_list)-1):
        if len(atoms_list) > 3:
            a3 = molecule.Atoms[atoms_list[ii - 3]]
        a2 = molecule.Atoms[atoms_list[ii - 2]]
        a1 = molecule.Atoms[atoms_list[ii - 1]]
        a0 = molecule.Atoms[atoms_list[ii]]
        
        delangle = line_to_line_angle(a0.coor,a1.coor,a2.coor,a1.coor) - ring_angle
        move_atoms = [a0.ID] + molecule.find_side_componend(a0.ID, a1.ID)
        tmp_coor = deepcopy(a0.coor)
        if atoms_list == [3, 4, 6, 7, 8, 9,] or a0.ID in [3,4]:
            for an in move_atoms:
                molecule.Atoms[an].coor = change_angle(a2.coor,a1.coor,tmp_coor,molecule.Atoms[an].coor,-1 * delangle)
            if  len(atoms_list) > 3 and ii > 1:
                tha = calculate_dihedral(a3.coor,a2.coor,a1.coor,a0.coor)
                deltha = 0 - tha
                move_atoms += [kk for kk in a1.connect if kk not in [a3.ID,a2.ID,a1.ID,a0.ID]]
                #move_atoms = list(set(move_atoms).difference(set(pre_used)).difference(set(atoms_list)))
                for an in move_atoms:
                    molecule.Atoms[an].coor = rotation_dihedral(a2.coor,a1.coor,molecule.Atoms[an].coor,deltha)
    if atoms_list == [3, 4, 6, 7, 8, 9,]:
        ring_connects_conformation_refine(molecule,atoms_list,ring_exist_atoms)

def ring_connects_conformation_refine(molecule,atoms_list,ring_exist_atoms):
    coors = [molecule.Atoms[ii].coor for ii in atoms_list]
    cog = find_center(coors, None, center_type="cog")
    v = cross_product_plane(coors[0],coors[1],coors[2])

    cog_v = [(cog[ii] - v[ii]) for ii in range(3)]
    print("atoms_list:", atoms_list)
    print("ring_exist_atoms:", ring_exist_atoms)
    for ii in range(1,len(atoms_list)):
        if atoms_list[ii] not in ring_exist_atoms:
            end_atom = molecule.Atoms[atoms_list[ii]]
            res_atoms = [molecule.Atoms[kk] for kk in end_atom.connect if kk not in atoms_list]
            if len(res_atoms) > 0:
                move_atoms = []
                for res in res_atoms:
                    move_atoms.append(res.ID)
                    move_atoms.extend(molecule.find_side_componend(res_atoms[0].ID, end_atom.ID))

                print("ii,move_atoms:",atoms_list[ii],move_atoms)
                end_v = [(end_atom.coor[ii] - v[ii]) for ii in range(3)]
                tha_end = calculate_dihedral(cog_v,end_v,end_atom.coor,res_atoms[0].coor)
                if end_atom.hybrid in plate_type:
                    deltha_end = 180 - tha_end
                    for an in move_atoms:
                        molecule.Atoms[an].coor = rotation_dihedral(end_v,end_atom.coor,molecule.Atoms[an].coor,deltha_end) 
                elif end_atom.hybrid not in line_type:
                    ra = rotation_dihedral(cog_v,cog,end_atom.coor,90.0)
                    rb = rotation_dihedral(cog_v,cog,end_atom.coor,45.0)
                    tha_cog = calculate_dihedral(ra,cog,end_atom.coor,res_atoms[0].coor)
                    tha_r = calculate_dihedral(ra,rb,end_atom.coor,res_atoms[0].coor)
                    deltha_cog = 90 - tha_cog
                    deltha_r = 120 - tha_r
                    deltha_end = 180 - tha_end
            
                    for an in move_atoms:
                        molecule.Atoms[an].coor = rotation_dihedral(end_v,end_atom.coor,molecule.Atoms[an].coor,deltha_end)
                    for an in move_atoms:
                        molecule.Atoms[an].coor = rotation_dihedral(cog,end_atom.coor,molecule.Atoms[an].coor,deltha_cog)
                    for an in move_atoms:
                        molecule.Atoms[an].coor = rotation_dihedral(rb,end_atom.coor,molecule.Atoms[an].coor,deltha_r)


def old_ring_connects_conformation_refine(molecule,atoms_list,ring_exist_atoms):

    #for ii in range(1,len(atoms_list)-1):
    for ii in range(1,len(atoms_list)):
        if atoms_list[ii] not in ring_exist_atoms:
            end_atom = molecule.Atoms[atoms_list[ii]]
            start_atom = molecule.Atoms[atoms_list[ii - 1]]
            next_atom = molecule.Atoms[atoms_list[ii - 2]]
            res_atoms = [molecule.Atoms[kk] for kk in end_atom.connect if kk not in atoms_list]
            if len(res_atoms) > 0:
                move_atoms = []
                for res in res_atoms:
                    move_atoms.append(res.ID)
                    move_atoms.extend(molecule.find_side_componend(res_atoms[0].ID, end_atom.ID))

                if end_atom.hybrid in plate_type:
                    tha = calculate_dihedral(next_atom.coor,start_atom.coor,end_atom.coor,res_atoms[0].coor)
                    deltha = 180 - tha if tha >= 0.0 else -180 - tha
                    for an in move_atoms:
                        molecule.Atoms[an].coor = rotation_dihedral(start_atom.coor,end_atom.coor,molecule.Atoms[an].coor,deltha)
                elif end_atom.hybrid not in line_type:
                    tha = calculate_dihedral(next_atom.coor,start_atom.coor,end_atom.coor,res_atoms[0].coor)
                    tha1 = calculate_dihedral(next_atom.coor,start_atom.coor,end_atom.coor,res_atoms[1].coor)
                    if abs(180 - abs(tha)) < 2.0:
                        tha = 180.0
                    if abs(180 - abs(tha1)) < 2.0:
                        tha1 = 180.0
                    if tha * tha1 <= 0.0:
                        deltha = 108 - tha if tha >= 0.0 else -108 - tha
                    else:
                        deltha = -108 - tha if tha >= 0.0 else 108 - tha
                    print(next_atom.ID,start_atom.ID,end_atom.ID,[res.ID for res in res_atoms],tha,tha1,deltha)
                    for an in move_atoms:
                        molecule.Atoms[an].coor = rotation_dihedral(start_atom.coor,end_atom.coor,molecule.Atoms[an].coor,deltha)
                    #node_atom = [molecule.Atoms[kk] for kk in start_atom.connect if kk not in atoms_list][0]
                    #beta = calculate_dihedral(node_atom.coor,start_atom.coor,end_atom.coor,res_atoms[0].coor)
                    #if abs(beta) <= 120.0:
                
def run_create(mtype,atoms,molecule,used):
    
    if mtype in ["ring","R"]:
        molecule.Atoms[atoms[-2]].hybrid = f"{molecule.Atoms[atoms[-2]].elem}{len(molecule.Atoms[atoms[-2]].connect)}"
        molecule.Atoms[atoms[-1]].hybrid = f"{molecule.Atoms[atoms[-1]].elem}{len(molecule.Atoms[atoms[-1]].connect)}"
        rs_connectivity = deepcopy(molecule.Atoms[atoms[-2]].connectivity)
        rs_bond_type = deepcopy(molecule.Atoms[atoms[-2]].bond_type)
        re_connectivity = deepcopy(molecule.Atoms[atoms[-1]].connectivity)
        re_bond_type = deepcopy(molecule.Atoms[atoms[-1]].bond_type)
        inds = molecule.Atoms[atoms[-2]].connectivity.index(atoms[-1])
        inde = molecule.Atoms[atoms[-1]].connectivity.index(atoms[-2])
        del molecule.Atoms[atoms[-2]].connectivity[inds]
        del molecule.Atoms[atoms[-2]].bond_type[inds]
        del molecule.Atoms[atoms[-1]].connectivity[inde]
        del molecule.Atoms[atoms[-1]].bond_type[inde]
        ring_exist_atoms = list(set(atoms).intersection(set(used)))
        pre_used = deepcopy(used)

    for ii in atoms:
        rr = molecule.Atoms[ii].connectivity
        if ii not in used:
            create_coordinates(
                                molecule.Atoms[ii],
                               None,
                               [molecule.Atoms[ii] for ii in molecule.Atoms[ii].connect],
                               )
            used.append(ii)
            used.extend(molecule.Atoms[ii].connect)
        else:
            exist_vertexs = list(set(molecule.Atoms[ii].connect).intersection(set(used)))
            predict_vertexs = list(set(molecule.Atoms[ii].connect).difference(set(used)))
            if len(predict_vertexs) != 0:
                if len(exist_vertexs) == 0:
                    create_coordinates(
                                        molecule.Atoms[ii],
                                       None,
                                       [molecule.Atoms[ii] for ii in predict_vertexs],
                                       )
                else:
                    create_coordinates(
                                       molecule.Atoms[ii],
                                       [molecule.Atoms[kk] for kk in exist_vertexs],
                                       [molecule.Atoms[kk] for kk in predict_vertexs],
                                       )
                    if mtype in ["chain","C","Chain"]:
                        if f"{molecule.Atoms[ii].elem}{len(molecule.Atoms[ii].connect)}" not in line_type:
                            refine_conformation(molecule,ii,predict_vertexs,exist_vertexs)
            used.append(ii)
            used.extend(predict_vertexs)
    if mtype in ["ring","R"]:
        #if atoms == [3, 4, 6, 7, 8, 9,]:
        ring_conformation(molecule,atoms,ring_exist_atoms,pre_used)
        molecule.Atoms[atoms[-2]].connectivity = rs_connectivity
        molecule.Atoms[atoms[-2]].bond_type = rs_bond_type
        molecule.Atoms[atoms[-1]].connectivity = re_connectivity
        molecule.Atoms[atoms[-1]].bond_type = re_bond_type
    return used


def creaate_ring(ring_atoms,molecule,used):
    for ii in ring_atoms:
        rr = molecule.Atoms.connectivity
        if ii not in used:
            create_coordinates(molecule,
                                molecule.Atoms[ii],
                               None,
                               [molecule.Atoms[ii] for ii in molecule.Atoms[ii].connect],
                               )
            used.append(ii)
            used.extend(molecule.Atoms[ii].connect)

def tri_ring(center,distance,vertexs,ring_atoms):
    pass
        
def tetra_ring():
    pass

def penta_ring():
    pass

def hexa_ring():
    pass

def hepta_ring():
    pass

def octa_ring():
    pass