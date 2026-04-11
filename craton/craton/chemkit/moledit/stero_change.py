from ...utils.geometry import *
from ...utils import logger

def calculate_structure(molecule,patoms):
    return calc_stru_para([molecule.Atoms[an].coordinates for an in patoms])

def change_structure(molecule,patoms,value,del_value=False,improper_flag=False):
    """
    依据某一拓扑结构，将分子的结构改变到某一数值
    拓扑由patoms决定：
        patoms是由于原子编号组成的数组，如[1,2],[5,8,12],[2,7,9,12]。
        改变方向为数组后面原子，及相连原子的的坐标
        len(patoms) == 2：沿某一键，拉伸距离
        len(patoms) == 3：沿某一角，弯曲角度
        len(patoms) == 4：沿某一二面角，旋转角度
    value 是改变的数值，距离单位A, 角度单位度
    del_value=Ture，表示改变的数值，False表示改变到的数值
    """
    _func = {
            2:change_bond,
            3:change_angle,
            4:change_dihedral,
             }
    if not del_value:
        curr_value = calculate_structure(molecule,patoms)
        del_value = value - curr_value
    else:
        del_value = value
    if len(patoms) == 2:
        change_atoms = [patoms[1]]
        change_atoms.extend(molecule.find_side_componend(patoms[1],patoms[0]))
    elif len(patoms) == 3:
        change_atoms = molecule.find_side_componend(patoms[1],patoms[0])
    elif len(patoms) == 4:
        if not improper_flag:
            change_atoms = molecule.find_side_componend(patoms[2],patoms[1])
        else:
            change_atoms = molecule.find_side_componend(patoms[3],patoms[2])
            change_atoms = [patoms[3]] + change_atoms
    else:
        logger.warning(f"lenght ({len(patoms)}) of patoms error: the length of patoms must be equal to 2,3,4")
        return
    
    for catom in change_atoms:
        molecule.Atoms[catom].coordinates = _func[len(patoms)](*[molecule.Atoms[an].coordinates for an in patoms],
                                                               molecule.Atoms[catom].coordinates,
                                                               del_value)
    return molecule