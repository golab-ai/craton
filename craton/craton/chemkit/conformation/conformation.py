import hashlib
import itertools
import numpy as np
from copy import deepcopy

from ..moledit import MolEdit as ME 
from ...chem.format._rdkit import RdkitMol
from ...utils import logger
from ...property import MolProperty as MP

__all__ = ["ConformType", "Conformation", "Scan"]

def get_conformation_RMSD(m1, m2):
    rdk = RdkitMol()
    return rdk._get_RMSD(m1.mol_script,m2.mol_script)

def is_in_same_ring(m, *ids):
    for k, v in m.ring_dict.items():
        if all(i in v[:-1] for i in ids):
            return True
            
def get_bond_angle_scan_term(molecule,ignore_ring=True,exists_type=None,inter_val=[-0.1,0.2,-5.0,5.0],add_n=1):
    scan_term = []
    for term in ["Bonds","Angles"]:
        items = getattr(molecule,term,[])
        for bb in items:
            ans = [bb.a1+add_n,bb.a2+add_n]
            if hasattr(bb,"a3"):
                ans.append(bb.a3+add_n)
            igr_flag = True
            if ignore_ring:
                if is_in_same_ring(molecule,ans[:-1]):
                    ignore_ring = False
            _flag = False
            if igr_flag:
                if exists_type is not None:
                    if len(set(bb.atom_type_names).intersection(set(exists_type))) == 0:
                        _flag = True
                        exists_type.append(bb.atom_type_name)
                else:
                    _flag = True
            if _flag:
                if term == "Bonds":
                    if inter_val[0] != 0.0:
                        nn = int(-1 * (bb.value - 0.5) / inter_val[0])
                        forward_n = nn if nn < 10 else 10
                        scan_term.append([ans,[forward_n,inter_val[0]]])
                    if inter_val[1] != 0.0:
                        scan_term.append([ans,[10,inter_val[1]]])
                elif term == "Angles":
                    if inter_val[2] != 0.0:
                        scan_term.append([ans,[6,inter_val[2]]])
                    if inter_val[3] != 0.0:
                        scan_term.append([ans,[6,inter_val[3]]])
    return scan_term
            

def bond_angle_extend_conformer(molecule,ignore_alkane=True):
    from ..structure.zmatrix import ZMatrix
    from ...chem.elements import Element
    import math,copy

    _alkane_types = ["h_1", "c_4", "c_4h", "c_4h2", "c_4h3", "c_4h4"]

    def is_atoms_alkane(m, *ids):
        return all(m.Atoms[i].atom_type_name in _alkane_types for i in ids)

    def is_in_same_ring(m, *ids):
        for k, v in m.ring_dict.items():
            if all(i in v[:-1] for i in ids):
                return True

        return False

    def determine_delta_bond(m, bond):
        if is_in_same_ring(m, bond.a1, bond.a2):
            return 0.05
        elif bond.type == "1":
            return 0.08
        elif bond.type == "3":
            return 0.05
        else:
            return 0.06

    def validate_pair_distance(mole, pos_array):
        for pair in getattr(mole, "Pair14", []) + getattr(mole, "Pair1n", []):
            elem1, elem2 = mole.Atoms[pair.a1].elem, mole.Atoms[pair.a2].elem
            vdw_diameter = Element.get(elem1).vdw_radius + Element.get(elem2).vdw_radius
            threshold = vdw_diameter / 2 ** (1 / 6)
            delta = pos_array[pair.a2] - pos_array[pair.a1]
            if delta.dot(delta) < threshold**2 and not is_in_same_ring(mole, pair.a1, pair.a2):
                return False

        return True
            

    zmat = ZMatrix()
    zmat.from_molecule(molecule)

    id_bonds = [
        i
        for i, b in enumerate(zmat.bonds)
        if b is not None and not (ignore_alkane and is_atoms_alkane(molecule, b.a1, b.a2))
    ]
    id_angles = [
        i
        for i, a in enumerate(zmat.angles)
        if a is not None
        and (a.type21 == "1" or a.type23 == "1")
        and not is_in_same_ring(molecule, a.a1, a.a3)
        and not (ignore_alkane and is_atoms_alkane(molecule, a.a1, a.a2, a.a3))
    ]

    _n_total = 0
    zmatrixes = {}

    half_len = math.ceil(len(id_bonds) / 2)
    for i, j in zip(id_bonds[:half_len], (id_bonds + [None])[half_len : 2 * half_len]):
        _n_total += 1
        zmatb = copy.deepcopy(zmat)
        zmatb.bonds[i].value += determine_delta_bond(molecule, zmatb.bonds[i])
        name = "b-%i" % (i + 1)
        if j is not None:
            zmatb.bonds[j].value -= determine_delta_bond(molecule, zmatb.bonds[j])
            name += "-%i" % (j + 1)
        if validate_pair_distance(molecule, zmatb.to_coordinates()):
            zmatrixes[name] = zmatb

    half_len = math.ceil(len(id_angles) / 2)
    for i, j in zip(id_angles[:half_len], (id_angles + [None])[half_len : 2 * half_len]):
        _n_total += 1
        zmata = copy.deepcopy(zmat)
        angle_i = zmata.angles[i]
        angle_i.value += 10 / 180 * math.pi
        angle_i.value_a += 10
        name = "a-%i" % (i + 1)
        if j is not None:
            angle_j = zmata.angles[j]
            angle_j.value -= 10 / 180 * math.pi
            angle_j.value_a -= 10
            name += "-%i" % (j + 1)
        if validate_pair_distance(molecule, zmata.to_coordinates()):
            zmatrixes[name] = zmata

    return zmatrixes


class ConformType:
    LOCAL_MINIMUM = "local minimum"
    STRETCH = "stretch"
    SCAN = "scan"
    CONSTRAINED_LOCAL_MINIMUM = "constrained local minimum"
    BARRIER = "barrier"
    ON_SLOPE = "on slope"
    OTHER = "other"
    SINGLE_POINT = "single point"
    OPTIMIZING = "optimizing"
    TORSION_SCAN_TYPES = [SCAN, CONSTRAINED_LOCAL_MINIMUM, BARRIER, ON_SLOPE, OTHER]

def conformation_id_hash(molecule):
    """
    Hash the object based on inchi_key, coordinates and constrain_term,
    so that can be queried efficiently in mongodb

    TODO Consider the translation and orientation of coordinates
    """
    try:
        string = molecule.inchi_key
    except:
        string = molecule.formula
    for coor in molecule.coordinates:
        for i in coor:
            string += "%.4f" % i
    if molecule.constrain_term is not None:
        for i in molecule.constrain_term:
            string += str(i)
    return hashlib.md5(string.encode()).hexdigest()

def ignore_alkane_torsion(molecule):
    scan_term = []
    for torsion in molecule.torsions:
        ats = [molecule.Atoms[an].atom_type_name for an in torsion]
        if not set(ats).issubset({"c_4", "c_4h", "c_4h2", "c_4h3", "c_4h4", "h_1"}):
            scan_term.append(torsion)
    molecule.torsions = scan_term

def get_scan_curve(molecules: list):
    """
    分析出沿某一自由度变化的势能曲线
    """
    def _base_dict(data_dict_arr):
        scan_curve = {}
        for data_dict in data_dict_arr:
            if "scan_term" in data_dict:
                if data_dict["mole_name"] not in scan_curve:
                    scan_curve[data_dict["mole_name"]] = {}
                scan_term_name = "-".join([str(aa) for aa in data_dict["constrain"][0][0]])
                if scan_term_name not in scan_curve[data_dict["mole_name"]]:
                    scan_curve[data_dict["mole_name"]][scan_term_name] = []
                scan_curve[data_dict["mole_name"]][scan_term_name].append(data_dict)
        
        for m_name,vv in scan_curve.items():
            for torsion,vvv in vv.items():
                scan_curve[m_name][torsion] = sorted(vvv, key=lambda m:m["constrain"][0][-1])
        return scan_curve

    def _base_molecule(molecules):
        scan_curve = {}
        for molecule in molecules:
            if hasattr(molecule,"constrain"):
                if molecule.inchi_key not in scan_curve:
                    scan_curve[molecule.inchi_key] = {}

                if molecule.constrain[0].name not in scan_curve[molecule.inchi_key]:
                    scan_curve[molecule.inchi_key][molecule.constrain[0].name] = []

                scan_curve[molecule.inchi_key][molecule.constrain[0].name].append(molecule)

        for m_name,vv in scan_curve.items():
            for torsion,vvv in vv.items():
                #_nn = [[ii,molecule.constrain[0].fix_value] for ii,molecule in enumerate(vvv) if 180 - abs(molecule.constrain[0].fix_value) < 1.0]
                #if len(_nn) > 1:
                #    if _nn[0][1] * _nn[1][1] > 0:
                #        vvv[_nn[1][0]].constrain[0].fix_value = vvv[_nn[1][0]].constrain[0].fix_value * -1
                scan_curve[m_name][torsion] = sorted(vvv, key=lambda m:m.energy)
        return scan_curve


    if isinstance(molecules[0],dict):
        return _base_dict(molecules)
    else:
        return _base_molecule(molecules)

def assign_scan_conf_type(scan_curve):
    """
    依据势能曲线，判断每个构象的类型，如局部最小点，能垒等
    """
    def _base_molecule(molecules,scan_curve):

        def assign_conf_type_run(molecules):
            for ii, curr in enumerate(molecules):
                curr_e = curr.energy
                ahead_e = molecules[ii - 1].energy
                behind_e = molecules[(ii + 1) - len(molecules)].energy
                if curr_e < ahead_e and curr_e < behind_e:
                    curr.conform_type = ConformType.CONSTRAINED_LOCAL_MINIMUM
                elif curr_e > ahead_e and curr_e > behind_e:
                    curr.conform_type = ConformType.BARRIER
                else:
                    curr.conform_type = ConformType.ON_SLOPE

        for terms in scan_curve.values():
            for molecules in terms.values():
                assign_conf_type_run(molecules)
                molecules[-1].conform_type = ConformType.OTHER

    def _base_dict(dict_arr):
        def assign_conf_type_run(dict_arr):
            for ii,curr in enumerate(dict_arr):
                curr_e = curr["energy"]
                ahead_e = dict_arr[ii - 1]["energy"]
                behind_e = dict_arr[(ii + 1) - len(dict_arr)]["energy"]

                if curr_e < ahead_e and curr_e < behind_e:
                    curr["conform_type"] = ConformType.CONSTRAINED_LOCAL_MINIMUM
                elif curr_e > ahead_e and curr_e > behind_e:
                    curr["conform_type"] = ConformType.BARRIER
                else:
                    curr["conform_type"] = ConformType.ON_SLOPE
            return dict_arr

        for terms in scan_curve.values():
            for dict_arr in terms.values():
                dict_arr = assign_conf_type_run(dict_arr)
                dict_arr[-1]["conform_type"] = ConformType.OTHER 

    for m_name,vv in scan_curve.items():
        for torsion,vvv in vv.items():
            if isinstance(vvv[0],dict):
                return _base_dict(scan_curve)
                break
            else:
                return _base_molecule(scan_curve)
                break

def local_minimum_pes(scan_curve):
    """
    找到一个分子local_minimum构象
    输入：
        scan_curve: List[Atom.No] 查找local_minimum的scan term
    输出：
        rlm_dict: Dict， 每个flexible torsion所包含的constrained local minimum值
    """
    def _base_molecule(scan_curve):

        rlm_dict = {}
        for m_name, terms in scan_curve.items():
            rlm_dict[m_name] = {}
            for term, molecules in terms.items():
                tmp = [molecule for molecule in molecules 
                       if molecule.conform_type in 
                       [ConformType.LOCAL_MINIMUM, ConformType.CONSTRAINED_LOCAL_MINIMUM]]

                tmp = sorted(tmp, key=lambda m:m.energy)
                rlm_dict[m_name][term] = [molecule.constrain[0].fix_value for molecule in tmp]
        return rlm_dict

    def _base_dict(scan_curve):
        rlm_dict = {}
        for m_name, terms in scan_curve.items():
            rlm_dict[m_name] = {}
            for term, data_dict_arr in terms.items():
                tmp = [data_dict for data_dict in data_dict_arr 
                       if data_dict["conform_type"] in 
                       [ConformType.LOCAL_MINIMUM, ConformType.CONSTRAINED_LOCAL_MINIMUM]]

                tmp = sorted(tmp, key=lambda m:m["energy"])
                rlm_dict[m_name][term] = [data_dict["constrain"][0][-1] for data_dict in tmp]
        return rlm_dict
    
    for m_name,vv in scan_curve.items():
        for torsion,vvv in vv.items():
            if isinstance(vvv[0],dict):
                return _base_dict(scan_curve)
            else:
                return _base_molecule(scan_curve)

def create_lm_by_combine_scan_rlm(molecule, rlm_dict_mol, n=64, create_constrain=False,idx=None):
    """
    输入：
        molecule: Molecule
        rlm_dict_mol: 从local_minimum_pes得到
        n: int, 生成local minimum结构不超过n
        create_constratin: 旋转的二面角是否变成constrain项
    输出：
        molecules: List[Molecule], 具有新坐标的molecule对象复制
    """
    molecule.update_topol_value()
    molecules, value_arr, rotation_atoms = [], [], []
    for aa, bb in rlm_dict_mol.items():
        value_arr.append(bb)
        atoms = [int(ss) for ss in aa.split("-")]
        rotation_atoms.append(atoms)
    for rr in list(itertools.product(*value_arr))[:n]:
        molecules.append(deepcopy(molecule))
        constrain_arr = []
        for i in range(len(rr)):
            if create_constrain:
                constrain_arr.append(
                    [
                        rotation_atoms[i][0],
                        rotation_atoms[i][1],
                        rotation_atoms[i][2],
                        rotation_atoms[i][3],
                        rr[i],
                    ]
                )

            ME._structure_change(molecules[-1],rotation_atoms[i],rr[i])
            
        if create_constrain:
            molecules[-1].create_constrain(constrain_arr)
    if idx is not None:
        return molecules, idx
    else:
        return molecules

def find_stablest_molecule(molecules):
    dict_tmp = {}
    for molecule in molecules:
        if molecule.inchi_key not in dict_tmp:
            dict_tmp[molecule.inchi_key] = []
        dict_tmp[molecule.inchi_key].append(molecule)
    try:
        return [sorted(vv,key=lambda m:m.energy)[0] for vv in dict_tmp.values()]
    except:
        return [vv[0] for vv in dict_tmp.values()]
    
def remove_similar_conformer(molecules,target_molecule=None):
    if target_molecule is not None:
        target_inertia = MP._inertia_calculate(target_molecule, ignore_hydrogen=True)[0]
    list_inertia = [MP._inertia_calculate(m, ignore_hydrogen=True)[0] for m in molecules]
    indx_keep = []

    for idx, inertia in enumerate(list_inertia):
        if target_molecule is not None:
            if np.allclose(inertia, target_inertia, rtol=0, atol=5.0):
                continue
        if any(np.allclose(inertia, list_inertia[ref], rtol=0, atol=5.0) for ref in indx_keep):
            continue
        indx_keep.append(idx)

    n_duplicate = len(molecules) - len(indx_keep)
    if n_duplicate:
        logger.debug(f"{n_duplicate}/{len(molecules)} duplicated conformations removed")

    return [molecules[ii] for ii in indx_keep]

