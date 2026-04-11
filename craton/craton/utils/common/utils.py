import os
import re
import sys
from collections.abc import Iterable
from typing import List
from copy import deepcopy


import numpy as np
import simplejson


#from .errors import OpenBabelError


def parse_string(str):
    if isinstance(str,list):
        tmp = []
        for rr in str:
            if rr[0] == "0":
                tmp.append(rr)
            else:
                try:
                    tmp.append(int(rr))
                except:
                    try:
                        tmp.append(float(rr))
                    except:
                        tmp.append(rr)
        return tmp
    else:
        if str[0] == "0":
            return str
        try:
            return int(str)
        except:
            try:
                return float(str)
            except:
                return str


def parse_inchi_key(string):
    pattern = r"[A-Z]{14}-[A-Z]{10}-[A-Z]"
    try:
        inchi = re.search(pattern, string).group()
        if not inchi:
            raise Exception(f"InChI key not found: {string}")
        return inchi
    except:
        #logger.warning("non inchi key string: %s" % string)
        return None


def greatest_common_divisor(numbers):
    """
    Calculate the greatest common divisor.

    Parameters
    ----------
    numbers : list of int

    Returns
    divisor : int
    """
    minimal = min(numbers)
    for i in range(minimal, 1, -1):
        flag = True
        for number in numbers:
            if number % i != 0:
                flag = False
        if flag:
            return i
    return 1


def random_string(length=8):
    """
    Generate a random string in specified length. The string contains only upper and lower case ascii letters.


    Parameters
    ----------
    length : int

    Returns
    -------
    string : str
    """
    import random
    import string

    return "".join(random.sample(string.ascii_letters, length))


def cd_or_create_and_cd(dir):
    """
    Go to the target directory. If not exist, create this directory and go to it.

    Parameters
    ----------
    dir : str
    """
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except:
            raise Exception("Cannot create directory: %s" % dir)

    try:
        os.chdir(dir)
    except:
        raise Exception("Cannot read directory: %s" % dir)


def create_mol_from_smiles(smiles: str, minimize=True, pdb_out=None, mol2_out=None, resname=None):
    """
    Create a openbabel molecule object from SMILES string.

    Parameters
    ----------
    smiles : str
    minimize : bool
    pdb_out : str, optional
    mol2_out : str, optional
    resname : str, optional

    Returns
    -------
    mol : pybel.Molecule
    """
    try:
        from openbabel import pybel
    except ImportError:
        raise ImportError("OpenBabel is required for parsing SMILES")

    try:
        py_mol = pybel.readstring("smi", smiles)
    except:
        raise OpenBabelError("Invalid SMILES")

    py_mol.addh()
    py_mol.make3D()
    if minimize:
        py_mol.localopt()

    if resname is not None:
        obmol = py_mol.OBMol
        res = obmol.GetResidue(0)
        if res is not None:
            res.SetName("UNL")

    if pdb_out is not None:
        py_mol.write("pdb", pdb_out, overwrite=True)
    if mol2_out is not None:
        py_mol.write("mol2", mol2_out, overwrite=True)
        if resname is not None:
            with open(mol2_out) as f:
                content = f.read()
            content = content.replace("UNL", resname[:3])
            with open(mol2_out, "w") as f:
                f.write(content)
    return py_mol


def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


class NumpyEncoder(simplejson.JSONEncoder):
    """Custom encoder for numpy data types"""

    def default(self, obj):
        if isinstance(
            obj,
            (
                np.int_,
                np.intc,
                np.intp,
                np.int8,
                np.int16,
                np.int32,
                np.int64,
                np.uint8,
                np.uint16,
                np.uint32,
                np.uint64,
            ),
        ):

            return int(obj)

        elif isinstance(obj, (np.float_, np.float16, np.float32, np.float64)):
            return float(obj)

        elif isinstance(obj, (np.complex_, np.complex64, np.complex128)):
            return {"real": obj.real, "imag": obj.imag}

        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()

        elif isinstance(obj, (np.bool_)):
            return bool(obj)

        elif isinstance(obj, (np.void)):
            return None

        return simplejson.JSONEncoder.default(self, obj)


class TopologyEncoder(simplejson.JSONEncoder):

    def default(self, obj):
        from ...chem.topology import Bond, Angle, Dihedral, Improper
        if isinstance(obj, Bond):
            return [obj.a1, obj.a2]
        elif isinstance(obj, Angle):
            return [obj.a1, obj.a2, obj.a3]
        elif isinstance(obj, Dihedral):
            return [obj.a1, obj.a2, obj.a3, obj.a4]
        elif isinstance(obj, Improper):
            return [obj.a1, obj.a2, obj.a3, obj.a4]
        elif isinstance(obj, set):
            return list(obj)
        return simplejson.JSONEncoder.default(self, obj)


def read_molecule_csv_file(f):
    with open(f) as inf:
        lines = inf.readlines()
    datas = []
    key_index = {term: ii for ii, term in enumerate(lines[0].strip().split(","))}
    for line in lines[1:]:
        if line[0] != "#":
            ss = line.strip().split(",")
            datas.append({kk: ss[vv] for kk, vv in key_index.items()})
    return datas


def add_two_dict(dict1, dict2):
    for key in dict1:
        if key not in dict2.keys():
            dict2[key] = dict1[key]
        else:
            if dict2[key] == "NONE":
                dict2[key] = dict1[key]
            else:
                if dict1[key] != "NONE":
                    for data in dict1[key]:
                        dict2[key].append(data)

    return dict2


def combine_arr(a: List[list]) -> List[list]:
    """Input a is a 2-D array. Merge all arrays which shares at least 1 common element
    e.g, [[1,2], [2,3], [4,5], [6,7], [8,9], [7,8]]
    returns [[1,2,3], [4,5], [6,7,8,9]]
    """
    b = len(a)
    for i in range(b):
        if a[i] == []:
            continue
        for j in range(b):
            if i == j:
                break
            if a[j] == []:
                continue
            merged = list(set(a[i] + a[j]))
            total_len = len(a[j]) + len(a[i])
            if len(merged) < total_len:  # there is intersection
                a[i], a[j] = merged, []
    return [ele for ele in a if ele != []]


def datas_to_bar(datas, start_point, end_point, interval):
    results = {vv + int(interval / 2): 0 for vv in range(start_point, end_point, interval)}
    for v in datas:
        key = (
            int(v / interval) * interval + int(interval / 2)
            if v >= 0
            else int(v / interval) * interval - int(interval / 2)
        )
        if key not in results:
            sys.exit("key error")
        results[key] += 1
    return results


def read_trosion_distribution_data(datas):
    def get_md_distribution(datas, keys):
        dist_tmp = {aa: [] for aa in keys}
        for item, value in datas.items():
            for rr in value:
                dist_tmp[rr["name"]].append(rr["value"])
        return dist_tmp

    tmp = {}

    for item, value in datas["energy"].items():
        angle = float(item)
        for rr in value:
            if rr["name"] not in tmp:
                tmp[rr["name"]] = {
                    "torsion_atoms": [int(aa) for aa in rr["name"].split("-")],
                    "potential_energy_surface": [],
                }
            tmp[rr["name"]]["potential_energy_surface"].append({"x": angle, "y": rr["value"], "conformation": "sdf"})
    for kk in ["rbfe", "rhfe"]:
        dist_tmp = get_md_distribution(datas[kk], list(tmp.keys()))
        for item, value in dist_tmp.items():
            tmp[item][f"{kk}_distribution"] = [
                {"x": aa, "y": bb} for aa, bb in datas_to_bar(value, -180, 180, 10).items()
            ]
    for rr in datas["rbfe"]["0"]:
        tmp[rr["name"]]["initial_position"] = (
            int(rr["value"] / 10) * 10 + 5 if rr["value"] >= 0 else int(rr["value"] / 10) * 10 - 5
        )
    for kk, vv in tmp.items():
        min_ener = min([v["y"] for v in vv["potential_energy_surface"]])
        for rr in vv["potential_energy_surface"]:
            rr["y"] = rr["y"] - min_ener

    return tmp


def read_md_json(f):
    datas = simplejson.loads(open(f).read())
    if isinstance(datas["conformation"], dict):
        for mol in datas["conformation"]:
            torsion_datas = {
                "molecule": datas["conformation"][mol],
                "items": [vv for kk, vv in read_trosion_distribution_data(datas["torsion"][mol]).items()],
            }
            interaction_datas = datas["interaction"][mol]["interaction"]
            energy_datas = {}
            for kk, vv in datas["interaction"][mol]["energy"].items():
                energy_datas[kk] = {
                    "total": {"mean": vv["total_mean"], "std": vv["total_std"]},
                    "Coul-SR": {"mean": vv["Coul-SR_mean"], "std": vv["Coul-SR_std"]},
                    "LJ-SR": {"mean": vv["LJ-SR_mean"], "std": vv["LJ-SR_std"]},
                }
            with open(f"{mol}_torsion.json", "w") as outf:
                outf.write(simplejson.dumps(torsion_datas))
            with open(f"{mol}_interaction.json", "w") as outf:
                outf.write(simplejson.dumps(interaction_datas))
            with open(f"{mol}_energy.json", "w") as outf:
                outf.write(simplejson.dumps(energy_datas))
    else:
        torsion_datas = {
            "molecule": datas["conformation"],
            "items": [vv for kk, vv in read_trosion_distribution_data(datas["torsion"])],
        }
        interaction_datas = datas["interaction"]["interaction"]
        energy_datas = {}
        for kk, vv in datas["interaction"]["energy"].items():
            energy_datas[kk] = {
                "total": {"mean": vv["total_mean"], "std": vv["total_std"]},
                "Coul-SR": {"mean": vv["Coul-SR_mean"], "std": vv["Coul-SR_std"]},
                "LJ-SR": {"mean": vv["LJ-SR_mean"], "std": vv["LJ-SR_std"]},
            }
        with open(f"torsion.json", "w") as outf:
            outf.write(simplejson.dumps(torsion_datas))
        with open(f"interaction.json", "w") as outf:
            outf.write(simplejson.dumps(interaction_datas))
        with open(f"energy.json", "w") as outf:
            outf.write(simplejson.dumps(energy_datas))

