import math
from urllib.parse import quote

import numpy as np

from ...chem.constants import CAL_TO_J, MOLAR_GAS_CONSTANT
#from compuchem.chemistry.format.mol2_file import Mol2File
#from compuchem.chemistry.format.mol_file import MolFile
#from compuchem.chemistry.structure import molecule_manager
#from compuchem.chemistry import MapMolecule as MM

try:
    import numexpr

    HAVE_NUMEXPR = True
except ImportError:
    HAVE_NUMEXPR = False


def free_energy_to_ic50(free_energy, temperature=310, unit="um"):
    # E = RTln(IC50)
    if unit == "um":
        factor = 10**6
    elif unit == "nm":
        factor = 10**9
    return math.exp(free_energy * 1000 * CAL_TO_J / MOLAR_GAS_CONSTANT / temperature) * factor


def ic50_to_free_energy(ic50, temperature=310, unit="um"):
    # to kcal/mol
    if unit == "um":
        factor = 10**6
    elif unit == "nm":
        factor = 10**9
    return MOLAR_GAS_CONSTANT * temperature * math.log(ic50 / factor) / (CAL_TO_J * 1000)


# def generate_ligand_from_sdf(sdffile):
#     with open(sdffile) as inf:
#         lines = inf.readlines()
#     script = []
#     for line in lines:
#         if line.strip() == "$$$$":
#             ligand_name = quote(script[0].strip()).replace("%", "_")
#             sdfobj = MolFile("normal")
#             sdfobj.read_file(script)
#             sdfm = sdfobj.export_moleobj()
#             sdfm = molecule_manager.MoleculeManager.prepare_molecule(sdfm)
#             sdfm = molecule_manager.MoleculeManager.assign_atom_type(sdfm, create_improper=True)
#             sdfm.update_mole_info()
#             sdfm.mole_name = ligand_name
#             for key in sdfobj.associated_data.keys():
#                 if "dg" in key:
#                     setattr(sdfm, "dg", float(sdfobj.associated_data[key]))
#                 if "ic50" in key:
#                     if "um" in key:
#                         setattr(sdfm, "ic50_um", float(sdfobj.associated_data[key]))
#                         setattr(sdfm, "dg", ic50_to_free_energy(float(sdfobj.associated_data[key]), unit="um"))
#                     else:  # default unit is nm
#                         setattr(sdfm, "ic50_nm", float(sdfobj.associated_data[key]))
#                         setattr(sdfm, "dg", ic50_to_free_energy(float(sdfobj.associated_data[key]), unit="nm"))
#             yield sdfm
#             script = []
#         else:
#             script.append(line)
#
def create_ligand_from_sdf(sdffile):
    sdfm_list = MM.molecule_create(sdffile)
    MM.molecule_structure(sdfm_list)
    MM.atom_type(sdfm_list,create_improper=True)
    for sdfm in sdfm_list:
        for key in sdfm.associated_data.keys():
            if "dg" in key:
                setattr(sdfm, "dg", float(sdfm.associated_data[key]))
            if "ic50" in key:
                if "um" in key:
                    setattr(sdfm, "ic50_um", float(sdfm.associated_data[key]))
                    setattr(sdfm, "dg", ic50_to_free_energy(float(sdfm.associated_data[key]), unit="um"))
                else:  # default unit is nm
                    setattr(sdfm, "ic50_nm", float(sdfm.associated_data[key]))
                    setattr(sdfm, "dg", ic50_to_free_energy(float(sdfm.associated_data[key]), unit="nm"))
    return sdfm_list

def old_create_ligand_from_sdf(sdffile):
    ligands_script = {}
    with open(sdffile) as inf:
        lines = inf.readlines()
    script = []
    for line in lines:
        if line.strip() == "$$$$":
            ln = quote(script[0].strip()).replace("%", "_")
            ligands_script[ln] = script
            script = []
        else:
            script.append(line)
    sdfm_list = []
    for ligand_name, sdf_script in ligands_script.items():
        sdfobj = MolFile("normal")
        sdfobj.read_file(sdf_script)
        sdfm = sdfobj.export_moleobj()
        #sdfm = molecule_manager.MoleculeManager.prepare_molecule(sdfm)
        #sdfm = molecule_manager.MoleculeManager.assign_atom_type(sdfm, create_improper=True)
        sdfm.update_mole_info()
        sdfm.mole_name = ligand_name
        sdfm_list.append(sdfm)
        for key in sdfobj.associated_data.keys():
            if "dg" in key:
                setattr(sdfm, "dg", float(sdfobj.associated_data[key]))
            if "ic50" in key:
                if "um" in key:
                    setattr(sdfm, "ic50_um", float(sdfobj.associated_data[key]))
                    setattr(sdfm, "dg", ic50_to_free_energy(float(sdfobj.associated_data[key]), unit="um"))
                else:  # default unit is nm
                    setattr(sdfm, "ic50_nm", float(sdfobj.associated_data[key]))
                    setattr(sdfm, "dg", ic50_to_free_energy(float(sdfobj.associated_data[key]), unit="nm"))
    return sdfm_list


def create_ligand_from_mol2(mol2_files):
    mol2_list = MM.molecule_create(mol2_files)
    MM.molecule_structure(mol2_list)
    for ii,mol2_file in enumerate(mol2_files):
        mol2_list[ii].mole_name = mol2_file.name.split(".")[0]
    #mol2_list = []
    #for mol2_file in mol2_files:
    #    with open(mol2_file) as f:
    #        mol2obj = Mol2File("normal")
    #        script = f.readlines()
    #        mol2obj.read_file(script)
    #        m = mol2obj.export_moleobj()
    #        m = molecule_manager.MoleculeManager.prepare_molecule(m)
    #        m.update_mole_info()
    #        m.mole_name = mol2_file.name.split(".")[0]
    #        mol2_list.append(m)
    return mol2_list


class LightEdge:
    def __init__(self, edge):
        self.structs = edge.structs
        self.name = edge.name
        self.node_names = [edge.nodes[0].name, edge.nodes[1].name]
        self.edge_cores = [edge.nodes[0].core, edge.nodes[1].core]
        self.atom_mapping = None
        self.similarity = None
        self.similarity_detail = None
        self.is_core_hopping = edge.is_core_hopping

        if edge.atom_mapping:
            n0, n1 = list(edge.atom_mapping)
            self.atom_mapping = {n0: edge.atom_mapping[n0], n1: edge.atom_mapping[n1]}
        else:
            self.atom_mapping = None


def logsumexp(a, axis=None, b=None, use_numexpr=True):
    """Compute the log of the sum of exponentials of input elements.

    Parameters
    ----------
    a : array_like
        Input array.
    axis : None or int, optional, default=None
        Axis or axes over which the sum is taken. By default `axis` is None,
        and all elements are summed.
    b : array-like, optional
        Scaling factor for exp(`a`) must be of the same shape as `a` or
        broadcastable to `a`.
    use_numexpr : bool, optional, default=True
        If True, use the numexpr library to speed up the calculation, which
        can give a 2-4X speedup when working with large arrays.

    Returns
    -------
    res : ndarray
        The result, ``log(sum(exp(a)))`` calculated in a numerically
        more stable way. If `b` is given then ``log(sum(b*exp(a)))``
        is returned.

    See Also
    --------
    numpy.logaddexp, numpy.logaddexp2, scipy.misc.logsumexp (soon to be replaced with  scipy.special.logsumexp)

    Notes
    -----
    This is based on scipy.misc.logsumexp but with optional numexpr
    support for improved performance.
    """

    a = np.asarray(a)

    a_max = np.amax(a, axis=axis, keepdims=True)

    if a_max.ndim > 0:
        a_max[~np.isfinite(a_max)] = 0
    elif not np.isfinite(a_max):
        a_max = 0

    if b is not None:
        b = np.asarray(b)
        if use_numexpr and HAVE_NUMEXPR:
            out = np.log(numexpr.evaluate("b * exp(a - a_max)").sum(axis))
        else:
            out = np.log(np.sum(b * np.exp(a - a_max), axis=axis))
    else:
        if use_numexpr and HAVE_NUMEXPR:
            out = np.log(numexpr.evaluate("exp(a - a_max)").sum(axis))
        else:
            out = np.log(np.sum(np.exp(a - a_max), axis=axis))

    a_max = np.squeeze(a_max, axis=axis)
    out += a_max

    return out
