import time
from collections import OrderedDict
from copy import deepcopy
from functools import partial
from itertools import combinations, repeat
from multiprocessing import Manager, Pool
from pathlib import Path

import psutil
from rdkit.Chem import AllChem, Draw

from ....utils import logger
####from compuchem.chemistry.format.call_rdkit import CallRdkit
from .fuzzy_mapping import get_matches
from .similarity import LightEdge ## , create_ligand_from_mol2, create_ligand_from_sdf
from ....utils.commons import parallel_run


_DEFAULT_NUM_PROCS = psutil.cpu_count(logical=False) - 1


def matcher(light_edge: LightEdge):
    #manager_dict
    struct0, struct1 = light_edge.structs
    fuzzy_params = {"timeout": 50}

    if sum(struct0.formal_charge) != sum(struct1.formal_charge):
        logger.info(f"struct {struct0.mole_name} and struct {struct1.mole_name} charge are different")
        logger.info(f"{struct0.mole_name}: {sum(struct0.formal_charge)}")
        logger.info(f"{struct1.mole_name}: {sum(struct1.formal_charge)}")
        fuzzy_params["allow_ring_breaking"] = False
    ####fuzzy_params["allow_ring_breaking"] = False
    allowed_atoms0 = light_edge.edge_cores[0]
    allowed_atoms1 = light_edge.edge_cores[1]
    #mapping_result = get_matches(
    #    deepcopy(struct0),
    #    deepcopy(struct1),
    #    allowed_atoms0=allowed_atoms0,
    #    allowed_atoms1=allowed_atoms1,
    #    **fuzzy_params,
    #)
    #print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    #print(mapping_result)
    #print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
    #try:
    mapping_result = get_matches(
        deepcopy(struct0),
        deepcopy(struct1),
        allowed_atoms0=allowed_atoms0,
        allowed_atoms1=allowed_atoms1,
        **fuzzy_params,
    )
    light_edge.atom_mapping = {
        struct0.mole_name: list(mapping_result.keys()),
        struct1.mole_name: list(mapping_result.values()),
    }
    #manager_dict[light_edge.name] = light_edge
    return light_edge

    #return manager_dict
    #except:
    #    logger.error(f"MCS Finding Failed, the error pair: {struct0.mole_name}, {struct1.mole_name}")


def assign_atom_mapping(edges, parallel=True):
    #logger.info("Start atom mapping ...")
    light_edges = [LightEdge(edge) for edge in edges]
    manager_dict = {}
    if parallel:
        this_edges = parallel_run(matcher,light_edges,keep_order=False)
        for this_edge in this_edges:
            manager_dict[this_edge.name] = this_edge
    else:
        for edge in edges:
            manager_dict[edge.name] = matcher(LightEdge(edge),manager_dict)
    for edge in edges:
        edge.atom_mapping = manager_dict[edge.name].atom_mapping

    #with Manager() as manager:
    #    manager_dict = manager.dict()
    #    with manager.Pool(processes=num_procs) as pool:
    #        pool.starmap(matcher, zip(light_edges, repeat(manager_dict, len(light_edges))))
    #        for edge in edges:
    #            edge.atom_mapping = manager_dict[edge.name].atom_mapping
    
    
    
    
    
    # for edge in edges:
    #     for key, value in edge.atom_mapping.items():
    #         logger.debug(f"{key}:{value}")


def _get_matches_manager(pair, mol_pair, **kwargs):
    return pair, get_matches(deepcopy(mol_pair[0]), deepcopy(mol_pair[1]), **kwargs)


# for debuging
def _generate_atom_mapping(input, pairs=None, num_procs=None, **kwargs):
    """Genearte atom mapping for given input file and pairs
    input: str or dict
        str, must be sdf file
        dict, mol_name: inf_mol
    """
    num_procs = num_procs or _DEFAULT_NUM_PROCS
    if isinstance(input, dict):
        mol_dict = input
    else:
        if str(input).endswith(".sdf"):
            mol_objs = create_ligand_from_sdf(input)
        else:
            mol2_files = list(Path(input).glob("*.mol2"))
            mol_objs = create_ligand_from_mol2(mol2_files)
        mol_dict = {mol.mole_name: mol for mol in mol_objs}

    if pairs is None:
        node_list = [mol_obj.mole_name for mol_obj in mol_objs]
        pairs = list(combinations(node_list, 2))
    else:
        pairs = [tuple(pair) for pair in pairs]
    mol_pairs = [(mol_dict[i], mol_dict[j]) for (i, j) in pairs]

    start = time.time()
    with Pool(processes=num_procs) as pool:
        multi_result = pool.starmap(partial(_get_matches_manager, **kwargs), zip(pairs, mol_pairs))

    # {(node1, node2): {1:2, 2:3, 3:4, ....}}
    map_result = {key: mapping_result for (key, mapping_result) in multi_result}
    logger.debug(f"Atom mmaping time: {time.time() - start:.3f}")
    return mol_dict, map_result


# def visualize_atom_mapping(mol_dict, map_result, pairs=None, filename="atom_mapping.svg"):
#     rdkmol_list = []
#     legends_list = []
#     highlight_atoms_list = []
#     if pairs is None:
#         pairs = combinations(mol_dict.keys(), 2)
#     for pair in pairs:
#         for i, node in enumerate(pair):
#             rdkobj = CallRdkit("3d")
#             rdkobj.import_moleobj(mol_dict[node])
#             rdkmh = rdkobj.rdkmh
#             AllChem.Compute2DCoords(rdkmh)
#             if i == 0:
#                 mcs = list(map_result[pair].keys())
#             else:
#                 mcs = list(map_result[pair].values())
#             rdkmol_list.append(rdkmh)
#             legends_list.append(node)
#             highlight_atoms_list.append(mcs)
#     molsPerRow = 4 if len(rdkmol_list) > 10 else 2
#     image = Draw.MolsToGridImage(
#         rdkmol_list,
#         molsPerRow=molsPerRow,
#         legends=legends_list,
#         highlightAtomLists=highlight_atoms_list,
#         subImgSize=[400, 400],
#         addAtomIndices=True,
#         useSVG=True,
#     )
#     with open(filename, "w") as f:
#         f.write(image)
#
#
# def visualize_atom_mapping_detail(mol_dict, map_result, pair, filename=None):
#     rdkobj = CallRdkit("3d")
#     rdkobj.import_moleobj(mol_dict[pair[0]])
#     rdkmh1 = rdkobj.rdkmh
#     rdkobj = CallRdkit("3d")
#     rdkobj.import_moleobj(mol_dict[pair[1]])
#     rdkmh2 = rdkobj.rdkmh
#     mapping_data = OrderedDict(sorted(map_result[pair].items()))
#
#     AllChem.Compute2DCoords(rdkmh1)
#     AllChem.Compute2DCoords(rdkmh2)
#
#     rdkmol_list = []
#     legends_list = []
#     highlight_atoms_list = []
#     for atom1, atom2 in mapping_data.items():
#         rdkmol_list.append(rdkmh1)
#         rdkmol_list.append(rdkmh2)
#         legends_list.append(str(atom1))
#         legends_list.append(str(atom2))
#         highlight_atoms_list.append([atom1])
#         highlight_atoms_list.append([atom2])
#
#     molsPerRow = 4 if len(rdkmol_list) > 10 else 2
#     image = Draw.MolsToGridImage(
#         rdkmol_list,
#         molsPerRow=molsPerRow,
#         legends=legends_list,
#         highlightAtomLists=highlight_atoms_list,
#         subImgSize=[400, 400],
#         addAtomIndices=True,
#         useSVG=True,
#     )
#     if filename is None:
#         filename = f"{pair[0]}-{pair[1]}.svg"
#     with open(filename, "w") as f:
#         f.write(image)
#
#
# def show_graph_atom_mapping(g):
#     mol_dict = {}
#     for name, node in g.nodes_dict.items():
#         mol_dict[name] = node.struct
#     map_result = {}
#     pairs = []
#     for edge in g.edges_iter():
#         pairs.append(edge.nodes_name)
#         key = tuple(edge.atom_mapping.keys())
#         value = dict(zip(*edge.atom_mapping.values()))
#         map_result[key] = value
#     visualize_atom_mapping(mol_dict, map_result, pairs)

def visualize_atom_mapping(mol1, mol2, matches):
    rdkobj = CallRdkit("3d")
    rdkobj.import_moleobj(mol1)
    rdkmh1 = rdkobj.rdkmh
    rdkobj = CallRdkit("3d")
    rdkobj.import_moleobj(mol2)
    rdkmh2 = rdkobj.rdkmh

    AllChem.Compute2DCoords(rdkmh1)
    AllChem.Compute2DCoords(rdkmh2)

    rdkmol_list = [rdkmh1, rdkmh2]
    legends_list = [mol1.name, mol2.name]
    highlight_atoms_list = [list(matches.keys()), list(matches.values())]

    molsPerRow = 2
    image = Draw.MolsToGridImage(
        rdkmol_list,
        molsPerRow=molsPerRow,
        legends=legends_list,
        highlightAtomLists=highlight_atoms_list,
        subImgSize=[400, 400],
        addAtomIndices=True,
        useSVG=True,
    )
    with open(f"{mol1.name}_to_{mol2.name}.svg", "w") as f:
        f.write(image)
