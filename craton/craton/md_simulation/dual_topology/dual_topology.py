import dataclasses
import json
import math
from copy import deepcopy
from itertools import product
from typing import Union

import simplejson
from joblib import Parallel, delayed
from openmm import unit
from simtk import openmm as mm

from ...utils import logger
###from compuchem.chemistry.constants.constants import FEPProperty
from ...chem.molecule import Molecule
##from compuchem.chemistry.molecule_edit.mole_analy import calc_stru_para
from ...software import gromacs
from ...chem.topology import Pair, Bond
from ...utils.common.utils import TopologyEncoder
####from compuchem.molecule_dynamics.algorithm.fep.mapping import visualize_atom_mapping
from .ring_breaking.ring_breaking import RingBreaking
#from ..mapping.algorithm.similarity import LightEdge
from .fep_ring_breaking import construct_broken_structure, nonbond_parameters, \
    bond_parameters, debug_info
####from compuchem.molecule_dynamics.run.restrain import Restraints
from ...utils.commons import parallel_run

EXCLUSION_DISTANCE = 5

def charge_abfe(molecule):
    if molecule.net_charge == 0:
        molecule.absolute_intra_flag = True
        return molecule
        #self.abfe_setting.intra_mol_flag = False  # if net change, we should separate the topology A and B
    for atom in molecule.Atoms:
        atom.ff_charge_m2 = 0.0
        atom.atom_type_name_m2 = "_D"
        atom.mass_m2 = atom.mass
        atom.parameter_m2 = [atom.parameter[0], 0.0]
    molecule.absolute_intra_flag = False
    return molecule

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
        self.is_charge_hopping = edge.is_charge_hopping

        if edge.atom_mapping:
            n0, n1 = list(edge.atom_mapping)
            self.atom_mapping = {n0: edge.atom_mapping[n0], n1: edge.atom_mapping[n1]}
        else:
            self.atom_mapping = None


class FEPTopology:
    def __init__(self, mol1: Molecule, mol2: Molecule, matches, is_core_hopping: bool = False, is_charge_hopping: bool = False):
        self.mol1 = mol1
        self.mol2 = mol2
        self._mol1 = deepcopy(mol1)
        self._mol2 = deepcopy(mol2)
        self._mol1_mapping = list(matches.keys())
        self._mol2_mapping = list(matches.values())
        self.matches = matches
        self.is_core_hopping = is_core_hopping
        self.is_charge_hopping = is_charge_hopping
        # visualize_atom_mapping(mol1, mol2, matches)

    def _update_topology_coordinates(self):
        system = mm.System()
        for atom in self._mol1.Atoms:
            system.addParticle(atom.mass)
        indexes1 = [i for i, atom in enumerate(self._mol1.Atoms) if not atom.atom_type_name == "_D"]
        indexes2 = [
            i
            for i, atom in enumerate(self._mol1.Atoms)
            if not getattr(atom, "atom_type_name_m2", atom.atom_type_name) == "_D"
        ]

        rmsd1 = mm.RMSDForce([a.coor for a in self._mol1.Atoms], indexes1)
        rmsd2 = mm.RMSDForce([a.coor for a in self._mol2.Atoms], indexes2)
        cv = mm.CustomCVForce("rmsd1^2+rmsd2^2")
        cv.addCollectiveVariable("rmsd1", rmsd1)
        cv.addCollectiveVariable("rmsd2", rmsd2)
        system.addForce(cv)
        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName("CPU")
        context = mm.Context(system, integrator, platform)
        context.setPositions([a.coor for a in self._mol1.Atoms])
        logger.debug(f"RMSD = {cv.getCollectiveVariableValues(context)}")
        mm.LocalEnergyMinimizer.minimize(context, tolerance=1e-10)
        logger.debug(f"RMSD = {cv.getCollectiveVariableValues(context)}")
        positions = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        for atom, position in zip(self._mol1.Atoms, positions):
            atom.coor = position.tolist()[:]

    def _update_atom_type_topology(self):
        """
        if atom in A not mapped:
            create a dummy atom for B
            map[not mapped A] = dummy atom in B
        if atom in B not mapped:
            create a dummy atom for A
            map[dummy atom in A] = not mapped B
        """

        self._mapping = {i: j for i, j in zip(self._mol1_mapping, self._mol2_mapping)}
        atom1_not_mapped = [atom for atom in self._mol1.Atoms if atom.No not in self._mol1_mapping]
        atom2_not_mapped = [atom for atom in self._mol2.Atoms if atom.No not in self._mol2_mapping]

        self._mol_A_dummy_atoms, self._mol_B_dummy_atoms = [], []
        for i, atom_in_A in enumerate(atom1_not_mapped, start=len(self._mol2.Atoms)):
            dummy_atom_for_mol_B = deepcopy(atom_in_A)
            dummy_atom_for_mol_B.connectivity = []
            dummy_atom_for_mol_B.No = i
            dummy_atom_for_mol_B.atom_type_name = "_D"
            dummy_atom_for_mol_B.elem = atom_in_A.elem
            dummy_atom_for_mol_B.formal_charge = 0
            dummy_atom_for_mol_B.ff_charge = 0.0
            dummy_atom_for_mol_B.parameter = [atom_in_A.parameter[0], 0.0]
            self._mol2.Atoms.append(dummy_atom_for_mol_B)
            self._mapping[atom_in_A.No] = i
            self._mol_B_dummy_atoms.append(dummy_atom_for_mol_B)

        for i, atom_in_B in enumerate(atom2_not_mapped, start=len(self._mol1.Atoms)):
            dummy_atom_for_mol_A = deepcopy(atom_in_B)
            dummy_atom_for_mol_A.connectivity = []
            dummy_atom_for_mol_A.No = i
            dummy_atom_for_mol_A.atom_type_name = "_D"
            dummy_atom_for_mol_A.elem = atom_in_B.elem
            dummy_atom_for_mol_A.formal_charge = 0
            dummy_atom_for_mol_A.ff_charge = 0.0
            dummy_atom_for_mol_A.parameter = [atom_in_B.parameter[0], 0.0]
            self._mol1.Atoms.append(dummy_atom_for_mol_A)
            self._mapping[i] = atom_in_B.No
            self._mol_A_dummy_atoms.append(dummy_atom_for_mol_A)

        for k, v in self._mapping.items():
            atom1 = self._mol1.Atoms[k]
            atom2 = self._mol2.Atoms[v]
            if any(
                    (
                            atom1.ff_charge != atom2.ff_charge,
                            atom1.atom_type_name != atom2.atom_type_name,
                            atom1.mass != atom2.mass,
                    )
            ):
                atom1.ff_charge_m2 = atom2.ff_charge
                atom1.atom_type_name_m2 = atom2.atom_type_name
                atom1.mass_m2 = atom2.mass
                atom1.parameter_m2 = atom2.parameter
        self._mol1.net_charge_m2 = self._mol2.net_charge
        self._mapping_reverse = {v: k for k, v in self._mapping.items()}

        for dummy_atom_for_mol_A in self._mol_A_dummy_atoms:
            self._mol1.Atoms[dummy_atom_for_mol_A.No].connect = [
                self._mapping_reverse[i] for i in dummy_atom_for_mol_A.connect
            ]
        for dummy_atom_for_mol_B in self._mol_B_dummy_atoms:
            self._mol2.Atoms[dummy_atom_for_mol_B.No].connect = [self._mapping[i] for i in dummy_atom_for_mol_B.connect]
        self._mol2.reorder_atoms({v: k for k, v in self._mapping.items()})
        self.idxes1_dummy = {atom.No for atom in self._mol_A_dummy_atoms}
        self.idxes2_dummy = {atom.No for atom in self._mol_B_dummy_atoms}

    def _make_restraints(self):
        pass
        # phi1 = calc_stru_para([self._mol1.Atoms[0].coor, self._mol1.Atoms[1].coor,
        #                       self._mol1.Atoms[2].coor, self._mol1.Atoms[3].coor])
        # phi2 = calc_stru_para([self._mol1.Atoms[52].coor, self._mol1.Atoms[51].coor,
        #                       self._mol1.Atoms[3].coor, self._mol1.Atoms[2].coor])
        # restraint = Restraints()
        # restraint.make_restraint(restrain_type="dihedral", restrain_func="default",
        #                             a1=1, a2=2, a3=3, a4=4, phi=phi1, dphi=10, fc=10.0)
        # restraint.make_restraint(restrain_type="dihedral", restrain_func="default",
        #                          a1=53, a2=52, a3=4, a4=3, phi=phi2, dphi=10, fc=10.0)
        # self._mol1.restraints = restraint


    def _update_connection_topology(self):
        def _equal_str_float(a: Union[str, float], b: Union[str, float], threshold=1e-6):
            if type(a) is str and type(b) is str:
                return a == b
            elif type(a) is float and type(b) is float:
                return math.isclose(a, b, abs_tol=threshold)
            elif (type(a) is float and type(b) is str) or (type(a) is str and type(b) is float):
                raise Exception(f"type of input is not the same: {type(a)}, {type(b)}")
            else:
                raise Exception("Not str or float: %s" % type(a))
        ######生成键、角、二面角等字典
        dicts_1 = {}
        dicts_2 = {}
        
        
        for attr in Molecule.attrs_topol:
            terms1 = getattr(self._mol1, attr, [])
            terms2 = getattr(self._mol2, attr, [])
            dicts_1 = {}
            dicts_2 = {}
            for term1 in terms1:
                dicts_1[term1.str] = term1
            for term2 in terms2:
                dicts_2[term2.str] = term2
            for term1 in terms1:
                #if term1 not in terms2:
                if term1.str not in dicts_2:
                    terms2.append(term1)
                    if attr == "Bonds":
                        self._mol2.Atoms[term1.a1].connectivity.append(term1.a2)
                        self._mol2.Atoms[term1.a1].bond_type.append(term1.get_type(self._mol1))
                        self._mol2.Atoms[term1.a1].bond_type_aromatic.append(term1.get_type_aromatic(self._mol1))
                        #self._mol2.Atoms[term1.a1].bond_type_conjugate.append(term1.get_type_aromatic(self._mol1))
                        #self._mol2.Atoms[term1.a1].connectivity_type.append(term1.get_type_aromatic(self._mol1))
                        self._mol2.Atoms[term1.a2].connectivity.append(term1.a1)
                        self._mol2.Atoms[term1.a2].bond_type.append(term1.get_type(self._mol1))
                        self._mol2.Atoms[term1.a2].bond_type_aromatic.append(term1.get_type_aromatic(self._mol1))
                        #self._mol2.Atoms[term1.a2].bond_type_conjugate.append(term1.get_type_aromatic(self._mol1))
                        #self._mol2.Atoms[term1.a2].connectivity_type.append(term1.get_type_aromatic(self._mol1))
                else:
                    term2 = dicts_2[term1.str]
                    if not all((_equal_str_float(para, term2.parameter[i]) for i, para in enumerate(term1.parameter))):
                        term1.use_parameter_m2 = "m1-m2"
                        term1.parameter_m2 = term2.parameter[:]
                    else:
                        term1.parameter_m2 = term1.parameter[:]
                
                
                #try:
                #    term2 = next(term for term in terms2 if term == term1)
                    # topologies in both molecules
                #    if not all((_equal_str_float(para, term2.parameter[i]) for i, para in enumerate(term1.parameter))):
                #        term1.use_parameter_m2 = "m1-m2"
                #        term1.parameter_m2 = term2.parameter[:]
                #    else:
                #        term1.parameter_m2 = term1.parameter[:]
                #except StopIteration:
                #    continue
            for term2 in terms2:  # topology only in second molecule
                #if term2 not in terms1:
                if term2.str not in dicts_1:
                    if attr == "Pair14" or attr == "Pair1n":  # for openmm molecule_dynamics
                        term2.parameter = [0.0, 0.0]
                        term2.charge_parameter = [0.0, 0.0]
                    term2.parameter_m2 = term2.parameter[:]
                    terms1.append(term2)
                    if attr == "Bonds":
                        self._mol1.Atoms[term2.a1].connectivity.append(term2.a2)
                        self._mol1.Atoms[term2.a1].bond_type.append(term2.get_type(self._mol2))
                        self._mol1.Atoms[term2.a1].bond_type_aromatic.append(term2.get_type_aromatic(self._mol2))
                        #self._mol2.Atoms[term2.a1].bond_type_conjugate.append(term1.get_type_aromatic(self._mol2))
                        #self._mol2.Atoms[term2.a1].connectivity_type.append(term1.get_type_aromatic(self._mol2))
                        self._mol1.Atoms[term2.a2].connectivity.append(term2.a1)
                        self._mol1.Atoms[term2.a2].bond_type.append(term2.get_type(self._mol2))
                        self._mol1.Atoms[term2.a2].bond_type_aromatic.append(term2.get_type_aromatic(self._mol2))
                        #self._mol2.Atoms[term2.a2].bond_type_conjugate.append(term1.get_type_aromatic(self._mol2))
                        #self._mol2.Atoms[term2.a2].connectivity_type.append(term1.get_type_aromatic(self._mol2))
    def _update_exclusion_list_topology(self):
        #logger.info(f"Determining exclusion list ...")
        excl_pairs = {(a1.No, a2.No) for a1, a2 in product(self._mol_A_dummy_atoms, self._mol_B_dummy_atoms)}
        if self.is_core_hopping:
            for pair in self._mol1.fep_exclusions_pairs:
                excl_pairs.add((pair.a1, pair.a2))
        self._mol1.rfe_exclusions = list(excl_pairs)

    def _update_ring_breaking_topology(self):
        ring_breaking = RingBreaking(self.mol1, self.mol2, self.matches)
        mol_a_bonds_formation = ring_breaking.wt_result.bond_formation
        mol_a_bonds_delete = ring_breaking.wt_result.bond_delete
        mol_b_bonds_formation = ring_breaking.mut_result.bond_formation
        mol_b_bonds_delete = ring_breaking.mut_result.bond_delete
        mol_b_bonds_formation = [(self._mapping_reverse[a1], self._mapping_reverse[a2]) for (a1, a2) in
                                 mol_b_bonds_formation]
        ic(mol_a_bonds_delete, mol_a_bonds_formation, mol_b_bonds_delete, mol_b_bonds_formation)
        mol_b_bonds_delete = [(self._mapping_reverse[a1], self._mapping_reverse[a2]) for (a1, a2) in mol_b_bonds_delete]
        mol_a_dummy, mol_b_dummy = construct_broken_structure(self._mol1, self._mol2, mol_a_bonds_formation,
                                                              mol_a_bonds_delete,
                                                              mol_b_bonds_formation, mol_b_bonds_delete)
        bond_para = bond_parameters(mol_a_dummy, mol_b_dummy, self.idxes1_dummy, self.idxes2_dummy)
        nonbond_para = nonbond_parameters(mol_a_dummy, mol_b_dummy, self.idxes1_dummy, self.idxes2_dummy)

        common_14_pair = nonbond_para.all_14_pair.difference(nonbond_para.a14_disappear).difference(
            nonbond_para.a14_grow)
        self._mol1.Bonds = bond_para.bonds_all
        self._mol1.Angles = bond_para.angles_all
        self._mol1.Dihedrals = bond_para.dihedrals_all
        self._mol1.Impropers = bond_para.impropers_all
        self._mol1.Pair14 = [Pair("LJ12_6", *pair) for pair in common_14_pair]
        self._mol1.fep_exclusions_pairs = [Pair("LJ12_6", *pair) for pair in nonbond_para.exclusion_pair]
        self._mol1.fep_extra_pairs = [Pair("grow", *pair) for pair in nonbond_para.a14_grow] + \
                                     [Pair("disappear", *pair) for pair in nonbond_para.a14_disappear]
        self._mol1.fep_extra_pairs_nb = [Pair("grow", *pair) for pair in nonbond_para.a1n_grow] + \
                                        [Pair("disappear", *pair) for pair in nonbond_para.a1n_disappear]
        gromacs.write_mole_itp(self._mol1, file_name="mol_final.itp")
        with open("debug_info.json", 'w') as f:
            simplejson.dump(dataclasses.asdict(debug_info), f, indent=4, cls=TopologyEncoder)

        ic(bond_para.bonds_grow, bond_para.bonds_disappear, bond_para.angles_disappear, bond_para.angles_grow,
           bond_para.dihedrals_grow, bond_para.dihedrals_disappear,bond_para.impropers_grow, bond_para.impropers_disappear,
           nonbond_para.a1n_disappear, nonbond_para.a1n_grow,
           nonbond_para.a14_disappear, nonbond_para.a14_grow)

    def dual_topology(self):
        self._update_atom_type_topology()
        self._update_connection_topology()
        if self.is_core_hopping:
            self._update_ring_breaking_topology()
        # self._update_topology_coordinates()
        self._make_restraints()
        self._update_exclusion_list_topology()
        self._mol1.atom_mapping = self._mapping
        self._mol1.left_molecule_name = self.mol1.mole_name
        self._mol1.right_molecule_name = self.mol2.mole_name
        if self.is_charge_hopping:
            self._mol1.dual_topology_type = "charge_hopping"
            self._mol1.dual_charge = 0
        elif self.is_core_hopping:
            self._mol1.dual_topology_type = "core_hopping"
        else:
            self._mol1.dual_topology_type = "r_group"
        self._mol1.mole_name = f"{self._mol1.mole_name}_to_{self._mol2.mole_name}"
        # self._mol1.draw(f"{self._mol1.name}_dual.svg")
        return self._mol1, self._mol2, self._mapping


def build_dual_topology_old(edge):
    mol1, mol2 = edge.structs
    atom_mapping = edge.atom_mapping
    mol_name1, mol_name2 = atom_mapping.keys()
    if mol1.name != mol_name1 and mol1.name == mol_name2:
        mol2_index, mol1_index = atom_mapping.values()
    else:
        mol1_index, mol2_index = atom_mapping.values()
    fep_topology = FEPTopology(mol1, mol2, {i: j for i, j in zip(mol1_index, mol2_index)}, edge.is_core_hopping)
    return fep_topology.dual_topology()


def dual_topology_assign_old(gg, num_procs=1):
    edges = list(gg.edges_iter())
    light_edges = [LightEdge(edge) for edge in edges]
    mol1,mol2 = edges[0].structs
    #topologies = Parallel(n_jobs=8)(delayed(build_dual_topology)(edge) for edge in light_edges)
    topologies = []
    for edge in light_edges:
        topologies.append(build_dual_topology_old(edge))
    for edge, (topology1, topology2, mapping) in zip(edges, topologies):
        if topology1.name != edge.name:
            raise RuntimeError("Dual topology name does not math the edge name")
        edge.set_data("topology", topology1)
        edge.set_data("topology2", topology2)
        edge.dual_atom_mapping = mapping
    return topologies

def assign_dual_topology_old(edges, num_procs=1):
    light_edges = [LightEdge(edge) for edge in edges]
    topologies = Parallel(n_jobs=8)(delayed(build_dual_topology_old)(edge) for edge in light_edges)
    for edge, (topology1, topology2, mapping) in zip(edges, topologies):
        if topology1.name != edge.name:
            raise RuntimeError("Dual topology name does not math the edge name")
        edge.set_data("topology", topology1)
        edge.set_data("topology2", topology2)
        edge.dual_atom_mapping = mapping

def build_dual_topology(edge,idx=None):
    mol1, mol2 = edge.structs
    atom_mapping = edge.atom_mapping
    mol_name1, mol_name2 = atom_mapping.keys()
    if mol1.name != mol_name1 and mol1.name == mol_name2:
        mol2_index, mol1_index = atom_mapping.values()
    else:
        mol1_index, mol2_index = atom_mapping.values()
    fep_topology = FEPTopology(mol1, mol2, {i: j for i, j in zip(mol1_index, mol2_index)}, is_core_hopping=edge.is_core_hopping,is_charge_hopping=edge.is_charge_hopping)
    #fep_topology = FEPTopology(mol1, mol2, {i: j for i, j in zip(mol1_index, mol2_index)}, is_core_hopping=False,is_charge_hopping=edge.is_charge_hopping)
    if idx is not None:
        return fep_topology.dual_topology(), idx
    return fep_topology.dual_topology()

def dual_topolgy_assign(gg,parallel=True):
    edges = list(gg.edges_iter())
    light_edges = [LightEdge(edge) for edge in edges]
    if parallel:
        dual_molecules = parallel_run(build_dual_topology,light_edges,keep_order=True,return_result=True)
    else:
        dual_molecules = []
        for edge in light_edges:
            dual_molecules.append(build_dual_topology(edge))
    return dual_molecules