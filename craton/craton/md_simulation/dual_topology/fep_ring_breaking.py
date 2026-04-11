import copy
from enum import Enum
from dataclasses import dataclass, field
from ...utils import logger


from ...chemkit import MolEdit as ME
from ...chem.topology import Pair

class FEPProperty(Enum):
    COMMON = 1
    GROW = 2
    DISAPPEAR = 3
    RESTRAIN = 4


@dataclass
class PairInteraction:
    all_14_pair: set
    a14_disappear: set
    a14_grow: set
    exclusion_pair: set
    a1n_disappear: set
    a1n_grow: set
    restrain_pair: set


@dataclass
class BondParaTransformation:
    bonds_all: set
    bonds_grow: set
    bonds_disappear: set
    angles_all: set
    angles_grow: set
    angles_disappear: set
    dihedrals_all: set
    dihedrals_grow: set
    dihedrals_disappear: set
    impropers_all: set
    impropers_grow: set
    impropers_disappear: set


@dataclass
class DebugInfo:
    mol_a_bonds: field(default_factory=list) = None
    mol_a_angles: field(default_factory=list) = None
    mol_a_dihedral: field(default_factory=list) = None
    mol_a_improper: field(default_factory=list) = None

    mol_b_bonds: field(default_factory=list) = None
    mol_b_angles: field(default_factory=list) = None
    mol_b_dihedral: field(default_factory=list) = None
    mol_b_improper: field(default_factory=list) = None

    mol_a_12: field(default_factory=list) = None
    mol_a_13: field(default_factory=list) = None
    mol_a_14: field(default_factory=list) = None
    mol_a_1n: field(default_factory=list) = None

    mol_b_12: field(default_factory=list) = None
    mol_b_13: field(default_factory=list) = None
    mol_b_14: field(default_factory=list) = None
    mol_b_1n: field(default_factory=list) = None

    all_14_pair: set = None
    a14_disappear: set = None
    a14_grow: set = None
    exclusion_pair: set = None
    a1n_disappear: set = None
    a1n_grow: set = None

    nonredundant_all_14_pair: set = None
    nonredundant_a14_disappear: set = None
    nonredundant_a14_grow: set = None
    nonredundant_exclusion_pair: set = None
    nonredundant_a1n_disappear: set = None
    nonredundant_a1n_grow: set = None

    bonds_all: set = None
    bonds_grow: set = None
    bonds_disappear: set = None
    angles_all: set = None
    angles_grow: set = None
    angles_disappear: set = None
    dihedrals_all: set = None
    dihedrals_grow: set = None
    dihedrals_disappear: set = None
    impropers_all: set = None
    impropers_grow: set = None
    impropers_disappear: set = None

debug_info = DebugInfo()


def set_debug_info(debug_info, mol_a_topo_with_dummy, mol_b_topo_with_dummy):
    debug_info.mol_a_bonds = mol_a_topo_with_dummy.Bonds
    debug_info.mol_a_angles = mol_a_topo_with_dummy.Angles
    debug_info.mol_a_dihedral = mol_a_topo_with_dummy.Dihedrals
    debug_info.mol_a_improper = mol_a_topo_with_dummy.Impropers
    debug_info.mol_b_bonds = mol_b_topo_with_dummy.Bonds
    debug_info.mol_b_angles = mol_b_topo_with_dummy.Angles
    debug_info.mol_b_dihedral = mol_b_topo_with_dummy.Dihedrals
    debug_info.mol_b_improper = mol_b_topo_with_dummy.Impropers


def construct_broken_structure(mol_a,
                               mol_b,
                               mol_a_bonds_formation,
                               mol_a_bonds_delete,
                               mol_b_bonds_formation,
                               mol_b_bonds_delete):
    """mol_a_bonds_formation: the bond need formation to reach mol_b; (exist in mol_b, need delete for top)
    mol_a_bonds_delete: the bond need to be deleted to reach mol_b: (exist in mol_a, need keep for top)
    mol_b_bonds_formation: the bond need formation to reach mol_a; (exist in mol_a, need keep for top)
    mol_b_bonds_delete: the bond need to be deleted to reach mol_a; (exist in mol_b, need delete for top)

    if mol_a, mol_b are a five-membered ring  and six-membered ring respectively, then the mol_a_bond_formation
    and mol_b_bond delete are emtpy, we should modify the topology according to the mol_b_bond_formation
    and mol_b_bond_delete.

    if mol_a, mol_b are a six-membered ring  and five-membered ring respectively, then the mol_b_bonds_formation and
    mol_b_bonds_delete are empty, we should modify the topology according to the mol_a_bonds_formation
    and mol_a_bonds_delete.
    """

    mol_a_topo_with_dummy = copy.deepcopy(mol_a)
    mol_b_topo_with_dummy = copy.deepcopy(mol_a)

    for bond in mol_a_bonds_formation + mol_b_bonds_delete:
        mol_a_topo_with_dummy.remove_bond(bond[0], bond[1])

    for bond in mol_a_bonds_delete + mol_b_bonds_formation:
        mol_b_topo_with_dummy.remove_bond(bond[0], bond[1])

    mol_a_topo_with_dummy.create_topols()
    mol_a_topo_with_dummy.create_improper(create_method="mix")

    mol_b_topo_with_dummy.create_topols()
    mol_b_topo_with_dummy.create_improper(create_method="mix")


    set_debug_info(debug_info, mol_a_topo_with_dummy, mol_b_topo_with_dummy)

    for ff_term in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
        terms = getattr(mol_a_topo_with_dummy, ff_term, None)
        if terms:
            mol_a_not_broken = {term: term.parameter for term in getattr(mol_a, ff_term)}
            for term in terms:
                term.parameter = mol_a_not_broken.get(term, "no_parameter")

        terms = getattr(mol_b_topo_with_dummy, ff_term, None)
        if terms:
            mol_b_not_broken = {term: term.parameter for term in getattr(mol_b, ff_term)}
            for term in terms:
                term.parameter = mol_b_not_broken.get(term, "no_parameter")

    return mol_a_topo_with_dummy, mol_b_topo_with_dummy


def remove_dummy_pair(pairs, dummy_a, dummy_b):
    temp_pair = set()
    for pair in pairs:
        if set(pair).intersection(dummy_a) and set(pair).intersection(dummy_b):
            temp_pair.add(pair)
    pairs.difference_update(temp_pair)


def remove_disappear_pair(pairs, dummy_a):
    temp_pair = set()
    for pair in pairs:
        if set(pair).intersection(dummy_a):
            temp_pair.add(pair)
    pairs.difference_update(temp_pair)
    return temp_pair


def nonbond_parameters(mol_a_with_dummy, mol_b_with_dummy, dummy_a, dummy_b):
    a12, a13, a14, a1n = mol_a_with_dummy.get_1234n_pairs()
    b12, b13, b14, b1n = mol_b_with_dummy.get_1234n_pairs()

    mol_a_with_dummy.Pair14 = [Pair("common", *pair) for pair in a14]
    mol_b_with_dummy.Pair14 = [Pair("common", *pair) for pair in a14]

    # gromacs.write_mole_itp(mol_a_with_dummy, file_name="mol_a.itp")
    # gromacs.write_mole_itp(mol_b_with_dummy, file_name="mol_b.itp")


    debug_info.mol_a_12, debug_info.mol_a_13, debug_info.mol_a_14, debug_info.mol_a_1n = a12.copy(), a13.copy(), a14.copy(), a1n.copy()
    debug_info.mol_b_12, debug_info.mol_b_13, debug_info.mol_b_14, debug_info.mol_b_1n = b12.copy(), b13.copy(), b14.copy(), b1n.copy()

    all_14_pair = set(a14).union(b14)
    a14_disappear = set.union(set(a14).intersection(b12), set(a14).intersection(b13), set(a14).intersection(b1n))
    a14_grow = set.union(set(b14).intersection(a12), set(b14).intersection(a13), set(b14).intersection(a1n))

    debug_info.all_14_pair, debug_info.a14_disappear, debug_info.a14_grow = all_14_pair.copy(), a14_disappear.copy(), a14_grow.copy()


    remove_dummy_pair(all_14_pair, dummy_a, dummy_b)
    remove_dummy_pair(a14_disappear, dummy_a, dummy_b)
    remove_dummy_pair(a14_grow, dummy_a, dummy_b)

    # a14_disappear
    # if the pair has any dummy_b atom(disappear), add them in pairs, let the program handle this automatically
    remove_disappear_pair(a14_disappear, dummy_b)
    # if the pair has any dummy_a atom(grow), just ignore them in the pairs
    a14_pair_contain_dummy_a_atom = remove_disappear_pair(a14_disappear, dummy_a)
    all_14_pair.difference_update(a14_pair_contain_dummy_a_atom)
    # if both of pair are real atom, we must need them to disappear

    # a14_grow
    # if the pair has any dummy_b atom(disappear), ignore them in the pairs
    removed_pair = remove_disappear_pair(a14_grow, dummy_b)
    all_14_pair.difference_update(removed_pair)
    # if the pair has any dummy_a atom(grow), add them in pairs, let the program handle this automatically
    remove_disappear_pair(a14_grow, dummy_a)
    # if both of them are real atom, we must need them to grow

    if len(a14_disappear) != 0 or len(a14_grow) != 0:
        logger.debug("a14_disappear, a14_grow: {}, {}".format(a14_disappear, a14_grow))

    a1n_disappear = set.union(set(a1n).intersection(b12), set(a1n).intersection(b13), set(a1n).intersection(b14))
    a1n_grow = set.union(set(b1n).intersection(a12), set(b1n).intersection(a13), set(b1n).intersection(a14))

    debug_info.a1n_grow, debug_info.a1n_disappear = a1n_grow.copy(), a1n_disappear.copy()

    # could keep the a1n disappear in the exclusion list
    exclusion_pair = a1n_disappear.copy().union(a1n_grow.copy())

    # if pair1n grow has any dummy_b atom(disappear), just ignore them
    remove_disappear_pair(a1n_grow, dummy_b)

    # if pair1n disappear has any dummy_a atom(grow), add them in exclusion list
    remove_disappear_pair(a1n_disappear, dummy_a)

    # the a1n_grow, and a1n_disappear should have zero length, if not have, logger an error
    if len(a1n_disappear) != 0 or len(a1n_grow) != 0:
        logger.debug(f"The a1n_disappear {a1n_disappear} or a1n_grow {a1n_grow} is not zero.")

    restrain_pair = set()
    b1n = set(b1n).difference(exclusion_pair)
    restrain_pair = set()
    for pair in b1n:
        temp_pair = list(pair)
        if temp_pair[0] in dummy_a and temp_pair[1] in dummy_a:
            distance = ME._structure_calculate([mol_a_with_dummy.Atoms[pair[0]].coor, mol_b_with_dummy.Atoms[pair[1]].coor])
            if distance < 3:
                restrain_pair.add(pair)
    exclusion_pair.update(restrain_pair)

    remove_dummy_pair(a1n_disappear, dummy_a, dummy_b)
    remove_dummy_pair(a1n_disappear, dummy_a, dummy_b)

    debug_info.nonredundant_all_14_pair, debug_info.nonredundant_a14_disappear, debug_info.nonredundant_a14_grow = all_14_pair.copy(), a14_disappear.copy(), a14_grow.copy()
    debug_info.nonredundant_exclusion_pair, debug_info.nonredundant_a1n_disappear, debug_info.nonredundant_a1n_grow = exclusion_pair.copy(), a1n_disappear.copy(), a1n_grow.copy()

    return PairInteraction(all_14_pair=all_14_pair, a14_disappear=a14_disappear, a14_grow=a14_grow,
                           exclusion_pair=exclusion_pair, a1n_disappear=a1n_disappear, a1n_grow=a1n_grow,
                           restrain_pair=restrain_pair)


def transformation_for_bond_topology(mol_a_with_dummy, mol_b_with_dummy, dummy_a, dummy_b):
    transformation_dict = {}
    for topo_type in ["Bonds", "Angles", "Dihedrals", "Impropers"]:
        terms_a = set(getattr(mol_a_with_dummy, topo_type, []))
        terms_b = set(getattr(mol_b_with_dummy, topo_type, []))
        remove_dummy_pair(terms_a, dummy_a, dummy_b)
        remove_dummy_pair(terms_b, dummy_a, dummy_b)
        transformation_dict[topo_type.lower() + "_all"] = terms_a.union(terms_b)
        transformation_dict[topo_type.lower() + "_grow"] = terms_b - terms_a
        transformation_dict[topo_type.lower() + "_disappear"] = terms_a - terms_b
        setattr(debug_info, topo_type.lower() + '_all', terms_a.union(terms_b))
        setattr(debug_info, topo_type.lower() + '_grow', terms_b - terms_a)
        setattr(debug_info, topo_type.lower() + '_disappear', terms_a - terms_b)
    return BondParaTransformation(**transformation_dict)


def bond_parameters(mol_a_with_dummy, mol_b_with_dummy, dummy_a, dummy_b):
    bond_topo = transformation_for_bond_topology(mol_a_with_dummy, mol_b_with_dummy, dummy_a, dummy_b)

    mol_a_bond = {bond: getattr(bond, 'parameter', 'undefined') for bond in mol_a_with_dummy.Bonds}
    mol_b_bond = {bond: getattr(bond, 'parameter', 'undefined') for bond in mol_b_with_dummy.Bonds}

    for bond in bond_topo.bonds_all:
        # For taylor expansion
        # https://www.symbolab.com/solver/taylor-series-calculator/taylor%20%5Cleft(1-e%5E%7B-x%7D%5Cright)%5E%7B2%7D?or=input
        # D(1-e^(-bx))^2 = D(bx)^2 - D(bx)^3 + ....
        # Db^2 = k
        if bond in bond_topo.bonds_disappear:
            bond.style = "morse"
            bond.fep_property = FEPProperty.DISAPPEAR
            harmonic_para = mol_a_bond[bond]
            bond.parameter = [harmonic_para[0], harmonic_para[1], 1]
            bond.parameter_m2 = [harmonic_para[0], 0, 1]
        elif bond in bond_topo.bonds_grow:
            bond.fep_property = FEPProperty.GROW
            bond.style = "morse"
            harmonic_para = mol_b_bond[bond]
            bond.parameter = [harmonic_para[0], 0, 1]
            bond.parameter_m2 = [harmonic_para[0], harmonic_para[1], 1]
        else:
            bond.parameter = mol_a_bond[bond][:]
            bond.fep_property = FEPProperty.COMMON
            bond.style = "harmonic"
            if bond.parameter == "no_parameter":
                bond.parameter = mol_b_bond[bond][:]
            bond.parameter_m2 = mol_b_bond[bond][:]
            if bond.parameter_m2 == "no_parameter":
                bond.parameter_m2 = mol_a_bond[bond][:]

    mol_a_angle = {angle: angle.parameter for angle in mol_a_with_dummy.Angles}
    mol_b_angle = {angle: angle.parameter for angle in mol_b_with_dummy.Angles}

    for angle in bond_topo.angles_all:
        angle.style = "harmonic"
        if angle in bond_topo.angles_disappear:
            angle.fep_property = FEPProperty.DISAPPEAR
            angle.parameter = mol_a_angle[angle][:]
            if angle.parameter == "no_parameter":
                angle.parameter = [0.0, 0.0]
            angle.parameter_m2 = [angle.parameter[0], 0.0]
        elif angle in bond_topo.angles_grow:
            angle.fep_property = FEPProperty.GROW
            angle.parameter_m2 = mol_b_angle[angle][:]
            if angle.parameter_m2 == "no_parameter":  # no parametermeters for angle grow
                angle.parameter_m2 = [0.0, 0.0]
            angle.parameter = [angle.parameter_m2[0], 0.0]
        else:
            angle.parameter = mol_a_angle[angle][:]
            angle.fep_property = FEPProperty.COMMON
            if angle.parameter == "no_parameter":
                angle.parameter = mol_b_angle[angle][:]
            angle.parameter_m2 = mol_b_angle[angle][:]
            if angle.parameter_m2 == "no_parameter":
                angle.parameter_m2 = mol_a_angle[angle][:]

    mol_a_dihedral = {dihedral: dihedral.parameter for dihedral in mol_a_with_dummy.Dihedrals}
    mol_b_dihedral = {dihedral: dihedral.parameter for dihedral in mol_b_with_dummy.Dihedrals}

    for dihedral in bond_topo.dihedrals_all:
        dihedral.style = "amber"
        if dihedral in bond_topo.dihedrals_disappear:
            dihedral.fep_property = FEPProperty.DISAPPEAR
            dihedral.parameter = mol_a_dihedral[dihedral][:]
            if dihedral.parameter == "no_parameter":
                dihedral.parameter = [0.0, 0.0, 0.0, 0.0]
            dihedral.parameter_m2 = [dihedral.parameter[0], 0.0, 0.0, 0.0]
        elif dihedral in bond_topo.dihedrals_grow:
            dihedral.fep_property = FEPProperty.GROW
            dihedral.parameter_m2 = mol_b_dihedral[dihedral][:]
            if dihedral.parameter_m2 == "no_parameter":
                dihedral.parameter_m2 = [0.0, 0.0, 0.0, 0.0]
            dihedral.parameter = [0.0, 0.0, 0.0, 0.0]
        else:
            dihedral.fep_property = FEPProperty.COMMON
            dihedral.parameter = mol_a_dihedral[dihedral][:]
            if dihedral.parameter == "no_parameter":
                dihedral.parameter = mol_b_dihedral[dihedral][:]
            dihedral.parameter_m2 = mol_b_dihedral[dihedral][:]
            if dihedral.parameter_m2 == "no_parameter":
                dihedral.parameter_m2 = mol_a_dihedral[dihedral][:]

    if hasattr(mol_a_with_dummy, "Impropers") or hasattr(mol_b_with_dummy, "Impropers"):
        mol_a_improper = {improper: improper.parameter for improper in mol_a_with_dummy.Impropers}
        mol_b_improper = {improper: improper.parameter for improper in mol_b_with_dummy.Impropers}
        for improper in bond_topo.impropers_all:
            improper.style = "amber"
            if improper in bond_topo.impropers_disappear:
                improper.fep_property = FEPProperty.DISAPPEAR
                improper.parameter = mol_a_improper[improper][:]
                if improper.parameter == "no_parameter":
                    improper.parameter = [0.0]
                improper.parameter_m2 = [0.0]
            elif improper in bond_topo.impropers_grow:
                improper.fep_property = FEPProperty.GROW
                improper.parameter_m2 = mol_b_improper[improper][:]
                if improper.parameter_m2 == "no_parameter":
                    improper.parameter_m2 = [0.0]
                improper.parameter = [0.0]
            else:
                improper.fep_property = FEPProperty.COMMON
                improper.parameter = mol_a_improper[improper][:]
                if improper.parameter == "no_parameter":
                    improper.parameter = mol_b_improper[improper][:]
                improper.parameter_m2 = mol_b_improper[improper][:]
                if improper.parameter_m2 == "no_parameter":
                    improper.parameter_m2 = mol_a_improper[improper][:]
    return bond_topo
