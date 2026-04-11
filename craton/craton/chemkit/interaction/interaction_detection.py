import itertools


from collections import defaultdict, namedtuple
from dataclasses import dataclass, fields
from typing import Dict, List, Tuple

import numpy as np


from .config import *
from ..structure.model import Site

from .util import euclidean3d, projection, vecangle

from .interaction_model import Interaction, AllInteraction


def hydrophobic_interactions(site_a: List[Site], site_b: List[Site]) -> List[Interaction]:
    pairings = []
    for an, bn in itertools.product(site_a, site_b):
        #a = an.atom
        #b = bn.atom
        d = np.linalg.norm(an.center - bn.center)
        if not MIN_DIST < d < HYDROPH_DIST_MAX:
            continue
        pairings.append(Interaction(type="Hydrophobic",donor=bn,acceptor=an,distance=d))
    return pairings


#def hbonds(acceptors: List[HBondAcceptor], donor_pairs: List[HBondDonor]) -> List[HBonds]:
def hbonds(acceptors: List[Site], donors: List[Site]) -> List[Interaction]:
    triplets = []
    for acc, don in itertools.product(acceptors, donors):
        dist_ad = np.linalg.norm(acc.center - don.center)
        vec1, vec2 = don.center - don.h_center, acc.center - don.h_center
        angle = vecangle(vec1, vec2)
        #if don.type == "regular":
        #    if not MIN_DIST < dist_ad < HBOND_DIST_MAX:
        #        continue
        #    if not angle > HBOND_DON_ANGLE_MIN:
        #        continue
        #else:
        #    if not MIN_DIST < dist_ad < WEAK_HBOND_DIST_MAX:
        #        continue
        #    if not angle > WEAK_HBOND_DON_ANGLE_MIN:
        #        continue
        if not MIN_DIST < dist_ad < HBOND_DIST_MAX:
            continue
        if not angle > HBOND_DON_ANGLE_MIN:
            continue
        triplets.append(Interaction(type="HBonds",acceptor=acc,donor=don,distance=dist_ad,angle=angle))
        #triplets.append(HBonds(d=don.d, h=don.h, a=acc.a, angle=angle, distance=dist_ad, type=don.type))
    return triplets

def saltbridge(pos_data: List[Site], neg_data: List[Site], prot_is_pos: bool) -> List[Interaction]:
    pairing = []
    for pos, neg in itertools.product(pos_data, neg_data):
        distance = np.linalg.norm(pos.center - neg.center)
        if not MIN_DIST < distance < SALTBRIDGE_DIST_MAX:
            continue
        pairing.append(Interaction(type="SaltBridge",acceptor=pos,donor=neg,distance=distance))
        #pairing.append(SaltBridge(pos, neg, distance, prot_is_pos))
    return pairing


def pistacking(ring_bs: List[Site], ring_lig: List[Site]) -> List[Interaction]:
    pairing = []
    for r, l in itertools.product(ring_bs, ring_lig):
        distance = euclidean3d(r.center, l.center)
        angle = vecangle(r.normal, l.normal)
        angle = min(angle, 180 - angle)

        proj1 = projection(l.normal, l.center, r.center)
        proj2 = projection(r.normal, r.center, l.center)
        if proj1 is not None and proj2 is not None:
            offset = min(euclidean3d(proj1, l.center), euclidean3d(proj2, r.center))
        else:
            continue
        if not MIN_DIST < distance < PISTACK_DIST_MAX:
            continue
        passed = False
        if 0 < angle < PISTACK_ANG_DEV and offset < PISTACK_OFFSET_MAX:
            ptype = "P"
            passed = True
        if 90 - PISTACK_ANG_DEV < angle < 90 + PISTACK_ANG_DEV and offset < PISTACK_OFFSET_MAX:
            ptype = "T"
            passed = True
        if passed:
            pairing.append(Interaction(type="PiStacking",
                                       acceptor=r,
                                       donor=l,
                                       distance=distance,
                                       angle=angle,
                                       offset=offset,
                                       subtype=ptype)
                            )
            #pairing.append(PiStacking(r, l, distance, angle, offset, ptype))
    return pairing


def pication(rings: List[Site], pos_charged: List[Site], prot_charged: bool) -> List[Interaction]:
    pairings = []
    for ring in rings:
        for p in pos_charged:
            d = euclidean3d(ring.center, p.center)
            proj = projection(ring.normal, ring.center, p.center)
            offset = euclidean3d(proj, ring.center)
            if not MIN_DIST < d < PICATION_DIST_MAX or not offset < PISTACK_OFFSET_MAX:
                continue
            pairings.append(Interaction(type="PiCation",acceptor=ring,donor=p,distance=d,offset=offset))
            #pairings.append(PiCation(ring, p, d, offset, prot_charged))
    return pairings


#def chpi(CH: List[HBondDonor], rings: List[Ring], ring_in_protein: bool) -> List[CHPi]:
def chpi(CH: List[Site], rings: List[Site], ring_in_protein: bool) -> List[Interaction]:
    pairings = []
    for donor, ring in itertools.product(CH, rings):
        if donor.type == "regular":
            continue
        proj = projection(ring.normal, ring.center, donor.h_center)
        d = euclidean3d(donor.h_center, proj)
        offset = euclidean3d(proj, ring.center)
        if d > CH_PI_DIST_MAX:
            continue
        else:
            if offset < CH_PI_CENTRAL_DIS:
                pairings.append(Interaction(
                                            type="chpi",
                                            acceptor=ring,
                                            donor=donor,
                                            distance=d,
                                            offset=offset,
                                            subtype="central",
                                            ring_in_protein=ring_in_protein
                                            )
                                            )
                #pairings.append(
                #    CHPi(
                #        ring=ring, ch=donor, distance=d, offset=offset, type="central", ring_in_protein=ring_in_protein
                #    )
                #)
            elif offset < CH_PI_PERIPHERAL_DIS:
                pairings.append(Interaction(
                                            type="chpi",
                                            acceptor=ring,
                                            donor=donor,
                                            distance=d,
                                            offset=offset,
                                            subtype="peripheral",
                                            ring_in_protein=ring_in_protein
                                            )
                                )
                #pairings.append(
                #    CHPi(
                ##        ring=ring,
                 #       ch=donor,
                 #       distance=d,
                 #       offset=offset,
                 #       type="peripheral",
                 #       ring_in_protein=ring_in_protein,
                 #   )
                #)
    return pairings


#def halogen(acceptor: List[HalAcceptor], donor: List[HalDonor]) -> List[Halogen]:
def halogen(acceptor: List[Site], donor: List[Site]) -> List[Interaction]:
    """Detect all halogen bonds of the type Y-O...X-C"""
    pairings = []
    for acc, don in itertools.product(acceptor, donor):
        #dist = euclidean3d(acc.o.coor, don.x.coor)
        dist = euclidean3d(acc.h_center, don.h_center)
        if not MIN_DIST < dist < HALOGEN_DIST_MAX:
            continue
        #vec1, vec2 = acc.y.coor - acc.o.coor, don.x.coor - acc.o.coor
        #vec3, vec4 = acc.o.coor - don.x.coor, don.c.coor - don.x.coor
        vec1, vec2 = acc.center - acc.h_center, don.h_center - acc.h_center
        vec3, vec4 = acc.h_center - don.h_center, don.center - don.h_center
        acc_angle, don_angle = vecangle(vec1, vec2), vecangle(vec3, vec4)
        if (
            not HALOGEN_ACC_ANGLE - HALOGEN_ANGLE_DEV
            < acc_angle
            < HALOGEN_ACC_ANGLE + HALOGEN_ANGLE_DEV
        ):
            continue
        if (
            not HALOGEN_DON_ANGLE - HALOGEN_ANGLE_DEV
            < don_angle
            < HALOGEN_DON_ANGLE + HALOGEN_ANGLE_DEV
        ):
            continue
        pairings.append(Interaction(type="Halogen",acceptor=acc,donor=don,distance=dist,angle=acc_angle,auxi_angle=don_angle))
        #pairings.append(Halogen(acc, don, dist, acc_angle, don_angle))
    return pairings


#def refine_hbonds_ldon(
#    all_hbonds: List[HBonds], salt_lneg: List[SaltBridge], salt_pneg: List[SaltBridge]
#) -> List[HBonds]:
def refine_hbonds_ldon(
    all_hbonds: List[Interaction], salt_lneg, salt_pneg
) -> List[Interaction]:
    """Refine selection of hydrogen bonds. Do not allow groups which already
    form salt bridges to form H-Bonds.postitive is acceptor, negative is donor
    atom_orig_idx is ID, hbond.d is donor.adj_atom, hbond.h is donor.atom, hbond.a is acceptor.atom"""
    removed_hbond = set()
    for i, hbond in enumerate(all_hbonds):
        if hbond.donor.adj_atom.ID in salt_pneg[1] and hbond.acceptor.atom.ID in salt_pneg[0]:
            removed_hbond.add(i)
        
        if hbond.donor.adj_atom.ID in salt_lneg[0] and hbond.acceptor.atom.ID in salt_lneg[1]:
            removed_hbond.add(i)

    # Allow only one hydrogen bond per donor, select interaction with larger donor angle
    second_set: Dict[int, Tuple[float, Interaction]] = {}
    hbls = [hbond for i, hbond in enumerate(all_hbonds) if i not in removed_hbond]
    for hbl in hbls:
        if hbl.donor.adj_atom.ID not in second_set:
            #second_set[hbl.donor.adj_atom.ID] = (hbl.angle, hbl)
            second_set[hbl.donor.adj_atom.ID] = (hbl.distance, hbl)
        else:
            #if second_set[hbl.donor.adj_atom.ID][0] < hbl.angle:
                #second_set[hbl.donor.adj_atom.ID] = (hbl.angle, hbl)
            if second_set[hbl.donor.adj_atom.ID][0] > hbl.distance:
                second_set[hbl.donor.adj_atom.ID] = (hbl.distance, hbl)
    return [hb[1] for hb in second_set.values()]


#def refine_hbonds_pdon(
#    all_hbonds: List[HBonds], salt_lneg: List[SaltBridge], salt_pneg: List[SaltBridge]
#) -> List[HBonds]:
def refine_hbonds_pdon(
    all_hbonds: List[Interaction], salt_lneg: List[Interaction], salt_pneg: List[Interaction]
) -> List[Interaction]:
    """Refine selection of hydrogen bonds. Do not allow groups which already form salt bridges to form H-Bonds with
    atoms of the same group.
    """
    removed_hbond = set()
    for i, hbond in enumerate(all_hbonds):
        if hbond.acceptor.atom.ID in salt_lneg[0] and hbond.donor.adj_atom.ID in salt_lneg[1]:
            removed_hbond.add(i)
        if hbond.acceptor.atom.ID in salt_pneg[1] and hbond.donor.adj_atom.ID in salt_pneg[0]:
            removed_hbond.add(i)
        

    # Allow only one hydrogen bond per donor, select interaction with larger donor angle
    second_set: Dict[int, Tuple[float, Interaction]] = {}
    hbps = [hbond for i, hbond in enumerate(all_hbonds) if i not in removed_hbond]
    for hbp in hbps:
        if hbp.donor.adj_atom.ID not in second_set:
            #second_set[hbp.donor.adj_atom.ID] = (hbp.angle, hbp)
            second_set[hbp.donor.adj_atom.ID] = (hbp.distance, hbp)
        else:
            #if second_set[hbp.donor.adj_atom.ID][0] < hbp.angle:
                #second_set[hbp.donor.adj_atom.ID] = (hbp.angle, hbp)
            if second_set[hbp.donor.adj_atom.ID][0] > hbp.distance:
                second_set[hbp.donor.adj_atom.ID] = (hbp.distance, hbp)
    return [hb[1] for hb in second_set.values()]


#def refine_pi_cation_laro(all_picat: List[PiCation], stacks: List[PiStacking]) -> List[PiCation]:
def refine_pi_cation_laro(all_picat: List[Interaction], stacks: List[Interaction]) -> List[Interaction]:
    """Just important for constellations with histidine involved. If the histidine ring is positioned in stacking
    position to an aromatic ring in the ligand, there is in most cases stacking and pi-cation interaction reported
    as histidine also carries a positive charge in the ring. For such cases, only report stacking.
    """
    i_set: List[Interaction] = []
    for picat in all_picat:
        exclude = False
        for stack in stacks:
            if stack.acceptor.atoms[0].residue == "HIS" and picat.acceptor.atoms == stack.donor.atoms:
                exclude = True
        if not exclude:
            i_set.append(picat)
    return i_set


#def refine_hydrophobic(all_h: List[Hydrophobic], pistacks: List[PiStacking], ch_pis: List[CHPi]) -> List[Hydrophobic]:
def refine_hydrophobic(all_h: List[Interaction], pistacks: List[Interaction], ch_pis: List[Interaction]) -> List[Interaction]:
    """Apply several rules to reduce the number of hydrophobic interactions."""
    exclude = set()
    #  1. Rings interacting via stacking can't have additional hydrophobic contacts between each other.
    #   or via ch-pi interaction
    for pistack, h in itertools.product(pistacks, all_h):
        h1, h2 = h.acceptor.atom.ID, h.donor.atom.ID
        brs, lrs = (
            [p1.ID for p1 in pistack.acceptor.atoms],
            [p2.ID for p2 in pistack.donor.atoms],
        )
        if h1 in brs and h2 in lrs:
            exclude.add((h1, h2))
    for ch_pi, hydrophobic in itertools.product(ch_pis, all_h):
        bs_atom, ligand_atom = hydrophobic.acceptor.atom.ID, hydrophobic.donor.atom.ID
        ch_pi_atom1, ch_pi_atoms = ch_pi.donor.atom.ID, [a.ID for a in ch_pi.acceptor.atoms]
        if ch_pi.ring_in_protein:
            if bs_atom in ch_pi_atoms and ligand_atom == ch_pi_atom1:
                exclude.add((bs_atom, ligand_atom))
        else:
            if bs_atom == ch_pi_atom1 and ligand_atom in ch_pi_atoms:
                exclude.add((bs_atom, ligand_atom))

    hydroph: List[Interaction] = [h for h in all_h if not (h.acceptor.atom.ID, h.donor.atom.ID) in exclude]

    sel2: Dict[Tuple[int, int], Interaction] = {}
    #  2. If a ligand atom interacts with several binding site atoms in the same residue,
    #  keep only the one with the closest distance
    for h in hydroph:
        if not (h.donor.atom.ID, h.acceptor.atom.ID) in sel2:
            sel2[(h.donor.atom.ID, h.acceptor.atom.ID)] = h
        else:
            if sel2[(h.donor.atom.ID, h.acceptor.atom.ID)].distance > h.distance:
                sel2[(h.donor.atom.ID, h.acceptor.atom.ID)] = h
    hydroph: List[Interaction] = [h for h in sel2.values()]

    #  3. If a protein atom interacts with several neighboring ligand atoms, just keep the one with the closest dist

    hydroph_final = []
    bsclust: Dict[int, List[Interaction]] = defaultdict(list)

    for h in hydroph:
        bsclust[h.acceptor.atom.ID].append(h)

    for _, hydros in bsclust.items():
        if len(hydros) == 1:
            hydroph_final.append(hydros[0])
        else:
            hydroph_final.append(sorted(hydros)[0])

    return hydroph_final


def interaction_detection(
    ligand_sites,
    binding_sites,
    coordinate = None,
    water_sites = None,
    metal=False,
):
    all_hbonds_ligand_donor = hbonds(binding_sites["hbond_acceptor"], ligand_sites["hbond_donor"])
    all_hbonds_protein_donor = hbonds(ligand_sites["hbond_acceptor"], binding_sites["hbond_donor"])

    protein_positive_charged = [charge for charge in binding_sites["charged"] if charge.subtype == "positive"]
    protein_negative_charged = [charge for charge in binding_sites["charged"] if charge.subtype == "negative"]
    ligand_negative_charged = [charge for charge in ligand_sites["charged"] if charge.subtype == "negative"]
    ligand_poisitve_charged = [charge for charge in ligand_sites["charged"] if charge.subtype == "positive"]
    saltbridge_ligand_negative = saltbridge(protein_positive_charged, ligand_negative_charged, True)
    saltbridge_protein_negative = saltbridge(ligand_poisitve_charged, protein_negative_charged, False)

    pneg_neg_atoms = [at.ID for salt in saltbridge_protein_negative for at in salt.donor.atoms]
    pneg_pos_atoms = [at.ID for salt in saltbridge_protein_negative for at in salt.acceptor.atoms]

    lneg_neg_atoms = [at.ID for salt in saltbridge_ligand_negative for at in salt.donor.atoms]
    lneg_pos_atoms = [at.ID for salt in saltbridge_ligand_negative for at in salt.acceptor.atoms]


    #hbonds_ldon = refine_hbonds_ldon(all_hbonds_ligand_donor, saltbridge_ligand_negative, saltbridge_protein_negative)
    #hbonds_pdon = refine_hbonds_pdon(all_hbonds_protein_donor, saltbridge_ligand_negative, saltbridge_protein_negative)

    hbonds_ldon = refine_hbonds_ldon(all_hbonds_ligand_donor, [lneg_neg_atoms,lneg_pos_atoms], [pneg_neg_atoms,pneg_pos_atoms])
    hbonds_pdon = refine_hbonds_pdon(all_hbonds_protein_donor, [lneg_neg_atoms,lneg_pos_atoms], [pneg_neg_atoms,pneg_pos_atoms])
    pi_stacking = pistacking(binding_sites["rings"], ligand_sites["rings"])
    
    all_pi_cation_laro = pication(ligand_sites["rings"], protein_positive_charged, True)
    pication_paro = pication(binding_sites["rings"], ligand_poisitve_charged, False)
    pication_laro = refine_pi_cation_laro(all_pi_cation_laro, pi_stacking)

    chpi_paro = chpi(ligand_sites["hbond_donor"], binding_sites["rings"], ring_in_protein=True)
    chpi_laro = chpi(binding_sites["hbond_donor"], ligand_sites["rings"], ring_in_protein=False)

    all_hydrophobic_contacts = hydrophobic_interactions(binding_sites["hydrophobic"], ligand_sites["hydrophobic"])
    hydrophobic_contacts = refine_hydrophobic(all_hydrophobic_contacts, pi_stacking, chpi_laro + chpi_paro)
    halogen_bonds = halogen(binding_sites["halogen_acceptor"], ligand_sites["halogen_donor"])

    if water_sites is not None:
        water_bridges = waterbridges(
            binding_sites["hbond_acceptor"],
            ligand_sites["hbond_acceptor"],
            binding_sites["hbond_donor"],
            ligand_sites["hbond_donor"],
            water_sites["binding_waters"],
        )
        water_bridges = refine_water_bridges(water_bridges, hbonds_ldon, hbonds_pdon)
    else:
        water_bridges = []
    # if metal:
    #     metal_complexation = metalcomplexation(
    #         self.ions,
    #         self.ligand.metal_binding,
    #         self.protein.metal_binding,
    #         self.water.water_metal_binding,
    #     )
    # else:
    #     metal_bindingl_bindingl_complexation = []

    return AllInteraction(
        saltbridge_lneg=saltbridge_ligand_negative,
        saltbridge_pneg=saltbridge_protein_negative,
        hbonds_ldon=hbonds_ldon,
        hbonds_pdon=hbonds_pdon,
        pi_stacking=pi_stacking,
        pication_paro=pication_paro,
        pication_laro=pication_laro,
        chpi_paro=chpi_paro,
        chpi_laro=chpi_laro,
        hydrophobic_contacts=hydrophobic_contacts,
        halogen_bonds=halogen_bonds,
        water_bridges=water_bridges,
        metal_complexation=[],
    )


#from .report import ReportMultipleInteraction
#from .visualize_by_plotly import fep_bar_plot_with_energy, fep_bar_plot

#def merge_interaction_dfs(df1: pd.DataFrame, df2: pd.DataFrame):
def merge_interaction_dfs(df1, df2):
    import pandas as pd
    index_all = df1.index.union(df2.index)
    interaction_all = df1.columns.union(df2.columns)

    for index in index_all:
        if index not in df1.index:
            interaction_dict = {interaction: 0 for interaction in interaction_all}
            df1 = df1._append(pd.Series(interaction_dict, name=index))
        if index not in df2.index:
            interaction_dict = {interaction: 0 for interaction in interaction_all}
            df2 = df2._append(pd.Series(interaction_dict, name=index))

    interaction_remain = interaction_all.difference(df1.columns)
    for interaction in interaction_remain:
        df1[interaction] = 0

    interaction_remain = interaction_all.difference(df2.columns)
    for interaction in interaction_remain:
        df2[interaction] = 0

    return df1.fillna(0).sort_index(
        key=lambda col: np.array([int(item.split("_")[0]) for item in col]), ascending=False
    ), df2.fillna(0).sort_index(key=lambda col: np.array([int(item.split("_")[0]) for item in col]), ascending=False)


def interaction_report(interactions_a,interactions_b,output_dir,calc_energy_flag=False):
    report1 = ReportMultipleInteraction(interactions_a)
    report2 = ReportMultipleInteraction(interactions_b)
    df1 = report1.info_df(show_figure=False)
    df2 = report2.info_df(show_figure=False)
    report1.info_detail(output_dir / "interaction_detail_1.txt")
    report2.info_detail(output_dir / "interaction_detail_2.txt")
    analyze_data = {}
    analyze_data['interaction_df1'], analyze_data['interaction_df2'] = merge_interaction_dfs(df1, df2)
    if calc_energy_flag:
        energy_analyze = EnergyFromGmxData(ana_data_handler.gro_a.parent, logs=output_dir)
        analyze_data['energy_df1'] = energy_analyze.binding_site_energy(
            ana_data_handler.traj_a, interaction_df=analyze_data['interaction_df1']
        )
        energy_analyze = EnergyFromGmxData(ana_data_handler.gro_b.parent)
        analyze_data['energy_df2'] = energy_analyze.binding_site_energy(
            ana_data_handler.traj_b, interaction_df=analyze_data['interaction_df2'], use_top_b=True
        )
        
    return analyze_data

def save_results(analyze_data,output_dir,calc_energy_flag=False):
    
    pair_name = output_dir.name
    dump_data = {}
    if calc_energy_flag:
        fep_bar_plot_with_energy(
            analyze_data['interaction_df1'],
            analyze_data['interaction_df2'],
            analyze_data['energy_df1'],
            analyze_data['energy_df2'],
            pair_name,
            output_dir / "interaction_result.html")
    else:
        fep_bar_plot(
            analyze_data['interaction_df1'], analyze_data['interaction_df2'],
            pair_name, output_dir / "interaction_result.html"
        )
    # fep_bar_plot(df1, df2, pair_name, self.output_dir / "interaction_result.html")
    if calc_energy_flag:
        dump_data["A"] = {
            "interaction": analyze_data['interaction_df1'].to_dict(orient="index"),
            "energy": analyze_data['energy_df1'].to_dict(orient="index"),
        }
        dump_data["B"] = {
            "interaction": analyze_data['interaction_df2'].to_dict(orient="index"),
            "energy": analyze_data['energy_df2'].to_dict(orient="index"),
        }
    else:
        dump_data["A"] = {"interaction": analyze_data['interaction_df1'].to_dict(orient="index")}
        dump_data["B"] = {"interaction": analyze_data['interaction_df2'].to_dict(orient="index")}