from collections import defaultdict, namedtuple
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Union

import numpy as np
import pandas as pd

from ...utils import logger
#from .find_atoms import Atom
from ...chem.atom import Atom
from .interaction_model import (
    AllInteraction,
    AllWaterLigandInteraction,
    WaterBridgeType,
)
from .visualize_by_plotly import (
    ligand_water_bar_plot,
    normal_md_bar_plot,
)
from .visualize_by_pymol import (
    CHPiInfoID,
    HalogenInfoID,
    HBondsInfoID,
    HydrophobicInfoID,
    MetalInfoID,
    PiCationInfoID,
    PiStackingInfoID,
    PymolInteractionData,
    SaltBridgeInfoID,
    SingleFrameViewer,
    WaterBridgeInfoID,
)


@dataclass
class InteractionDual:
    type: str
    ligand: List[Atom]
    protein: List[Atom]

    def __str__(self):
        protein_str, ligand_str = [], []
        for atom in self.protein:
            try:
                protein_str.append(f"{atom.residue}-{atom.residue_ID}-{atom.name}-{atom.ID}")
            except:
                protein_str.append(f"{atom.name}-{atom.ID}")
        for atom in self.ligand:
            ligand_str.append(f"{atom.name}-{atom.ID}")
        return f"{self.type}\n   {','.join(protein_str)}\n   {','.join(ligand_str)}\n"

    @property
    def ligand_atom_id(self):
        return [atom.ID for atom in self.ligand]

    @property
    def protein_atom_id(self):
        return [atom.ID for atom in self.protein]


def merge_interaction_dfs(df1: pd.DataFrame, df2: pd.DataFrame):
    index_all = df1.index.union(df2.index)
    interaction_all = df1.columns.union(df2.columns)

    for index in index_all:
        if index not in df1.index:
            interaction_dict = {interaction: 0 for interaction in interaction_all}
            df1 = pd.concat([df1, pd.DataFrame([interaction_dict], index=[index])])
        if index not in df2.index:
            interaction_dict = {interaction: 0 for interaction in interaction_all}
            df2 = pd.concat([df2, pd.DataFrame([interaction_dict], index=[index])])

    interaction_remain = interaction_all.difference(df1.columns)
    for interaction in interaction_remain:
        df1[interaction] = 0

    interaction_remain = interaction_all.difference(df2.columns)
    for interaction in interaction_remain:
        df2[interaction] = 0

    return df1.fillna(0).sort_index(
        #key=lambda col: np.array([int(item.split("_")[0]) for item in col]), 
        ascending=False
    ), df2.fillna(0).sort_index(
                                #key=lambda col: np.array([int(item.split("_")[0]) for item in col]), 
                                ascending=False)


class ReportInteractions:
    def __init__(self, interactions: AllInteraction, protein_name="protein", ligand_name="LIG"):
        self.interactions = interactions
        self.protein_name = protein_name
        self.ligand_name = ligand_name
    
    __ATTRS = [["hbonds_ldon","hbond"],["hbonds_pdon","hbond"],
               ["chpi_laro","chpi"],["chpi_paro","chpi"],
               ["halogen_bonds","halogen"],
               ["pi_stacking","pistacking"],
               ["pication_laro","pication"],["pication_paro","pication"],
               ["saltbridge_lneg","saltbridge"],["saltbridge_pneg","saltbridge"],
               ["hydrophobic_contacts","hydrophobic"],
               ["water_bridges","waterbridges"]]
        
    __INFO = {
        "protein":
            {
                "saltbridge":{"saltbridge_lneg":["acceptor","atoms",0],"saltbridge_pneg":["donor","atoms",0]},
                "hbond":{"hbonds_ldon":["acceptor","atom",-1],"hbonds_pdon":["donor","adj_atom",-1]},
                "pistacking":{"pi_stacking":["acceptor","atoms",0]},
                "pication":{"pication_laro":["donor","atoms",0],"pication_paro":["acceptor","atoms",0]},
                "weakhbond":{},
                "chpi":{"chpi_laro":["donor","adj_atom",-1],"chpi_paro":["acceptor","atoms",0]},
                "hydrophobic":{"hydrophobic_contacts":["acceptor","atom",-1]},
                "halogen":{"halogen_bonds":["acceptor","atom",-1]},
                "waterbridges":{}
            },
        "ligand":
            {
               "saltbridge":{"saltbridge_lneg":["donor","atoms",0],"saltbridge_pneg":["acceptor","atoms",0]},
                "hbond":{"hbonds_ldon":["donor","adj_atom",-1],"hbonds_pdon":["acceptor","atom",-1]},
                "pistacking":{"pi_stacking":["donor","atoms",0]},
                "pication":{"pication_laro":["acceptor","atoms",0],"pication_paro":["donor","atoms",0]},
                "weakhbond":{},
                "chpi":{"chpi_laro":["acceptor","atoms",0],"chpi_paro":["donor","adj_atom",-1]},
                "hydrophobic":{"hydrophobic_contacts":["donor","atom",-1]},
                "halogen":{"halogen_bonds":["donor","atom",-1]}, 
                "waterbridges":{}
            }
    }

    def atom_info(self):
        #ligand_atoms: Dict[str, List[Atom]] = defaultdict(list)
        total_atoms = {}
        for mm,dd in self.__INFO.items():
            total_atoms[mm] = defaultdict(list)
        
            for attr,ee in dd.items():
                total_atoms[mm][attr] = []
                for kk,vv in ee.items():
                    terms = getattr(self.interactions,kk,[])
                    for term in terms:
                        _tmp_atom = getattr(getattr(term,vv[0]),vv[1])
                        if vv[2] != -1:
                            this_atom = _tmp_atom[vv[2]]
                        else:
                            this_atom = _tmp_atom
                        total_atoms[mm][attr].append(this_atom)
        return total_atoms["ligand"], total_atoms["protein"]

    def residue_info(self):
        ligand_atoms,protein_atoms = self.atom_info()
        interaction_dict: Dict[str, Set[str]] = defaultdict(set)
        for wbridge in self.interactions.water_bridges:
            if wbridge.type == WaterBridgeType.pro_acc_lig_acc:
                protein_atoms["waterbridge"].append(wbridge.hbond_pair[0].a)
            elif wbridge.type == WaterBridgeType.pro_don_lig_don:
                protein_atoms["waterbridge"].append(wbridge.hbond_pair[0].d)
            elif wbridge.type == WaterBridgeType.pro_don_lig_acc:
                protein_atoms["waterbridge"].append(wbridge.hbond_pair[0].d)
            else:  # wbridge.type == WaterBridgeType.pro_acc_lig_don:
                protein_atoms["waterbridge"].append(wbridge.hbond_pair[0].a)
        if len(protein_atoms):
            for key, value in protein_atoms.items():
                for atom in value:
                    try:
                        interaction_dict[key].add(f"{atom.residue_ID}_{atom.residue}")
                    except:
                        interaction_dict[key].add(f"{atom.ID}_{atom.name}")
        return interaction_dict


    def interaction_general_info(self):  # pdb by residue index, ligand by atom index
        residue_info = self.residue_info()



    def residue_ligand_atom_idx_detail(self):  # for pymol
        #####代码还未改正
        
        
        # Hydrophobic Contacts
        # Contains IDs of contributing binding site, ligand atoms and the pairings
        pymol_data = PymolInteractionData(ligand_info=self.ligand_name, protein_name=self.protein_name)

        pymol_data.hydrophobic_contacts = [
            HydrophobicInfoID(
                protein_atom=hydrophobic.acceptor.atom.ID, ligand_atom=hydrophobic.donor.ID
            )
            for hydrophobic in self.interactions.hydrophobic_contacts
        ]

        pymol_data.hydrogen_bonds = [
            HBondsInfoID(hbond.d.atom_orig_idx, hbond.a.atom_orig_idx, False, hbond.type)
            for hbond in self.interactions.hbonds_ldon
        ] + [
            HBondsInfoID(hbond.d.atom_orig_idx, hbond.a.atom_orig_idx, True, hbond.type)
            for hbond in self.interactions.hbonds_pdon
        ]

        pymol_data.halogen_bonds = [
            HalogenInfoID(don_id=h.don.x.atom_orig_idx, acc_id=h.acc.o.atom_orig_idx)
            for h in self.interactions.halogen_bonds
        ]

        # Pistacking
        pymol_data.pi_stacking = [
            PiStackingInfoID(
                protein_ring_atoms=[atom.atom_orig_idx for atom in pistack.protein.atoms],
                protein_ring_center=list(pistack.protein.center),
                ligand_ring_atoms=[atom.atom_orig_idx for atom in pistack.ligand.atoms],
                ligand_ring_center=list(pistack.ligand.center),
                type=pistack.type,
            )
            for pistack in self.interactions.pi_stacking
        ]

        # Pi-cation interactions
        pymol_data.pi_cation = [
            PiCationInfoID(
                ring_center=list(picat.ring.center),
                charge_center=list(picat.charge.center),
                ring_atoms=[atom.atom_orig_idx for atom in picat.ring.atoms],
                charge_atoms=[atom.atom_orig_idx for atom in picat.charge.atoms],
                protein_charged=picat.prot_charged,
            )
            for picat in self.interactions.pication_paro + self.interactions.pication_laro
        ]

        pymol_data.ch_pi = [
            CHPiInfoID(
                ring_center=list(ch_pi.ring.center),
                don_id=ch_pi.ch.d.atom_orig_idx,
                ring_atoms=[atom.atom_orig_idx for atom in ch_pi.ring.atoms],
                ring_in_protein=ch_pi.ring_in_protein,
            )
            for ch_pi in self.interactions.chpi_laro + self.interactions.chpi_paro
        ]

        pymol_data.salt_bridges = [
            SaltBridgeInfoID(
                positive_atoms=[atom.atom_orig_idx for atom in sbridge.positive.atoms],
                negative_atoms=[atom.atom_orig_idx for atom in sbridge.negative.atoms],
                positive_center=list(sbridge.positive.center),
                negative_center=list(sbridge.negative.center),
                protein_is_positive=sbridge.prot_is_pos,
            )
            for sbridge in self.interactions.saltbridge_lneg + self.interactions.saltbridge_pneg
        ]

        for wbridge in self.interactions.water_bridges:
            if wbridge.type == WaterBridgeType.pro_acc_lig_acc:
                protein = wbridge.hbond_pair[0].a.atom_orig_idx
                ligand = wbridge.hbond_pair[1].a.atom_orig_idx
            elif wbridge.type == WaterBridgeType.pro_don_lig_don:
                protein = wbridge.hbond_pair[0].d.atom_orig_idx
                ligand = wbridge.hbond_pair[1].d.atom_orig_idx
            elif wbridge.type == WaterBridgeType.pro_don_lig_acc:
                protein = wbridge.hbond_pair[0].d.atom_orig_idx
                ligand = wbridge.hbond_pair[1].a.atom_orig_idx
            else:  # wbridge.type == WaterBridgeType.pro_acc_lig_don:
                protein = wbridge.hbond_pair[0].a.atom_orig_idx
                ligand = wbridge.hbond_pair[1].d.atom_orig_idx
            water = wbridge.water.o.atom_orig_idx
            pymol_data.water_bridges.append(
                WaterBridgeInfoID(protein_atom_id=protein, ligand_atom_id=ligand, water_id=water)
            )

        # metal binding
        pymol_data.metal_complexes = [
            MetalInfoID(
                metal_id=metal_interaction.metal.atom_orig_idx,
                target_id=metal_interaction.target.atom_orig_idx,
                location=metal_interaction.location,
            )
            for metal_interaction in self.interactions.metal_complexation
        ]

        return pymol_data

    def get_atoms(self, hydrophobic=True):
        atoms = {"ligand":[],"protein":[]}
        for kk,vv in atoms.items():
            
            for attr in self.__ATTRS:
                _attr1 = attr[1]
                _attr = attr[0]
                _tmp1 = getattr(self.interactions,_attr)
                for tmp in _tmp1:
                    _label = self.__INFO[kk][_attr1][_attr]
                    _tmp2 = getattr(tmp,_label[0])
                    _tmp3 = getattr(_tmp2,_label[1])
                    if isinstance(_tmp3,list):
                        vv.extend(_tmp3)
                    else:
                        vv.append(_tmp3)
        return atoms["ligand"], atoms["protein"]
            
    def residue_ligand_atom_idx(self):
        res_str = ""
        for attr in self.__ATTRS:
            _attr1 = attr[1]
            _attr = attr[0]
            _tmp1 = getattr(self.interactions,_attr)
            for tmp in _tmp1:
                _label1 = self.__INFO["ligand"][_attr1][_attr]
                _label2 = self.__INFO["protein"][_attr1][_attr]
                _tmp21 = getattr(tmp,_label1[0])
                _tmp31 = getattr(_tmp21,_label1[1])
                _tmp22 = getattr(tmp,_label2[0])
                _tmp32 = getattr(_tmp22,_label2[1])
                if not isinstance(_tmp31,list):
                    _tmp31 = [_tmp31]
                if not isinstance(_tmp32,list):
                    _tmp32 = [_tmp32]
                
                res_str += str(InteractionDual(type=_attr1,ligand=_tmp31,protein=_tmp32))
                
        return res_str 


    def interaction_to_dict(self):
        result = defaultdict(list)
        

    def interation_to_dict(self):  # for database
        #####代码还没有改正
        
        
        result = defaultdict(list)

        # Hydrogen bond interaction
        for hbond in self.interactions.hbonds_ldon:
            hbond_dict = {
                "donor": hbond.d.to_dict(),
                "acceptor": hbond.a.to_dict(),
                "hydrogen": hbond.h.to_dict(),
                "angle": hbond.angle,
                "distance": hbond.distance,
                "protein_is_donor": False,
            }
            result["hbond"].append(hbond_dict)

        for hbond in self.interactions.hbonds_pdon:
            hbond_dict = {
                "donor": hbond.d.to_dict(),
                "acceptor": hbond.a.to_dict(),
                "hydrogen": hbond.h.to_dict(),
                "angle": hbond.angle,
                "distance": hbond.distance,
                "protein_is_donor": True,
            }
            result["hbond"].append(hbond_dict)

        for sb in self.interactions.saltbridge_lneg + self.interactions.saltbridge_pneg:
            sb_dict = {}
            pos_charged_property = {
                "atoms": [atom.to_dict() for atom in sb.positive.atoms],
                "center": sb.positive.center.tolist(),
                "group": sb.positive.group,
            }
            neg_charged_property = {
                "atoms": [atom.to_dict() for atom in sb.negative.atoms],
                "center": sb.negative.center.tolist(),
                "group": sb.negative.group,
            }
            sb_dict["positive"] = pos_charged_property
            sb_dict["negative"] = neg_charged_property
            sb_dict["distance"] = sb.distance
            sb_dict["protein_is_positive"] = sb.prot_is_pos
            result["salt_bridge"].append(sb_dict)

        for hydro in self.interactions.hydrophobic_contacts:
            hydro_dict = {
                "protein": hydro.bsatom.to_dict(),
                "ligand": hydro.bsatom.to_dict(),
                "distance": hydro.distance,
            }
            result["hydrophobic"].append(hydro_dict)

        for pipi in self.interactions.pi_stacking:
            pipi_dict = {}
            protein_property = {
                "atoms": [atom.to_dict() for atom in pipi.protein.atoms],
                "normal": pipi.protein.normal.tolist(),
                "center": pipi.protein.center.tolist(),
            }

            ligand_property = {
                "atoms": [atom.to_dict() for atom in pipi.ligand.atoms],
                "normal": pipi.ligand.normal.tolist(),
                "center": pipi.ligand.center.tolist(),
            }
            pipi_dict["protein"] = protein_property
            pipi_dict["ligand"] = ligand_property
            pipi_dict["angle"] = pipi.angle
            pipi_dict["distance"] = pipi.distance
            pipi_dict["offset"] = pipi.offset
            pipi_dict["type"] = pipi.type
            result["pi_stacking"].append(pipi_dict)

        for pi in self.interactions.pication_laro:
            pi_dict = {}
            charged_property = {
                "atoms": [atom.to_dict() for atom in pi.charge.atoms],
                "center": pi.charge.center.tolist(),
                "group": pi.charge.group,
            }
            ring_property = {
                "atoms": [atom.to_dict() for atom in pi.ring.atoms],
                "normal": pi.ring.normal.tolist(),
                "center": pi.ring.center.tolist(),
            }
            pi_dict["charged"] = charged_property
            pi_dict["ring"] = ring_property
            pi_dict["distance"] = pi.distance
            pi_dict["offset"] = pi.offset
            pi_dict["protein_charged"] = pi.prot_charged
            result["pi_cation"].append(pi_dict)

        # Halogen
        for hal in self.interactions.halogen_bonds:
            hal_dict = {
                "acceptor": [hal.acc.o.to_dict(), hal.acc.y.to_dict()],
                "donor": [hal.don.x.to_dict(), hal.don.c.to_dict()],
                "acceptor_angle": hal.acc_angle,
                "donor_angle": hal.don_angle,
                "distance": hal.distance,
            }
            result["halogen"].append(hal_dict)

        # WaterBridge
        for wb in self.interactions.water_bridges:
            wb_dict = {
                "hbond1": {
                    "donor": wb.hbond_pair[0].d.to_dict(),
                    "acceptor": wb.hbond_pair[0].a.to_dict(),
                    "hydrogen": wb.hbond_pair[0].h.to_dict(),
                    "angle": wb.hbond_pair[0].angle,
                    "distance": wb.hbond_pair[0].distance,
                },
                "hbond2": {
                    "donor": wb.hbond_pair[1].d.to_dict(),
                    "acceptor": wb.hbond_pair[1].a.to_dict(),
                    "hydrogen": wb.hbond_pair[1].h.to_dict(),
                    "angle": wb.hbond_pair[1].angle,
                    "distance": wb.hbond_pair[1].distance,
                },
                "type": wb.type.value,
            }
            result["water_bridge"].append(wb_dict)
        return result
        # for metal in self.int.metal_complexation:
        #     metal.complex_num

    def to_molstar_data(self, total_chain=2):
        #####代码还未改正
        """
        Unknown = 0,
        Ionic = 1,
        CationPi = 2,
        PiStacking = 3,
        HydrogenBond = 4,
        HalogenBond = 5,
        Hydrophobic = 6,
        MetalCoordination = 7,
        WeakHydrogenBond = 8,
        """

        def coor_to_dict(coor):
            return {"x": coor[0], "y": coor[1], "z": coor[2]}

        interaction_tuple = namedtuple("interation_tuple", ["chainId", "type", "locationProtein", "locationLigand"])
        interContacts = []
        for hbond in self.interactions.hbonds_ldon:
            interContacts.append(
                interaction_tuple(
                    chainId=hbond.a.chain,
                    type=4,
                    locationProtein=coor_to_dict(hbond.a.coor),
                    locationLigand=coor_to_dict(hbond.d.coor),
                )
            )

        for hbond in self.interactions.hbonds_pdon:
            interContacts.append(
                interaction_tuple(
                    chainId=hbond.d.chain,
                    type=4,
                    locationProtein=coor_to_dict(hbond.d.coor),
                    locationLigand=coor_to_dict(hbond.a.coor),
                )
            )
        for saltbridge in self.interactions.saltbridge_lneg:
            interContacts.append(
                interaction_tuple(
                    chainId=saltbridge.positive.atoms[0].chain,
                    type=1,
                    locationProtein=coor_to_dict(saltbridge.positive.center),
                    locationLigand=coor_to_dict(saltbridge.negative.center),
                )
            )
        for saltbridge in self.interactions.saltbridge_pneg:
            interContacts.append(
                interaction_tuple(
                    chainId=saltbridge.negative.atoms[0].chain,
                    type=1,
                    locationProtein=coor_to_dict(saltbridge.negative.center),
                    locationLigand=coor_to_dict(saltbridge.positive.center),
                )
            )
        for pi in self.interactions.pi_stacking:
            interContacts.append(
                interaction_tuple(
                    chainId=pi.protein.atoms[0].chain,
                    type=3,
                    locationProtein=coor_to_dict(pi.protein.center),
                    locationLigand=coor_to_dict(pi.ligand.center),
                )
            )
        for pi in self.interactions.pication_laro:
            interContacts.append(
                interaction_tuple(
                    chainId=pi.charge.atoms[0].chain,
                    type=2,
                    locationProtein=coor_to_dict(pi.charge.center),
                    locationLigand=coor_to_dict(pi.ring.center),
                )
            )
        for pi in self.interactions.pication_paro:
            interContacts.append(
                interaction_tuple(
                    chainId=pi.ring.atoms[0].chain,
                    type=2,
                    locationProtein=coor_to_dict(pi.ring.center),
                    locationLigand=coor_to_dict(pi.charge.center),
                )
            )
        for contact in self.interactions.hydrophobic_contacts:
            interContacts.append(
                interaction_tuple(
                    chainId=contact.bsatom.chain,
                    type=6,
                    locationProtein=coor_to_dict(contact.bsatom.coor),
                    locationLigand=coor_to_dict(contact.ligatom.coor),
                )
            )
        for halogen in self.interactions.halogen_bonds:
            interContacts.append(
                interaction_tuple(
                    chainId=halogen.acc.o.chain,
                    type=5,
                    locationProtein=coor_to_dict(halogen.acc.o.coor),
                    locationLigand=coor_to_dict(halogen.don.x.coor),
                )
            )

        result = {
            "interContacts": [contact._asdict() for contact in interContacts],
            "proteinData": {"chainTotalCount": total_chain},
        }
        return result

    def visualize_by_pymol(self, *source_file, output_dir="."):
        #####代码还未改正
        try:
            import pymol
            from pymol import cmd
        except ImportError:
            logger.info("Pymol is not installed, the visualization function won't work")
            return
        pymol_data = self.residue_ligand_atom_idx_detail()
        SingleFrameViewer(*source_file, report=pymol_data).save_pymol_presentation(Path(output_dir))


    def old_atom_info(self):
        ligand_atoms: Dict[str, List[Atom]] = defaultdict(list)
        ligand_atoms["saltbridge"] = [
            salt_bridge.donor.atoms[0] for salt_bridge in self.interactions.saltbridge_lneg
        ] + [salt_bridge.acceptor.atoms[0] for salt_bridge in self.interactions.saltbridge_pneg]
        ligand_atoms["hbond"] = [hbond.donor.adj_atom for hbond in self.interactions.hbonds_ldon ] + [
            hbond.acceptor.atom for hbond in self.interactions.hbonds_pdon
        ]
        ligand_atoms["pistacking"] = [pi.donor.atoms[0] for pi in self.interactions.pi_stacking]
        ligand_atoms["pication"] = [pi.acceptor.atoms[0] for pi in self.interactions.pication_laro] + [
            pi.donor.atoms[0] for pi in self.interactions.pication_paro
        ]
        #ligand_atoms["weakhbond"] = [hbond.d for hbond in self.interactions.hbonds_ldon if hbond.type == "weak"] + [
        #    hbond.a for hbond in self.interactions.hbonds_pdon if hbond.type == "weak"
        #]
        ligand_atoms["weakhbond"] = []
        ligand_atoms["chpi"] = [chpi.donor.adj_atom for chpi in self.interactions.chpi_paro] + [
            chpi.acceptor.atoms[0] for chpi in self.interactions.chpi_laro
        ]
        ligand_atoms["hydrophobic"] = [contact.donor.atom for contact in self.interactions.hydrophobic_contacts]
        ligand_atoms["halogen"] = [halogen.donor.atom for halogen in self.interactions.halogen_bonds]

        protein_atoms: Dict[str, List[Atom]] = defaultdict(list)
        protein_atoms["saltbridge"] = [
            salt_bridge.acceptor.atoms[0] for salt_bridge in self.interactions.saltbridge_lneg
        ] + [salt_bridge.donor.atoms[0] for salt_bridge in self.interactions.saltbridge_pneg]
        protein_atoms["hbond"] = [hbond.acceptor.atom for hbond in self.interactions.hbonds_ldon] + [
            hbond.donor.adj_atom for hbond in self.interactions.hbonds_pdon
        ]
        protein_atoms["pistacking"] = [pi.acceptor.atoms[0] for pi in self.interactions.pi_stacking]
        protein_atoms["pication"] = [pi.donor.atoms[0] for pi in self.interactions.pication_laro] + [
            pi.acceptor.atoms[0] for pi in self.interactions.pication_paro
        ]
        #protein_atoms["weakhbond"] = [hbond.a for hbond in self.interactions.hbonds_ldon if hbond.type == "weak"] + [
        #    hbond.d for hbond in self.interactions.hbonds_pdon if hbond.type == "weak"
        #]
        protein_atoms["weakhbond"] = []
        protein_atoms["chpi"] = [chpi.acceptor.atoms[0] for chpi in self.interactions.chpi_paro] + [
            chpi.donor.adj_atom for chpi in self.interactions.chpi_laro
        ]
        protein_atoms["hydrophobic"] = [contact.acceptor.atom for contact in self.interactions.hydrophobic_contacts]
        protein_atoms["halogen"] = [halogen.acceptor.atom for halogen in self.interactions.halogen_bonds]
        return ligand_atoms, protein_atoms

    def old_residue_info(self):  # for traj analyze (bar plot)
        atoms: Dict[str, List[Atom]] = defaultdict(list)
        interaction_dict: Dict[str, Set[str]] = defaultdict(set)
        atoms["saltbridge"] = [salt_bridge.acceptor.atoms[0] for salt_bridge in self.interactions.saltbridge_lneg] + [
            salt_bridge.donor.atoms[0] for salt_bridge in self.interactions.saltbridge_pneg
        ]
        atoms["hbond"] = [hbond.acceptor.atom for hbond in self.interactions.hbonds_ldon] + [
            hbond.donor.adj_atom for hbond in self.interactions.hbonds_pdon
        ]
        atoms["pistacking"] = [pi.acceptor.atoms[0] for pi in self.interactions.pi_stacking]
        atoms["pication"] = [pi.donor.atoms[0] for pi in self.interactions.pication_laro] + [
            pi.acceptor.atoms[0] for pi in self.interactions.pication_paro
        ]
        #atoms["weakhbond"] = [hbond.a for hbond in self.interactions.hbonds_ldon if hbond.type == "weak"] + [
        #    hbond.d for hbond in self.interactions.hbonds_pdon if hbond.type == "weak"
        #]
        atoms["weakhbond"] = []
        atoms["chpi"] = [chpi.acceptor.atoms[0] for chpi in self.interactions.chpi_paro] + [
            chpi.donor.adj_atom for chpi in self.interactions.chpi_laro
        ]
        
        atoms["hydrophobic"] = [contact.acceptor.atom for contact in self.interactions.hydrophobic_contacts]
        atoms["halogen"] = [halogen.acceptor.atom for halogen in self.interactions.halogen_bonds]
        for wbridge in self.interactions.water_bridges:
            if wbridge.type == WaterBridgeType.pro_acc_lig_acc:
                atoms["waterbridge"].append(wbridge.hbond_pair[0].a)
            elif wbridge.type == WaterBridgeType.pro_don_lig_don:
                atoms["waterbridge"].append(wbridge.hbond_pair[0].d)
            elif wbridge.type == WaterBridgeType.pro_don_lig_acc:
                atoms["waterbridge"].append(wbridge.hbond_pair[0].d)
            else:  # wbridge.type == WaterBridgeType.pro_acc_lig_don:
                atoms["waterbridge"].append(wbridge.hbond_pair[0].a)
        if len(atoms):
            for key, value in atoms.items():
                for atom in value:
                    interaction_dict[key].add(f"{atom.residue_ID}_{atom.residue}")
        return interaction_dict

    def old_get_atoms(self, hydrophobic=True):    
        ligand_atoms = set()
        protein_atoms = set()
        for hb in self.interactions.hbonds_ldon:
            ligand_atoms.add(hb.d)
            protein_atoms.add(hb.a)
        for hb in self.interactions.hbonds_pdon:
            ligand_atoms.add(hb.a)
            protein_atoms.add(hb.d)
        for hb in self.interactions.halogen_bonds:
            ligand_atoms.add(hb.don.x)
            protein_atoms.add(hb.acc.o)
        for saltbridge in self.interactions.saltbridge_lneg:
            for atom in saltbridge.negative.atoms:
                ligand_atoms.add(atom)
            for atom in saltbridge.positive.atoms:
                protein_atoms.add(atom)
        for saltbridge in self.interactions.saltbridge_pneg:
            for atom in saltbridge.positive.atoms:
                ligand_atoms.add(atom)
            for atom in saltbridge.negative.atoms:
                protein_atoms.add(atom)
        if hydrophobic:
            for hydro in self.interactions.hydrophobic_contacts:
                ligand_atoms.add(hydro.ligatom)
                protein_atoms.add(hydro.bsatom)
        return ligand_atoms, protein_atoms

    def old_residue_ligand_atom_idx(self):  # for ai fep
        res_str = ""
        
        for hb in self.interactions.hbonds_ldon:
            if hb.type == "regular":
                hbond = InteractionDual(type="hbond", ligand=[hb.d], protein=[hb.a])
            else:
                hbond = InteractionDual(type="weakhbond", ligand=[hb.d], protein=[hb.a])
            res_str += str(hbond)

        for hb in self.interactions.hbonds_pdon:
            if hb.type == "regular":
                hbond = InteractionDual(type="hbond", ligand=[hb.a], protein=[hb.d])
            else:
                hbond = InteractionDual(type="weakhbond", ligand=[hb.d], protein=[hb.a])
            res_str += str(hbond)

        for chpi in self.interactions.chpi_laro:
            res_str += str(InteractionDual(type="chpi", ligand=chpi.ring.atoms, protein=[chpi.ch.d]))
        for chpi in self.interactions.chpi_paro:
            res_str += str(InteractionDual(type="chpi", ligand=[chpi.ch.d], protein=chpi.ring.atoms))

        for hb in self.interactions.halogen_bonds:
            res_str += str(InteractionDual(type="Halogen", ligand=[hb.don.x], protein=[hb.acc.o]))
        for pistack in self.interactions.pi_stacking:
            res_str += str(
                InteractionDual(type="pistacking", ligand=pistack.ligand.atoms, protein=pistack.protein.atoms)
            )
        for pication in self.interactions.pication_laro:
            res_str += str(InteractionDual(type="pication", ligand=pication.ring.atoms, protein=pication.charge.atoms))
        for pication in self.interactions.pication_paro:
            res_str += str(InteractionDual(type="pication", ligand=pication.charge.atoms, protein=pication.ring.atoms))

        for saltbrige in self.interactions.saltbridge_lneg:
            res_str += str(
                InteractionDual(
                    type="saltbridge",
                    ligand=saltbrige.negative.atoms,
                    protein=saltbrige.positive.atoms,
                )
            )
        for saltbridge in self.interactions.saltbridge_pneg:
            res_str += str(
                InteractionDual(
                    type="saltbridge",
                    ligand=saltbridge.positive.atoms,
                    protein=saltbridge.negative.atoms,
                )
            )
        for hydro in self.interactions.hydrophobic_contacts:
            res_str += str(InteractionDual(type="hydro", ligand=[hydro.ligatom], protein=[hydro.bsatom]))

        for wbridge in self.interactions.water_bridges:

            water_residue_idx = wbridge.water.o.residue_idx
            if wbridge.type == WaterBridgeType.pro_acc_lig_acc:
                protein = wbridge.hbond_pair[0].a
                ligand = wbridge.hbond_pair[1].a
            elif wbridge.type == WaterBridgeType.pro_don_lig_don:
                protein = wbridge.hbond_pair[0].d
                ligand = wbridge.hbond_pair[1].d
            elif wbridge.type == WaterBridgeType.pro_don_lig_acc:
                protein = wbridge.hbond_pair[0].d
                ligand = wbridge.hbond_pair[1].a
            else:  # wbridge.type == WaterBridgeType.pro_acc_lig_don:
                protein = wbridge.hbond_pair[0].a
                ligand = wbridge.hbond_pair[1].d
            res_str += str(InteractionDual(type="waterbridge", ligand=[ligand], protein=[protein]))
            res_str += f"\nWater\n{water_residue_idx}"
        return res_str

class ReportMultipleInteraction:
    def __init__(self, interactions=None):
        if interactions is None:
            self.interactions: List[Union[AllInteraction, AllWaterLigandInteraction]] = []
        else:
            self.interactions = interactions

    def add(self, interaction: AllInteraction):
        self.interactions.append(interaction)

    def info_df(self, show_figure=True, figure_name="interaction.html"):
        if isinstance(self.interactions[0], AllInteraction):
            return self._residue_info_df(show_figure=show_figure, figure_name=figure_name)
        else:  #
            return self._atom_info_df(show_figure=show_figure, figure_name=figure_name)

    def info_detail(self, filename):
        if isinstance(self.interactions[0], AllInteraction):
            return self._residue_info_detail(filename)
        else:
            return self._atom_info_detail(filename)

    def _residue_info_df(self, show_figure=True, figure_name="interaction.html"):
        if len(self.interactions) == 0:
            logger.error("analyze analyze interaction first!")
            return None

        residue_count = defaultdict(lambda: defaultdict(int))
        effective_interaction = set()

        for interaction in self.interactions:
            result = ReportInteractions(interaction).residue_info()
            for key, value in result.items():
                effective_interaction.add(key)
                for residue in value:
                    residue_count[residue][key] += 1

        df = pd.DataFrame(columns=["residue"] + list(effective_interaction))
        for residue in residue_count.keys():
            dict_value = residue_count[residue]
            pd_row = dict()
            for interaction in effective_interaction:
                pd_row["residue"] = residue
                pd_row[interaction] = dict_value.get(interaction, 0) / len(self.interactions)
            df = pd.concat((df, pd.DataFrame.from_records([pd_row])), ignore_index=True)

        df.sort_values(
            by="residue",
            #key=lambda col: np.array([int(item.split("_")[0]) for item in col]),
            inplace=True,
            ascending=False,
        )
        df.set_index("residue", inplace=True)

        if show_figure:
            normal_md_bar_plot(df, figure_name)
        return df

    def _atom_info_df(self, show_figure=True, figure_name="interaction.html"):
        atom_count = defaultdict(lambda: defaultdict(int))
        for interaction in self.interactions:
            for hbond in interaction.hbonds_ldon:
                atom_count[f"{hbond.d.atom_orig_idx}_{hbond.d.atom_name}"]["donor"] += 1
            for hbond in interaction.hbonds_pdon:
                atom_count[f"{hbond.a.atom_orig_idx}_{hbond.a.atom_name}"]["acceptor"] += 1
        df = pd.DataFrame.from_dict(atom_count).T
        df.fillna(0, inplace=True)
        df /= len(self.interactions)
        df.sort_index(
            key=lambda col: np.array([int(item.split("_")[0]) for item in col]), inplace=True, ascending=False
        )
        if show_figure:
            ligand_water_bar_plot(df, figure_name)
        return df

    def _residue_info_detail(self, filename="detail.txt"):
        with open(filename, "w") as f:
            for i, interaction in enumerate(self.interactions):
                result = ReportInteractions(interaction)
                f.write(f"Frame {i}\n")
                f.write(result.residue_ligand_atom_idx())
                f.write("\n")

    def _atom_info_detail(self, filename="detail.txt"):
        with open(filename, "w") as f:
            for i, interaction in enumerate(self.interactions):
                f.write(f"Frame {i}\n")
                for hbond in interaction.hbonds_ldon:
                    f.write(str(InteractionDual("hbond_ligand_donor", ligand=[hbond.d], protein=[hbond.a])))
                for hbond in interaction.hbonds_pdon:
                    f.write(str(InteractionDual("hbond_water_donor", ligand=[hbond.a], protein=[hbond.d])))
                f.write("\n")
