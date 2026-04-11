import numpy as np
from dataclasses import dataclass, fields
from typing import Dict, List, Optional
from textwrap import indent

from ..structure.model import Site
from .config import BS_DIST
from ...utils.geometry import find_center, calc_radius
from ..structure.model import Model
from enum import Enum


@dataclass
class Interaction:
    type: str = None
    subtype: str = None
    donor: Site = None
    acceptor: Site = None
    distance: float = 0.0
    angle: float = None
    offset: float = None
    auxi_angle: float = None
    ring_in_protein: bool = False

    def __lt__(self, other):
        return self.distance < other.distance


@dataclass(repr=False)
class AllInteraction:
    saltbridge_lneg: List[Interaction]
    saltbridge_pneg: List[Interaction]
    hbonds_ldon: List[Interaction]
    hbonds_pdon: List[Interaction]
    pi_stacking: List[Interaction]
    halogen_bonds: List[Interaction]
    pication_paro: List[Interaction]
    pication_laro: List[Interaction]
    chpi_paro: List[Interaction]
    chpi_laro: List[Interaction]
    hydrophobic_contacts: List[Interaction]
    water_bridges: List[Interaction]
    metal_complexation: List[Interaction]

    def __repr__(self):
        str = ""
        for field in fields(self):
            values = getattr(self, field.name)
            if values:
                str += field.name + "\n"
                for value in values:
                    str += indent(repr(value), "    ") + "\n"
        return str
    
    @property
    def __interaction__(self):
        interaction = []
        for field in fields(self):
            values = getattr(self, field.name)
            if values:
                interaction.extend(values)
        return interaction


@dataclass(repr=False)
class AllWaterLigandInteraction:
    hbonds_ldon: List[Interaction]
    hbonds_wdon: List[Interaction]

    def __repr__(self):
        str = ""
        for field in fields(self):
            values = getattr(self, field.name)
            if values:
                str += field.name + "\n"
                for value in values:
                    str += indent(repr(value), "    ") + "\n"
        return str


class WaterBridgeType(Enum):
    pro_don_lig_don = "protein_donor_ligand_donor"  # protein donor and ligand donor
    pro_don_lig_acc = "protein_donor_ligand_acceptor"
    pro_acc_lig_don = "protein_acceptor_ligand_donor"
    pro_acc_lig_acc = "protein_acceptor_ligand_acceptor"
    


class InteractionSite:
    def __init__(self,molecule,probe,coords=None):
        #self.molecule = molecule
        #self.probe = probe
        #self.coords = coords
        pass

    @property
    def centroids(self):
        return self.find_center(self.molecule)

    @property
    def radius(self):
        return calc_radius(self.molecule)

    def find_center(self,molecule):
        return find_center([atom.coor for atom in molecule.Atoms])

    def calc_radius(self,molecule):
        return calc_radius([atom.coor for atom in molecule.Atoms])

    def find_binding_residue(self,molecule,probe):
        binding_site_residues = []
        probe_centroids = np.array(self.find_center(probe))
        probe_radius = self.calc_radius(probe)
        for residue in molecule.Groups:
            centroids = np.array(find_center([molecule.Atoms[an].coor for an in residue.atoms]))
            if np.linalg.norm(centroids - probe_centroids) < BS_DIST + probe_radius:
                binding_site_residues.append(residue)
        atoms = [an for group in binding_site_residues for an in group.atoms]
        
        return binding_site_residues,atoms

    def get_interaction_site(self,molecule,atoms=None):
        mole_model = Model(molecule,atoms=atoms)
        return mole_model.run()
    


