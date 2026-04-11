import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple

try:
    import pymol
    from pymol import cmd
except ImportError:
    pymol = None
    cmd = None
    # logger.info("pymol has not installed, the visualize cannot work")


@dataclass
class HBondsInfoID:
    donor: int
    acceptor: int
    protein_is_donor: bool
    type: str


@dataclass
class HydrophobicInfoID:
    protein_atom: int
    ligand_atom: int


@dataclass
class HalogenInfoID:
    don_id: int
    acc_id: int


@dataclass
class PiStackingInfoID:
    protein_ring_atoms: List[int]
    protein_ring_center: List[float]
    ligand_ring_atoms: List[int]
    ligand_ring_center: List[float]
    type: str


@dataclass
class PiCationInfoID:
    ring_center: List[float]
    charge_center: List[float]
    ring_atoms: List[int]
    charge_atoms: List[int]
    protein_charged: bool


@dataclass
class CHPiInfoID:
    ring_center: List[float]
    don_id: int
    ring_atoms: List[int]
    ring_in_protein: bool


@dataclass
class SaltBridgeInfoID:
    positive_atoms: List[int]
    negative_atoms: List[int]
    positive_center: List[float]
    negative_center: List[float]
    protein_is_positive: bool


@dataclass
class WaterBridgeInfoID:
    protein_atom_id: int
    ligand_atom_id: int
    water_id: int


@dataclass
class MetalInfoID:
    metal_id: int
    target_id: int
    location: str


@dataclass  # python 3.8 does not support slots
class PymolInteractionData:
    ligand_info: str
    protein_name: str
    hydrogen_bonds: List[HBondsInfoID] = field(default_factory=list)
    hydrophobic_contacts: List[HydrophobicInfoID] = field(default_factory=list)
    halogen_bonds: List[HalogenInfoID] = field(default_factory=list)
    pi_stacking: List[PiStackingInfoID] = field(default_factory=list)
    pi_cation: List[PiCationInfoID] = field(default_factory=list)
    ch_pi: List[CHPiInfoID] = field(default_factory=list)
    salt_bridges: List[SaltBridgeInfoID] = field(default_factory=list)
    water_bridges: List[WaterBridgeInfoID] = field(default_factory=list)
    metal_complexes: List[MetalInfoID] = field(default_factory=list)


class SingleFrameViewer:
    def __init__(self, *source_file, report: PymolInteractionData):
        self.source_file = source_file
        self.protein_name = report.protein_name
        self.ligand_name = report.ligand_info
        self.report = report

    def save_pymol_presentation(self, output_directory=Path("")):
        vis = PyMOLVisualizer(self.report)
        pymol.finish_launching(args=["pymol", "-pcq", "-K"])
        pymol.cmd.reinitialize()
        vis.set_initial_representations()
        for i, file in enumerate(self.source_file):
            cmd.load(file)
            if i == 1:
                cmd.select("LIG", Path(file).name)
        if len(self.source_file) != 2:
            cmd.select(self.ligand_name, f"resname {self.ligand_name}")
        current_name = cmd.get_object_list(selection="(all)")[0]
        cmd.set_name(current_name, self.protein_name)
        cmd.hide("everything", "all")
        # if self.report.ligand_info[1]:
        #     cmd.select(self.ligand_name, f"resname {self.ligand_name} and chain {self.report.ligand_info[1]}")
        # else:
        cmd.show("sticks", self.ligand_name)
        cmd.color("myblue")
        cmd.color("myorange", self.ligand_name)
        cmd.util.cnc("all")

        cmd.deselect()
        vis.make_initial_selections()
        vis.show_hydrophobic()  # Hydrophobic Contacts
        vis.show_hbonds()
        vis.show_halogen()
        vis.show_stacking()
        vis.show_cationpi()
        vis.show_chpi()
        vis.show_sbridges()
        vis.show_wbridges()
        vis.show_metal()

        vis.refinements()
        vis.zoom_to_ligand()
        vis.selections_cleanup()
        vis.selections_group()
        vis.additional_cleanup()
        vis.save_session(output_directory)


class PyMOLVisualizer:
    def __init__(self, interaction_data: PymolInteractionData):
        if interaction_data is not None:
            self.protein_name = interaction_data.protein_name
            self.ligand_name = interaction_data.ligand_info
            self.interaction_data = interaction_data

    def set_initial_representations(self):
        """General settings for PyMOL"""
        self.standard_settings()
        cmd.set("dash_gap", 0)  # Show not dashes
        cmd.set("ray_shadow", 0)  # Turn on ray shadows for clearer ray-traced images
        cmd.set("cartoon_color", "mylightblue")

        # Set clipping planes for full view
        cmd.clip("far", -1000)
        cmd.clip("near", 1000)

    @staticmethod
    def make_initial_selections():
        """Make empty selections for structures and interactions"""
        for group in [
            "Hydrophobic-P",
            "Hydrophobic-L",
            "HBondDonor-P",
            "HBondDonor-L",
            "HBondAcceptor-P",
            "HBondAcceptor-L",
            "WaterBridge-P",
            "WaterBridge-L",
            "HalogenAccept",
            "HalogenDonor",
            "Water",
            "MetalIons",
            "StackRings-P",
            "PosCharge-P",
            "PosCharge-L",
            "NegCharge-P",
            "NegCharge-L",
            "PiCatRing-P",
            "StackRings-L",
            "PiCatRing-L",
            "CHPiRing-L",
            "CHPiRing-P",
            "Metal-M",
            "Metal-P",
            "Metal-W",
            "Metal-L",
            "Unpaired-HBA",
            "Unpaired-HBD",
            "Unpaired-HAL",
            "Unpaired-RINGS",
        ]:
            cmd.select(group, "None")

    def standard_settings(self):
        """Sets up standard settings for a nice visualization."""
        cmd.set("bg_rgb", [1.0, 1.0, 1.0])  # White background
        cmd.set("depth_cue", 0)  # Turn off depth cueing (no fog)
        cmd.set("cartoon_side_chain_helper", 1)  # Improve combined visualization of sticks and cartoon
        cmd.set("cartoon_fancy_helices", 1)  # Nicer visualization of helices (using tapered ends)
        cmd.set("transparency_mode", 1)  # Turn on multilayer transparency
        cmd.set("dash_radius", 0.05)
        self.set_custom_colorset()

    @staticmethod
    def set_custom_colorset():
        """Defines a colorset with matching colors. Provided by Joachim."""
        cmd.set_color("myorange", "[253, 174, 97]")
        cmd.set_color("mygreen", "[171, 221, 164]")
        cmd.set_color("myred", "[215, 25, 28]")
        cmd.set_color("myblue", "[43, 131, 186]")
        cmd.set_color("mylightblue", "[158, 202, 225]")
        cmd.set_color("mylightgreen", "[229, 245, 224]")

    @staticmethod
    def select_by_ids(sel_name, id_list, selection_exists=False, chunksize=20, restrict=None):
        """Selection with a large number of ids concatenated into a selection
        list can cause buffer overflow in PyMOL. This function takes a selection
        name and and list of IDs (list of integers) as input and makes a careful
        step-by-step selection (packages of 20 by default)"""
        id_list = list(set(id_list))  # Remove duplicates
        if not selection_exists:
            cmd.select(sel_name, "None")  # Empty selection first
        id_chunks = [id_list[i : i + chunksize] for i in range(0, len(id_list), chunksize)]
        for id_chunk in id_chunks:
            cmd.select(sel_name, "%s or (id %s)" % (sel_name, "+".join(map(str, id_chunk))))
        if restrict is not None:
            cmd.select(sel_name, "%s and %s" % (sel_name, restrict))

    @staticmethod
    def object_exists(object_name):
        """Checks if an object exists in the open PyMOL session."""
        return object_name in cmd.get_names("objects")

    def show_hydrophobic(self):
        """Visualizes hydrophobic contacts."""
        hydrophobics = self.interaction_data.hydrophobic_contacts
        protein_atoms = [h.protein_atom for h in hydrophobics]
        ligand_atoms = [h.ligand_atom for h in hydrophobics]

        if hydrophobics:
            self.select_by_ids("Hydrophobic-P", protein_atoms, restrict=self.protein_name)
            self.select_by_ids("Hydrophobic-L", ligand_atoms, restrict=self.ligand_name)

            for item in hydrophobics:
                cmd.select("tmp_bs", "id %i & %s" % (item.protein_atom, self.protein_name))
                cmd.select("tmp_lig", "id %i & %s" % (item.ligand_atom, self.ligand_name))
                cmd.distance("Hydrophobic", "tmp_bs", "tmp_lig")
            if self.object_exists("Hydrophobic"):
                cmd.set("dash_gap", 0.5, "Hydrophobic")
                cmd.set("dash_color", "grey50", "Hydrophobic")
        else:
            cmd.select("Hydrophobic-P", "None")

    def show_hbonds(self):
        """Visualizes hydrogen bonds."""
        hydrogen_bonds = self.interaction_data.hydrogen_bonds
        (
            protein_atom_donor,
            ligand_atom_donor,
            protein_donor_type,
            protein_atom_acceptor,
            ligand_atom_acceptor,
            ligand_donor_type,
        ) = ([], [], [], [], [], [])
        for bond in hydrogen_bonds:
            if bond.protein_is_donor:
                protein_atom_donor.append(bond.donor)
                ligand_atom_acceptor.append(bond.acceptor)
                protein_donor_type.append(bond.type)
            else:
                ligand_atom_donor.append(bond.donor)
                protein_atom_acceptor.append(bond.acceptor)
                ligand_donor_type.append(bond.type)

        if hydrogen_bonds:
            self.select_by_ids("HBondDonor-P", protein_atom_donor, restrict=self.protein_name)
            self.select_by_ids("HBondDonor-L", ligand_atom_donor, restrict=self.ligand_name)
            self.select_by_ids("HBondAcceptor-P", protein_atom_acceptor, restrict=self.protein_name)
            self.select_by_ids("HBondAcceptor-L", ligand_atom_acceptor, restrict=self.ligand_name)
        for protein_atom, ligand_atom, bond_type in itertools.chain(
            zip(protein_atom_donor, ligand_atom_acceptor, protein_donor_type),
            zip(protein_atom_acceptor, ligand_atom_donor, ligand_donor_type),
        ):
            cmd.select("tmp_bs", "id %i & %s" % (protein_atom, self.protein_name))
            cmd.select("tmp_lig", "id %i & %s" % (ligand_atom, self.ligand_name))
            if bond_type == "weak":
                cmd.distance("WeakHBonds", "tmp_bs", "tmp_lig")
            else:
                cmd.distance("HBonds", "tmp_bs", "tmp_lig")
        if self.object_exists("HBonds"):
            cmd.set("dash_color", "blue", "HBonds")
        elif self.object_exists("WeakHBonds"):
            cmd.set("dash_color", "lightblue", "WeakHBonds")
            cmd.set("dash_gap", 0.5, "WeakHBonds")
            cmd.set("dash_length", 0.6, "WeakHBonds")

    def show_halogen(self):
        """Visualize halogen bonds."""
        halogen = self.interaction_data.halogen_bonds
        all_don_x, all_acc_o = [], []
        for h in halogen:
            all_don_x.append(h.don_id)
            all_acc_o.append(h.acc_id)
            cmd.select("tmp_bs", "id %i & %s" % (h.acc_id, self.protein_name))
            cmd.select("tmp_lig", "id %i & %s" % (h.don_id, self.ligand_name))

            cmd.distance("HalogenBonds", "tmp_bs", "tmp_lig")
        if not len(all_acc_o) == 0:
            self.select_by_ids("HalogenAccept", all_acc_o, restrict=self.protein_name)
            self.select_by_ids("HalogenDonor", all_don_x, restrict=self.ligand_name)
        if self.object_exists("HalogenBonds"):
            cmd.set("dash_color", "greencyan", "HalogenBonds")

    def show_stacking(self):
        """Visualize pi-stacking interactions."""
        stacks = self.interaction_data.pi_stacking
        for i, stack in enumerate(stacks):
            pires_ids = "+".join(map(str, stack.protein_ring_atoms))
            piligand_atom_ids = "+".join(map(str, stack.ligand_ring_atoms))
            cmd.select("StackRings-P", "StackRings-P or (id %s & %s)" % (pires_ids, self.protein_name))
            cmd.select("StackRings-L", "StackRings-L or (id %s & %s)" % (piligand_atom_ids, self.ligand_name))
            cmd.select("StackRings-P", "byres StackRings-P")
            cmd.show("sticks", "StackRings-P")

            cmd.pseudoatom("ps-pistack-1-%i" % i, pos=stack.protein_ring_center)
            cmd.pseudoatom("ps-pistack-2-%i" % i, pos=stack.ligand_ring_center)
            cmd.pseudoatom("Centroids-P", pos=stack.protein_ring_center)
            cmd.pseudoatom("Centroids-L", pos=stack.ligand_ring_center)

            if stack.type == "P":
                cmd.distance("PiStackingP", "ps-pistack-1-%i" % i, "ps-pistack-2-%i" % i)
            if stack.type == "T":
                cmd.distance("PiStackingT", "ps-pistack-1-%i" % i, "ps-pistack-2-%i" % i)
        if self.object_exists("PiStackingP"):
            cmd.set("dash_color", "green", "PiStackingP")
            cmd.set("dash_gap", 0.3, "PiStackingP")
            cmd.set("dash_length", 0.6, "PiStackingP")
        if self.object_exists("PiStackingT"):
            cmd.set("dash_color", "smudge", "PiStackingT")
            cmd.set("dash_gap", 0.3, "PiStackingT")
            cmd.set("dash_length", 0.6, "PiStackingT")

    def show_cationpi(self):
        """Visualize cation-pi interactions."""
        for i, p in enumerate(self.interaction_data.pi_cation):
            cmd.pseudoatom("ps-picat-1-%i" % i, pos=p.ring_center)
            cmd.pseudoatom("ps-picat-2-%i" % i, pos=p.charge_center)
            if p.protein_charged:
                cmd.pseudoatom("Chargecenter-P", pos=p.charge_center)
                cmd.pseudoatom("Centroids-L", pos=p.ring_center)
                piligand_atom_ids = "+".join(map(str, p.ring_atoms))
                cmd.select("PiCatRing-L", "PiCatRing-L or (id %s & %s)" % (piligand_atom_ids, self.ligand_name))
                for a in p.charge_atoms:
                    cmd.select("PosCharge-P", "PosCharge-P or (id %i & %s)" % (a, self.protein_name))
            else:
                cmd.pseudoatom("Chargecenter-L", pos=p.charge_center)
                cmd.pseudoatom("Centroids-P", pos=p.ring_center)
                pires_ids = "+".join(map(str, p.ring_atoms))
                cmd.select("PiCatRing-P", "PiCatRing-P or (id %s & %s)" % (pires_ids, self.protein_name))
                for a in p.charge_atoms:
                    cmd.select("PosCharge-L", "PosCharge-L or (id %i & %s)" % (a, self.ligand_name))
            cmd.distance("PiCation", "ps-picat-1-%i" % i, "ps-picat-2-%i" % i)
        if self.object_exists("PiCation"):
            cmd.set("dash_color", "orange", "PiCation")
            cmd.set("dash_gap", 0.3, "PiCation")
            cmd.set("dash_length", 0.6, "PiCation")

    def show_chpi(self):
        for i, p in enumerate(self.interaction_data.ch_pi):
            cmd.pseudoatom("ps-chpi-1-%i" % i, pos=p.ring_center)
            if not p.ring_in_protein:
                cmd.pseudoatom("Centroids-L", pos=p.ring_center)
                piligand_atom_ids = "+".join(map(str, p.ring_atoms))
                cmd.select("CHPiRing-L", "CHPiRing-L or (id %s & %s)" % (piligand_atom_ids, self.ligand_name))
                cmd.select("tmp_donor", "id %i & %s" % (p.don_id, self.ligand_name))
            else:
                cmd.pseudoatom("Centroids-P", pos=p.ring_center)
                pires_ids = "+".join(map(str, p.ring_atoms))
                cmd.select("CHPiRing-P", "CHPiRing-P or (id %s & %s)" % (pires_ids, self.protein_name))
                cmd.select("tmp_donor", "id %i & %s" % (p.don_id, self.ligand_name))
            cmd.distance("CHPi", "ps-chpi-1-%i" % i, "tmp_donor")
        if self.object_exists("CHPi"):
            cmd.set("dash_color", "lightorange", "CHPi")
            cmd.set("dash_gap", 0.3, "CHPi")
            cmd.set("dash_length", 0.6, "CHPi")

    def show_sbridges(self):
        """Visualize salt bridges."""
        for i, saltb in enumerate(self.interaction_data.salt_bridges):
            if saltb.protein_is_positive:
                for patom in saltb.positive_atoms:
                    cmd.select("PosCharge-P", "PosCharge-P or (id %i & %s)" % (patom, self.protein_name))
                for latom in saltb.negative_atoms:
                    cmd.select("NegCharge-L", "NegCharge-L or (id %i & %s)" % (latom, self.ligand_name))
                for sbgroup in [
                    ["ps-sbl-1-%i" % i, "Chargecenter-P", saltb.positive_center],
                    ["ps-sbl-2-%i" % i, "Chargecenter-L", saltb.negative_center],
                ]:
                    cmd.pseudoatom(sbgroup[0], pos=sbgroup[2])
                    cmd.pseudoatom(sbgroup[1], pos=sbgroup[2])
                cmd.distance("Saltbridges", "ps-sbl-1-%i" % i, "ps-sbl-2-%i" % i)
            else:
                for patom in saltb.negative_atoms:
                    cmd.select("NegCharge-P", "NegCharge-P or (id %i & %s)" % (patom, self.protein_name))
                for latom in saltb.positive_atoms:
                    cmd.select("PosCharge-L", "PosCharge-L or (id %i & %s)" % (latom, self.ligand_name))
                for sbgroup in [
                    ["ps-sbp-1-%i" % i, "Chargecenter-P", saltb.negative_center],
                    ["ps-sbp-2-%i" % i, "Chargecenter-L", saltb.positive_center],
                ]:
                    cmd.pseudoatom(sbgroup[0], pos=sbgroup[2])
                    cmd.pseudoatom(sbgroup[1], pos=sbgroup[2])
                cmd.distance("Saltbridges", "ps-sbp-1-%i" % i, "ps-sbp-2-%i" % i)

        if self.object_exists("Saltbridges"):
            cmd.set("dash_color", "yellow", "Saltbridges")
            cmd.set("dash_gap", 0.5, "Saltbridges")

    def show_wbridges(self):
        """Visualize water bridges."""
        for bridge in self.interaction_data.water_bridges:
            cmd.select("WaterBridge-L", "WaterBridge-L or (id %i & %s)" % (bridge.ligand_atom_id, self.ligand_name))
            cmd.select("WaterBridge-P", "WaterBridge-P or (id %i & %s)" % (bridge.protein_atom_id, self.protein_name))
            cmd.select("tmp_don", "id %i & %s" % (bridge.ligand_atom_id, self.ligand_name))
            cmd.select("tmp_acc", "id %i & %s" % (bridge.protein_atom_id, self.protein_name))
            cmd.select("Water", "Water or (id %i)" % bridge.water_id)
            cmd.select("tmp_water", "id %i" % bridge.water_id)
            cmd.distance("WaterBridges", "tmp_acc", "tmp_water")
            cmd.distance("WaterBridges", "tmp_don", "tmp_water")
        if self.object_exists("WaterBridges"):
            cmd.set("dash_color", "lightblue", "WaterBridges")
        cmd.delete("tmp_water or tmp_acc or tmp_don")
        cmd.color("lightblue", "Water")
        cmd.show("spheres", "Water")

    def show_metal(self):
        """Visualize metal coordination."""
        metal_complexes = self.interaction_data.metal_complexes
        metal_ids = [metal_complex.metal_id for metal_complex in metal_complexes]
        if not len(metal_complexes) == 0:
            self.select_by_ids("Metal-M", metal_ids)
            for metal_complex in metal_complexes:
                cmd.select("tmp_m", "id %i" % metal_complex.metal_id)
                cmd.select("tmp_t", "id %i" % metal_complex.target_id)
                if metal_complex.location == "water":
                    cmd.select("Metal-W", "Metal-W or id %s" % metal_complex.target_id)
                if metal_complex.location.startswith("protein"):
                    cmd.select("tmp_t", "tmp_t & %s" % self.protein_name)
                    cmd.select(
                        "Metal-P",
                        "Metal-P or (id %s & %s)" % (metal_complex.target_id, self.protein_name),
                    )
                if metal_complex.location == "ligand":
                    cmd.select("tmp_t", "tmp_t & %s" % self.ligand_name)
                    cmd.select(
                        "Metal-L",
                        "Metal-L or (id %s & %s)" % (metal_complex.target_id, self.ligand_name),
                    )
                cmd.distance("MetalComplexes", "tmp_m", "tmp_t")
                cmd.delete("tmp_m or tmp_t")
        if self.object_exists("MetalComplexes"):
            cmd.set("dash_color", "violetpurple", "MetalComplexes")
            cmd.set("dash_gap", 0.5, "MetalComplexes")
            # Show water molecules for metal complexes
            cmd.show("spheres", "Metal-W")
            cmd.color("lightblue", "Metal-W")

    def selections_group(self):
        """Group all selections"""
        cmd.group("Structures", f"{self.protein_name} {self.ligand_name} {self.protein_name}Cartoon")
        cmd.group(
            "Interactions",
            "Hydrophobic HBonds WeakHBonds HalogenBonds WaterBridges PiCation PiStackingP PiStackingT "
            "WaterBridges CHPi "
            "Saltbridges MetalComplexes",
        )
        cmd.group("Atoms", "")
        cmd.group(
            "Atoms.Protein",
            "Hydrophobic-P HBondAcceptor-P HBondDonor-P HalogenAccept Centroids-P PiCatRing-P "
            "WaterBridge-P CHPiRing-P "
            "StackRings-P PosCharge-P NegCharge-P AllBSRes Chargecenter-P  Metal-P",
        )
        cmd.group(
            "Atoms.Ligand",
            "Hydrophobic-L HBondAcceptor-L HBondDonor-L HalogenDonor Centroids-L NegCharge-L "
            "WaterBridge-L "
            "PosCharge-L NegCharge-L ChargeCenter-L StackRings-L PiCatRing-L Metal-L Metal-M "
            "Unpaired-HBA Unpaired-HBD Unpaired-HAL Unpaired-RINGS",
        )
        cmd.group("Atoms.Other", "Water Metal-W")
        cmd.order("*", "y")

    def additional_cleanup(self):
        """Cleanup of various representations"""

        cmd.remove('not alt ""+A')  # Remove alternate conformations
        cmd.hide("labels", "Interactions")  # Hide labels of lines
        cmd.disable("%sCartoon" % self.protein_name)
        cmd.hide("everything", "hydrogens")

    def zoom_to_ligand(self):
        """Zoom in too ligand and its interactions."""
        cmd.center(self.ligand_name)
        cmd.orient(self.ligand_name)
        cmd.turn("x", 110)  # If the ligand is aligned with the longest axis, aromatic rings are hidden
        if "AllBSRes" in cmd.get_names("selections"):
            cmd.zoom("%s or AllBSRes" % self.ligand_name, 3)
        else:
            if self.object_exists(self.ligand_name):
                cmd.zoom(self.ligand_name, 3)
        cmd.origin(self.ligand_name)

    def save_session(self, output_directory):
        filename = "%s_%s" % (self.protein_name.upper(), self.ligand_name)
        cmd.save(output_directory / f"{filename}.pse")

    def adapt_for_intra(self):
        """Adapt visualization for intra-protein interactions"""

    def refinements(self):
        """Refinements for the visualization"""

        # Show sticks for all residues interacing with the ligand
        cmd.select(
            "AllBSRes",
            "byres (Hydrophobic-P or HBondDonor-P or HBondAcceptor-P or PosCharge-P or NegCharge-P or "
            "WaterBridge-p or StackRings-P or PiCatRing-P or HalogenAcc or Metal-P or Water or CHPiRing-P)",
        )
        cmd.show("sticks", "AllBSRes")
        # Show spheres for the ring centroids
        cmd.hide("everything", "centroids*")
        cmd.show("nb_spheres", "centroids*")
        # Show spheres for centers of charge
        if self.object_exists("Chargecenter-P") or self.object_exists("Chargecenter-L"):
            cmd.hide("nonbonded", "chargecenter*")
            cmd.show("spheres", "chargecenter*")
            cmd.set("sphere_scale", 0.4, "chargecenter*")
            cmd.color("yellow", "chargecenter*")

        cmd.set("valence", 1)  # Show bond valency (e.g. double bonds)
        # Optional cartoon representation of the protein

        cmd.copy(f"{self.protein_name}Cartoon", self.protein_name)
        cmd.show("cartoon", f"{self.protein_name}Cartoon")
        cmd.show("sticks", f"{self.protein_name}Cartoon")
        cmd.set("stick_transparency", 1, f"{self.protein_name}Cartoon")

        # Resize water molecules. Sometimes they are not heteroatoms HOH, but part of the protein
        cmd.set("sphere_scale", 0.2, "resn HOH or Water")  # Needs to be done here because of the copy made
        cmd.set("sphere_transparency", 0.4, "!(resn HOH or Water)")

        if "Centroids*" in cmd.get_names("selections"):
            cmd.color("grey80", "Centroids*")
        cmd.hide("spheres", "%sCartoon" % self.protein_name)
        cmd.hide("cartoon", "%sCartoon and resn DA+DG+DC+DU+DT+A+G+C+U+T" % self.protein_name)  # Hide DNA/RNA Cartoon
        if self.ligand_name == "SF4":  # Special case for iron-sulfur clusters, can't be visualized with sticks
            cmd.show("spheres", "%s" % self.ligand_name)

        cmd.hide("everything", "resn HOH &!Water")  # Hide all non-interacting water molecules
        cmd.hide(
            "sticks", "%s and !%s and !AllBSRes" % (self.protein_name, self.ligand_name)
        )  # Hide all non-interacting residues

    def selections_cleanup(self):
        """Cleans up non-used selections"""

        selections = cmd.get_names("selections")
        for selection in selections:
            try:
                empty = len(cmd.get_model(selection).atom) == 0
            except:
                empty = True
            if empty:
                cmd.delete(selection)
        cmd.deselect()
        cmd.delete("tmp*")
        cmd.delete("ps-*")

    def set_fancy_ray(self):
        """Give the molecule a flat, modern look."""
        cmd.set("light_count", 6)
        cmd.set("spec_count", 1.5)
        cmd.set("shininess", 4)
        cmd.set("specular", 0.3)
        cmd.set("reflect", 1.6)
        cmd.set("ambient", 0)
        cmd.set("direct", 0)
        cmd.set("ray_shadow", 0)  # Gives the molecules a flat, modern look
        cmd.set("ambient_occlusion_mode", 1)
        cmd.set("ray_opaque_background", 0)  # Transparent background
