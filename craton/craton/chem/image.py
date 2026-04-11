#!/usr/bin/env python
import os
from typing import List

from .format._rdkit import RdkitMol
from ..utils import logger
from .. import CRATON_CONFIGURE
tmp_path = CRATON_CONFIGURE["path"]["tmp"]

from rdkit import Chem
from rdkit.Chem import AllChem, Draw


def visualize_atom_mapping(molecule1, molecule2, matches):
    rdkmh1 = Chem.MolFromMolBlock(molecule1.mol_script, removeHs=False)
    rdkmh2 = Chem.MolFromMolBlock(molecule2.mol_script, removeHs=False)
    

    AllChem.Compute2DCoords(rdkmh1)
    AllChem.Compute2DCoords(rdkmh2)

    rdkmol_list = [rdkmh1, rdkmh2]
    legends_list = [molecule1.name, molecule2.name]
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
    with open(f"{molecule1.name}_to_{molecule2.name}.svg", "w") as f:
        f.write(image)

class MoleculeImage:
    
    def __init__(self, molecule, save_file=True, target=None,fpath=None,extra=None,TD_flag=False,remove_H_flag=True):
        self.rdk = None
        self.TD_flag = TD_flag
        self.remove_H_flag = remove_H_flag
        if molecule is not None:
            self.molecule = molecule
            self._generate_rdk(molecule,TD_flag=TD_flag)
        if self.rdk is None:
            logger.warning("Molecule is None: a Molecule object is need")
            return
        if fpath is None:
            self.path = tmp_path
        else:
            self.path = fpath
        if target is None:
            self.target = ["molecule"]
        elif target == "all":
            self.target = list(self.__FUNC__.keys())
        elif target == "fragment":
            self.target = list(self.__frag_label.keys())
        else:
            if isinstance(target,list):
                self.target = target
            else:
                self.target = [target]

        self.save_file = save_file
        self.extra = extra

        self.imgs = []
        for attr in self.target:
            self.imgs.extend(self._run(attr))
        
    __frag_label = {
            "EF": "elem_frag",
            "SF": "seco_frag",
            "CSF": "csf_frag",
            "RSF": "rsf_frag",
            "TF": "tf_frag",
            "RF": "ring_frag",
            "SSSRF": "sssr_frag",
            "SaF": "scaffold_frag",
            "SkF": "sketch_frag",
        }

    def _generate_rdk(self,molecule,TD_flag=False):
        self.rdk = RdkitMol()
        self.rdk._convert(molecule.mol_script,normalization=False,draw_image=True,TD_flag=TD_flag)
        self.molecule_name = molecule.mole_name

    def _run(self,attr):
        if attr not in self.__FUNC__:
            if len(getattr(self.molecule, "Atoms", [])) > 0 and hasattr(self.molecule.Atoms[0],attr):
                extra = attr
                attr = "atom"
            else:
                extra = attr
                attr = "molecule"
        else:
            extra = self.extra

        fnames, highlighttype, highlightarrs, title = self.__FUNC__[attr](self,attr,extra=extra)
        imgs = []
        if not isinstance(fnames,list):
            fnames = [fnames]
            highlightarrs = [highlightarrs]
        if not self.save_file:
            fnames = [None for ii in range(len(fnames))]

        for ii,fname in enumerate(fnames):
            imgs.append(self.rdk.draw_figure(fname, highlighttype=highlighttype, highlightarrs=highlightarrs[ii],mole_name=title,removeH_flag=self.remove_H_flag))

        return imgs

    def _gen_atom_property_image(self,attr,extra=None):
        if attr == "atom":
            attr = extra
        ats = []
        for atom in self.molecule.Atoms:
            if hasattr(atom,attr) and atom.element != "H":
                if isinstance(getattr(atom,attr),float):
                    ats.append(str(round(getattr(atom,attr,None),3)))
                else:
                    ats.append(str(getattr(atom,attr,None)))
        mole_name = self.molecule_name.split("/")[-1]
        fname = os.path.join(self.path, f"{mole_name}_{attr}.svg")
        return fname,"label",ats, self.molecule_name

    def _gen_molecule_image(self,attr,extra=None):
        mole_name = self.molecule_name.split("/")[-1]
        fname = os.path.join(self.path, f"{mole_name}.png")
        title = self.molecule_name
        if extra is not None:
            value = getattr(self.molecule, extra, None)
            if value is None:
                logger.warning(f"molecule property not found for image title: {extra}")
            else:
                title += f"\n {str(value)}"
        return fname, None, None , title

    def _gen_molecule_ring_image(self,attr,extra=None):
    
        def craton_ids_to_rdkit_nohighlight(molecule, craton_atom_ids):
            id_to_rdk = {}
            rdk_i = 0
            for i, atom in enumerate(molecule.Atoms):
                if getattr(atom, "elem", None) == "H":
                    continue
                id_to_rdk[i] = rdk_i
                rdk_i += 1

            out = []
            for aid in craton_atom_ids:
                if aid not in id_to_rdk:
                    continue  # 氢或异常
                idx = id_to_rdk[aid]
                out.append(idx)
            return out
        fnames = []
        atoms = []
        n = 0
        for aa, bb in self.molecule.ring_dict.items():
            n += 1
            mole_name = self.molecule_name.split("/")[-1]
            this_fname = os.path.join(self.path, f"{mole_name}_ring_{n}.png")
            fnames.append(this_fname)
            # atoms.append(bb[:-1])
            atoms.append(craton_ids_to_rdkit_nohighlight(self.molecule, bb[:-1]))
 
        return fnames, "atoms", atoms, self.molecule_name

    def _gen_molecule_break_bond_image(self,attr,extra=None):
        arrs = []
        for a in self.molecule.Atoms:
            if not hasattr(a, "break_bond"):
                continue
            for i in range(len(a.connect)):
                if i < len(a.break_bond) and a.break_bond[i]:
                    if [a.No, a.connect[i]] not in arrs and [a.connect[i], a.No] not in arrs:
                        arrs.append([a.No, a.connect[i]])
        bonds = []
        for aa in arrs:
            bond = self.rdk.rdkm.GetBondBetweenAtoms(aa[0], aa[1])
            if bond is not None:
                bonds.append(bond.GetIdx())
        mole_name = self.molecule_name.split("/")[-1]
        fname = os.path.join(self.path, f"{mole_name}_break_bond.png")
        return fname, "bonds", bonds, self.molecule_name

    def _gen_torsion_image(self,attr,extra=None):
        fnames = []
        atoms = []
        bonds = []
        n_rdk_atoms = self.rdk.rdkm.GetNumAtoms()

        for tt in getattr(self.molecule, "torsions", []):
            if len(tt) < 4:
                continue
            # Keep only atom ids that exist in the H-removed rdkit molecule used for drawing.
            atom_ids = [ii for ii in tt if isinstance(ii, int) and 0 <= ii < n_rdk_atoms]
            if len(atom_ids) >= 4:
                atoms.append(atom_ids[:4])
            a1, a2 = int(tt[1]), int(tt[2])
            if not (0 <= a1 < n_rdk_atoms and 0 <= a2 < n_rdk_atoms):
                logger.warning(f"skip torsion with out-of-range bond atoms: {tt}")
                continue
            bond = self.rdk.rdkm.GetBondBetweenAtoms(a1, a2)
            if bond is None:
                logger.warning(f"skip torsion without rdkit bond: {tt}")
                continue
            bonds.append(bond.GetIdx())
        if extra is None:
            mole_name = self.molecule_name.split("/")[-1]
            fname = os.path.join(self.path, f"{mole_name}_torsion_scan.png")
            return fname, "bonds", bonds, self.molecule_name
        elif isinstance(extra,list):
            mole_name = self.molecule_name.split("/")[-1]
            fname = os.path.join(self.path, f"{mole_name}_torsion_{'_'.join([str(an) for an in extra])}.png")
            if len(extra) == 4:
                return fname, "atoms", extra, self.molecule_name
            elif len(extra) == 2:
                bond = self.rdk.rdkm.GetBondBetweenAtoms(extra[0], extra[1])
                if bond is None:
                    logger.warning(f"skip torsion bond highlight for invalid bond atoms: {extra}")
                    return fname, "bonds", [], self.molecule_name
                return fname, "bonds", [bond.GetIdx()], self.molecule_name
        else:
            for ii,bond in enumerate(bonds):
                mole_name = self.molecule_name.split("/")[-1]
                fnames.append(os.path.join(self.path, f"{mole_name}_torsion_scan_{ii}.png"))
            if extra == "bond":
                bonds_arr = [[bond] for bond in bonds]
                return fnames, "bonds", bonds_arr, self.molecule_name
            elif extra == "atom":
                return fnames, "atoms", atoms, self.molecule_name

    def _gen_fragmentation_image(self,attr,extra=None):
        attr = self.__frag_label[attr]
        fnames = []
        atoms =[]
        if hasattr(self.molecule,attr):
            for index, value in getattr(self.molecule,attr).items():
                mole_name = self.molecule_name.split("/")[-1]
                fnames.append( os.path.join(self.path, f"{mole_name}_{attr}_{index}.png"))
                atoms.append([an for an in value["components"] if self.molecule.Atoms[an].elem != "H"])
        return fnames, "atoms", atoms, self.molecule_name
    

    def _gen_normal_property_image(self,attr,extra=None):
        if extra is None:
            extra = {"type":"atoms","highlights":[],"fname_pre":"0"}
        stype = extra["type"]
        highlights = extra["highlights"]
        fname_pre = extra["fname_pre"]
        mole_name = self.molecule_name.split("/")[-1]
        fname = os.path.join(self.path,f"{mole_name}_{fname_pre}.png")
        return fname, stype, highlights, self.molecule_name


    __FUNC__ = {
        "normal": _gen_normal_property_image,
        "ff_charge": _gen_atom_property_image,
        "esp_charge": _gen_atom_property_image,
        "mulliken_charge": _gen_atom_property_image,
        "atom_type_name": _gen_atom_property_image,
        "hybrid": _gen_atom_property_image,
        "atom": _gen_atom_property_image,
        "molecule": _gen_molecule_image,
        "break_bond": _gen_molecule_break_bond_image,
        "ring": _gen_molecule_ring_image, 
        "torsion": _gen_torsion_image,
        "EF": _gen_fragmentation_image,
        "SF": _gen_fragmentation_image,
        "CSF": _gen_fragmentation_image,
        "RSF": _gen_fragmentation_image,
        "TF": _gen_fragmentation_image,
        "RF": _gen_fragmentation_image,
        "SSSR": _gen_fragmentation_image,
        "SSSRF": _gen_fragmentation_image,
        "SaF": _gen_fragmentation_image,
        "SkF": _gen_fragmentation_image,
    }
