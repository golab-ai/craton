import itertools
import os
import sys
from copy import deepcopy

import itertools
from ..chem.molecule import Molecule
from ..chem import FormatMolecule as FM
from ..chemkit import Structure as Stru


config = {"reference":"molecules[0]",
          "ref_rule":{
              "def_bonds":[[19,18],[4,3]],
              "def_atoms":[],
          },
          #"frag":{
          #    "frag_type":"elem_frag",
          #    "ring_number":2,
          #    "ring_size":[6,6],
          #    "ring_property":["ar1","ar4"]},
          "frag":"frag.csv",
          "replace_rule":[
              {
                  "element":["C"],
                  "ring_property":["ar4"],
                  "position":[["local:CF","m"]]
              },
              {
                  "element":["C"],
                  "ring_property":["ar1"]
              }
          ],
          "spatial_rule":[[0,1,"G"]],
          "output_directory":"./enumeration"
         }

class Enumeration:
    def __init__(self,config=None) -> None:
        self.ref_molecule_all = config["reference"]
        self.ref_rule = config["ref_rule"]
        self.frag_rule = config["frag"]
        self.replace_rule = config["replace_rule"]
        self.spatial_rule = config["spatial_rule"]
        self.output_directory = config["output_directory"]

    def get_ref_anchor(self):
        self.ref_anchor = []
        self.ref_move_atoms = []
        if self.ref_rule["def_atoms"]:
            pass
        elif self.ref_rule["def_bonds"]:
            tmps = []
            for bond in self.ref_rule["def_bonds"]:
                self.ref_anchor.append(bond[0])
                tmp = self.ref_molecule_all.find_side_componend(bond[1],bond[0])
                tmp.append(bond[1])
                tmps.append(tmp)
            _tmp = set(tmps[0])
            for tmp in tmps[1:]:
                _tmp = _tmp.intersection(set(tmp))
            self.ref_move_atoms = list(_tmp)
        Atoms = []
        _labels = {}
        ii = -1
        for atom in self.ref_molecule_all.Atoms:
            if atom.ID not in self.ref_move_atoms:
               ii += 1
               _labels[atom.ID] = ii
               Atoms.append(atom)
        for atom in Atoms:
            connect = []
            bond_type =[]
            for ii,an in enumerate(atom.connectivity):
                if an not in self.ref_move_atoms:
                    connect.append(_labels[an])
                    bond_type.append(atom.bond_type[ii])
            atom.connectivity = connect
            atom.bond_type = bond_type
            atom.ID = _labels[atom.ID] 
        self.ref_anchor = [_labels[an] for an in self.ref_anchor]
        self.ref_molecule = Molecule("molecule")
        self.ref_molecule.Atoms = Atoms
    
    def get_frags(self):
        if isinstance(self.frag_rule,dict):
            inputs = None
            selector = {}
            for kk,vv in self.frag_rule.items():
                if kk in ["ring_size","ring_property"]:
                    selector[kk] = {"$in":[list(rr) for rr in itertools.permutations(vv)]}
                else:
                    selector[kk] = vv
            extra_var = {"selector":selector}
        else:
            inputs = self.frag_rule
            extra_var=None
        self.frag_molecules = FM._parse(inputs,extra_var=extra_var)
        self.frag_molecules = Stru._basic_structure_analyze(self.frag_molecules,ignore_existing=True,parallel=True)
        FM._convert(self.frag_molecules,otype="png",opath=f"{self.output_directory}/frag")
        FM._convert(self.frag_molecules,otype="csv",ofilename="frag",opath=self.output_directory)
    
    def parse_replace_rule(self):
        _label_spatial = {"s":0,"o":1,"m":2,"p":3,"A":1,"B":2,"G":3,"D":4,"E":5,"Z":6}
        replace_rule = []
        spatial_rule = []
        for rule in self.replace_rule:
            tmp = {}
            for rr,vv in rule.items():
                if rr == "position":
                    tmp["position"] = [[rrr[0],_label_spatial[rrr[1]]]for rrr in vv]
                else:
                    ss = rr.split(":")
                    kk = ss[0] 
                    flag  = ss[1] if len(ss) > 1 else None
                    if flag is None:
                        tmp[kk] = [vv,True]
                    elif flag == "not":
                        tmp[kk] = [vv, False]
            replace_rule.append(tmp)
        self.replace_rule = replace_rule
        
        for rule in self.spatial_rule:
            spatial_rule.append([rule[0],rule[1],_label_spatial[rule[2]]])
        self.spatial_rule = spatial_rule
    
    def match_rule_for_frag_atom(self,rule,atom,molecule):
        for kk,vv in rule.items():
            if kk in ["ring_size","ring_property"]:
                if (len(set(getattr(atom,f"has_{kk}",[])).intersection(set(vv[0]))) > 0) != vv[1]:
                    return False
            elif kk in ["position"]:
                for rr in vv:
                    kkk = rr[0].split(":")[0]
                    vvv = rr[0].split(":")[1]
                    _tmp = [an.ID for an in molecule.Atoms if getattr(an,kkk,None) == vvv]
                    if len(_tmp) > 0:
                        
                        nns = [molecule.calc_bond_distance(atom.ID,aa) - 1 for aa in _tmp]
                        if min(nns) != rr[1]:
                            return False
                    else:
                        return False 
                    
            else:
                if (getattr(atom,kk,None) in vv[0]) != vv[1]:
                    return False
        return True
        
    def get_frag_anchor(self):
        self.cut_atoms = {}
        for molecule in self.frag_molecules:
            _total_tmp = []
            flag = True
            for rule in self.replace_rule:
                tmp = []
                for atom in molecule.Atoms:
                    connect_elems = [molecule.Atoms[an].elem for an in atom.connectivity]
                    H_ndx = connect_elems.index("H") if "H" in connect_elems else -1
                    if H_ndx != -1:
                        this_flag = self.match_rule_for_frag_atom(rule,atom,molecule)
                        #if self.match_rule_for_frag_atom(rule,atom):
                        if this_flag:
                            tmp.append([atom.ID,atom.connectivity[H_ndx]])
                if len(tmp) == 0:
                    flag = False
                    break
                else:
                    _total_tmp.append(tmp)
            if flag:
                self.cut_atoms[molecule.mole_name] = _total_tmp
        self.frag_anchor = []
        for molecule in self.frag_molecules:
            tmp = []
            if molecule.mole_name not in self.cut_atoms:
                continue
            for arr in itertools.product(*self.cut_atoms[molecule.mole_name]):
                flag = True
                for rr in self.spatial_rule:
                    nk = molecule.calc_bond_distance(arr[rr[0]][0],arr[rr[1]][0]) - 1
                    if nk != rr[2]:
                        flag = False
                        break
                if flag:
                    tmp.append(arr)
            if tmp:
                self.frag_anchor.append([molecule,tmp])
                
    def get_new_molecule(self):
        self.target_molecules = []
        nn = len(self.ref_molecule.Atoms)
        for rr in self.frag_anchor:
            molecule_template = rr[0]
            frag_anchors = rr[1]
            for anchors in frag_anchors:
                del_atoms = [an[1] for an in anchors]
                anchor_atoms = [an[0] for an in anchors]
                frag_molecule = deepcopy(molecule_template)
                this_molecule = deepcopy(self.ref_molecule)
                Atoms = []
                ii = nn - 1
                _label_ndx = {}
                
                for atom in frag_molecule.Atoms:
                    if atom.ID not in del_atoms:
                        ii += 1
                        _label_ndx[atom.ID] = ii
                        Atoms.append(atom)
                
                for atom in Atoms:
                    connect = []
                    bond_type =[]
                    for jj,an in enumerate(atom.connectivity):
                        if an not in del_atoms:
                            connect.append(_label_ndx[an])
                            bond_type.append(atom.bond_type[jj])
                    
                    anchor_ndx = anchor_atoms.index(atom.ID) if atom.ID in anchor_atoms else -1
                    if anchor_ndx != -1:
                        
                        ref_anchor_an = self.ref_anchor[anchor_ndx]
                        connect.append(ref_anchor_an)
                        bond_type.append("1")
                        this_molecule.Atoms[ref_anchor_an].connectivity.append(_label_ndx[atom.ID])
                        this_molecule.Atoms[ref_anchor_an].bond_type.append("1")
                        
                    atom.connectivity = connect
                    atom.bond_type = bond_type
                    atom.ID = _label_ndx[atom.ID]
                this_molecule.Atoms.extend(Atoms)
                this_molecule.create_topols()
                self.target_molecules.append(this_molecule)

        

    def get_smiles(self):
        self.smiles = [mol.smiles for mol in self.target_molecules]
        molecules = FM._parse(self.smiles)
        
        FM._convert(molecules,otype="png",opath=self.output_directory)
        FM._convert(molecules,otype="csv",ofilename="new_molecules",opath=self.output_directory)
        return self.smiles
                        
    def run(self):
        self.parse_replace_rule()
        self.get_ref_anchor()
        self.get_frags()
        self.get_frag_anchor()
        self.get_new_molecule()
        return self.get_smiles()

