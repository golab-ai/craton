import math
import os
import numpy as np
from .packmol import Packmol
from copy import deepcopy

__attrs =  [ "counter_ion","alchemical_ion","virtual","solute","solvent",]

def arrangement(molecules,distance=2.0,direction=2):
    max_coor = max([atom.coor[direction] for atom in molecules[0].Atoms])
    point_coor = max_coor + distance
    min_coor = min([atom.coor[direction] for atom in molecules[1].Atoms])
    transfer = point_coor - min_coor
    molecule_1 = deepcopy(molecules[1])
    for atom in molecule_1.Atoms:
        atom.coor[direction] = atom.coor[direction] + transfer
    return molecules[0],molecule_1

    



class SlabStyle:

    def __init__(self,style,molecules,box,fn,output_dir,chunks=None) -> None:
        self.style = style
        self.molecules = molecules
        self.box = box
        self.fn = fn
        self.output_dir = output_dir
        self.chunks = chunks

    def __call__(self):
        func = self.__FUNC[self.style]
        return func(self)

    def solution(self,):
        sections = []
        if self.chunks is not None:  
            for rr in self.chunks:  
                box_position = ["inside","box"] + self.box
                sections.append([rr[4],rr[1],box_position,{"center":rr[5]}])

        for rr in self.molecules:
            box_position = ["inside","box"] + [ff + 1.0 for ff in self.box[:3]] + [ff - 1.0 for ff in self.box[3:]]
            
            sections.append([rr[4],rr[1],box_position,None])

        Packmol.packmol(sections,self.fn,self.output_dir,pbc=self.box)
        
    def layer(self):
        packmol_input = f"{self.output_dir}/{self.fn}.inp"
        packmol_output = f"{self.fn}.xyz"
        text = "seed 0\n"
        text += "tolerance 2.0\n"
        text += "filetype xyz\n\n"
        text += f"pbc {' '.join([str(round(ll,3)) for ll in self.box])}\n"
        text += f"output {packmol_output}\n\n"
    
        if self.chunks is not None:
            for rr in self.chunks:
                text += f"structure {rr[4]}\n"
                text += f"  number {rr[1]}\n"
                text += f"  center\n"
                text = " ".join([str(round(ss,3)) for ss in rr[5]["center"]] + ["0.0","0.0","0.0"])
                text += f"  fixed "
                

        sections = []
        if self.chunks is not None:  
            for rr in self.chunks:  
                box_position = ["inside","box"] + self.box
                sections.append([rr[4],rr[1],box_position,{"center":rr[5]}])
        
        for rr in self.molecules:
                box_position = ["inside","box"]
                sections.append([rrr[0],rrr[1],["inside box"]+rrr[2],rrr[3],rrr[4]])
        Packmol.packmol(sections,self.fn,self.output_dir,pbc=self.box)
            

    def _monolayer(self):
        sections = []
        if self.chunks is not None:
            for rr in self.chunks:
                box_position = ["inside","box"] + self.box
                sections.append([rr[4],rr[1],box_position,{"center":rr[5]}])
        for rr in self.molecules:
            pass

    def _micelle(self):
        pass

    def _vesicle(self):
        pass

    def _sphere(self):
        pass

    def _dimer(self):
        pass

    

    __FUNC = {
        "solution":solution,
        "layer": layer,
    }