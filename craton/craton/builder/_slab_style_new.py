import math
import os
import numpy as np
from .packmol_new import Packmol
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

    def __init__(self,style,molecules,box,fn,output_dir,chunks=None,method=None) -> None:
        self.style = style
        self.molecules = molecules
        self.box = box
        self.fn = fn
        self.output_dir = output_dir
        self.chunks = chunks
        self.method = method

    def __call__(self):
        if self.method == "penetrate":
            text,inf,outf, dir = self.penetrate()
        else:
            func = self.__FUNC[self.style]
            text0, inf, outf, dir = self.head_text()
        
            text1 = func(self)
            text = text0 + text1
        
        return Packmol.packmol(text,inf,f"{dir}/{outf}",dir,self.style), outf

    def head_text(self):
        text = ""
        packmol_input = f"{self.output_dir}/{self.fn}.inp"
        packmol_output = f"{self.fn}.xyz"
        text += "seed 0\n"
        text += "tolerance 2.0\n"
        text += "filetype xyz\n\n"
        text += f"pbc {' '.join([str(round(ll,3)) for ll in self.box])}\n"
        text += f"output {packmol_output}\n\n"
        
        if self.chunks is not None:  
            for rr in self.chunks:  
                text += f"structure {rr[4]}\n"
                text += f"  number {rr[1]}\n"
                text += f"  center\n"
                ss = " ".join([str(round(ss,3)) for ss in rr[5]] + ["0.0","0.0","0.0"])
                text += f"    fixed {ss}\n"
                # box = self.box
                # box_position = [ff + 1.0 for ff in box[:3]] + [ff - 1.0 for ff in box[3:]]
                # text_box = " ".join([str(round(ss,3)) for ss in box_position])
                # text += f"  inside box {text_box}\n"
                text += "end structure\n\n"
        return text,packmol_input,packmol_output,self.output_dir
        
    def _run_solution(self,rr,box_flag=False):
        if box_flag:
            box = rr[-2]
        else:
            box = self.box
        text = f"structure {rr[4]}\n"
        text += f"  number {rr[1]}\n"
        box_position = [ff + 1.0 for ff in box[:3]] + [ff - 1.0 for ff in box[3:]]
        text_box = " ".join([str(round(ss,3)) for ss in box_position])
        text += f"  inside box {text_box}\n"
        text += f"end structure\n\n"
        return text

    def solution(self):
        text = ""
        for rr in self.molecules:
            text += self._run_solution(rr)
        return text

    def old_solution(self,):
        text = ""
        for rr in self.molecules:
            text += f"structure {rr[4]}\n"
            text += f"  number {rr[1]}\n"
            box_position = [ff + 1.0 for ff in self.box[:3]] + [ff - 1.0 for ff in self.box[3:]]
            text_box = " ".join([str(round(ss,3)) for ss in box_position])
            text += f"  inside box {text_box}\n"
            text += f"end structure\n\n"
        return text

    def _run_layer(self,rr,box_flag=False):
        if box_flag:
            box = rr[-2]
        else:
            box = self.box
        text = f"structure {rr[4]}\n"
        text += f"  number {rr[1]}\n"
        box_position = [ff + 1.0 for ff in box[:3]] + [ff - 1.0 for ff in box[3:]]
        text_box = " ".join([str(round(ss,3)) for ss in box_position])
        text += f"  inside box {text_box}\n"
        sn = rr[9]
        _tmp = ["0.0" if ii != sn else "1.0" for ii in range(3)]
            
        text += f"  atoms {' '.join([str(an+1) for an in rr[5]])}\n"
        if rr[10] == -1:
            ln = box[sn] + rr[7] + 1.0
            tmp = _tmp + [str(round(ln,3))]
            text += f"    below plane {' '.join(tmp)}\n"
        else:
            ln = box[sn+3] - rr[7] - 1.0
            tmp = _tmp + [str(round(ln,3))]
            text += f"    above plane {' '.join(tmp)}\n"
        text += f"  end atoms\n"
            
        text += f"  atoms {' '.join([str(an+1) for an in rr[6]])}\n"
        if rr[10] == -1:
            ln = box[sn+3] - rr[8] - 1.0
            tmp = _tmp + [str(round(ln,3))]
            text += f"    above plane {' '.join(tmp)}\n"
        else:
            ln = box[sn] + rr[8] + 1.0
            tmp = _tmp + [str(round(ln,3))]
            text += f"    below plane {' '.join(tmp)}\n"
        text += f"  end atoms\n"
            
        text += f"end structure\n\n"
        return text

    def layer(self):
        text = ""
        for rr in self.molecules:
            text += self._run_layer(rr)
        return text

    def old_layer(self):
        text = ""
        for rr in self.molecules:
            text += f"structure {rr[4]}\n"
            text += f"  number {rr[1]}\n"
            box_position = [ff + 1.0 for ff in self.box[:3]] + [ff - 1.0 for ff in self.box[3:]]
            text_box = " ".join([str(round(ss,3)) for ss in box_position])
            text += f"  inside box {text_box}\n"
            sn = rr[9]
            _tmp = ["0.0" if ii != sn else "1.0" for ii in range(3)]
            
            text += f"  atoms {' '.join([str(an+1) for an in rr[5]])}\n"
            if rr[10] == -1:
                ln = self.box[sn] + rr[7] + 1.0
                tmp = _tmp + [str(round(ln,3))]
                text += f"    below plane {' '.join(tmp)}\n"
            else:
                ln = self.box[sn+3] - rr[7] - 1.0
                tmp = _tmp + [str(round(ln,3))]
                text += f"    above plane {' '.join(tmp)}\n"
            text += f"  end atoms\n"
            
            text += f"  atoms {' '.join([str(an+1) for an in rr[6]])}\n"
            if rr[10] == -1:
                ln = self.box[sn+3] - rr[8] - 1.0
                tmp = _tmp + [str(round(ln,3))]
                text += f"    above plane {' '.join(tmp)}\n"
            else:
                ln = self.box[sn] + rr[8] + 1.0
                tmp = _tmp + [str(round(ln,3))]
                text += f"    below plane {' '.join(tmp)}\n"
            text += f"  end atoms\n"
            
            text += f"end structure\n\n"
        return text
            
        Packmol.packmol(text,packmol_input,packmol_output,self.output_dir)

    def penetrate(self):
        text = ""
        packmol_input = f"{self.output_dir}/{self.fn}.inp"
        packmol_output = f"{self.fn}.xyz"
        text += "seed 0\n"
        text += "tolerance 2.0\n"
        text += "filetype xyz\n\n"
        text += f"pbc {' '.join([str(round(ll,3)) for ll in self.box])}\n"
        text += f"output {packmol_output}\n\n"
        rr = self.molecules[0]
        text += f"structure {rr[4]}\n"
        text += f"  number {rr[1]}\n"
        ss = " ".join([str(round(ss,3)) for ss in rr[5]])
        text += f"  inside box {ss}\n"
        text += "end structure\n\n"
        
        for rr in self.molecules:
            if rr[-1] == "layer":
                text += self._run_layer(rr,box_flag=True)
            elif rr[-1] == "solution":
                text += self._run_solution(rr,box_flag=True)
        return text,packmol_input,packmol_output,self.output_dir
    
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
        "penetrate": penetrate,
    }