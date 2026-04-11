from pathlib import Path
import os,sys,json
import numpy as np
import math
from ..utils import logger
from .packmol import Packmol
from copy import deepcopy
#from compuchem.chemistry.format import CallRdkit

from .. import CRATON_CONFIGURE
from ..chem.molecule import Molecule
from ..chem.chemsystem import System
from ..chem import FormatMolecule as FM
from ..property.unit import Unit
from ..property.expt_equ import require_expt_data
from ..force_field import MolForceField as MFF #checkout_force_field

from ._slab_style import SlabStyle
from ._slab_style import arrangement


DEFAULT_SURFACTANT_DENSITY = 0.5
DEFAULT_SOLVENT = ["water",1.0]
CHUNK_GAP = 10.0

molecule_templeate = CRATON_CONFIGURE["path"]["molecule"]
JOB_DIR = os.getcwd()
#ion_templeate = CRATON_CONFIGURE["path"]["ion"]

class Builder:
    """
    slab类型：solution, membrane, surface
    slab组成：chunk,solute,solvent,counter ion, alchemical ion, virtual
    chunk: 
        biomacromolecule, 
        oligomer (protein + small molecule), 
        aggregate(micelle), 
        substrate(surface), 
        nano-structure(carbon nanotube)
    注意：
        builder模块只生成在一定区间内，满足一定条件的随机结构。如果是有序结构，请通过其他方式构建
    slab:[chunk,solvent,solute,type=solution,box_size=None,]
    solvent: [molecule,number,density or concentration,density or concentration unit, xyz_file]
    chunk: [component[molecule,number],chunk_layout,chunk_type,solvent,solute,type=solution,layout,box_size,]]
    """

    def __init__(self,configure,style="",counter_ion_flag=True):
        self.style = style
        self.configure = configure
        self.counter_ion_flag = counter_ion_flag
        self._counter_anion = self.configure["_anion"][0] if "_anion" in self.configure else "Cl-"
        self._counter_cation = self.configure["_cation"][0] if "_cation" in self.configure else "Na+"

        self._alion = self.configure["_alion"][0] if "_alion" in self.configure else "ANA"
        if "output_dir" in self.configure:
            if self.configure["output_dir"] is not None:
                Path(self.configure["output_dir"]).mkdir(exist_ok=True)
        
        self._initi_config()

    __attrs =  [
                    "chunk","counter_ion","alchemical_ion","virtual","solute","solvent","component",
                    "box_size","solution_density","type","layout","extra_range","chunk_type","chunk_layout",
                  ]

    __default = {
        "chunk":None,
        "counter_ion":None,
        "alchemical_ion":None,
        "virtual":None,
        "solute":None,
        "solvent":None,
        "component":None,
        "box_size":None,
        "solution_density":None,
        "type":"solution",
        "layout":None,
        "extra_range":5.0,
        "chunk_type":"merge",
        "chunk_layout":"center",
    }

    def _write_xyz_file(self,molecules, xyz_file):
        num_atoms = sum([len(molecule.Atoms) for molecule in molecules])
        with open(xyz_file, "w") as f:
            f.write(f"{num_atoms}\n")
            f.write("cpy create xyz file \n")
            for molecule in molecules:
                for atom in molecule.Atoms:
                    f.write(f"{atom.elem} {atom.coor[0]:.3f} {atom.coor[1]:.3f} {atom.coor[2]:.3f}\n")

    def _check_molecule(self,molecule):
        if molecule == "vaccum":
            return None,None
        template_file = f"{molecule_templeate}/{molecule}.mtx"
        molecule = molecule if isinstance(molecule,Molecule) \
                    else FM._parse(template_file)[0] if os.path.exists(template_file) \
                    else FM._parse(molecule)[0]
        if molecule.style not in ["protein"]:
            xyz_file = f"{molecule.mole_name}"
        else:
            xyz_file = molecule.mole_name
        self._write_xyz_file([molecule],f"{self.configure['output_dir']}/{xyz_file}.xyz")
        return molecule,f"{xyz_file}.xyz"

    def _check_molecule_arrs(self,slab):
        def _check_value(arr):
            _tmp = [None for _ in arr]
            for ii,rr in enumerate(arr):
                if isinstance(rr,int):
                    _tmp[ii] = 1
                elif isinstance(rr,float):
                    _tmp[ii] = 2
                elif isinstance(rr,str):
                    _tmp[ii] = 3
            return _tmp

        for typ in self.__attrs[1:7]:
            if slab[typ] is not None:
                if np.array(slab[typ]).ndim == 1:
                    slab[typ] = [slab[typ]]
                for ii,rr in enumerate(slab[typ]):
                    if len(rr) == 1:
                        slab[typ][ii] = [rr[0],None,None,None,None]
                    elif len(rr) > 1:
                        _tmp = [rr[0],None,None,None,None]
                        _tmp_ = _check_value(rr[1:])
                        
                        for jj,vv in enumerate(_tmp_):
                            _tmp[vv] = rr[jj+1]
                        if _tmp[2] is not None and _tmp[3] is None:
                            _tmp[3] = "g/ml" if typ == "solvent" else "%"
                        
                        slab[typ][ii] = _tmp
                    molecule = deepcopy(slab[typ][ii][0])
                    slab[typ][ii][0],slab[typ][ii][4]= self._check_molecule(molecule)

    def _check_attrs(self,slab):
        for typ in self.__attrs[:12]:
            if typ not in slab:
                slab[typ] = self.__default[typ]

        if slab["type"] == "solution":
            if slab["solvent"] is None:
                slab["solvent"] = [DEFAULT_SOLVENT]

        if slab["solvent"] == "vacuum" or slab["solvent"][0] == "vacuum":
            slab["solvent"] = None
        if slab["solvent"] is not None:
            if not isinstance(slab["solvent"][0],Molecule):
                if slab["solvent"][0][0] == "vacuum":
                    slab["solvent"] = None
        
        self._check_molecule_arrs(slab)

        if slab["chunk"] is not None:
            for chunk in slab["chunk"]:
                for typ in self.__attrs[1:]:
                    if typ not in chunk:
                        chunk[typ] = self.__default[typ]
                if any(chunk[typ] is not None for typ in self.__attrs[1:6]) and chunk["layout"] is None:
                    chunk["layout"] = "outside"
                self._check_molecule_arrs(chunk)

    def _initi_config(self,):
        #### 输出文件夹
        if "output_dir" not in self.configure:
            self.configure["output_dir"] = "./"
        else:    
            if self.configure["output_dir"] is None:
                self.configure["outptu_dir"] = "./"

        #### 温度和压强
        if "temperature" not in self.configure:
            self.configure["temperature"] = 298.15
        if "pressure" not in self.configure:
            self.configure["pressure"] = 1.0
            
        #### 盒子大小   
        if "box_size" not in self.configure:
            self.configure["box_size"] = None

        if self.configure["box_size"] is None:
            if "length" in self.configure:
                self.configure["box_size"] = [self.configure["length"] for _ in range(3)]
            if "box" in self.configure:
                self.configure["box_size"] = self.configure["box"]

        #### slab 方向
        _label = {"x":0,"y":1,"z":2}

        if "slab_direction" not in self.configure:
                self.configure["slab_direction"] = 2
        else:
            if not isinstance(self.configure["slab_direction"],int):
                self.configure["slab_direction"] = _label[self.configure["slab_direction"]]

        #### 生成一个slab
        if "slab" not in self.configure:
            self.configure["slab"] = [{typ:self.configure[typ] if typ in self.configure else None for typ in self.__attrs}]
            if self.configure["slab"][0]["type"] is None:
                self.configure["slab"][0]["type"] = "solution"
            if "protein" in self.configure:
                self.configure["slab"][0]["chunk"] = [{
                    "component":[[self.configure["protein"],1]],
                    "extra_range":20.0,
                    }]
            if "complex" in self.configure:
                self.configure["slab"][0]["chunk"] = [{
                    "component":[[molecule,1] for molecule in self.configure["complex"]],
                    "extra_range":20.0
                    }]

        ### slab处理
        for slab in self.configure["slab"]:
            self._check_attrs(slab)
            if slab["chunk"] is None:
                if slab["solute"] is not None:
                    if len(slab["solute"]) == 1 and slab["solute"][0][1] == 1:
                        _tmp_ = {}
                        for typ in self.__attrs[1:]:
                            _tmp_[typ] = self.__default[typ]
                        slab["chunk"] = [_tmp_]
                        slab["chunk"][0]["component"] = slab["solute"]
                        slab["solute"] = None   

    def _get_net_charge(self,molecules):
        net_charge = int(round(sum([rr[0].net_charge*rr[1] for rr in molecules])))
        return net_charge
        #if net_charge >= 0:
        #    return net_charge
        #else:
        #    return net_charge * -1

    def _get_excess_charge(self,molecules):
        excess_charge = 0.0000
        for rr in molecules:
            if hasattr(rr[0],"excess_charge_flag"):
                m_excess_charge = (sum([round(atom.ff_charge,4) for atom in rr[0].Atoms]) - rr[0].net_charge) * rr[1]
                excess_charge += m_excess_charge
        return excess_charge

    def _get_box_volume(self,box):
        if len(box) == 3:
            volume = np.prod(box)
        elif len(box) == 6:
            volume = np.prod([box[3]-box[0],box[4]-box[1],box[5]-box[2]])
        return volume
    
    def get_box_size_molecule_size(self,*molecules,cube_flag=False):
        molecule_size = []
        length = [10.0,10.0,10.0]
        for molecule in molecules:
            if isinstance(molecule,list):
                Atoms = []
                for sub_molecule in molecule:
                    Atoms += sub_molecule.Atoms
            else:
                Atoms = molecule.Atoms
            min_coor = [
                        min(Atoms[0].coordinate[0],Atoms[1].coordinate[0]),
                        min(Atoms[0].coordinate[1],Atoms[1].coordinate[1]),
                        min(Atoms[0].coordinate[2],Atoms[1].coordinate[2]),
                        ]

            max_coor = [
                        max(Atoms[0].coordinate[0],Atoms[1].coordinate[0]),
                        max(Atoms[0].coordinate[1],Atoms[1].coordinate[1]),
                        max(Atoms[0].coordinate[2],Atoms[1].coordinate[2]),
                        ]
            for atom in Atoms[2:]:
                min_coor = [
                            min(min_coor[0],atom.coordinate[0]),
                            min(min_coor[1],atom.coordinate[1]),
                            min(min_coor[2],atom.coordinate[2]),
                            ]
                max_coor = [
                            max(max_coor[0],atom.coordinate[0]),
                            max(max_coor[1],atom.coordinate[1]),
                            max(max_coor[2],atom.coordinate[2]),
                            ]
            molecule_size.append([max_coor[ii] - min_coor[ii] for ii in range(3)])
            length = [max(length[ii],molecule_size[-1][ii]) for ii in range(3)]


        if cube_flag:
            le = max(length)
            length = [le,le,le]
        
        return [le for le in length]

    def _get_box_size_based_concentration(self,solvent,others=None):
        others_mass = 0.0
        if others is not None:
            others_mass = sum([rr[0].mass*rr[1] for rr in others])
        length = math.ceil(((solvent[1]*solvent[0].mass + others_mass)/solvent[2]/0.602)**(1.0/3.0))
        return [length,length,length]

    def _get_box_size_based_solution(self,ref_solutes,ref_box_size,others=None):
        solutes = deepcopy(ref_solutes)
        box_size = [30.0,30.0,30.0] if ref_box_size is None else [max(rr,30.0) for rr in ref_box_size]
        
        box_sizes = []
        for solute in solutes:
            if solute[1] is None:
                solute[1] = 1
            if solute[2] is not None:
                box_sizes.append(self._get_box_size_based_concentration(solute,others=others))
            else:
                box_sizes.append(box_size)
        max_l = max([ln[0] for ln in box_sizes])
                
        if max_l ** 3 >= np.prod(box_size):
            box_size = [max(max_l, box_size[ii]) for ii in range(3)]

        return box_size
    
    def get_box_size_molecule_position(self,molecules,cube_flag=False):
        max_x, max_y, max_z = -9e10, -9e10, -9e10
        min_x, min_y, min_z = 9e10, 9e10, 9e10
        for molecule in molecules:
            for atom in molecule.Atoms:
                min_x = min(min_x, atom.coor[0])
                min_y = min(min_y, atom.coor[1])
                min_z = min(min_z, atom.coor[2])

                max_x = max(max_x, atom.coor[0])
                max_y = max(max_y, atom.coor[1])
                max_z = max(max_z, atom.coor[2])

        box_size = np.array([max_x, max_y, max_z]) - np.array([min_x, min_y, min_z])
        if cube_flag:
            le = max(box_size)
            box_size = [le,le,le]
        return box_size 

    def _get_density(self,molecule):
        try:
            density = require_expt_data(molecule,"density",temperatures=self.configure["temperature"],pressures=self.configure["pressure"])
        except:
            density = None
        if density is not None:
            return density
        else:
            return 0.75, "g/ml"

    def setting_concentration(self,rr,box_size=None,others=None,nn=None,solution_density=None,solvent_density=None,solvent_flag=False):
        if solvent_flag:
            if rr[1] is not None:
                return
            if rr[2] is None:
                if box_size is not None and rr[1] is not None:
                    solute_mass = 0.0
                    if others is not None:
                        solute_mass = sum([rr[0].mass*rr[1] for rr in others])
                    rr[2] = (rr[1]*rr[0].mass + solute_mass)/0.602/np.prod(box_size)
                    rr[3] = "g/ml"
                else:
                    vvs = self._get_density(rr[0])
                    
                    rr[2] = vvs[0] / nn
                    rr[3] = vvs[1]
            if rr[2] is not None:
                rr[2] = Unit(rr[2],rr[3],"g/ml",extra={"mass":rr[0].mass})()
                rr[3] = "g/ml"    
        else:
            if rr[2] is not None:
                if solution_density is not None:
                    rr[2] = Unit(rr[2],rr[3],"g/ml",extra={"mass":rr[0].mass,"density":solution_density})()
                    rr[3] = "g/ml"
                else:
                    rr[2] = Unit(rr[2],rr[3],"g/ml",extra={"mass":rr[0].mass,"solvent_density":solvent_density})()
                    rr[3] = "g/ml"

    def slab_setting_concentration(self,slab):
        solvent_density = None
        if slab["solvent"] is not None:
            nn = len(slab["solvent"])
            for rr in slab["solvent"]:
                others = None
                if "chunk" in slab and slab["chunk"] is not None:
                    others =[[rr[0],rr[1]] for chunk in slab["chunk"] for rr in chunk["molecules"]] 
                self.setting_concentration(rr,box_size=slab["box_size"],others=others,nn=nn,solvent_flag=True)
            solvent_density = sum([rr[2] for rr in slab["solvent"] if rr[2] is not None])
        
        for attr in self.__attrs[1:5]:
            if slab[attr] is not None:
                for rr in slab[attr]:
                    self.setting_concentration(rr,solution_density=slab["solution_density"],solvent_density=solvent_density)
            
    def get_number_base_on_concentration(self, box_size, concentration,mass,others=None):  # 
        volume = self._get_box_volume(box_size)

        solute_mass = 0.0
        if others is not None:
            solute_mass = sum([rr[0].mass*rr[1] for rr in others])
        
        solvent_number = int((concentration * volume * 0.602 - solute_mass) / mass)
        return solvent_number

    def get_number_of_layer(self,box, concentration,direction):
        area = np.prod([ln for ii,ln in enumerate(box) if ii != direction ])
        return int(area * concentration)
        
    def get_molecule_number(self,rr,box,others=None,slab_type="solution",direction=2):
        if slab_type != "layer":
            if rr[0] is not None:
                if rr[2] is None:
                    if rr[1] is None:
                        rr[1] = 1
                else:
                    rr[1] = self.get_number_base_on_concentration(box,rr[2],rr[0].mass,others=others)
        else:
            if rr[1] is None:
                rr[1] = self.get_number_of_layer(box,rr[2],direction)

    def _chunk_run(self,ii,jj,chunk):
        #
        if chunk["chunk_type"] == "merge":
            chunk["box_size"] = self.get_box_size_molecule_position([rr[0] for rr in chunk["component"]])
            chunk["xyz_file"] = f"chunk-{ii}-{jj}.xyz"
            chunk["chunk_xyz_file"] = chunk["xyz_file"]
            chunk["molecules"] = [[rr[0],1] for rr in  chunk["component"]]
            self._write_xyz_file([rr[0] for rr in chunk["component"]],f"{self.configure['output_dir']}/{chunk['chunk_xyz_file']}")
        elif chunk["chunk_type"] == "arrangment":
            new_molecules = arrangement([rr[0] for rr in chunk["component"]],distance=chunk["distance"],direction=chunk["direction"])
            if new_molecules[0].mole_name == new_molecules[1].mole_name:
                chunk["molecules"] = [[new_molecules[0],2]]
            else:
                chunk["molecules"] = [[molecule,1] for molecule in new_molecules]
            chunk["box_size"] = self.get_box_size_molecule_position(new_molecules)
            chunk["xyz_file"] = f"chunk-{ii}-{jj}.xyz"
            chunk["chunk_xyz_file"] = chunk["xyz_file"]
            self._write_xyz_file([new_molecules],f"{self.configure['output_dir']}/{chunk['chunk_xyz_file']}")

        elif chunk["chunk_type"] in ["micelle"]:
            pass

        box_flag = True
        if any(chunk[typ] is not None for typ in self.__attrs[1:6]):
            if chunk["layout"] == "outside":
                box_flag = False
                extra_range = 0.0 if chunk["extra_range"] is None else chunk["extra_range"]
                chunk["box_size"] = [extra_range + rr for rr in chunk["box_size"]]
                if len(chunk["box_size"]) == 6:
                    box = chunk["box_size"]
                else:
                    box = [0.0,0.0,0.0]
                    box.extend(chunk["box_size"])
                self.slab_setting_concentration(chunk)
                for attr in self.__attrs[1:6]:
                    if chunk[attr] is not None:
                        for rr in chunk[attr]:
                            self.get_molecule_number(rr,box,others=chunk["molecules"]) 
                chunk["xyz_file"] = f"chunk-extend-{ii}-{jj}.xyz"
                
                #os.chdir(self.configure["output_dir"])
                ss = []
                for typ in self.__attrs[1:6]:
                    if chunk[typ] is not None:
                        ss.extend(chunk[typ])
                SlabStyle("solution",
                        [rr for typ in self.__attrs[1:6] if chunk[typ] is not None for rr in chunk[typ]],
                        box,
                        f"chunk-extend-{ii}-{jj}",
                        self.configure["output_dir"],
                        chunks=[[None,1,None,None,chunk["chunk_xyz_file"],
                                 [box[0] / 2.0 + box[3] / 2.0, box[1] / 2.0 + box[4] / 2.0, box[2] / 2.0 + box[5] / 2.0]]]
                        )()
                #os.chdir(JOB_DIR)
                chunk["molecules"].extend([[rr[0],rr[1]] for typ in self.__attrs[1:6] if chunk[typ] is not None for rr in chunk[typ]])
            else:
                pass

        if box_flag:
            if chunk["extra_range"] is not None:
                chunk["box_size"] = [rr + chunk["extra_range"] for rr in chunk["box_size"]]

    def _chunk_arrangement(self,nn,ss):
        if ss == "solution":
            n = math.ceil(nn**(1.0/3.0))
            if n**3 - nn - n**2 >= 0:
                if n**3 - nn - n**2 - n*(n-1) >= 0:
                    extend_arr = [n,n-1,n-1]
                else:
                    extend_arr = [n,n,n-1]
            else:
                extend_arr = [n,n,n]
        elif ss == "membrane":
            pass
        
        return extend_arr

    def _total_chunk_run(self,ii,chunks,slab_type="solution"):    
        for jj,chunk in enumerate(chunks):
            self._chunk_run(ii,jj,chunk)
        chunk_n = jj + 1
        chunk_box = [0.0,0.0,0.0]
        for chunk in chunks:
            chunk_box = [max(chunk["box_size"][ii],chunk_box[ii]) for ii in range(3)]
        extend_arr = self._chunk_arrangement(chunk_n,slab_type)
        if chunk_n > 1:
            chunk_box = [(chunk_box[ii]+CHUNK_GAP)*extend_arr[ii] for ii in range(3)]
        return chunk_box, extend_arr 

    def get_layer_box_size(self,slab):
        nn = slab["slab_direction"]
        mm = len(slab["solvent"])
        setting_box = [30.0,30.0,30.0]
        if slab["box_size"] is not None:
            setting_box = slab["box_size"]
        lengths = []
        for rr in slab["solvent"]:
            phorbic_coors = [rr[0].Atoms[an].coordinates[nn] for an in rr[0].phobic_terminal]
            phorbic_coors = sorted(phorbic_coors,key=lambda x:abs(x),reverse=True)
            philic_coors = [rr[0].Atoms[an].coordinates[nn] for an in rr[0].philic_terminal]
            philic_coors = sorted(philic_coors,key=lambda x:abs(x),reverse=True)
            lengths.append(abs(philic_coors[0] - phorbic_coors[0]))
            if rr[1] is None:
                if rr[2] is not None:
                    rr[2] = rr[2]/mm
                else:
                    rr[2] = DEFAULT_SURFACTANT_DENSITY/mm
                rr[3] = "N/nm^2"
        length = max(lengths)
        setting_box[nn] = length
        return setting_box
        
    def _slab_box_size(self,ii,slab):
        # 根据chunk的设置box大小
        if slab["type"] == "solution":
            chunk_box = [0.0,0.0,0.0]
            if slab["chunk"] is not None:
                chunk_box,extend_arr = self._total_chunk_run(ii,slab["chunk"]) 
                slab["chunk_accum"] = extend_arr
            slab["molecules"] = []
            if slab["chunk"] is not None:
                for chunk in slab["chunk"]:
                    slab["molecules"].extend(chunk["molecules"])
            self.slab_setting_concentration(slab)
            # 根据solvent的设置box大小
            solvent_box = [0.0,0.0,0.0]
            if slab["solvent"] is not None:
                solvent_box = self._get_box_size_based_solution(slab["solvent"],slab["box_size"],others=slab["molecules"])

            # 根据solute的设置box大小
            solute_box = [0.0,0.0,0.0]        
            if slab["solute"] is not None:
                solute_box = self._get_box_size_based_solution(slab["solute"],slab["box_size"],others=slab["molecules"])

            # 设置的box大小
            setting_box = [30.0,30.0,30.0]
            if slab["box_size"] is not None:
                setting_box = slab["box_size"]
            box_sizes = sorted([setting_box,solvent_box,solute_box,chunk_box],key=lambda x:np.prod(x))
            box_size = [max(box_sizes[-1][ii],chunk_box[ii]) for ii in range(3)]
            if "slab_deepth" in slab and slab["slab_deepth"] is not None:
                box_size = [ln if ii != slab["slab_direction"] else slab["slab_deepth"] for ln in box_size]
            return box_size
        elif slab["type"] == "layer":
            return self.get_layer_box_size(slab)

    def _get_chunk_center(self,box,arr):
        chunk_centers = []
        for ii in range(1,arr[0]*2,2):
            for jj in range(1,arr[1]*2,2):
                for kk in range(1,arr[2]*2,2):
                    chunk_centers.append([box[3]*ii/(arr[0]*2),box[4]*jj/(arr[1]*2),box[5]*kk/(arr[2]*2)])
        return chunk_centers

    def _get_alchem_ion(self,molecules):
        alchem_charge = 0.0
        for rr in molecules:
            _net_charge_m1 = sum([getattr(atom,"ff_charge") for atom in rr[0].Atoms])
            _net_charge_m2 = sum([getattr(atom,"ff_charge_m2",atom.ff_charge) for atom in rr[0].Atoms])
            if abs(_net_charge_m1 - _net_charge_m2) > 0.0001:
                alchem_charge += (_net_charge_m1 * 1 + _net_charge_m2 * -1) * rr[1]

        if abs(alchem_charge) > 0.001: #### != 0.0:
            logger.info("the charge is not zero between A and B ligand: %.4f" %alchem_charge)
            setattr(self._alion.Atoms[0],"ff_charge_m2",alchem_charge)
            #setattr(self._alion.Atoms[0],"atom_type_name_m2","_D")
            setattr(self._alion.Atoms[0],"atom_type_name_m2",self._alion.Atoms[0].atom_type_name)
            setattr(self._alion.Atoms[0],"mass_m2",self._alion.Atoms[0].mass)
            setattr(self._alion.Atoms[0],"parameter_m2",self._alion.Atoms[0].parameter)
            alion_m,alion_xyz = self._check_molecule(self._alion)
            return [[alion_m,1,None,None,alion_xyz]]
        else:
            return None

    def _get_counter_ion(self,molecules):
        _net_charge = self._get_net_charge(molecules)
        if _net_charge != 0:
            if _net_charge < 0:
                ion = self._counter_cation
                #ion = "Na+"
            else:
                ion = self._counter_anion 
                #ion = "Cl-"
            ion_m,ion_xyz = self._check_molecule(ion)
            return [[ion_m,abs(_net_charge),None,None,ion_xyz]]
        else:
            return None

    def _slab_number(self,slab_boxs):
        pass


    def _slab_run(self,ii,slab):
        for rr in slab["molecules"]:
            self.molecules.append(rr[0])
            self.molecule_number.append(rr[1])
        _molecules = deepcopy(slab["molecules"])
        for attr in self.__attrs[1:6]:
            if slab[attr] is not None:
                for rr in slab[attr]:
                    self.get_molecule_number(rr,slab["box_size"],others=slab["molecules"],slab_type=slab["type"],direction=slab["slab_direction"]) 
                    _molecules.append([rr[0],rr[1]])
                    #self.molecules.append(rr[0])
                    #self.molecule_number.append(rr[1])
        this_counter_ion = None
        if self.counter_ion_flag:
            this_counter_ion = self._get_counter_ion(_molecules)
        if this_counter_ion is not None:
            slab["counter_ion"] = this_counter_ion
            
        this_alchem_ion = None
        if self.counter_ion_flag:    
            this_alchem_ion = self._get_alchem_ion(_molecules)
        if this_alchem_ion is not None:
            slab["alchemical_ion"] = this_alchem_ion
        #_net_charge = self._get_net_charge(_molecules)



        if slab["type"] == "solution":
            sections = []
            chunks=[]
            if slab["chunk"] is not None:
                chunk_centers = self._get_chunk_center(slab["box_size"],slab["chunk_accum"])
                for jj ,chunk in enumerate(slab["chunk"]):
                    chunks.append([None,1,None,None,chunk["xyz_file"],chunk_centers[jj]])

            for attr in self.__attrs[1:6]:
                if slab[attr] is not None:
                    for rr in slab[attr]:
                        if rr[0] is not None:
                            #slab["molecules"].append([rr[0],rr[1]])
                            sections.append(rr)
                            self.molecules.append(rr[0])
                            self.molecule_number.append(rr[1])

            slab["xyz_file"] = f"slab_{ii}.xyz"
            SlabStyle("solution",sections,slab["box_size"],f"slab_{ii}",self.configure["output_dir"],chunks=chunks)()

    def _make_system(self,system_info):
        sm = System()
        sm.set_info(system_info)
        for attr in ["name","style","simulation_type","output_dir","md_setting","env_setting","mdengine"]:
            if attr in self.configure:
                setattr(sm,attr,self.configure[attr])
        if "this_ff" in self.configure:
            this_ff = self.configure["this_ff"]
            setattr(sm,"ff",MFF.checkout_force_field(sm.molecules,this_ff))
        return sm

    def run(self):
        self.molecules = []
        self.molecule_number = []

        sd = self.configure["slab_direction"]
        slab_boxs = []
        for ii,slab in enumerate(self.configure["slab"]):
            slab_boxs.append(self._slab_box_size(ii,slab))
        max_box = [max([rr[0] for rr in slab_boxs]),max([rr[1] for rr in slab_boxs]),max([rr[2] for rr in slab_boxs])]
        total_boxs = []
        box = [0.0,0.0,0.0]
        box.extend([slab_boxs[0][jj] if jj == sd else max_box[jj] for jj in range(3)])
        total_boxs.append(box)
        self.configure["slab"][0]["box_size"] = box
        self._slab_run(0,self.configure["slab"][0])

        for ii,sbox in enumerate(slab_boxs[1:]):
            pre_ln = self.configure["slab"][ii]["box_size"][sd + 3]
            box_0 = [pre_ln + 1.0 if jj == sd else 0.0 for jj in range(3)]
            box_0 = box_0 + [box_0[jj] + sbox[jj] if jj == sd else max_box[jj] for jj in range(3)]
            total_boxs.append(box_0)

            box = [0.0,0.0,0.0]
            box.extend([sbox[jj] if jj == sd else max_box[jj] for jj in range(3)])

            self.configure["slab"][ii+1]["box_size"] = box
            self._slab_run(ii+1,self.configure["slab"][ii+1])

        if len(self.configure["slab"]) == 1:
            #molecules = self.configure["slab"][0]["molecules"]
            coordinates = Packmol._get_coordinates(self.configure["slab"][0]["xyz_file"],self.configure["output_dir"])
            lattics = self.configure["slab"][0]["box_size"][3:]
        else:
            all_box = total_boxs[0][:3]
            all_box.extend(total_boxs[-1][3:])
            #molecules = []
            
            total_slabs = []
            for ii,slab in enumerate(self.configure["slab"]):
                #molecules.extend(slab["molecules"])
                sbox = total_boxs[ii]
                center = [sbox[0] / 2.0 + sbox[3] / 2.0,sbox[1] / 2.0 + sbox[4] / 2.0,sbox[2] / 2.0 + sbox[5] / 2.0]
                total_slabs.append([None,1,None,None,slab["xyz_file"],center])

            SlabStyle("solution",[],all_box,"total",self.configure["output_dir"],chunks=total_slabs)()
            
            coordinates = Packmol._get_coordinates("total.xyz",self.configure["output_dir"])
            lattics = all_box
        system_info = {"lattics": lattics, 
                       "molecules":self.molecules, 
                       "molecule_number":self.molecule_number,
                       "coordinates": coordinates}
        syses = self._make_system(system_info)
        return syses
        #return self._make_system(system_info)
    