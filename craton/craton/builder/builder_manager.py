from typing import Any
from copy import deepcopy
from .builder_new import Builder
from pathlib import Path
from ..utils.commons import parallel_run
from ..utils import logger
import time

def BuildManager(molecules,config,parallel=True):
    OBJ = {
        "vacuum": SolutionBuilder,
        "solution": SolutionBuilder,
        "liquid": SolventBuilder,
        "complex": ComplexBuilder,
        "protein": SolutionBuilder,
        "rbfe": ComplexSolutionBuilder,
        "rhfe": SolutionBuilder,
        "rlogp":TwoSolutionBuilder,
        "rlogs": LogsBuilder,
        "mem-rbfe": ComplexSolutionBuilder,
        "cov-rbfe": ComplexSolutionBuilder,
        "mutation": SolutionBuilder,
        "pep-rbfe": ComplexSolutionBuilder,
        "rna-rbfe": ComplexSolutionBuilder,
        "abfe": ComplexSolutionBuilder,
        "ahfe": SolutionBuilder,
        "hfe": SolutionBuilder,
        "alogp": TwoSolutionBuilder,
        "bilayer": BiLayerBuilder,
        "biomembrane": BioMembraneBuilder,
        "mem-protein": ProteinBiLayerBuilder,
    }
    simulation_type = config["MDSetting"]["simulation_type"]
    RUN = OBJ[simulation_type](molecules,config,parallel=parallel)
    return RUN()

class NormalBuilder:
    def __init__(self,molecules,config,style=None,parallel=True):
        self.config = config
        if style is None:
            self.style = self.config["MDSetting"]["simulation_type"]
        self.molecules = molecules
        self.parallel = parallel
        self.building_setting = self.config["BuildingSetting"]
        self.bu_config = {
                    "md_setting": self.config["MDSetting"]["md"],
                    "env_setting": self.config["EnvironmentSetting"],
                    "mdengine": self.config["MDSetting"]["mdengine"],
                    "_anion":self.molecules["anion"],
                    "_cation":self.molecules["cation"],
                    "_alion":self.molecules["alion"],
                    "box_size": self.building_setting["box_size"]
                    #"box_size":[153.0,153.0,153.0]
                    }
        if "this_ff" in self.config["ForceFieldSetting"]:
            self.bu_config["this_ff"] = self.config["ForceFieldSetting"]["this_ff"]
        self.ref_density = self.__check_parameter_building_setting(self.building_setting["density"],"density")
        self.ref_concentration = self.__check_parameter_building_setting(self.building_setting["concentration"],"concentration")
        self.ref_molecules_numbers = self.__check_molecules_numbers_setting(
            self.building_setting.get("molecules_numbers")
        )
        self.counter_ion_flag = self.building_setting["counter_ion"]

    def __check_molecules_numbers_setting(self, value):
        if value is None:
            return None
        if isinstance(value, int):
            if value <= 0:
                logger.warning("Invalid molecules_numbers setting: must be positive")
                return None
            return [value]
        if isinstance(value, list):
            if not value:
                logger.warning("Invalid molecules_numbers setting: empty list")
                return None
            if not all(isinstance(rr, int) and rr > 0 for rr in value):
                logger.warning("Invalid molecules_numbers setting: must be a list of positive integers")
                return None
            return value
        logger.warning("Invalid molecules_numbers setting due to incorrect variable type")
        return None

    def _resolve_molecules_number(self, index, n_total):
        ref = self.ref_molecules_numbers
        if ref is None:
            return None
        if len(ref) == 1 or len(ref) != n_total:
            return int(ref[0])
        return int(ref[index])

    def __check_parameter_building_setting(self,value,pp_name):
        if value is not None:
            if not isinstance(value,list):
                if isinstance(value,float) or isinstance(value,int):
                    return [value]
                else:
                    logger.warning(f"Invalid f{pp_name} setting due to incorrect variable type")
                    return None
            else:
                if all([isinstance(rr,float) or isinstance(rr,int) for rr in value]):
                    return value
                else:
                    logger.warning(f"Invalid f{pp_name} setting due to incorrect variable type")
                    return None
        return None

    __slab_attrs =  {
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
                    "extra_range":None,
                    "chunk_type":None,
                    "chunk_layout":None,
    }
    
    __chunk_attrs = {
                    'counter_ion': None, 
                    'alchemical_ion': None, 
                    'virtual': None, 
                    'solute': None, 
                    'solvent': None, 
                    'component': None, 
                    'box_size': None, 
                    'solution_density': None, 
                    'type': 'solution', 
                    'layout': None, 
                    'extra_range': 5.0, 
                    'chunk_type': 'merge', 
                    'chunk_layout': 'center'
    }

    def _run_solution_type_(self,solutes,output_dir="./",solvents=None,simulation_type=None,idx=None):
        if not isinstance(solutes,list):
            solutes = [solutes]
        if solvents is not None:
            if not isinstance(solvents,list):
                solvents = [solvents]
        
        
        bu_config = deepcopy(self.bu_config)
        bu_config["output_dir"] = output_dir 
        bu_config["simulation_type"] = simulation_type

        if solvents == "vacuum":
            bu_config["solvent"] = "vacuum"
        else:
            if self.ref_density is None:
                bu_config["solvent"] = [solvents[0],getattr(solvents[0],"density",0.75)] if len(solvents) == 1 \
                                    else [[solvent,getattr(solvent,"density",0.75)] for solvent in solvents if solvent is not None]

            else:
                ref_density = deepcopy(self.ref_density)
                if len(solvents) != len(ref_density):
                    ref_density = [ref_density[0] for __ in solvents]
                bu_config["solvent"] = [solvents[0],ref_density[0]] if len(solvents) == 1 \
                                    else [[solvent,density] for solvent,density in zip(solvents,ref_density) if solvent is not None]
              
        if solutes[0] is not None:
            if self.ref_concentration is None:
                bu_config["solute"] = [solutes[0],getattr(solutes[0],"concentration",1)] if len(solutes) == 1 \
                                   else [[solute,getattr(solute,"concentration",1)] for solute in solutes if solute is not None]
            else:
                ref_concentration = deepcopy(self.ref_concentration)
                if len(solutes) != len(ref_concentration):
                    ref_concentration = [ref_concentration[0] for __ in solutes]
                bu_config["solute"] = [solutes[0],ref_concentration[0]] if len(solutes) == 1 \
                                    else [[solute,concentration] for solute,concentration in zip(solutes,ref_concentration) if solute is not None]

        bu_config["name"] = solutes[0].mole_name if solutes[0] is not None else solvents[0].mole_name
        bu_config["style"] = simulation_type
        if idx is not None:
            return bu_config, idx
        else:
            return bu_config
        
    def _run_complex_type_(self,ligand, output_dir="./",solvents=None,protein=None,simulation_type=None,coligands=None):
        if not isinstance(ligand,list):
            ligand = [ligand]
        
        bu_config = deepcopy(self.bu_config)
        bu_config["output_dir"] = output_dir
        bu_config["complex"] = ligand
        bu_config["simulation_type"] = simulation_type

        if coligands is not None:
            bu_config["complex"].extend(coligands)
        bu_config["complex"].append(protein)

        if solvents == "vacuum":
            bu_config["solvent"] = "vacuum"
        else:
            if self.ref_density is None:
                bu_config["solvent"] = [solvents[0],getattr(solvents[0],"density",0.75)] if len(solvents) == 1 \
                                    else [[solvent,getattr(solvent,"density",0.75)] for solvent in solvents if solvent is not None]
                #bu_config["solvent"] = [solvents[0],45000]
            else:
                ref_density = deepcopy(self.ref_density)
                if len(solvents) != len(ref_density):
                    ref_density = [ref_density[0] for __ in solvents]
                bu_config["solvent"] = [solvents[0],ref_density[0]] if len(solvents) == 1 \
                                    else [[solvent,density] for solvent,density in zip(solvents,ref_density) if solvent is not None]
        bu_config["name"] = ligand[0].mole_name
        bu_config["style"] = simulation_type
        return bu_config
    
    def create_build_config(self):
        pass

    def run_single_builder(self,bu_config,idx=None):
        bu = Builder(bu_config,counter_ion_flag = self.counter_ion_flag)
        sm = bu.run()

        #sm.name = f"{ligand.mole_name}-{coligand.mole_name}"
        #sm.style = simulation_type
        #sm.output_dir = output_dir
        if idx is None:
            return sm
        else:
            return sm, idx
        
    def run_builder(self):
        if self.parallel:
            return parallel_run(self.run_single_builder,self.bu_configs,keep_order=True,return_result=True)
        else:
            sms = []
            for bu_config in self.bu_configs:
                sms.append(self.run_single_builder(bu_config))
            return sms

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        
        self.create_build_config()
        return self.run_builder()

class MemberRbfeBuilder(NormalBuilder):
    def __init__(self,molecules,configure, style="",parallel=True):
        super().__init__(molecules,configure,style,parallel)

    def _run_membrane_type_(self,surfactant,output_dir="./",solvents=None,solutes=None,simulation_type=None,idx=None):
        if not isinstance(surfactant,list):
            surfactant = [surfactant]
        
        simulation_type = self.config["MDSetting"]["simulation_type"]
        bu_config["name"] = surfactant[0].mole_name
        bu_config["style"] = simulation_type
        bu_config = deepcopy(self.bu_config)
        bu_config["output_dir"] = output_dir
        bu_config["simulation_type"] = simulation_type
        bu_config["slab_direction"] = self.building_setting["slab_direction"]
        bu_config["slab_deepth"] = self.building_setting["slab_deepth"]
        bu_config["membrane_density"] = self.building_setting["membrane_density"]
        
        solvents = self.molecules["solvent"]
        solutes = self.molecules["solute"]
        
        if solutes is None:
            solutes = [None for __ in bu_config["slab_deepth"]]
        bu_config["slab"] = []
        
        tmp = {"chunk":None,"counter_ion":None,"alchemical_ion":None,"virtual":None,"solute":None,"solvent":None,"component":None,
                    "box_size":None,"solution_density":None,"type":None,"layout":None,"extra_range":None,"chunk_type":None,"chunk_layout":None,}


        if self.ref_density is None:
            ref_density = [None for __ in solvents]
        else:
            ref_density = self.ref_density
        
        membrane_type = "monolayer" if simulation_type in [] else "bilayer"
        
        tmp0 = {"solute":solutes[0],"solvent":[solvents[0],ref_density[0]],
                    "box_size":[None,None,bu_config["slab_deepth"][0]],"solution_density":None,"type":"solution",}
        
        tmp1 = {"solute":None,"solvent":[surfactant,bu_config["membrane_density"]],
                    "box_size":[None,None,None],"solution_density":None,"type":membrane_type,}
        
        tmp2 = {"solute":solutes[1],"solvent":[solvents[1],ref_density[1]],
                    "box_size":[None,None,bu_config["slab_deepth"][1]],"solution_density":None,"type":"solution",}
        bu_config["slab"].append(tmp0)
        bu_config["slab"].append(tmp1)
        bu_config["slab"].append(tmp2)
        
        if simulation_type in []:
        
            tmp3 = {"solute":None,"solvent":[surfactant,bu_config["membrane_density"]],
                    "box_size":[None,None,None],"solution_density":None,"type":"membrane",}
        
            tmp4 = {"solute":solutes[2],"solvent":[solvents[2],ref_density[2]],
                    "box_size":[None,None,bu_config["slab_deepth"][2]],"solution_density":None,"type":membrane_type,}
        
            bu_config["slab"].append(tmp3)
            bu_config["slab"].append(tmp4)
        
        return bu_config

    def create_build_config(self):
        """
        build box for vacuum, ahfe, rhfe, solution,protein, mutation calculation. they have one solution task
        """
        simulation_type = self.config["MDSetting"]["simulation_type"]
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        surfactants = []
        for attr in ["surfactants","molecules"]:
            if self.molecules[attr] is not None:
                surfactants.extend(self.molecules[attr])
        solvents = self.molecules["solvent"]
            
        
        if self.parallel:
            args = [{"output_dir":f'{parent_dir}/{molecule.mole_name}',"solvents":solvents,"simulation_type":simulation_type} for molecule in surfactants]
            self.bu_configs = parallel_run(self._run_membrane_type_,surfactants,kwds=args,keep_order=True,return_result=True)
        else:
            self.bu_configs = []
            for molecule in surfactants:
                self.bu_configs.append(
                    self._run_membrane_type_(
                    molecule,
                    output_dir=f'{parent_dir}/{molecule.mole_name}',
                    solvents=solvents,
                    simulation_type=simulation_type
                    )
                )

class ProteinBiLayerBuilder(NormalBuilder):
    def __init__(self,molecules,configure, style="",parallel=True):
        super().__init__(molecules,configure,style,parallel)

    def _run_membrane_type_(self,ligand,idx=None):
        protein = self.molecules["protein"][0]
        coligand = self.molecules["coligands"]
        surfactant = self.molecules["surfactant"]
        if not isinstance(surfactant,list):
            surfactant = [surfactant]
        
        box_size = self.building_setting["box_size"]
        
        solvents = self.molecules["solvent"]
        solutes = self.molecules["solute"]
        simulation_type = self.config["MDSetting"]["simulation_type"]

        complex = []
        if ligand is not None:
            complex.append(ligand)
        if coligand is not None:
            complex.extend(coligand)
        complex.append(protein)

        bu_config = deepcopy(self.bu_config)
        if ligand is not None:
            bu_config["name"] = ligand.mole_name
            output_dir = f'{self.config["EnvironmentSetting"]["output_directory"]}/{ligand.mole_name}'
        else:
            bu_config["name"] = protein.mole_name
            output_dir = f'{self.config["EnvironmentSetting"]["output_directory"]}/{protein.mole_name}'
        bu_config["style"] = simulation_type
        bu_config["output_dir"] = output_dir
        bu_config["simulation_type"] = simulation_type
        bu_config["slab_direction"] = self.building_setting["slab_direction"]
        bu_config["slab_deepth"] = self.building_setting["slab_deepth"]
        
        membrane_site = []
        for rr in self.building_setting["membrane_site"]:
            tmp = []
            ss = rr.split("-")
            aa = int(ss[0])
            bb = int(ss[1])
            for group in protein.Groups:
                try:
                    if int(group.group_idx) >= aa and int(group.group_idx) <= bb:
                        tmp.append(group.group_str)
                except:
                    pass
            membrane_site.append(tmp) 
        bu_config["membrane_site"] = membrane_site
        bu_config["membrane_site_aided"] = self.building_setting["membrane_site_aided"]
        bu_config["membrane_site_rotation"] = self.building_setting["membrane_site_rotation"]
        bu_config["complex"] = complex
        bu_config["chunk"] = [{
                    "component":[[molecule,1] for molecule in bu_config["complex"]],
                    "extra_range":20.0
                    }]
        bu_config["type"] = "penetrate"
        bu_config["box_size"] = box_size
        
        if solutes is None:
            solutes = [None for __ in range(len(bu_config["membrane_site"]) + 1)]
        
        bu_config["solutes"] = solutes
        
        if not isinstance(self.building_setting["membrane_density"],list):
            _t_ = self.building_setting["membrane_density"]
            bu_config["membrane_density"] = [_t_ for __ in surfactant]
        else:
            bu_config["membrane_density"] = self.building_setting["membrane_density"]
        
        
        if self.ref_density is None:
            ref_density = [None for _ in solvents]
        else:
            ref_density = self.ref_density
        bu_config["solvent"] = ["vacuum" if solvent == "vacuum" else [solvent,getattr(solvent,"density",0.75)] \
                                if ref_density[ii] is None else [solvent,ref_density[ii]] for ii, solvent in enumerate(solvents)]
        
        
        bu_config["surfactant"] = [[sur,bu_config["membrane_density"][ii]] for ii,sur in enumerate(surfactant)]
        
        total_bu_config = deepcopy(bu_config)
        total_bu_config["slab"] = [bu_config]
        return total_bu_config

    def create_build_config(self):
        """
        build box for vacuum, ahfe, rhfe, solution,protein, mutation calculation. they have one solution task
        """
        if self.molecules["ligands"] is None:
            self.molecules["ligands"] = [None]
        
        if self.parallel:
            self.bu_configs = parallel_run(self._run_membrane_type_,self.molecules["ligands"],keep_order=True,return_result=True)
        else:
            self.bu_configs = []
            for molecule in self.molecules["ligands"]:
                self.bu_configs.append(
                    self._run_membrane_type_(
                    molecule,)
                )



class BioMembraneBuilder(NormalBuilder):
    def __init__(self,molecules,configure, style="",parallel=True):
        super().__init__(molecules,configure,style,parallel)

    def _run_membrane_type_(self,surfactant,idx=None):
        if not isinstance(surfactant,list):
            surfactant = [surfactant]
        
        box_size = self.building_setting["box_size"]
        
        solvents = self.molecules["solvent"]
        solutes = self.molecules["solute"]
        simulation_type = self.config["MDSetting"]["simulation_type"]
        output_dir = f'{self.config["EnvironmentSetting"]["output_directory"]}/{surfactant[0].mole_name}'

        bu_config = deepcopy(self.bu_config)
        bu_config["name"] = surfactant[0].mole_name
        bu_config["style"] = simulation_type
        bu_config["output_dir"] = output_dir
        bu_config["simulation_type"] = simulation_type
        bu_config["slab_direction"] = self.building_setting["slab_direction"]
        bu_config["slab_deepth"] = self.building_setting["slab_deepth"]
        if not isinstance(self.building_setting["membrane_density"],list):
            _t_ = self.building_setting["membrane_density"]
            bu_config["membrane_density"] = [_t_ for __ in surfactant]
        else:
            bu_config["membrane_density"] = self.building_setting["membrane_density"]
        
        
        
        if solutes is None:
            solutes = [None for __ in bu_config["slab_deepth"]]
        bu_config["slab"] = []


        if self.ref_density is None:
            ref_density = [None for _ in solvents]
        else:
            ref_density = self.ref_density
        
        
        tmp0 = {
                    "solute":solutes[0],
                    "solvent": "vacuum" if solvents[0] == "vacuum" else [solvents[0],getattr(solvents[0],"density",0.75)] if ref_density[0] is None else [solvents[0],ref_density[0]],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"solution",
                    "slab_direction":bu_config["slab_direction"],
                    "slab_deepth":bu_config["slab_deepth"][0],
                }
        
        tmp1 = {
                    "solute":None,
                    "solvent":[[sur,bu_config["membrane_density"][ii]] for ii,sur in enumerate(surfactant)],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"layer",
                    "slab_direction":bu_config["slab_direction"],
                    "philic_direction": -1,
                }
        
        tmp2 = {
                    "solute":None,
                    "solvent":[[sur,bu_config["membrane_density"][ii]] for ii,sur in enumerate(surfactant)],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"layer",
                    "slab_direction":bu_config["slab_direction"],
                    "philic_direction": 1,
                }
        
        tmp3 = {
                    "solute":solutes[1],
                    "solvent": "vacuum" if solvents[1] == "vacuum" else [solvents[1],getattr(solvents[1],"density",0.75)] if ref_density[1] is None else [solvents[1],ref_density[1]],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"solution",
                    "slab_direction":bu_config["slab_direction"],
                    "slab_deepth":bu_config["slab_deepth"][1],
                }
        
        bu_config["slab"].append(tmp0)
        bu_config["slab"].append(tmp1)
        bu_config["slab"].append(tmp2)
        bu_config["slab"].append(tmp3)

        return bu_config

    def create_build_config(self):
        """
        build box for vacuum, ahfe, rhfe, solution,protein, mutation calculation. they have one solution task
        """
            
        if self.parallel:
            self.bu_configs = parallel_run(self._run_membrane_type_,self.molecules["surfactant"],keep_order=True,return_result=True)
        else:
            self.bu_configs = []
            for molecule in self.molecules["surfactant"]:
                self.bu_configs.append(
                    self._run_membrane_type_(
                    molecule,)
                )

class BiLayerBuilder(NormalBuilder):
    def __init__(self,molecules,configure, style="",parallel=True):
        super().__init__(molecules,configure,style,parallel)

    def _run_membrane_type_(self,surfactant,idx=None):
        if not isinstance(surfactant,list):
            surfactant = [surfactant]
        
        solvents = self.molecules["solvent"]
        if len(solvents) == 1:
            solvents = ["vacuum",solvents[0],"vacuum"]
        
        
        solutes = self.molecules["solute"]
        simulation_type = self.config["MDSetting"]["simulation_type"]
        output_dir = f'{self.config["EnvironmentSetting"]["output_directory"]}/{surfactant[0].mole_name}'
        
        bu_config = deepcopy(self.bu_config)
        bu_config["name"] = surfactant[0].mole_name
        bu_config["style"] = simulation_type
        bu_config["output_dir"] = output_dir
        bu_config["simulation_type"] = simulation_type
        bu_config["slab_direction"] = self.building_setting["slab_direction"]
        bu_config["slab_deepth"] = self.building_setting["slab_deepth"]
        if not isinstance(self.building_setting["membrane_density"],list):
            _t_ = self.building_setting["membrane_density"]
            bu_config["membrane_density"] = [_t_ for __ in surfactant]
        else:
            bu_config["membrane_density"] = self.building_setting["membrane_density"]
        
        box_size = self.building_setting["box_size"]
        

        
        if solutes is None:
            solutes = [None for __ in bu_config["slab_deepth"]]
        bu_config["slab"] = []


        if self.ref_density is None:
            ref_density = [None for _ in solvents]
        else:
            ref_density = self.ref_density
        
        tmp0 = {
                    "solute":solutes[0],
                    "solvent": "vacuum" if solvents[0] == "vacuum" else [solvents[0],getattr(solvents[0],"density",0.75)] if ref_density[0] is None else [solvents[0],ref_density[0]],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"solution",
                    "slab_direction":bu_config["slab_direction"],
                    "slab_deepth":bu_config["slab_deepth"][0],
                }
        
        tmp1 = {
                    "solute":None,
                    "solvent":[[sur,bu_config["membrane_density"][ii]] for ii,sur in enumerate(surfactant)],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"layer",
                    "slab_direction":bu_config["slab_direction"],
                    "philic_direction": 1,
                }
        
        tmp2 = {
                    "solute":solutes[1],
                    "solvent": "vacuum" if solvents[1] == "vacuum" else [solvents[1],getattr(solvents[1],"density",0.75)] if ref_density[1] is None else [solvents[1],ref_density[1]],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"solution",
                    "slab_direction":bu_config["slab_direction"],
                    "slab_deepth":bu_config["slab_deepth"][1],
                }
        
        tmp3 = {
                    "solute":None,
                    "solvent":[[sur,bu_config["membrane_density"][ii]] for ii,sur in enumerate(surfactant)],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"layer",
                    "slab_direction":bu_config["slab_direction"],
                    "philic_direction": -1,
                }
        
        tmp4 = {
                    "solute":solutes[2],
                    "solvent": "vacuum" if solvents[2] == "vacuum" else [solvents[2],getattr(solvents[2],"density",0.75)] if ref_density[2] is None else [solvents[2],ref_density[2]],
                    "box_size":box_size,
                    "solution_density":None,
                    "type":"solution",
                    "slab_direction":bu_config["slab_direction"],
                    "slab_deepth":bu_config["slab_deepth"][2],
                }
        
        bu_config["slab"].append(tmp0)
        bu_config["slab"].append(tmp1)
        bu_config["slab"].append(tmp2)
        bu_config["slab"].append(tmp3)
        bu_config["slab"].append(tmp4)

        return bu_config

    def create_build_config(self):
        """
        build box for vacuum, ahfe, rhfe, solution,protein, mutation calculation. they have one solution task
        """
            
        if self.parallel:
            self.bu_configs = parallel_run(self._run_membrane_type_,self.molecules["surfactant"],keep_order=True,return_result=True)
        else:
            self.bu_configs = []
            for molecule in self.molecules["surfactant"]:
                self.bu_configs.append(
                    self._run_membrane_type_(
                    molecule,)
                )

class SolutionBuilder(NormalBuilder):
    def __init__(self, molecules,configure, style="",parallel=True):
        super().__init__(molecules,configure, style,parallel)
    
    def create_build_config(self):
        """
        build box for vacuum, ahfe, rhfe, solution,protein, mutation calculation. they have one solution task
        """
        simulation_type = self.config["MDSetting"]["simulation_type"]
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        molecules = []
        for attr in ["molecules","ligands"]:
            if self.molecules[attr] is not None:
                molecules.extend(self.molecules[attr])

        if simulation_type == "vacuum":
            solvents = "vacuum"
        else:
            solvents = self.molecules["solvent"]
        
        if self.parallel:
            args = [{"output_dir":f'{parent_dir}/{molecule.mole_name}',"solvents":solvents,"simulation_type":simulation_type} for molecule in molecules]
            self.bu_configs = parallel_run(self._run_solution_type_,molecules,kwds=args,keep_order=True,return_result=True)
        else:
            self.bu_configs = []
            for molecule in molecules:
                self.bu_configs.append(
                    self._run_solution_type_(
                    molecule,
                    output_dir=f'{parent_dir}/{molecule.mole_name}',
                    solvents=solvents,
                    simulation_type=simulation_type
                    )
                )

class TwoSolutionBuilder(NormalBuilder):
    def __init__(self, molecules, configure, style="", parallel=True):
        super().__init__(molecules,configure, style, parallel)
    
    def create_build_config(self):
        """
        build box for rlogp/d alogp/d calculation. they have one solution task
        """
        simulation_type = self.config["MDSetting"]["simulation_type"]
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        molecules = []
        for attr in ["molecules","ligands"]:
            if self.molecules[attr] is not None:
                molecules.extend(self.molecules[attr])
        
        self.bu_configs = []
        
        for molecule in molecules:
            Path(f'{parent_dir}/{molecule.mole_name}').mkdir(exist_ok=True)
            self.bu_configs.append(
                self._run_solution_type_(
                    molecule,
                    output_dir=f'{parent_dir}/{molecule.mole_name}/water',
                    solvents=self.molecules["solvent"][0],
                    simulation_type=simulation_type
                )
            )
            self.bu_configs.append(
                self._run_solution_type_(
                    molecule,
                    output_dir=f'{parent_dir}/{molecule.mole_name}/oil',
                    solvents=self.molecules["solvent"][1],
                    simulation_type=simulation_type
                )
            )
        
class SolventBuilder(NormalBuilder):
    def __init__(self, molecules, configure, style="", parallel=True):
        super().__init__(molecules,configure, style,parallel)

    def _run_solution_type_(self, solutes, output_dir="./", solvents=None, simulation_type=None, idx=None, molecule_index=0, n_molecules=1):
        bu_config = super()._run_solution_type_(
            solutes, output_dir=output_dir, solvents=solvents,
            simulation_type=simulation_type, idx=idx,
        )
        if self.ref_molecules_numbers is not None and solvents not in (None, "vacuum"):
            if not isinstance(solvents, list):
                solvents = [solvents]
            bu_config["solvent"] = [
                solvents[0], self._resolve_molecules_number(molecule_index, n_molecules)
            ]
        return bu_config
    
    def create_build_config(self):
        """
        build box for pure liquid calculation. they have one solution task
        """
        simulation_type = self.config["MDSetting"]["simulation_type"]
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        molecules = []
        for attr in ["molecules","ligands"]:
            if self.molecules[attr] is not None:
                molecules.extend(self.molecules[attr])
            
        self.bu_configs = []
        n_molecules = len(molecules)
        for molecule_index, molecule in enumerate(molecules):
            self.bu_configs.append(
                self._run_solution_type_(
                    None,
                    solvents=molecule,
                    output_dir=f'{parent_dir}/{molecule.mole_name}',
                    simulation_type=simulation_type,
                    molecule_index=molecule_index,
                    n_molecules=n_molecules,
                )
            )

class LogsBuilder(NormalBuilder):
    def __init__(self, molecules,configure, style="", parallel=True):
        super().__init__(molecules, configure, style,parallel)
    
    def create_build_config(self):
        """
        build box for rlogs calculation. they have two solution task
        """
        simulation_type = self.config["MDSetting"]["simulation_type"]
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        molecules = []
        for attr in ["molecules","ligands"]:
            if self.molecules[attr] is not None:
                molecules.extend(self.molecules[attr])
            
        self.bu_configs = []
        for molecule in molecules:
            Path(f'{parent_dir}/{molecule.mole_name}').mkdir(exist_ok=True)
            self.bu_configs.append(
                self._run_solution_type_(
                    None,
                    output_dir=f'{parent_dir}/{molecule.mole_name}/amorphous',
                    solvents=molecule,
                    simulation_type=simulation_type
                )
            )
            self.bu_configs.append(
                self._run_solution_type_(
                    molecule,
                    output_dir=f'{parent_dir}/{molecule.mole_name}/solution',
                    solvents=self.molecules["solvent"],
                    simulation_type=simulation_type
                )
            )

class ComplexBuilder(NormalBuilder):
    def __init__(self, molecules, config=None, style=None, parallel=True):
        super().__init__(molecules, config, style, parallel)
        
    def create_build_config(self):
        """
        build box for complex calculation. they have two solution task
        """
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        ligands = self.molecules["ligands"]
        
        self.bu_configs = []
        for ligand in ligands:
            Path(f'{parent_dir}/{ligand.mole_name}').mkdir(exist_ok=True)
            self.bu_configs.append(
                self._run_complex_type_(
                    ligand,
                    output_dir = f'{parent_dir}/{ligand.mole_name}',
                    solvents = self.molecules["solvent"],
                    protein = self.molecules["protein"][0],
                    simulation_type = self.config["MDSetting"]["simulation_type"],
                    coligands = self.molecules["coligands"],
                )
            )
        
class ComplexSolutionBuilder(NormalBuilder):
    def __init__(self, molecules, config=None, style=None, parallel=True):
        super().__init__(molecules, config, style, parallel)
        
    def create_build_config(self):
        """
        build box for complex calculation. they have two solution task
        """
        parent_dir = self.config["EnvironmentSetting"]["output_directory"]
        ligands = self.molecules["ligands"]
        
        self.bu_configs = []
        for ligand in ligands:
            Path(f'{parent_dir}/{ligand.mole_name}').mkdir(exist_ok=True)
            self.bu_configs.append(
                self._run_complex_type_(
                    ligand,
                    output_dir = f'{parent_dir}/{ligand.mole_name}/bfe',
                    solvents=self.molecules["solvent"],
                    protein=self.molecules["protein"][0],
                    simulation_type=self.config["MDSetting"]["simulation_type"],
                    coligands=self.molecules["coligands"],
                )
            )
            
            self.bu_configs.append(
                self._run_solution_type_(
                    ligand,
                    solvents=self.molecules["solvent"],
                    simulation_type=self.config["MDSetting"]["simulation_type"],
                    output_dir=f'{parent_dir}/{ligand.mole_name}/hfe'
                )
            )

