from typing import Any
import os
import time
from copy import deepcopy
from craton import molxpert as MX
from pathlib import Path

import sys

JOB_TYPE = {
                "RFE":["rbfe","rhfe","rlogp","rlogs","mem-rbfe","cov-rbfe","pep-rbfe","rna-rbfe"],
                "AFE":["abfe","ahfe"],
                "box":["solution","liquid"]
                  }


def Simulation_run(config):
    OBJ = {
        "vacuum": NormalSimulation,
        "solution": NormalSimulation,
        "liquid": NormalSimulation,
        "complex": NormalSimulation,
        "protein": NormalSimulation,
        "rbfe": RFESimulation,
        "rhfe": RFESimulation,
        "rlogp":RFESimulation,
        "rlogs": RFESimulation,
        "mem-rbfe": RFESimulation,
        "cov-rbfe": RFESimulation,
        "mutation": ProteinRMutationSimulation,
        "pep-rbfe": PepRbfeSimulation,
        "rna-rbfe": PepRbfeSimulation,
        "abfe": AFESimulation,
        "ahfe": AFESimulation,
        "hfe": AFESimulation,
        "alogp": AFESimulation,
        "bilayer": NormalSimulation,
        "biomembrane": NormalSimulation,
        "mem-protein": NormalSimulation,
    }
    simulation_type = config["simulation_type"]
    RUN = OBJ[simulation_type](config)
    RUN()


class MDSimulation:
    """
    MD simulation 分为通用和固定流程两大类型。
    固定流程设定“simulation_type"参数：
        rbfe:
        rhfe:
        abfe:
        ahfe:
        complex:
        protein:
        vaccum:
        solution:
        alogp:
        rlogp:
        alogs:
        rlogs:

    所有流程分为：
        1.  parse input,
        2.  molecule create, 
        3.  force field, 
        4.  docking, optional
        5.  alignment, optional
        6.1 initial pair network, optional
        7.  atom mapping, optional
        8.  silimarlity, optional
        6.2 final pair network, optional
        9. topol integral, optional
        10. build box, 
        11. property, optional
        12. md input file
        13. bash script
        14. run md
        15. analyze data
    """
    
    def __init__(self, config):
        self.user_config = config
        self.molecules = {"protein":None,
                     "chunk":None,
                     "solute":None,
                     "solvent":None,
                     "ligands":None,
                     "coligands":None,
                     "molecules":None,
                     "anion":None,
                     "cation":None,
                     "alion":None,
                     "surfactant":None,
                     }
        self.__JOB_TYPE = {
                "RFE":["rbfe","rhfe","rlogp","rlogs","mem-rbfe","cov-rbfe","pep-rbfe","rna-rbfe"],
                "AFE":["abfe","ahfe"],
                "box":["solution","liquid"]
                  }
        self._default_gmx_res = {
                            "solute":"SUT",
                            "solvent":"SOL",
                            "ligands":"LIG",
                            "coligands":"CLI",
                            "molecules":"MOL",
                            "anion":"CLA",
                            "cation":"NAN",
                            "alion":"ANA",
                            "surfactant":"SUR",
                            }
        
    def update_config(self):
        self.config = MX.update_configure(self.user_config)
        Path(self.config['EnvironmentSetting']['output_directory']).mkdir(exist_ok=True)
        self.molecule_paths = [self.config['path']['molecule']] + [f"{self.config['path']['molecule']}/{dd}" 
                              for dd in os.listdir(self.config["path"]["molecule"]) 
                              if os.path.isdir(f"{self.config['path']['molecule']}/{dd}" )]
        
    def molecule_create(self):
        for key in self.molecules.keys():
            molecule_input = self.config["MoleculeFileSetting"][key]
            if molecule_input is not None:
                
                molecules = MX.molecule_create(molecule_input,template_path = self.molecule_paths, show_figure=False,parallel=False)
                molecules = MX.molecule_structure(molecules)
                for molecule in molecules:
                    if molecule.style not in ["pdb","protein","template","dna","rna","DNA","RNA","Protein"]:
                        if not hasattr(molecule,"gmx_residue_name"):
                            molecule.gmx_residue_name = self._default_gmx_res[key] if key in self._default_gmx_res else molecule.mole_name[:2]
                self.molecules[key] = molecules  
        
    def get_force_field(self):
        _this_ff = MX.force_field_read(self.config["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"],use_scalevdw=self.config["ForceFieldSetting"]["use_scalevdw"])
        self.this_ff = {"general":_this_ff["general"]}
        if "pairwise" in _this_ff:
            self.this_ff["pairwise"] = _this_ff["pairwise"]
        if "scalevdw" in _this_ff:
            self.this_ff["scalevdw"] = _this_ff["scalevdw"]
        self.config["ForceFieldSetting"]["this_ff"] = self.this_ff
        for key,molecules in self.molecules.items():
            if molecules is not None:
                _molecules_ = MX.atom_type(molecules)
                self.molecules[key] = MX.grasp_force_field(
                                                    #self.molecules[key],
                                                   _molecules_,
                                                   atom_type_file=self.config["ForceFieldSetting"]["DEFAULT_TYPING_FILE"],
                                                   force_field_file=self.config["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"],
                                                   empi_ff_flag=True,
                                                   charge_method=self.config["ForceFieldSetting"]["charge_method"],
                                                   use_scalevdw=self.config["ForceFieldSetting"]["use_scalevdw"],
                                                   #reassign_atom_type=True,
                                                   return_ff=False,
                                                   parallel=False,
                                                   )
        
    def molecule_docking(self):
        pass
    def molecule_alignment(self):
        pass
    def pair_network_init(self):
        pass
    def atom_mapping(self):
        pass
    def molecule_similiarity(self):
        pass
    def pair_network_final(self):
        pass
    def molecule_topol_integral(self):
        pass
    def force_field_postprepare(self):
        for typ,molecules in self.molecules.items():
            if molecules is not None:
                for molecule in molecules:
                    for atom in molecule.Atoms:
                        if hasattr(atom,"parameter") and not hasattr(atom,"_ff_parameter"):
                            atom._ff_parameter = {"name":atom.atom_type_used_name,"mass":atom.mass,
                                              "pstyle":atom.pstyle,"fix_parameter":[],"parameter":atom.parameter,"ptag":"null","pscore":"nan","pcount":"nan"}
                            if hasattr(atom,"parameter_m2"):
                                atom._ff_parameter_m2 = {"name":atom.atom_type_name_m2,
                                                     "pstyle":atom.pstyle,"mass":atom.mass_m2,"fix_parameter":[],
                                                     "parameter":atom.parameter_m2,"ptag":"null","pscore":"nan","pcount":"nan"}

                    for term in ["Bonds","Angles","Dihedrals","Impropers"]:
                        for item in getattr(molecule,term,[]):
                            if hasattr(item,"parameter") and not hasattr(item,"_ff_parameter"):
                                item._ff_parameter = {"name":item.atom_type_used_name,
                                                  "pstyle":item.pstyle,"fix_parameter":[],"parameter":item.parameter,"ptag":"dummy","pscore":"nan","pcount":"nan"}
    
    def build_box(self):
        self.sms = MX.builder(self.molecules,config=self.config)
    
    def property_setting(self):
        pass
    def save_job_info(self):
        MX.write_job_infos(self.sms) 
    def md_input_files(self):
        MX.write_md_input_files(self.sms)
    def bash_script(self):
        MX.write_bash_files(self.sms)
    def md_run(self):
        pass
    def analyze_result(self):
        pass
    
    def __call__(self):
        self.update_config()
        self.molecule_create()
        self.get_force_field()
        self.molecule_docking()
        self.molecule_alignment()
        self.pair_network_init()
        self.atom_mapping()
        self.molecule_similiarity()
        self.pair_network_final()
        self.molecule_topol_integral()
        self.force_field_postprepare()
        self.build_box()
        self.property_setting()
        self.save_job_info()
        self.md_input_files()
        self.bash_script()
        self.md_run()
        self.analyze_result()

class PepRbfeSimulation(MDSimulation):
    def __init__(self, config):
        super().__init__(config)

    def molecule_alignment(self):
        if self.config["AlignmentSetting"]["mutation"] is not None:
            changes = self.config["AlignmentSetting"]["mutation"]
        else:
            changes = []

        if self.config["AlignmentSetting"]["sequences"] is not None:
            changes.extend(MX.protein_sequence_mutation(self.molecules["ligands"][0],self.config["AlignmentSetting"]["sequences"]))    
        pre_pep = self.molecules["ligands"][0]
        pre_name = pre_pep.mole_name
        for ii,change in enumerate(changes):
            this_protein = MX.protein_process(pre_pep,change)
            this_protein.mole_name = f"{pre_name}_{ii}"
            self.molecules["ligands"].append(this_protein)               
        #self.molecules["ligands"] = MX.molecule_structure(self.molecules["ligands"])
            #self.molecules["ligands"].append(MX.protein_process(self.molecules["ligands"][0],change))
        self._before_integral_ligands = {molecule.mole_name: deepcopy(molecule) for molecule in self.molecules["ligands"]}

    def old_molecule_alignment(self):
        changes = []
        if self.config["AlignmentSetting"]["mutation"] is not None:
            mutations = self.config["AlignmentSetting"]["mutation"]
            if not isinstance(mutations[0],list):
                changes.append(mutations + ["mutation"])
            else:
                for mutation in mutations:
                    changes.append(mutation + ["mutation"])
        if self.config["AlignmentSetting"]["modify"] is not None:
            modifies = self.config["AlignmentSetting"]["modify"]
            if not isinstance(modifies[0],list):
                changes.append(modifies + ["modify"])
            else:
                for modify in modifies:
                    changes.append(modify + ["modify"])

        for change in changes:
            self.molecules["ligands"].append(MX.protein_process(self.molecules["ligands"][0],change))
        self._before_integral_ligands = {molecule.mole_name: deepcopy(molecule) for molecule in self.molecules["ligands"]}

    def pair_network_init(self):
        self.ggs = []
        for molecule in self.molecules["ligands"][1:]:
            self.ggs.append([self.molecules["ligands"][0],molecule])
    
    def atom_mapping(self):
        for gg in self.ggs:
            atom_mapping = MX.protein_atom_mapping(gg[0],gg[1])
            gg.append(atom_mapping)
    
    def molecule_topol_integral(self):
        topologys = MX.dual_topology(
            self.config["MDSetting"]["simulation_type"],
            self.ggs,
            output_directory=self.config["EnvironmentSetting"]["output_directory"],
            parallel=True,
        )
        self.molecules["ligands"] = [topology[0] for topology in topologys]
    
    def property_setting(self):
        for sm in self.sms:
            fep_type = sm.molecules[0].dual_topology_type
            sm.md_setting["free_energy_auixed"]["lambdas"] = MX.get_fep_lambda(
                sm.md_setting["free_energy"],
                fep_type=fep_type,
                mixed_lambda=sm.md_setting["free_energy_auixed"]["mixed_lambda"],
                is_relative=sm.md_setting["free_energy_auixed"]["is_relative"]
            )    

class ProteinRMutationSimulation(PepRbfeSimulation):
    def __init__(self, config):
        super().__init__(config)
    
    def molecule_alignment(self):
        changes = []
        if self.config["AlignmentSetting"]["mutation"] is not None:
            mutations = self.config["AlignmentSetting"]["mutation"]
            if not isinstance(mutations[0],list):
                changes.append(mutations + ["mutation"])
            else:
                for mutation in mutations:
                    changes.append(mutation + ["mutation"])
        if self.config["AlignmentSetting"]["modify"] is not None:
            modifies = self.config["AlignmentSetting"]["modify"]
            if not isinstance(modifies[0],list):
                changes.append(modifies + ["modify"])
            else:
                for modify in modifies:
                    changes.append(modify + ["modify"])
        

        self.molecules["ligands"] = self.molecules["protein"]
        for change in changes:
            self.molecules["ligands"].append(MX.protein_process(self.molecules["ligands"][0],change))
        self._before_integral_ligands = {molecule.mole_name: deepcopy(molecule) for molecule in self.molecules["ligands"]}

class ProteinRBFESimulation(PepRbfeSimulation):
    def __init__(self, config):
        super().__init__(config)
    
    def molecule_alignment(self):
        changes = []
        if self.config["AlignmentSetting"]["mutation"] is not None:
            mutations = self.config["AlignmentSetting"]["mutation"]
            if not isinstance(mutations[0],list):
                changes.append(mutations + ["mutation"])
            else:
                for mutation in mutations:
                    changes.append(mutation + ["mutation"])
        if self.config["AlignmentSetting"]["modify"] is not None:
            modifies = self.config["AlignmentSetting"]["modify"]
            if not isinstance(modifies[0],list):
                changes.append(modifies + ["modify"])
            else:
                for modify in modifies:
                    changes.append(modify + ["modify"])
        

        self.molecules["ligands"] = self.molecules["protein"]
        self.molecules["protein"] = [self.molecules["protein"]]
        for change in changes:
            self.molecules["protein"].append(MX.protein_process(self.molecules["protein"][0],change))
        self._before_integral_ligands = {molecule.mole_name: deepcopy(molecule) for molecule in self.molecules["protein"]}

class RFESimulation(MDSimulation):
    def __init__(self, config):
        super().__init__(config)
    
    def molecule_alignment(self):
        if self.config["AlignmentSetting"]["align_flag"]:
            if self.rbfe_fep_setting.bias_nodes:
                align_node = self.rbfe_fep_setting.bias_nodes[0]
            else:
                align_node = self.rbfe_fep_setting.align_index
            ligand_align = LigandAlign(
                self.molecules["ligands"],
                output_file="_aligned.sdf",
                reference_ligand=align_node,
                flexible_align=True,
                align_backend=self.rbfe_fep_setting.align_backend,
            )
            ligand_align.align()

    def pair_network_init(self):
        self.gg = MX.init_pair_network(self.molecules["ligands"],
                              topology=self.config["PairNetworkSetting"]["topology"],
                              user_pair_list=self.config["PairNetworkSetting"]["user_pair_list"],
                              bias_nodes=self.config["PairNetworkSetting"]["bias_nodes"],
                              core=self.config["PairNetworkSetting"]["core"],
                              nbunch=self.config["PairNetworkSetting"]["nbunch"])
        
    def atom_mapping(self):
        MX.atom_mapping(self.gg)
        
    def molecule_similiarity(self):
        MX.molecule_similiarity(self.gg)

    def pair_network_final(self):
        MX.final_pair_nework(self.gg,
                        #topolgy=self.config["PairNetworkSetting"]["topology"],
                       nbunch=self.config["PairNetworkSetting"]["nbunch"],
                       bias_nodes=self.config["PairNetworkSetting"]["bias_nodes"],
                       )
        self._before_integral_ligands = {molecule.mole_name: deepcopy(molecule) for molecule in self.molecules["ligands"]}
                    ####MX.molecule_show(molecule,attrs=None,save_file=True,opath=None,extra=None,show_image=False,TD_flag=False)

    def molecule_topol_integral(self):
        topologys = MX.dual_topology("rfe",self.gg,
                                                    output_directory=self.config["EnvironmentSetting"]["output_directory"],
                                                    parallel=True,
                                                    )
        if not topologys:
            raise ValueError(
                "rbfe: no FEP pairs produced. Ensure --ligands SDF has at least 2 molecules and the pair network "
            )
        self.molecules["ligands"] = [topology[0] for topology in topologys]
        
    def property_setting(self):
        for sm in self.sms:
            fep_type = sm.molecules[0].dual_topology_type
            sm.md_setting["free_energy_auixed"]["lambdas"] = MX.get_fep_lambda(sm.md_setting["free_energy"],
                                                                fep_type=fep_type,
                                                                mixed_lambda=sm.md_setting["free_energy_auixed"]["mixed_lambda"],
                                                                is_relative=sm.md_setting["free_energy_auixed"]["is_relative"])
            #if sm.simulation_type in JOB_TYPE:
                #sm.info_dir = Path(sm.output_dir).parent
            sm.atom_mapping = [
                                self._before_integral_ligands[sm.molecules[0].left_molecule_name],
                                self._before_integral_ligands[sm.molecules[0].right_molecule_name],
                                sm.molecules[0].atom_mapping,
                                ]
            
            m1_nonH = {}
            m2_nonH = {}
            for ii,atom in enumerate([atom for atom in sm.atom_mapping[0].Atoms if atom.elem != "H"]):
                m1_nonH[atom.ID] = ii
            for ii,atom in enumerate([atom for atom in sm.atom_mapping[1].Atoms if atom.elem != "H"]):
                m2_nonH[atom.ID] = ii    
            
            atom_mapping = {}
            atom_mapping_nonH = {}
            
            
            nn1 = len(sm.atom_mapping[0].Atoms)
            nn2 = len(sm.atom_mapping[1].Atoms)
            for an,bn in sm.molecules[0].atom_mapping.items():
                if an < nn1 and bn < nn2:
                    atom_mapping[an] = bn
                    if sm.atom_mapping[0].Atoms[an].elem != "H" and sm.atom_mapping[1].Atoms[bn].elem != "H":
                        atom_mapping_nonH[m1_nonH[an]] = m2_nonH[bn]
            sm.atom_mapping.append(atom_mapping)
            sm.atom_mapping.append(atom_mapping_nonH)

class AFESimulation(MDSimulation):
    def __init__(self, config):
        super().__init__(config)
        
    def molecule_topol_integral(self):
        self.molecules["ligands"] = MX.dual_topology("afe",self.molecules["ligands"],
                                                    output_directory=self.config["EnvironmentSetting"]["output_directory"],
                                                    )
    
    def property_setting(self):
        for sm in self.sms:
            if self.config["MDSetting"]["simulation_type"] in ["abfe","ahfe","hfe"]:
                sm.md_setting["free_energy_auixed"]["absolute_intra_flag"] = sm.molecules[0].absolute_intra_flag
            sm.md_setting["free_energy_auixed"]["lambdas"] = MX.get_fep_lambda(
                                                                                    sm.md_setting["free_energy"],
                                                                                    fep_type=f"afe",
                                                                                    mixed_lambda=sm.md_setting["free_energy_auixed"]["mixed_lambda"],
                                                                                    is_relative=sm.md_setting["free_energy_auixed"]["is_relative"]
                                                                                    )
                    
            if hasattr(sm,"protein_force_field"):
                sm.intermolecular_interaction = MX.get_intermolecule_interaction(sm)
            #if Path(sm.output_dir).name == "abfe":
            #    sm.info_dir = Path(sm.output_dir).parent

class NormalSimulation(MDSimulation):
    def __init__(self, config):
        super().__init__(config)

    def property_setting(self):
        self.sms = MX.build_property_system(self.sms)

