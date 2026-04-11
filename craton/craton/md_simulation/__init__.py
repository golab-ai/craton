import os

from .md import write_system,write_bash,write_job_info
from .mapping.pairnetwork import PairNetwork
from .fep.lambda_schedule import LambdaSchedule
from .fep.abfe import get_intermolecular_restrain
from .dual_topology.ring_breaking.match import get_num_soft_bonds
from .dual_topology.dual_topology import dual_topolgy_assign,charge_abfe, FEPTopology
from pathlib import Path
from ..chem import FormatMolecule as FM

from .md_analyze.gmx_md_check import CheckGmxRun
from .md_analyze.gmx_analyze import GmxAnalyze

class MDSimulation:
    def __init__(self) -> None:
        pass

    @staticmethod
    def write_input_file(systems,parallel=True):
        write_system(systems,parallel=parallel)

    @staticmethod
    def write_bash_file(systems,parallel=True):
        write_bash(systems,parallel=parallel)

    @staticmethod
    def write_info(systems,parallel=True):
        write_job_info(systems,parallel=parallel)
        
class FEPTool:
    def __init__(self) -> None:
        pass
    
    @staticmethod
    def pair_network_init(ligands,
                          topology="normal",
                          user_pair_list=None,
                          bias_nodes=None,
                          core=None,
                          nbunch=None):
        return PairNetwork.create_graph_from_molecules(ligands,
                                                   topology=topology,
                                                   user_pair_list=user_pair_list,
                                                   bias_nodes=bias_nodes,
                                                   core=core,
                                                   nbunch=nbunch)
    #######mapping
    @staticmethod
    def atom_mapping_calculate(gg):
        PairNetwork.calculate_atom_mapping(gg)

    @staticmethod
    def molecule_similiarity_calculate(gg):
        PairNetwork.calculate_similarity(gg)

    @staticmethod
    def pair_network_final(gg,topology="normal",nbunch=None,bias_nodes=None):
        if  topology== "normal":
            PairNetwork.reduce_normal_graph(gg,nbunch=nbunch,bias_nodes=bias_nodes)

    @staticmethod
    def graph_attributes_report(gg,topology="normal"):
        PairNetwork.report_graph_attributes(gg,topology=topology)
    
    ##########fep    
    @staticmethod
    def get_lambda_schedule(fep_setting,fep_type="r_group",mixed_lambda=False,is_relative=False):
        if fep_type == "r_group":
            return LambdaSchedule(fep_setting=fep_setting,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()
        elif fep_type == "charge_hopping":
            return LambdaSchedule(fep_setting=fep_setting, is_charge_hopping=True,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()
        elif fep_type == "core_hopping":
            return LambdaSchedule(fep_setting=fep_setting, is_core_hopping=True,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()
        elif fep_type == "couple":
            return LambdaSchedule.hfe_lambda()
        else:
            return LambdaSchedule(fep_setting=fep_setting,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()

    @staticmethod
    def abfe_intermolecule(system):
        return get_intermolecular_restrain(system)
    
    ########dual_topology
    @staticmethod
    def num_soft_bonds(wt, mut, wt_core, mut_core):
        return get_num_soft_bonds(wt,mut,wt_core,mut_core)
    
    @staticmethod
    def assign_dual_topology(fep_type, gg,output_directory=".",parallel=True):
        if fep_type == "afe":
            molecules = []
            for m in gg:
                molecules.append(charge_abfe(m))
            return molecules
        elif fep_type in ["mutation","pep-rbfe","rna-rbfe"]:
            topologys = []
            for g in gg:
                FEP_topoly = FEPTopology(*g)
                topology = FEP_topoly.dual_topology()
                Path(output_directory + f"/{topology[0].mole_name}").mkdir(parents=True, exist_ok=True)
                FM._convert(topology[0],otype="mol",ofilename="left",opath=f"{output_directory}/{topology[0].mole_name}")
                FM._convert(topology[1],otype="mol",ofilename="right",opath=f"{output_directory}/{topology[0].mole_name}")
                topologys.append(topology)
            
            return topologys
        else:
            topologys = dual_topolgy_assign(gg,parallel=parallel)
            for topology in topologys:
                Path(output_directory + f"/{topology[0].mole_name}").mkdir(parents=True, exist_ok=True)
                FM._convert(topology[0],otype="mol",ofilename="left",opath=f"{output_directory}/{topology[0].mole_name}")
                FM._convert(topology[1],otype="mol",ofilename="right",opath=f"{output_directory}/{topology[0].mole_name}")

            return topologys
        
class MDAnalyze:
    def __init__(self) -> None:
        pass

    @staticmethod
    def check_md(directory,batchfile="batch_0.txt",mdengine="gmx",parallel=True):
        if mdengine in ["gmx","gromacs"]:
            checker = CheckGmxRun(directory,batchfile=batchfile,parallel=parallel)
            results = checker._get_info()
            for kk in ["done_jobs","wait_jobs","run_jobs","error_jobs"]:
                vv = results[kk]
                if len(vv) > 0:
                    print("%s:" %kk)
                    print("%s" %"\n".join([rr[0] for rr in vv]))
                    print("####################\n\n")
            return results

    @staticmethod
    def gmx_analyzer(ptype,args):
        ga = GmxAnalyze()
        ga.run(ptype,args)