from .ring_breaking.match import get_num_soft_bonds
from .dual_topology import dual_topolgy_assign,charge_abfe, FEPTopology
from pathlib import Path
from ...chem import FormatMolecule as FM

def num_soft_bonds(wt, mut, wt_core, mut_core):
    return get_num_soft_bonds(wt,mut,wt_core,mut_core)

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