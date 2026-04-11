import copy
import math
import statistics
from copy import deepcopy
import itertools
from typing import Iterable
import statistics

from ..utils import logger

from .atom_type import DecoratedAtomType
from ..chemkit.conformation.conformation import ConformType

#from . import AtomType
from ..utils.commons import parallel_run

from .force_field import ForceField, _topol_to_ff_term, get_hybrid_froce_field
from .charge_parameters import get_nonbinc_charge


class ForceFieldManager:
    def __init__(self):
        pass

    #######爬参数######################################
    @staticmethod
    def assign_single_force_field(molecule,ffjson=None,empi_ffjson=None,this_terms=None,ignore_existing=False,idx=None):
        if not ignore_existing and "force field" in molecule.steps:
            pass
        else:
            if ffjson is not None:
                ForceField.assign_para(molecule,ffjson,empi_ff=empi_ffjson,this_terms=this_terms)
        if idx is not None:
            return molecule,idx
        else:
            return molecule

    @staticmethod
    def assign_force_field_parameter(molecules,ffjson,empi_ffjson=None,this_terms=None,ignore_existing=False,parallel=True):
        ####assign force field######################################
        if this_terms is None:
            #this_terms = {molecule.mole_name:None for molecule in molecules}
            this_terms = [None for _ in molecules]
        
        if parallel:
            new_molecules = parallel_run(ForceFieldManager.assign_single_force_field,molecules,
                                     kwds=[{"this_terms":this_terms[jj],"ffjson":ffjson,"empi_ffjson":empi_ffjson,
                                            "ignore_existing":ignore_existing} for jj,molecule in enumerate(molecules)],keep_order=True)
            
        else:
            new_molecules = []
            for jj,molecule in enumerate(molecules):
                new_molecules.append(ForceFieldManager.assign_single_force_field(molecule,ffjson=ffjson,empi_ffjson=empi_ffjson,this_terms=this_terms[jj],ignore_existing=ignore_existing))
        #for molecule in new_molecules:
        #    molecule.steps.append("force field")
        ###########################################################
        return new_molecules

    @staticmethod
    def get_force_field(molecules,
                           atom_type_file,
                           force_field_file,
                           this_terms=None,
                           reassign_atom_type=False,
                           compensating_at_file=None,
                           compensating_ff_file=None,
                           charge_method=None,
                           charge_ff = None,
                           ignore_existing=False,
                           return_ff=False,
                           empi_ff = None,
                           use_scalevdw=True,
                           parallel=True,):
        flag = False
        this_empi_ff = empi_ff
        if compensating_ff_file is not None and compensating_ff_file != force_field_file:
            flag = True
            this_empi_ff = None
        
        arrs = ForceFieldManager.assign_force_field(molecules,
                                                    atom_type_file,
                                                    force_field_file,
                                                    this_terms=this_terms,
                                                    reassign_atom_type=reassign_atom_type,
                                                    charge_method=charge_method,
                                                    charge_ff = charge_ff,
                                                    ignore_existing=ignore_existing,
                                                    return_ff=return_ff,
                                                    empi_ff = this_empi_ff,
                                                    use_scalevdw=use_scalevdw,
                                                    parallel=parallel,
                                                    )
        
        if return_ff:
            molecules = arrs[0]
            ff_json = arrs[1]
        else:
            molecules = arrs

        if flag:
            compen_molecules = ForceFieldManager.assign_force_field(deepcopy(molecules),
                                                    compensating_at_file,
                                                    compensating_ff_file,
                                                    this_terms=[molecule.loss_para_items for molecule in molecules],####
                                                    reassign_atom_type=True,
                                                    charge_method=charge_method,
                                                    charge_ff = charge_ff,
                                                    ignore_existing=True,
                                                    return_ff=False,
                                                    empi_ff = empi_ff,
                                                    use_scalevdw=use_scalevdw,
                                                    parallel=parallel,
                                                    )
            for ii,molecule in enumerate(molecules):
                for term in ["Atoms","Bonds","Angles","Dihedrals","Impropers","Pair12","Pair13","Pair14","Pair1n"]:
                    if term in molecule.loss_para_items:
                        for jj in molecule.loss_para_items[term]:
                            item = getattr(molecule,term)[jj]
                            compen_item = getattr(compen_molecules[ii],term)[jj]
                            item.parameter = compen_item.parameter
                            item._ff_parameter = compen_item._ff_parameter
                            if item.parameter is not None:
                                item.ptag = "transfer"
                                item._ff_parameter["ptag"] = "transfer"
                            item.pscore = compen_item.pscore
                            item.pcount = compen_item.pcount
                if "binc" in molecule.loss_para_items:
                    for an,vv in molecule.loss_para_items["binc"].items():
                        item = molecule.Atoms[an]
                        compen_item = compen_molecules[ii].Atoms[an]
                        item.point_charge = compen_item.point_charge
                        for nn in vv:
                            item.binc_parameter[nn] = compen_item.binc_parameter[nn]
                            item.binc_tag[nn] = "transfer"
                            item.binc_score[nn] = compen_item.binc_score[nn]
                            item.binc_style[nn] = compen_item.binc_style[nn]
                            item.binc_count[nn] = compen_item.binc_count[nn]
                            item._ff_binc_parameter[nn] = compen_item._ff_binc_parameter[nn]
                            item._ff_binc_parameter[nn]["ptag"] = "transfer"
                if "pair_charge" in molecule.loss_para_items:
                    for term,vv in molecule.loss_para_items["pair_charge"].items():
                        for jj in vv:
                            item = getattr(molecule,term)[jj]
                            compen_item = getattr(compen_molecules[ii],term)[jj]
                            item.charge_parameter = compen_item.charge_parameter
                                
                
        if return_ff:
            return molecules, ff_json
        else:
            return molecules                


    @staticmethod
    def assign_force_field(molecules,
                           atom_type_file,
                           force_field_file,
                           this_terms=None,
                           reassign_atom_type=True,
                           charge_method=None,
                           charge_ff = None,
                           ignore_existing=False,
                           return_ff=False,
                           empi_ff = None,
                           use_scalevdw=True,
                           parallel=True,
                           ):
        """
        爬参数，可能会缺失参数
        """        
        ####assign atom type########################################
        from . import AtomType
        if reassign_atom_type:
            AT = AtomType(atf=atom_type_file)
            molecules = AT._atom_type(molecules,parallel=parallel,this_terms=this_terms,ignore_existing=True)
        ############################################################
        ####读取力场参数#############################################
        
        ffjson = ForceField.read_files(force_field_file,use_scalevdw)
        if empi_ff is not None:
            empi_ffjson = ForceField.read_files(empi_ff,False)
        else:
            empi_ffjson = None
        
        ####生成 非常规charge参数。只有binc，atc 参数在力场文件中
        if charge_ff is not None:
            ffjson["general"]["qmodel"] = "manual"
            ffjson["manual"] = charge_ff 
        
        else:
            if charge_method is not None and charge_method not in ["binc","atc"]:

                ffjson["general"]["qmodel"] = charge_method
                ffjson[charge_method] = get_nonbinc_charge(molecules,charge_method,ignore_existing=ignore_existing)
            
        ############################################################
        for ii,sp in enumerate(ffjson["general"]["special_bond"]):
            if sp in [None,"None"]:
                _pair_ = f"Pair1{ii+2}"
                for molecule in molecules:
                    if hasattr(molecule,_pair_):
                        delattr(molecule,_pair_)
            
        molecules = ForceFieldManager.assign_force_field_parameter(molecules,ffjson,empi_ffjson=empi_ffjson,this_terms=this_terms,ignore_existing=ignore_existing,parallel=parallel)

        if not return_ff:
            return molecules
        else:
            if use_scalevdw:
                ffjson["general"]["use_scalevdw"] = True
            return molecules, ffjson 
    ################################################################################################


    @staticmethod
    def old_assign_force_field(molecules,
                           atom_type_file,
                           force_field_file,
                           this_terms=None,
                           reassign_atom_type=True,
                           charge_method=None,
                           charge_ff = None,
                           ignore_existing=False,
                           return_ff=False,
                           empi_ff = None,
                           use_scalevdw=True,
                           parallel=True,
                           ):
        """
        爬参数，可能会缺失参数
        """        
        ####assign atom type########################################
        from . import AtomType
        if reassign_atom_type:
            AT = AtomType(atf=atom_type_file)
            molecules = AT._atom_type(molecules,parallel=parallel,ignore_existing=True)
        ############################################################
        ####读取力场参数#############################################
        
        ffjson = ForceField.read_files(force_field_file,use_scalevdw)
        if empi_ff is not None:
            empi_ffjson = ForceField.read_files(empi_ff,False)
        else:
            empi_ffjson = None
        
        ####生成 非常规charge参数。只有binc，atc 参数在力场文件中
        if charge_ff is not None:
            ffjson["general"]["qmodel"] = "manual"
            ffjson["manual"] = charge_ff 
        
        else:
            if charge_method is not None and charge_method not in ["binc","atc"]:

                ffjson["general"]["qmodel"] = charge_method
                ffjson[charge_method] = get_nonbinc_charge(molecules,charge_method,ignore_existing=ignore_existing)
            
        ############################################################
        for ii,sp in enumerate(ffjson["general"]["special_bond"]):
            if sp in [None,"None"]:
                _pair_ = f"Pair1{ii+2}"
                for molecule in molecules:
                    if hasattr(molecule,_pair_):
                        delattr(molecule,_pair_)
            
        molecules = ForceFieldManager.assign_force_field_parameter(molecules,ffjson,empi_ffjson=empi_ffjson,this_terms=this_terms,ignore_existing=ignore_existing,parallel=parallel)

        if not return_ff:
            return molecules
        else:
            if use_scalevdw:
                ffjson["general"]["use_scalevdw"] = True
            return molecules, ffjson 
    ################################################################################################

    #########从分子中提取参数##########
    @staticmethod
    def checkout_force_field(molecules,ffjson,isfittig_flag=False):
        loss_parameters = []
        empi_parameters = []
        this_ff = {"general":ffjson["general"]}

        # TODO
        #####isfitting_flag 是做什么的？？？？？？？
        if isfittig_flag:
            name_attrs = "atom_type_used_names"
            name_attr = "atom_type_used_name"
        else:
            #name_attrs = "atom_type_names"
            #name_attr = "atom_type_name"
            name_attrs = "atom_type_used_names"
            name_attr = "atom_type_used_name"
        #############
        
        for molecule in molecules:
            #if molecule.style == "protein":
            #    continue
            for term,item in _topol_to_ff_term.items():
                tts = getattr(molecule,term,[])
                if len(tts) > 0 and item not in this_ff:
                        this_ff[item] = {}
                for tt in tts:
                    _name = list(set(getattr(tt,name_attrs)).intersection(this_ff[item].keys()))
                    if len(_name) > 0:
                        if getattr(tt,name_attr) != _name[0]:
                            setattr(tt,name_attr,_name[0])
                            tt._ff_parameter["name"] = _name[0]
                    else:
                        if item not in ["pairwise","pairwise12","pairwise13","pairwise14"]:
                            this_ff[item][getattr(tt,name_attr)] = tt._ff_parameter
                            tt._ff_parameter["name"] = getattr(tt,name_attr)
                            if tt._ff_parameter["ptag"] == "null":
                                loss_parameters.append(f"{item}:{tt.atom_type_used_name}")
                            elif tt._ff_parameter["ptag"] == "empi":
                                empi_parameters.append(f"{item}:{tt.atom_type_used_name}")
                if term == "Atoms":
                    for atom in molecule.Atoms:
                        if hasattr(atom,"atom_type_name_m2"):
                            if atom.atom_type_name_m2 not in this_ff[item].keys():
                                this_ff[item][atom.atom_type_name_m2] ={
                                    "name": atom.atom_type_name_m2,
                                    "pstyle":atom.pstyle,
                                    "fix_parameter":[],
                                    "parameter":atom.parameter_m2,
                                    "fit_parameter":[],
                                    "mass":atom.mass_m2,
                                    "ptag":"fit",
                                    "pscore":"nan",
                                    "pcount":"nan"
                                }

                if term == "binc":
                    for atom in molecule.Atoms:
                        binc_tts = getattr(atom,"_ff_binc_parameter",[])
                        if len(binc_tts) > 0 and item not in this_ff:
                            this_ff[item] = {}
                        for tt in binc_tts:
                            if tt["name"] not in this_ff[item].keys():
                                this_ff[item][tt["name"]] = tt
                                if tt["ptag"] == "null":
                                    loss_parameters.append(f"{item}:{tt['name']}")
                                if tt["ptag"] == "empi":
                                    empi_parameters.append(f"{item}:{tt['name']}")
        
        all_ats = list(this_ff["atomtype"].keys())
        if "pairwise" in ffjson:
            if "pairwise" not in this_ff:
                this_ff["pairwise"] = {}
            for pairterm in itertools.permutations(all_ats,2):
                if pairterm in ffjson["pairwise"]:
                    this_ff["pairwise"][pairterm] = deepcopy(ffjson["pairwise"][pairterm])
        
        if this_ff["general"]["use_scalevdw"] and "scalevdw" in ffjson:
            parent_term = list(set(all_ats).intersection(set(ffjson["scalevdw"])))
            if len(parent_term) > 0:
                if "pairwise" not in this_ff:
                    this_ff["pairwise"] = {}
                for parent_at in parent_term:
                    
                    para_parent_at = this_ff["atomtype"][parent_at]
                    for at in all_ats:
                        if at in ffjson["scalevdw"][parent_at]:
                            if f"{parent_at}${at}" not in this_ff["pairwise"] and f"{at}${parent_at}" not in this_ff["pairwise"]:
                                scale_factor = ffjson["scalevdw"][parent_at][at]["parameter"][0]
                                para_at = this_ff["atomtype"][at]
                                sigma_min = 0.5 * (para_parent_at["parameter"][0] + para_at["parameter"][0])
                                epsilon_min = math.pow(para_parent_at["parameter"][1] * para_at["parameter"][1], 0.5)
                                sigma = sigma_min / math.pow(scale_factor,1.0 / 6.0)
                                epsilon = epsilon_min * math.pow(scale_factor,2.0)
                                parameter = deepcopy(para_at)
                                parameter["name"] = f"{parent_at}${at}"
                                parameter["parameter"] = [sigma,epsilon]
                                this_ff["pairwise"][parameter["name"]] = parameter
                                #sigma_1 = o_2w_para[0]
                                #epsilon_1 = o_2w_para[1]
                                #sigma_2 = aa.parameter[0]
                                #epsilon_2 = aa.parameter[1]
                                #sigma_mix = 0.5 * (sigma_1 + sigma_2)
                                #epsilon_mix = math.pow(epsilon_1 * epsilon_2, 0.5)
                                #sigma = sigma_mix / math.pow(aa.parameter[2], 1.0 / 6.0)
                                #epsilon = epsilon_mix * math.pow(aa.parameter[2], 2.0)


        if len(loss_parameters) > 0:
            logger.error("Missing parameters which will bring some heavy error in other jobs: %s" %loss_parameters)
            this_ff["_loss_parameters"] = loss_parameters
        if len(empi_parameters) > 0:
            logger.warning("Please check parameters, because some parameter from empirical force field: %s" %empi_parameters)
            this_ff["_empi_parameters"] = empi_parameters

        return this_ff

    ####拟定要拟合的参数##############################################################################
    ###目前仅针对bond,angle是harmonic形式，dihedral是fourier形式，improper是fourier_2n形式
    ###其他函数形式会有很大的不同，勿用
    @staticmethod
    def get_haved_scan_para(molecules):
        """
        确定有scan数据可以进行拟合的dihedral参数。
        如果没有相应的scan数据可以拟合，则该dihedral参数仅仅用hessian数据拟合，结果不是很可靠
        """
        haved_scan_para = []
        _s_at = ["h_1","h_1=","h_1=2"]
        for m in molecules:
            if hasattr(m, "constrain"):
                for rr in m.constrain:
                    tmp = []
                    for dihe in m.Dihedrals:
                        if (dihe.a2 == rr.a2 and dihe.a3 == rr.a3) or (dihe.a3 == rr.a2 and dihe.a2 == rr.a3):
                            if dihe.a1_atom_type_used not in _s_at and dihe.a4_atom_type_used not in _s_at:
                                tmp.append(dihe.atom_type_used_names)
                    for rr in tmp:
                        haved_scan_para.append(rr[0])
                        haved_scan_para.append(rr[1])

        return list(dict.fromkeys(haved_scan_para).keys())

    @staticmethod
    def get_guess_bond_angle_equ_value(molecules):
        opt_values = {}
        for molecule in molecules:
            for bond in getattr(molecule, "Bonds", []):
                if bond.atom_type_used_name not in opt_values:
                    opt_values[bond.atom_type_used_name] = []
                opt_values[bond.atom_type_used_name].append(bond.value)
            for angle in getattr(molecule, "Angles", []):
                name = angle.atom_type_used_name
                if name not in opt_values:
                    opt_values[name] = []
                opt_values[name].append(angle.value)
        opt_values = {k: statistics.mean(v) for k, v in opt_values.items()}
        return opt_values

    @staticmethod
    def _get_used_infos(molecules,flag_fitting_scan_parameter,preprocessing_fitting_parameter):
        if flag_fitting_scan_parameter and molecules is not None: 
            scan_molecules = [
                              molecule 
                              for molecule in molecules 
                              if molecule.conform_type in ConformType.TORSION_SCAN_TYPES
                             ]
        else:
            scan_molecules = None

        if scan_molecules is not None:
            # 提取有scan数据的二面角参数
            haved_scan_para = ForceFieldManager.get_haved_scan_para(scan_molecules)
        else:
            haved_scan_para = None

        if molecules is not None and preprocessing_fitting_parameter:
            optimized_molecules=[
                                  molecule 
                                  for molecule in molecules 
                                  if molecule.conform_type in ConformType.TORSION_SCAN_TYPES + [ConformType.LOCAL_MINIMUM]
                                 ]
        else:
            optimized_molecules = None

        if preprocessing_fitting_parameter and optimized_molecules is not None:
            opt_values = ForceFieldManager.get_guess_bond_angle_equ_value(optimized_molecules)
        else:
            opt_values = {}
        return haved_scan_para,opt_values

    @staticmethod
    def normal_fitting_parameter(this_ff,terms,fix_tag):
        for term, items in this_ff.items():
            if term not in terms:
                continue
            for name,item in items.items():
                if item["ptag"] in fix_tag:
                    continue

                ss = name.split("$")

                # validate binc parameters between same atoms
                if term == "binc" and ss[0] == ss[1]:
                    item["ptag"] = "Fit"
                    continue

                item["ptag"] = "isfitting"
                ####dihedralterm的偏移角固定，不拟合，即1，3是0; 2，4是180
                if term == "dihedralterm":
                    item["isfitting"] = [0,2,4,6]
                    continue

                item["isfitting"] = [ii for ii in range(len(item["parameter"])) if ii not in item["fix_parameter"]]
        return this_ff

    @staticmethod
    def preprocessing_fitting_parameter(this_ff,terms,opt_values,atom_type_file,haved_scan_para):
        for term,items in this_ff.items():
            if term not in terms:
                continue
            if term in ["bondterm","angleterm","improperterm"]:
                this_ff[term] = ForceFieldManager.preprocessing_non_dihedral_fitting_parameter(term,items,opt_values)
            elif term in ["dihedralterm"]:
                this_ff[term] = ForceFieldManager.preprocessing_dihedral_fitting_parameter(items,atom_type_file,haved_scan_para)
        return this_ff

    @staticmethod
    def preprocessing_dihedral_fitting_parameter(items,atom_type_file,haved_scan_para,):
        from . import AtomType
        if atom_type_file is None:
            typer = AtomType()
        else:
            typer = AtomType(atf=atom_type_file)
        AtomTypeDecoChar = typer.at.at_types["root"]["AtomTypeDecoChar"]
        at_types = typer.at.at_types

        for name,item in items.items():
            if item["ptag"] == "isfitting":
                ss = name.split("$")
                if haved_scan_para and item["name"] not in haved_scan_para:
                    item["isfitting"] = []
                    continue
                # zero torsion if central atoms are degenerated
                #if AtomTypeDecoChar.DEG in (ss[1], ss[2]):
                if AtomTypeDecoChar["DEG"] in (ss[1], ss[2]):
                    item["parameter"] = [0.0, 0.0, 0.0, 180.0, 0.0, 0.0, 0.0, 180.0]
                    item["isfitting"] = []
                    continue
                # we don't like too large torsion parameters and k4 is zero normally
                for ii in [0,2,4]:
                    item["parameter"][ii] = max(-6.0, min(item["parameter"][ii], 6.0))
                #item["para"][:3] = [max(-6.0, min(v, 6.0)) for v in item["para"][:3]]
                item["parameter"][6] = 0
                try:
                    item["isfitting"].remove(6)
                except ValueError:
                    pass

                deco2 = DecoratedAtomType(ss[1])
                deco3 = DecoratedAtomType(ss[2])
                # identify special torsion
                _is_degenerated = ss[0] == AtomTypeDecoChar["DEG"]
                _is_h = "h_1" in (ss[0], ss[-1])
                # currently all ring torsion are degenerated
                # but it's not guaranteed in the future
                _is_endo = deco2.tag_endo and deco3.tag_endo
                # pre-defined planer torsion
                if typer.at and at_types["root"]["planer_dihedral"]:
                    deco_types = [i for i in at_types["root"]["planer_dihedral"] if isinstance(i, DecoratedAtomType)]
                    deco_pairs = [i for i in at_types["root"]["planer_dihedral"] if isinstance(i, Iterable)]
                    try:
                        if any(deco2.is_subset_of(deco) for deco in deco_types) and any(
                            deco3.is_subset_of(deco) for deco in deco_types
                        ):
                            raise StopIteration
                        if any(DecoratedAtomType.is_pair_subset_of((deco2, deco3), pair) for pair in deco_pairs):
                            raise StopIteration
                    except StopIteration:
                        item["parameter"] = [0.0, 0.0, 8.9, 180.0, 0.0, 0.0, 0.0, 180.0]
                        item["isfitting"] = []
                        continue
                # hybridization
                ep2 = typer.at.get_electron_pairs(ss[1])
                ep3 = typer.at.get_electron_pairs(ss[2])
                # sp2-sp2
                if {ep2, ep3} == {3}:
                    if _is_endo and _is_degenerated:
                        item["parameter"] = [0.0, 0.0, 1.8, 180.0, 0.0, 0.0, 0.0, 180.0]
                        item["isfitting"] = []
                    elif _is_h or _is_endo or _is_degenerated:
                        item["parameter"] = [0.0, 0.0, 0.0, 180.0, 0.0, 0.0, 0.0, 180.0]
                        item["isfitting"] = [2]
                    else:
                        item["parameter"] = [0.0, 0.0, 0.0, 180.0, 0.0, 0.0, 0.0, 180.0]
                        item["isfitting"] = [0, 2]
                    continue
                # sp3-sp3
                if {ep2, ep3} == {4}:
                    if _is_h or _is_endo or _is_degenerated:
                        item["parameter"] = [0.0, 0.0, 0.0, 180.0, 0.0, 0.0, 0.0, 180.0]
                        item["isfitting"] = [4]
                    continue
                # sp2-sp3
                if {ep2, ep3} == {3, 4}:
                    if _is_endo and _is_degenerated:
                        item["parameter"] = [0.0, 0.0, 0.0, 180.0, 0.0, 0.0, 0.0, 180.0]
                        item["isfitting"] = []
                    continue
        return items

    @staticmethod
    def preprocessing_non_dihedral_fitting_parameter(typ,items,opt_values):
        if typ == "bondterm":
            for name,item in items.items():
                if item["ptag"] == "isfitting":
                    ss = name.split("$")
                    try:
                        item["isfitting"].remove(0)
                    except ValueError:
                        pass
                    else:
                        if name in opt_values:
                            item["parameter"][0] = opt_values[name]
                # using large force for PO4 bond
                    if ss[0][:4] == "p_4=" or ss[1][:4] == "p_4=":
                        item["parameter"][1] = 300.0
        # freeze angle value
        if typ == "angleterm":
            for name,item in items.items():
                if item["ptag"] == "isfitting":
                    ss = name.split("$")
                    try:
                        item["isfitting"].remove(0)
                    except ValueError:
                        pass
                    else:
                        if name in opt_values:
                            item["parameter"][0] = opt_values[name]
                    # using large force for OH, NH, PO4 angle
                    if "h_1o" in (ss[0], ss[2]):
                        item["parameter"][1] = 100.0
                    elif "h_1n" in (ss[0], ss[2]):
                        item["parameter"][1] = 100.0
                    elif ss[1][:4] == "p_4=":
                        item["parameter"][1] = 125.0
        # freeze improper parameters
        if typ == "improperterm":
            for name,item in items.items():
                if item["ptag"] == "isfitting":
                    ss = name.split("$")
                    if ss[0][:3] == "c_3" or ss[0][3:4] == "a":
                        item["parameter"][0] = 6.0
                    item["isfitting"] = []
        return items

    @staticmethod
    def get_fitting_parameters(
        this_ff,
        fix_tag=["V", "Fit","amber99sb"],
        terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
        preprocessing_fitting_parameter=True,
        molecules=None,
        flag_fitting_scan_parameter=False,
        atom_type_file=None,
    ):
        """
        确定需要拟合的参数，并修改某些参数的初始值
        this_ff will be modified. `isfitting` attribute will be added for fitting terms
        输入：
            this_ff: Dict, 字典形式的力场参数
            fitting_scan_para_flag: True or False, 是否只拟合有scan数据的二面角参数
            moles： List[Molecule],如果fitting_scan_para_flag=True，则需要
            fix_flag: List[str], 参数如果是这些类型的tag标记，则不会被拟合，会被继承
            terms: List[str], 需要拟合的拓扑结构项
            remove_ring_para: True or False, 是否改变环的参数
            modify_dihe_para: True or False, 是否修改二面角参数，以使参数能在更大数值范围内变动
            fix_hydrogen_dihe: True or False, 是否固定含氢的二面角
            fitting_fourier_n: None or List[int]。确定那些傅立叶展开形式的二面角参数需要拟合。如果为None则所有参数都拟合
        输出：
            this_ff: Dict, 字典形式的力场参数，还有isfitting标记
        """

        haved_scan_para,opt_values = ForceFieldManager._get_used_infos(molecules,
                                                                       flag_fitting_scan_parameter,
                                                                       preprocessing_fitting_parameter)
        this_ff = ForceFieldManager.normal_fitting_parameter(this_ff,terms,fix_tag)
        if preprocessing_fitting_parameter:
            this_ff = ForceFieldManager.preprocessing_fitting_parameter(this_ff,terms,opt_values,atom_type_file,haved_scan_para)
        return this_ff
    ###########################################################################

    @staticmethod
    def get_special_infos(molecules,this_ff):
        torsion_infos = {}
        _torsion_para_ = []
        _fitting_para_ = [para for para,vv in this_ff["dihedralterm"].items() if vv["ptag"] == "isfitting"]
        exist_torsion = []
        for molecule in molecules:
            if hasattr(molecule,"constrain"):
                if molecule.mole_name not in torsion_infos:
                    torsion_infos[molecule.mole_name] = {}
                for cons in molecule.constrain:
                    atoms = cons.atoms
                    if len(atoms) == 4:
                        if f"{molecule.mole_name}_{cons.name}" in exist_torsion:
                            continue
                        exist_torsion.append(f"{molecule.mole_name}_{cons.name}")
                        tmp = [[],[],[]]
                        for dihe in molecule.Dihedrals:
                            if (dihe.a2 == atoms[1] and dihe.a3 == atoms[2]) or (dihe.a2 == atoms[2] and dihe.a3 == atoms[1]):
                                tmp[0].append(dihe.atoms)
                                tmp[1].append(dihe.atom_type_used_name)
                                if dihe.atom_type_used_name in _fitting_para_:
                                    tmp[2].append(1)
                                    _torsion_para_.append(dihe.atom_type_used_name)
                                else:
                                    tmp[2].append(0)
                        
                        torsion_infos[molecule.mole_name][cons.name] = tmp

        this_ff["torsion_infos"] = torsion_infos
        this_ff["torsion_para"] = _torsion_para_
                

    @staticmethod
    def reset_term_to_loss(
        this_ff,
        term=["atomtype", "pair14term", "pair1nterm", "pair12term", "pair13term"],
    ):
        """
        把某些拓扑结构设定为未找到力场参数的状态
        """
        tmp1_ff = {}
        tmp2_ff = {}
        for aa, bb in this_ff.items():
            if aa in ["qmodel", "special_bond"]:
                tmp2_ff[aa] = bb
            else:
                if aa in term:
                    if aa not in tmp1_ff.keys():
                        tmp1_ff[aa] = []
                    for aaa, bbb in bb.items():
                        tmp1_ff[aa].append(aaa)
                else:
                    if aa not in tmp2_ff.keys():
                        tmp2_ff[aa] = {}
                    for aaa, bbb in bb.items():
                        tmp2_ff[aa][aaa] = bbb
        return tmp2_ff, tmp1_ff
    ###################################################################################################

    ####写力场参数文件##################################################################################
    @staticmethod
    def write_term(term_name, jdict, outf):
        """
        被self.transfer所调用
        """
        for k, v in jdict[term_name].items():
            terms = k.split("$")
            temp_str = f"{term_name} "
            for t in terms:
                temp_str += "%8s " % t
            temp_str += "%s " % v["pstyle"]
            if term_name == "atomtype":
                #if "mass" in v:
                #    mass = v["mass"]
                #else:
                #    elem = k.capitalize().split("_")[0]
                #    mass = Element.get(elem).mass
                temp_str += "%8.4f " % v["mass"]

            for i in range(len(v["parameter"])):
                p = v["parameter"][i]
                #if term_name == "angleterm" and i == 0:
                #if term_name in ["angleterm","dihedralterm"]:
                #    if _degree_to_radian(term_name,v["style"],i,0):
                #        p = p * 180 / math.pi

                if v["ptag"] == "isfitting" and "isfitting" in v.keys() and i not in v["isfitting"]:
                    temp_str += "%8.4f* " % p
                else:
                    temp_str += "%8.4f  " % p
            temp_str += "%4s " % v["ptag"]
            temp_str += f"{v['pscore']} nan\n"

            outf.write(temp_str)

        outf.write("\n")

    @staticmethod
    def transfer(jdict, output_log_path):
        """
        把字典形式的力场参数，写成一个文本格式的力场参数文件
        输入：
            jdict: Dict，以字典形式记录的力场参数文件
            output_log_path: 要生成的力场参数的文件名
        """
        # with open(input_log_path, 'r') as inf:
        #     jdict = json.load(inf)

        with open(output_log_path, "w") as outf:
            keys = jdict.keys()
            # write combination_rule
            if "pairterm" in keys:
                for k, v in jdict["pairterm"].items():
                    outf.write("combination_rule %s\n" % v["combination_rule"])
                    break
            elif "pair14term" in keys:
                for k, v in jdict["pair14term"].items():
                    outf.write("combination_rule %s\n" % v["combination_rule"])
                    break
            elif "pair1nterm" in keys:
                for k, v in jdict["pair1nterm"].items():
                    outf.write("combination_rule %s\n" % v["combination_rule"])
                    break
            else:
                outf.write("combination_rule LB\n")


            # write special_bond
            if "special_bond" in jdict["general"]:
                outf.write("special_bond")
                ###for lst in jdict["general"]["special_bond"]:
                outf.write(" " + " ".join(map(str, jdict["general"]["special_bond"])))
                outf.write("\n")

            # write qmodel
            if "qmodel" in jdict["general"]:
                qmodel = jdict["general"]["qmodel"]
                if qmodel not in ("binc", "atc"):
                    qmodel = "None"
                outf.write("qmodel %s\n" % qmodel)

            # write equ_table
            if "equ_table" in jdict["general"]:
                outf.write("equ_table %s\n" % jdict["general"]["equ_table"])

            outf.write("\n")

            terms_to_write = [
                "atomtype",
                "bondterm",
                "angleterm",
                "dihedralterm",
                "improperterm",
                "binc",
            ]

            for term in terms_to_write:
                if term in jdict.keys():
                    ForceFieldManager.write_term(term, jdict, outf)

    @staticmethod
    def combine_ff(ff_ref, ff_new, prior_tags=None):
        """
        合并两个字典形式的力场参数
        In case a term exists in both ff_ref and ff_new, the term with tag in prior_tags will be kept.
        If prior_tags is None, the term from ff_new will always overwrite the one from ff_ref.
        """
        ff = copy.deepcopy(ff_ref)
        for typ, terms in ff_new.items():
            if typ not in ["torsion_infos","torsion_para","general"]:
                if type(terms) is not dict:
                    continue
                if typ not in ff.keys():
                    ff[typ] = {}
                for name, term in terms.items():
                    if name in ff[typ] and prior_tags is not None and term["ptag"] not in prior_tags:
                        continue
                    ff[typ][name] = term

        return ff

    @staticmethod
    def combine_ff_file(ff_file_ref, ff_file_new,used_vdw,out_path):
        """
        将两个文本形式的力场文件，合并成一个文件
        第二个力场中标记为 V, Fit, isfitting, pseudofitting 的条目会覆盖第一个力场中相应的条目
        同时第二个力场中标记为 isfitting 的条目会被标识为 Fit
        """
        ffobj = ForceField()
        ff_ref = ffobj.read_files(ff_file_ref,used_vdw)
        ff_new = ffobj.read_files(ff_file_new,used_vdw)

        for typ, terms in ff_new.items():
            if typ in ["general","torsion_infos","torsion_para"]:
                continue
            for kk, vv in terms.items():
                if vv["ptag"] == "isfitting":
                    vv["ptag"] = "Fit"
        
        ff = ForceFieldManager.combine_ff(ff_ref, ff_new, prior_tags=["V", "Fit", "pseudofitting"])
        ForceFieldManager.transfer(ff, out_path)

    ###########some tools##########
    @staticmethod
    def _convert_force_field_format_to_json(fff):
        return ForceField.read_files(fff,use_scalevdw=True)

    @staticmethod
    def _get_empirical_force_field(atf,fff1,fff2):
        get_hybrid_froce_field(atf,fff1,fff2)

    @staticmethod
    def read_force_field_file(fff,use_scalevdw):

        return ForceField.read_files(fff,use_scalevdw)
    ############################################################################################################
    
    @staticmethod
    def grasp_protein_amber_ff(fff,protein):
        pass