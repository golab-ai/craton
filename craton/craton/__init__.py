import os,sys,time
from pathlib import Path
import inspect

import matplotlib
import psutil
import yaml
from copy import deepcopy
from IPython.display import display,Image,SVG

matplotlib.use("agg")  # disable X server

CRATON_DIR = Path(__file__).parent
ROOT_DIR = CRATON_DIR.parent
f = open(f"{ROOT_DIR}/configure/configure.yaml")
CRATON_CONFIGURE = yaml.safe_load(f.read())

DATABASE_EXISTS = Path(f"{ROOT_DIR}/database").exists()
if DATABASE_EXISTS:
    dbf = open(f"{ROOT_DIR}/database/configure.yaml")
    DB_CONFIGURE = yaml.safe_load(dbf.read())
    CRATON_CONFIGURE.update(DB_CONFIGURE)

for key,ff in CRATON_CONFIGURE["ForceFieldSetting"].items():
    if key.find("FILE") != -1:
        CRATON_CONFIGURE["ForceFieldSetting"][key] = f"{ROOT_DIR}/{ff}"

for key,pp in CRATON_CONFIGURE["path"].items():
    CRATON_CONFIGURE["path"][key] = f"{ROOT_DIR}/{pp}"

if not Path(CRATON_CONFIGURE["path"]["tmp"]).exists():
    Path(CRATON_CONFIGURE["path"]["tmp"]).mkdir(exist_ok=True)

CRATON_CONFIGURE["EnvironmentSetting"]["NUM_PROCS"] = psutil.cpu_count(logical=False)

####run in jupyter notebook or not
if 'ipykernel' in sys.modules:
    NOTEBOOK_FLAG = True
else:
    NOTEBOOK_FLAG = False

def check_file_path(config):
    for key,ff in config["ForceFieldSetting"].items():
        if key.find("FILE") != -1:
            if ff[0] == "." or ff[0] == "/":
                pass
            else:
                config["ForceFieldSetting"][key] = f"{ROOT_DIR}/{ff}"

def update_configure(section,user_config):
    for kk,vv in section.items():
        if isinstance(vv,dict):
            if kk in ["selector","document"]:
                if kk in user_config:
                    section[kk] = user_config[kk]
            else:
                update_configure(vv,user_config)
        else:
            if kk in user_config:
                if user_config[kk] is not None:
                    section[kk] = user_config[kk]
    return section

def combine_dicts(dict1,dict2):
    for kk,vv in dict2.items():
        if isinstance(vv,dict):
            if kk not in dict1:
                dict1[kk] = {}
            combine_dicts(dict1[kk],vv)
        else:
            dict1[kk] = vv
    return dict1

def update_md_configure(CRATON_CONFIGURE,user_config):
    
    _tmp = user_config["simulation_type"] if "simulation_type" in user_config else None
    md_type = _tmp if _tmp is not None else CRATON_CONFIGURE["MDSetting"]["simulation_type"]
    if md_type in ["rbfe","abfe","ahfe","rhfe","hfe","mutation","rlogs","rlogp","alogp","mem-rbfe","cov-rbfe","pep-rbfe"]:
        _md_type_ = md_type
        if md_type in ["ahfe","rhfe"]:
            _md_type_ = "hfe"
        if md_type in ["mutation","rlogs","rlogp","mem-rbfe","cov-rbfe","pep-rbfe"]:
            _md_type_ = "rbfe"
        if md_type in ["alogp"]:
            _md_type_ = "abfe"
        md_type_section = f"{_md_type_.upper()}FEPMDSetting"
        CRATON_CONFIGURE["MDSetting"] = combine_dicts(CRATON_CONFIGURE["MDSetting"],CRATON_CONFIGURE["FEPMDSetting"])
        CRATON_CONFIGURE["MDSetting"] = combine_dicts(CRATON_CONFIGURE["MDSetting"],CRATON_CONFIGURE[md_type_section])
    if md_type in ["pull"]:
        CRATON_CONFIGURE["MDSetting"] = combine_dicts(CRATON_CONFIGURE["MDSetting"],CRATON_CONFIGURE["PullMDSetting"])
    CRATON_CONFIGURE["MDSetting"] = update_configure(CRATON_CONFIGURE["MDSetting"],user_config)
    return CRATON_CONFIGURE["MDSetting"]

def check_md_setting(config):
    simulation_time = None if "simulation_time" not in config["MDSetting"] else config["MDSetting"]["simulation_time"]
    #logger.info(f"the run simulation time is {simulation_time} ns")
    if simulation_time:
        timestep = config["MDSetting"]["md"]["timestep"]
        config["MDSetting"]["md"]["nsteps"][-1] = int( float(simulation_time) * 1e3 / timestep[-1])
    if config["MDSetting"]["debug"]:
        config["MDSetting"]["md"]["nsteps"] = [1000 for _ in config["MDSetting"]["md"]["jobs"]] 
    return config

def total_update_configure(CRATON_CONFIGURE,usr_config):
    for kk,vv in CRATON_CONFIGURE.items():
        if kk.find("MDSetting") == -1:
            CRATON_CONFIGURE[kk] = update_configure(vv,usr_config)
    CRATON_CONFIGURE["MDSetting"] = update_md_configure(deepcopy(CRATON_CONFIGURE),usr_config)
    CRATON_CONFIGURE = check_md_setting(CRATON_CONFIGURE)
    check_file_path(CRATON_CONFIGURE)
    return CRATON_CONFIGURE

def expand_path(file_path: str):
    if (path := Path.cwd() / file_path).exists():
        return str(path)
    if file_path.startswith("."):
        return os.path.join(os.getcwd(), file_path)
    elif file_path.startswith("/"):
        return file_path
    else:
        return os.path.join(ROOT_DIR, file_path)

class MolXpert:
    """
    all avialiable function for craton
    chem:
        molecule_create: 生成分子对象，支持的文件格式包括 smiles, mol, mol2, sdf, pdb, 
        format_convert: 转换分子的文件格式，例如从smiles到sdf格式，或从sdf到smiles
    """
    def __init__(self):
        pass

    def update_configure(self,usr_config):
        """
        更新configure参数
        默认参数见 craton.configure.configure.yaml
        
        :param usr_config: 用户自定义的参数
        """
        self.CRATON_CONFIGURE = total_update_configure(CRATON_CONFIGURE,usr_config)
        return total_update_configure(CRATON_CONFIGURE,usr_config)

    def help(cls):
        """
        默认包含的功能函数
        
        """
        members = inspect.getmembers(cls)
    
        # 属性（排除函数与特殊方法）
        attrs = [name for name, obj in members 
             if not (inspect.isroutine(obj) or name.startswith('__'))]
    
        # 方法
        methods = [name for name, obj in members 
                   if inspect.isfunction(obj) or inspect.ismethod(obj)]
        print(f"\n⚙️ Methods ({len(methods)}):")
        for name in methods:
            sig = str(inspect.signature(getattr(cls, name)))
            print(f"  • {name}{sig}")

    ###############chemical structure#########################################################
    ### 读输入文件，生成分子对象
    @staticmethod
    def molecule_create(input_files,extra_var=None,parallel=True,show_figure=True,opath=None,template_path=None):
        """
        读入smiles/mol/sdf/pdb等分子文件，生成craton内部的Molecule对象
        通常情况下，craton包所有的功能都需要最先运行该函数
        
        :param input_files: 输入的分子文件格式，可以是smiles, mol/sdf/pdb等文件，也可是一个smiles的数组，或一个包含所有分子文件的目录路径
        :param extra_var: 额外的参数
        :param parallel: 是否进行并行操作
        :param show_figure: 显示分子的结构图（是jupyter环境下适用）
        :param opath: (分子图片存储的路径)
        :param template_path: 临时文件的路径

        :return [Molecules]: 包括分子对象的数组（注意：如果只输入一个分子的smiles或文件，返回的仍然是Molecule对象的数组） 
        """

        from .chem import FormatMolecule as FM
        from .chem.molecule import Molecule

        if template_path is None:
            template_path = [CRATON_CONFIGURE['path']['molecule']] + [f"{CRATON_CONFIGURE['path']['molecule']}/{dd}" 
                              for dd in os.listdir(CRATON_CONFIGURE["path"]["molecule"]) 
                              if os.path.isdir(f"{CRATON_CONFIGURE['path']['molecule']}/{dd}" )]
        if not isinstance(input_files,list):
            input_files = [input_files]
        input_files_tmp = []
        for ff in input_files:
            if isinstance(ff,Molecule):
                input_files_tmp.append(ff)
            elif isinstance(ff,dict):
                input_files_tmp.append(ff)
            
            elif os.path.isdir(ff):
                input_files_tmp.extend([f"{ff}/{fff}" for fff in os.listdir(ff)])
            else:
                for pp in template_path:
                    if os.path.isfile(f"{pp}/{ff}.mtx"):
                        ff = f"{pp}/{ff}.mtx"
                        break
                input_files_tmp.append(ff)
        molecules = FM._parse(input_files_tmp,extra_var=extra_var,parallel=parallel)

        if NOTEBOOK_FLAG and show_figure:
            _tmp = []
            nn = 0
            for molecule in molecules:
                if hasattr(molecule,"Bonds"):
                    #if molecule.inchi_key not in _tmp:
                    nn += 1
                    #MolXpert.molecule_show(molecule,attrs=["ID"],opath=opath)
                    MolXpert.molecule_show(molecule,opath=opath)
                    _tmp.append(molecule.mole_name)
                    if nn >= 100:
                        break
        return molecules
    
    #### 文件格式转换
    @staticmethod
    def format_convert(molecules,otype=None,ofilename=None,opath=None,extra_var=None):
        """
        转换分子文件格式
        如将smiles转成sdf/mol/mol2等文件格式
        
        :param molecules: 分子对象的数组
        :param otype: 要转换的文件格式，可以是smiles, inchi key, mol/mol2/sdf/pdb等文件格式，或分子2D图，或一个csv文件
        :param ofilename: 新文件的名称（如果是None,将以分子的mole_nmae进行命名）
        :param opath: 新文件存储的路径
        :param extra_var: 额外的参数
        """

        from .chem import FormatMolecule as FM
        from .chem.molecule import Molecule

        if otype == "json":
            _attrs = {"inchi_key":"","inchi":"","smiles":"","mass":0.0,"formula":"",
                "elements":[],"connectivity":[],"bond_type":[],"formal_charge":[],"coordinates":[],
                "energy":None,"force":None,"dipole":None,"esp_charge":None,
                }
            datas = []
            for molecule in molecules:
                op = {attr:getattr(molecule,attr,value) for attr,value in _attrs.items()}
                if hasattr(molecule,"constrain"):
                    op["constrain"] = molecule.constrain[0].atoms
                    op["constrain_value"] = molecule.constrain[0].fix_value
                    constrain_type = "Bonds" if len(op["constrain"]) == 2 else "Angles" if len(op["constrain"]) == 3 else "Dihedrals"
                    op["constrain_type"] = None
                    for bb in getattr(molecule,constrain_type,[]):
                        if bb.atoms == op["constrain"] or bb.atoms[::] == op["constrain"]:
                            if hasattr(bb,"a1_atom_type"):
                                op["constrain_type"] = bb.atom_type_name
                            break

                datas.append(op)
            return datas

        if not isinstance(molecules,list):
            molecules = [molecules]
        if opath is not None:
            os.makedirs(opath,exist_ok=True)
        target_molecule = []
        for molecule in molecules:
            if not isinstance(molecule,Molecule):
                target_molecule.extend(MolXpert.molecule_create(molecule,extra_var=extra_var))
            else:
                target_molecule.append(molecule)
        texts = FM._convert(target_molecule,otype=otype,ofilename=ofilename,opath=opath,extra_var=extra_var)
        return texts

    @staticmethod
    def pdbqt_file(molecule,extra_var=None):
        from .chem import pdbqt_file
        return pdbqt_file(molecule,extra_var=extra_var)

    ### 分子拓扑结构
    @staticmethod
    def molecule_structure(molecules,ignore_existing=False,parallel=True):
        """
        分子拓扑结构分析，如查找环、环的大小、环的芳香性、是否共轭，以及其他结构相关的属性
        通常在生成分子对象后，需要进行该操作，以支持后续更多函数的功能
        
        :param molecules: 分子对象或分子对象的数组
        :param ignore_existing: 如果分子进行过该操作，是否忽略
        :param parallel: 是否进行并行操作

        :return [Molecules] :
        """

        from .chemkit import Structure as Stru
        return Stru._basic_structure_analyze(molecules,ignore_existing=ignore_existing,parallel=parallel)

    @staticmethod
    def molecule_prepare(molecules,ph_min=7.4,ph_max=7.4):
        from dimorphite_dl import protonate_smiles
        prepare_smiles=[]
        for molecule in molecules:
            prepare_smiles.append(".".join(protonate_smiles(molecule.smiles,ph_min=ph_min,ph_max=ph_max)))
        
        return MolXpert.molecule_create(prepare_smiles)


    ### 手性原子
    @staticmethod
    def molecule_chiral(molecules,parallel=True):
        from .chemkit import Structure as Stru
        return Stru._get_chiral_atom(molecules,parallel=parallel)

    ### 可旋转键
    @staticmethod
    def molecule_torsion(molecules,parallel=True):
        from .chemkit import Structure as Stru

        return Stru._get_flexible_torsion(molecules,parallel=parallel)

    ### 原子杂化类型
    @staticmethod
    def molecule_hybrid(molecules,parallel=True):
        from .chemkit import Structure as Stru

        return Stru._assign_hybrid(molecules,parallel=parallel)

    ### 潜在的作用位点
    @staticmethod
    def molecule_model(molecules):
        from .chemkit import Structure as Stru
        
        return Stru._model_atom(molecules)

    ### 链
    @staticmethod
    def chain_search(molecules,parallel=True):
        from .chemkit import Structure as Stru 

        return Stru._get_chain(molecules,parallel=parallel)

    ### 更新键长、键角、二面角等数值
    @staticmethod
    def update_structure_topol(molecules,parallel=True):
        from .chemkit import Structure as Stru

        return Stru._update_topol_values(molecules,parallel=parallel)

    ### 管能团
    @staticmethod
    def molecule_function_group(molecules,parallel=True):
        from .chemkit import MolScalpel as MolS
        return MolS._function_group(molecules,parallel=parallel)

    ### 管能团描述
    @staticmethod
    def function_group(molecules,parallel=True):
        from .chemkit import MolScalpel as Mols
        return Mols._function_group_info(molecules,parallel=parallel)

    ### 存储管能团描述
    @staticmethod
    def save_molecule_info(molecules,fpath=".",parallel=True):
        from .chemkit import save_txt
        save_txt(molecules,fpath=fpath,parallel=parallel)

    ### 分子碎片化
    @staticmethod
    def fragmentation(molecules,frag_types=None,frag_info=False,parallel=True):
        """
        将分子进行碎片化，根据输入可以得到一级、二级、三级等碎片
        
        :param molecules: 分子对象
        :param frag_types: 碎片类型
        :param frag_info: 是否保留碎片信息
        :param parallel: 是否进行并行操作

        :return : 分子对象
        """
        from .chemkit import MolScalpel as MolS
        
        return MolS._fragmentation(molecules,frag_types=frag_types,frag_info=frag_info,parallel=parallel)
    
    ### 分子粗粒化
    @staticmethod
    def atom_cluster(molecules,frag_types=None,frag_info=False,parallel=True):
        """
        将分子折成粗粒化的bead。bead为一种特殊的碎片
        进行粗粒化MD计算或进行粗粒化力场开发时，需要进行该操作
        
        :param molecules: 分子对象
        :param frag_types: 碎片的类型
        :param frag_info: 是否保留碎片的信息
        :param parallel: 是否并行操作
        """
        from .chemkit import MolScalpel as MolS
        
        return MolS._atom_cluster(molecules,parallel=parallel)
    
    ### 合并氢原子电荷
    @staticmethod
    def combine_hydrogen_charge_(molecule):
        """
        把非极性氢原子的电荷合并到与其相连的重原子上
        当分子拓扑结构去氢的时候，为保证电荷为电中性，需要进行该操作
        
        molecule: 分子对象
        """
        for term in ["esp_charge","mulliken_charge","apt_charge"]:
            if hasattr(molecule.Atoms[0],term):
                new_term = f"{term}_CH"
                for atom in molecule.Atoms:
                    if atom.elem != "H":
                        setattr(atom,new_term,getattr(atom,term))
                        for an in atom.connect:
                            if molecule.Atoms[an].elem == "H":
                                setattr(atom,new_term,getattr(atom,new_term) + getattr(molecule.Atoms[an],term))
        return molecule
    
    @staticmethod
    def molecule_show(molecule,attrs=None,save_file=True,opath=None,extra=None,show_image=False,TD_flag=False,return_flag=False):
        from .chem import _show_molecule
        if extra is not None and extra == "combine_H":
            new_molecule = MolXpert.combine_hydrogen_charge_(deepcopy(molecule))
            extra = None
            attrs = [attr if attr not in ["esp_charge","mulliken_charge","apt_charge"] else f"{attr}_CH" for attr in attrs]
        else:
            new_molecule = molecule
        molshow = _show_molecule(new_molecule,attrs=attrs,extra=extra,save_file=save_file,opath=opath,TD_flag=TD_flag)
        #if show_image:
        #    return molshow.imgs
        if NOTEBOOK_FLAG:
            imgs = [ img[0] if img[1] == "img" else SVG(data=img[0]) for img in molshow.imgs]
            display(*imgs)
        if return_flag:
            return molshow.imgs

    @staticmethod
    def figure_show(datas,figure_type,args=None,show_image=False):
        from .utils import _show_figure

        file_name = _show_figure(datas,figure_type,args=args)
        
        #if show_image:
        #    return file_name  
        if NOTEBOOK_FLAG:
            display(Image(filename=file_name))  
    
    ### 生成constrain
    @staticmethod
    def constrain_create(molecule,atoms,value=None):
        if value is None:
            value = MolXpert.structure_calculate(molecule,atoms)
        molecule.create_constrain(atoms+[value])

    ### molecule edit
    @staticmethod
    def structure_calculate(molecule,patoms) :
        """
        structure_calculate 的 Docstring
        
        :param molecule: 说明
        :param patoms: 说明
        """
        from .chemkit import MolEdit as ME

        return ME._structure_calculate(molecule,patoms)
    
    @staticmethod
    def structure_change(molecule,patoms,value,del_value=False,improper_flag=False):
        """
        对某个分子进行bond stretching, angle blending, diherdal rotation等结构微调操作

        :param molecule: 分子对象
        :param patoms: 要操作的原子编号，
                       如果两个原子，进行bond stretching编辑；
                       如果三个原子，进行angle blending编辑；
                       如果四个原子，进行diherdal rotation编辑
        :param value: 变化的数值，bond stretching范围不限，但推荐不过2.5，angle blending范围是0-180，diherdal rotation范围-180.0 - 180.0
        :param del_value: True, value表示在当前值基础改变的数值；False表示目标变化的绝对值
        :param improper_flag: 是否忽略
        """
        from .chemkit import MolEdit as ME
        return ME._structure_change(molecule,patoms,float(value),del_value=del_value,improper_flag=improper_flag)

    @staticmethod
    def molecule_topolgy_update(molecules1,molecules2,match_key=None):
        """
        molecule_topolgy_update 的 Docstring
        
        :param molecules1: 说明
        :param molecules2: 说明
        :param match_key: 说明
        """
        from .chemkit import MolEdit as ME

        ME._molecule_topolgy_update(molecules1,molecules2,match_key=match_key)

    @staticmethod
    def conformer_expand(molecules1,molecules2,attrs=["coordinates"],match_key="inchi_key"):
        """
        conformer_expand 的 Docstring
        
        :param molecules1: 说明
        :param molecules2: 说明
        :param attrs: 说明
        :param match_key: 说明
        """
        from .chemkit import MolEdit as ME

        return ME._conformer_expand(molecules1,molecules2,attrs=attrs,match_key=match_key)
    
    ### conformation
    @staticmethod
    def conformer_RMSD(molecule1,molecule2):
        from .chemkit import MolConformer as MConf #_conformer_RMSD

        return MConf._conformer_RMSD(molecule1,molecule2)
    
    ################################################################################

    ###############生物大分子#########################################################
    ### 生物大分子拓扑结构
    @staticmethod
    def protein_structure(protein):
        from .chemkit import Structure as Stru
        return Stru._protein_structure(protein)

    ### 蛋白结构准备
    @staticmethod
    def protein_prepare(protein):
        from .chemkit import Protein
        return Protein._pdb_prepare(protein)
    
    ### 蛋白mutation, modify等
    @staticmethod
    def protein_process(protein,arg):
        from .chemkit import Protein
        return Protein._pdb_process(protein,arg)

    ### 生物大分子原子匹配
    @staticmethod
    def protein_atom_mapping(protein1,protein2):
        from .chemkit import Protein
        return Protein._pdb_atom_mapping(protein1,protein2)
    
    ### 蛋白序列突变
    @staticmethod
    def protein_sequence_mutation(protein,sequences):
        from .chemkit import Protein
        return Protein._sequence_create_mutation(protein,sequences)
    
    ################################################################################

    ###############分子力场#########################################################
    @staticmethod
    def atom_type(molecules,atf=None,
                  atoms_arr=None,
                   assign_atom_type_flag=True,
                   check_at_types=False,
                   ignore_existing=False,
                   ignore_ff_existing=False,
                   this_terms=None,
                   parallel=True):
        """
        atom_type 的 Docstring
        
        :param molecules: 说明
        :param atf: 说明
        :param atoms_arr: 说明
        :param assign_atom_type_flag: 说明
        :param check_at_types: 说明
        :param ignore_existing: 说明
        :param this_terms: 说明
        :param parallel: 说明
        """
        from .force_field import AtomType


        AT = AtomType(atf=atf,check_at_types=check_at_types)
        return AT._atom_type(molecules,
                             atoms_arr=atoms_arr,
                             assign_atom_type_flag=assign_atom_type_flag,
                             this_terms=this_terms,
                             ignore_existing=ignore_existing,parallel=parallel)

    @staticmethod
    def _convert_atom_type_file_to_json_(atf):
        from .force_field import AtomType

        AtomType._convert_atom_type_file_to_json(atf)

    @staticmethod
    def _convert_force_field_file_to_json_(fff):
        from .force_field import MolForceField as MFF
        MFF._convert_force_field_file_to_json(fff)
    ###################
    
    ### force field ###
    @staticmethod
    def assign_force_field(molecules,this_ff,this_terms=None,parallel=True):
        from .force_field import MolForceField as MFF
        return MFF.assign_force_field_parameter(molecules,this_ff,this_terms=this_terms,parallel=parallel)

    @staticmethod
    def get_force_field(molecules,
                           atom_type_file=None,
                           force_field_file=None,
                           this_terms=None,
                           reassign_atom_type=True,
                           compensating_at_file=None,
                           compensating_ff_file=None,
                           charge_method=None,
                           charge_ff = None,
                           ignore_existing=False,
                           return_ff=False,
                           empi_ff_flag = True,
                           use_scalevdw=True,
                           parallel=True,
                           ):
        from .force_field import MolForceField as MFF

        return MFF.get_force_field(molecules,
                           atom_type_file,
                           force_field_file,
                           this_terms=this_terms,
                           reassign_atom_type=reassign_atom_type,
                           compensating_at_file=compensating_at_file,
                           compensating_ff_file=compensating_ff_file,
                           charge_method=charge_method,
                           charge_ff = charge_ff,
                           ignore_existing=ignore_existing,
                           return_ff=return_ff,
                           empi_ff_flag = empi_ff_flag,
                           use_scalevdw=use_scalevdw,
                           parallel=parallel,
                               )

    @staticmethod
    def grasp_force_field(molecules,
                           atom_type_file=None,
                           force_field_file=None,
                           reassign_atom_type=False,
                           charge_method=None,
                           charge_ff=None,
                           ignore_existing=False,
                           empi_ff_flag = True,
                           use_scalevdw=True,
                           return_ff=False,
                           parallel=True,
                           ):
        from .force_field import MolForceField as MFF

        return MFF.grasp_force_field(
                          molecules,
                          atom_type_file,
                          force_field_file,
                          reassign_atom_type=reassign_atom_type,
                          charge_method=charge_method,
                          charge_ff=charge_ff,
                          ignore_existing=ignore_existing,
                          empi_ff_flag=empi_ff_flag,
                          use_scalevdw=use_scalevdw,
                          return_ff=return_ff,
                          parallel=parallel
                          )
    
    @staticmethod
    def force_field_checkout(molecules,ffjson):
        from .force_field import MolForceField as MFF

        return MFF.checkout_force_field(molecules,ffjson)

    @staticmethod
    def get_fitting_parameters(this_ff,
            fix_tag=["V", "Fit"],
            terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
            preprocessing_fitting_parameter=True,
            molecules=None,
            flag_fitting_scan_parameter=False,
            atom_type_file=None,
            ):
        
        from .force_field import MolForceField as MFF
        return MFF.get_fitting_parameters(
                                        this_ff,
                                        fix_tag=fix_tag,
                                        terms=terms,
                                        preprocessing_fitting_parameter=preprocessing_fitting_parameter,
                                        molecules=molecules,
                                        flag_fitting_scan_parameter=flag_fitting_scan_parameter,
                                        atom_type_file=atom_type_file,
                                    )

    @staticmethod
    def _combine_ff(ff_ref, ff_new, prior_tags=None):
        from .force_field import MolForceField as MFF

        return MFF.combine_ff(ff_ref,ff_new,prior_tags=prior_tags)

    @staticmethod
    def _combine_ff_file(ff_file_ref, ff_file_new, used_vdw,out_path):
        from .force_field import MolForceField as MFF

        MFF.combine_ff_file(ff_file_ref, ff_file_new, used_vdw,out_path)

    @staticmethod
    def _write_ff_file(this_ff,opath):
        from .force_field import MolForceField as MFF
        
        MFF.write_ff_file(this_ff,opath) 

    @staticmethod
    def empirical_force_field_file(atf,fff1,fff2):
        from .force_field import MolForceField as MFF
        MFF.create_empirical_force_field(atf,fff1,fff2)

    @staticmethod
    def force_field_read(fff,use_scalevdw=True):
        from .force_field import MolForceField as MFF
        return MFF.read_force_field(fff,use_scalevdw)
    
    ### fitting ###
    @staticmethod
    def _fitting_validation(
                    molecules,
                    this_ff,
                    output_dir="./",
                    optimizer="openmm",
                    hessian_flag=False,
                    fitting_info=None,
                    init_this_ff=None,
                ):
        from .force_field import ForceFieldFitting as FFF

        return FFF._validation(molecules,
                                this_ff,
                                output_dir=output_dir,
                                optimizer=optimizer,
                                hessian_flag=hessian_flag,
                                fitting_info=fitting_info,
                                init_this_ff=init_this_ff
                                )
    
    @staticmethod
    def _binc_fitting(molecules, this_ff, target= "esp"):
        from .force_field import ForceFieldFitting as FFF

        FFF._binc_fitting(molecules,this_ff,target=target)

    @staticmethod
    def _intra_fitting(this_ff,
                        molecules,
                        fitting_terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
                        target_prop = ["energy", "force", "hessian", "penalty_torsion"],
                        torsion_constraint_step=None,
                        optimizer="openmm"
                        ):
        
        from .force_field import ForceFieldFitting as FFF
        return FFF._intra_fitting(this_ff,
                        molecules,
                        fitting_terms=fitting_terms,
                        target_prop = target_prop,
                        torsion_constraint_step=torsion_constraint_step,
                        optimizer=optimizer,
        )

    #### analyze
    @staticmethod
    def _mm_qm_analyze(molecules,
                       results_path="./",
                       force_field=None,
                       atom_type_file=None,
                       optimizer="openmm",
                       done_fitting=None,
                       init_this_ff=None,
                       validation_terms=None,
                       optimize_flag = True
                       ):
        from .force_field import FittingAnlayze as FA
        if validation_terms is None:
            validation_terms = ["energy", "pes", "esp_charge", "Bonds", "Angles", "Dihedrals", "rmsd", "Pair1n","hessian","freq"]
        return FA.analyze_qm_mm(
                                molecules,
                                results_path=results_path,
                                force_field=force_field,
                                atom_type_file=atom_type_file,
                                optimizer=optimizer,
                                done_fitting=done_fitting,
                                init_this_ff=init_this_ff,
                                validation_terms=validation_terms,
                                optimize_flag=optimize_flag
                            )

    @staticmethod
    def _mm_qm_result_show(data, param_val_data=None,save_path="./"):
        from .force_field import FittingAnlayze as FA # show_figure_fitting, show_figure_parameter
        prop_figure = []
        para_figure = []
        prop_figure = FA.show_figure_fitting(data,save_path=save_path)
        if param_val_data is not None:
            para_figure = FA.show_figure_parameter(param_val_data,save_path=save_path)
        if NOTEBOOK_FLAG:
            for fn in prop_figure:
                display(Image(filename=fn))
            for fn in para_figure:
                display(Image(filename=fn))
    ###########################################################################
    
    ##########dock#############################################################
    @staticmethod
    def pocket_analyze(protein,
                       csv_file_flag=False,
                       cavity_file_flag=False,
                       no_dispaly=True,
                       step=0.6,
                       probe_in=1.4,
                       probe_out=1.4,
                       removal_distance=2.4,
                       volume_cutoff=5,
                       include_depth=True,
                       include_hydropathy=False,
                       verbose=False,
                       output_directory="."):
        from .dock import get_pocket
        from .chem.molecule import Molecule
        if isinstance(protein,Molecule):
            MolXpert.format_convert(protein,otype="pdb",ofilename=f"{protein.mole_name}",opath=output_directory)
            protein_file = f"{output_directory}/{protein.mole_name}"
        else:
            protein_file = protein
        pre_file_name = protein_file[:-4]
        if csv_file_flag:
            csv_file = f"{output_directory}/{pre_file_name}_pocket.csv"
        if cavity_file_flag:
            cavity_file =f"{output_directory}/{pre_file_name}_pocket.pdb"

        return get_pocket(protein_file,
                          csv_file=csv_file,
                          cavity_file=cavity_file,
                          no_dispaly=no_dispaly,
                          step=step,
                          probe_in=probe_in,
                          probe_out=probe_out,
                          removal_distance=removal_distance,
                          volume_cutoff=volume_cutoff,
                          include_depth=include_depth,
                          include_hydropathy=include_hydropathy,
                          verbose=verbose)

    @staticmethod  
    def molecule_docking(protein,ligands,center,box_size,output_directory=".",parallel=True):
        from .dock import Dock
        docs = Dock(protein,ligands,center,box_size,output_directory=output_directory, parallel=parallel)
        return docs.docking()
    

    ###########################################################################

    ###############QM计算#######################################################
    @staticmethod
    def _ignore_alkane_torsion_(molecules):
        from .chemkit import MolConformer as MConf

        return MConf._ignore_alkane_torsion_(molecules)

    @staticmethod
    def scan_curve(molecules):
        from .chemkit import MolConformer as MConf

        return MConf._scan_curve(molecules)

    @staticmethod
    def scan_curve_data(scan_curve):
        from .chemkit import MolConformer as MConf
        
        return MConf._scan_curve_data(scan_curve)

    @staticmethod
    def scan_conf_type_(scan_curve):
        from .chemkit import MolConformer as MConf

        MConf._scan_conf_type(scan_curve)

    @staticmethod
    def pes_local_minimum(scan_curve):
        from .chemkit import MolConformer as MConf

        return MConf._pes_local_minimum(scan_curve)

    @staticmethod
    def lm_by_combine_scan_curve(molecules,rlm_dicts,n=64, create_constrain=False,parallel=True):
        from .chemkit import MolConformer as MConf

        return MConf._lm_by_combine_scan_curve(molecules,rlm_dicts,n=n,create_constrain=create_constrain,parallel=parallel)

    @staticmethod
    def find_stablest_conformer(molecules):
        from .chemkit import MolConformer as MConf

        return MConf._find_stablest_molecule(molecules)
    
    @staticmethod
    def remove_similar_conformer(molecules,target_molecule=None):
        from .chemkit import MolConformer as MConf

        return MConf._remove_similar_conformer(molecules,target_molecule=target_molecule)
    
    @staticmethod
    def Q8_bond_angle_scan(molecule,inter_val=[0.1,0.2,5.0,5.0],ignore_ring=True,exists_type=None):
        from .chemkit import MolConformer as MConf
        return MConf._get_bond_angle_scan_term(molecule,inter_val=inter_val,ignore_ring=ignore_ring,exists_type=exists_type)

    @staticmethod
    def Q2_bond_angle_conformer_(molecules,ignore_alkane=True):
        from .chemkit import MolConformer as MConf

        tmp_dict = {}
        for molecule in molecules:
            tmp_dict[molecule.inchi_key] = MConf._bond_angle_extend_conformer(molecule,ignore_alkane=ignore_alkane)
        return tmp_dict
    
    @staticmethod
    def qm_input_file(molecules,qmpara=None,step=None,local_path="./",indexs=None,zmatrixs=None,fpath_pre="",parallel=True):
        from .qm_calculation import QMCalc as QC
        if qmpara is None:
            if step is not None:
                if step in ["Q0","Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q100"]:
                    qmpara = CRATON_CONFIGURE[f"{step}QMSetting"]
                    #qm_setting_file = CRATON_CONFIGURE["setting_file"][f"{step}_SETTING_FILE"]

        if not isinstance(molecules,list):
            molecules = [molecules]
        if indexs is not None:
            if not isinstance(indexs,list):
                indexs = [indexs]
        if zmatrixs is not None:
            if not isinstance(zmatrixs,list):
                zmatrixs = [zmatrixs]
        
        QC.qm_input_file(molecules, qmpara=qmpara,step=step,local_path=local_path,indexs=indexs,zmatrixs=zmatrixs,fpath_pre=fpath_pre,parallel=parallel)

    ###########################################################################

    ###############MM计算#########################################################
    ### calculator
    @staticmethod
    def energy(molecules,prop="energy",parallel=True):
        from .mm_calculator import Calculator as Calc
        return Calc._energy(molecules,prop=prop,parallel=parallel)


    @staticmethod
    def _optimize(molecules, optimizer="openmm", all_torsion_constraint = 0.0, write_mol=None,):
        from .mm_calculator import Calculator as Calc
        return Calc._optimize(molecules,optimizer=optimizer,all_torsion_constraint=all_torsion_constraint,write_mol=write_mol)

    @staticmethod
    def _torsion_scan(molecules,scan_interval=30,parallel=True):
        from .mm_calculator import Calculator as Calc
        return Calc._torsion_scan(molecules,scan_interval=scan_interval,parallel=parallel)

    @staticmethod
    def _am1bcc_charge(molecules):
        from .mm_calculator import Calculator as Calc

        for molecule in molecules:
            Calc.am1bcc_charge(molecule)
        return molecules

    @staticmethod
    def molecule_volume_surface(molecules,parallel=True):
        from .mm_calculator import MolGeo as MG
        return MG(molecules,parallel=parallel).volume_surface()
    
    @staticmethod
    def molecule_multipole_moment(molecules,parallel=True):
        from .mm_calculator import Calculator as Calc
        for molecule in molecules:
            Calc.am1bcc_charge(molecule)
        from .mm_calculator import MolGeo as MG
        return MG(molecules,parallel=parallel).multipole_moments()
    
    @staticmethod
    def molecule_inertia(molecules,parallel=True):
        from .mm_calculator import MolGeo as MG
        return MG(molecules,parallel=parallel).moment_of_inertia()
    
    @staticmethod
    def molecule_center(molecules,parallel=True):
        from .mm_calculator import MolGeo as MG
        return MG(molecules,parallel=parallel).molecule_center()

    ###########################################################################


    @staticmethod
    def insert_to_db(data_type,molecules,config=None):
        if DATABASE_EXISTS:
            from ..database import DataDB as DDB
            DDB.db_insert(molecules,data_type=data_type,config=config)
        else:
            sys.exit("数据库模块不存在")
    
    @staticmethod
    def get_from_db(config=None):
        if DATABASE_EXISTS:
            from ..database import DataDB as DDB
            return DDB.db_get(config=config)
        else:
            sys.exit("数据库模块不存在")
    
    @staticmethod
    def conformation_id(molecules):
        from .chemkit.conformation import conformation_id_hash
        for molecule in molecules:
            hash_id = conformation_id_hash(molecule)
            molecule.confID = hash_id

    ###builder
    @staticmethod
    def _builder(configure):
        from .builder.builder import Builder
        builder = Builder(configure)
        return builder.run()
    
    @staticmethod
    def old_builder(molecules,config=None,style=None,parallel=True):
        from .builder import builder_manager

        return builder_manager(molecules,config=config,style=style,parallel=parallel)

    @staticmethod
    def builder(molecules,config=None,style=None,parallel=True):
        from .builder import builder_manager

        return builder_manager(molecules,config=config,style=style,parallel=parallel)

    @staticmethod
    def build_property_system(sms):
        from .builder import make_propety_system
        return make_propety_system(sms)
    ### atom mapping for rfe calculation
    @staticmethod
    def init_pair_network(ligands,
                      topology="normal",
                      user_pair_list=None,
                      bias_nodes=None,
                      core=None,
                      nbunch=None):
        from .md_simulation import FEPTool as FEPT #pair_network_init

        return FEPT.pair_network_init(ligands,
                                 topology=topology,
                                 user_pair_list=user_pair_list,
                                 bias_nodes=bias_nodes,
                                 core=core,
                                 nbunch=nbunch
                                 )
    
    @staticmethod
    def final_pair_nework(gg,nbunch=None,bias_nodes=None):
        from .md_simulation import FEPTool as FEPT #pair_network_final

        FEPT.pair_network_final(gg,nbunch=nbunch,bias_nodes=bias_nodes)
    
    @staticmethod
    def atom_mapping(gg):
        from .md_simulation import FEPTool as FEPT #atom_mapping_calculate

        FEPT.atom_mapping_calculate(gg)
    
    @staticmethod
    def molecule_similiarity(gg):
        from .md_simulation import FEPTool as FEPT #molecule_similiarity_calculate
        
        FEPT.molecule_similiarity_calculate(gg)
    
    @staticmethod
    def graph_attributes(gg,topology="normal"):
        from .md_simulation import FEPTool as FEPT #graph_attributes_report

        FEPT.graph_attributes_report(gg,topology=topology)

    @staticmethod
    def dual_topology(fep_type,gg,output_directory=".",parallel=True):
        from .md_simulation import FEPTool as FEPT #assign_dual_topology

        return FEPT.assign_dual_topology(fep_type,gg,output_directory=output_directory,parallel=parallel)
    
    @staticmethod
    def get_fep_lambda(fep_setting,fep_type="r_group",mixed_lambda=False,is_relative=False):
        from .md_simulation import FEPTool as FEPT #get_lambda_schedule
        return FEPT.get_lambda_schedule(fep_setting,fep_type=fep_type,mixed_lambda=mixed_lambda,is_relative=is_relative)

    @staticmethod
    def get_intermolecule_interaction(system):
        from .md_simulation import FEPTool as FEPT #abfe_intermolecule
        return FEPT.abfe_intermolecule(system)

    @staticmethod
    def write_md_input_files(systems,parallel=True):
        from .md_simulation import MDSimulation as MDS #write_input_file
        MDS.write_input_file(systems,parallel=parallel)

    @staticmethod
    def write_bash_files(systems,parallel=True):
        from .md_simulation import MDSimulation as MDS #write_bash_file
        MDS.write_bash_file(systems,parallel=parallel)

    @staticmethod
    def write_job_infos(systems,parallel=False):
        from .md_simulation import MDSimulation as MDS #write_info
        MDS.write_info(systems,parallel=parallel)

    @staticmethod
    def mdrun_check(directory,mdengine="gmx",batchfile="batch_0.txt",parallel=True):
        from .md_simulation import MDAnalyze as MDA #check_md
        return MDA.check_md(directory,mdengine="gmx",batchfile=batchfile,parallel=parallel)

    #@staticmethod
    #def analyze_rbfe(input_dir,output_dir="md_results",exp_file=None, pka_file=None, two_stages=False,parallel=True):
    #    from .md_analyze import rbfe_analyze
    #    rbfe_analyze(input_dir, output_dir=output_dir, exp_file=exp_file, pka_file=pka_file, two_stages=two_stages,parallel=parallel)

    @staticmethod
    def analyzer_gmx(ptype,args):
        from .md_simulation import MDAnalyze as MDA #gmx_analyzer

        MDA.gmx_analyzer(ptype,args)

    # these are tools
    @staticmethod
    def ic50_and_free_energy(stype,value,temp, unit):
        from .fep.utils import ic50_to_free_energy, free_energy_to_ic50
        if stype in ["ic50_to_g","ic50tog","ic502g"]:
            return ic50_to_free_energy(value,temp,unit)
        elif stype in ["g2ic50","g_to_ic50","gtoic50"]:
            return free_energy_to_ic50(value,temp,unit)
        else:
            return None

    @staticmethod
    def get_property(inputs,props,molecule_type=None,sources=None,temperatures=None,pressures=None,condinations=None):
        from .property import MolProperty as MP
        return MP.thermodyna_property(inputs,props,molecule_type=molecule_type,sources=sources,temperatures=temperatures,pressures=pressures,condinations=condinations)
    
    @staticmethod
    def find_interaction_model(molecule,probe):
        from .chemkit import InteractionModel as IM #get_interaction_model

        return IM.get_interaction_model(molecule,probe)
    
    @staticmethod
    def assign_AA_atom_name(molecules):
        from .chemkit.biomacromolecule import Protein
        return Protein._assign_AA_atom_name(molecules)
    
    @staticmethod
    def AA_registered(molecules):
        from .chemkit.biomacromolecule import Protein
        Protein._register_non_AA(molecules)
    
    @staticmethod
    def create_AA_template(molecules):
        from .chemkit.biomacromolecule import Protein
        return Protein._create_AA_template(molecules)
    
    @staticmethod
    def create_peptide(n,left_cap="ACE",right_cap="NME",terminal_flag=True,templates=None):
        from .molgen import create_aminoacid
        return create_aminoacid(n,left_cap=left_cap,right_cap=right_cap,terminal_flag=terminal_flag,templates=templates)
    
    @staticmethod
    def create_dnarna(n,templates=None):
        from .molgen import create_rnadna
        return create_rnadna(n,templates=templates)
    
    #####utils
    @staticmethod
    def _aminoacid_json(ff):
        from .software import AmberUtils #aminoacid_json
        AmberUtils.aminoacid_json(ff)
        
    @staticmethod
    def _amberff_to_ff(atf,nonbf,bondf):
        from .software import AmberUtils
        AmberUtils.amberff_to_ff(atf,nonbf,bondf)
        
    @staticmethod
    def analyze_chem_space(inf,output_dir="."):
        from .chemkit import chem_space_analyze
        chem_space_analyze(inf,output_dir=output_dir)
    
    @staticmethod
    def get_all_property(inf,temperature=298.15,output_dir=".",outfn="property"):
        from .property import MolProperty as MP
        results,error = MP.get_all_property_at_temperature(inf,temperature,output_dir=output_dir,outfn=outfn)
        
    @staticmethod
    def split_train_test(inf,prop,style,value=None,output_dir=".",test_flag=True):
        from .property import MolProperty as MP
        MP.split_train_set(inf,prop,style,value=value,output_dir=output_dir,test_flag=test_flag)
        
    @staticmethod
    def get_pubchem(input,input_type="smiles",output_type="name",file_flag=False):
        from .chemkit import get_moleinfo_from_pubchem
        return get_moleinfo_from_pubchem(input,input_type=input_type,output_type=output_type,file_flag=file_flag)

    @staticmethod
    def pubchem_info(strs,typ="smiles",print_flag=True,opath="."):
        from .utils import _pubchem_info
        return _pubchem_info(strs,typ,print_flag=print_flag,opath=opath)
    
    @staticmethod
    def pdb_file(info_file=None,pdb_id=None,output_dir=".",output_format="pdb"):
        from .utils import _get_pdb_file
        return _get_pdb_file(info_file=info_file,pdb_id=pdb_id,output_dir=output_dir,output_format=output_format)
    
    @staticmethod
    def uniport_info(target=None,uniprot_id=None,output_directory="."):
        from .utils import _get_uniport
        return _get_uniport(target=target,uniprot_id=uniprot_id,output_directory=output_directory)

molxpert = MolXpert()
