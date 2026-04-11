import os
from copy import deepcopy
from pathlib import Path
from craton import molxpert as MX
from craton import CRATON_CONFIGURE
from craton.craton.utils import logger
from craton.craton.chem.molecule import Molecule
MAX_TASK_NUM = 10
####from craton.craton.force_field import default_complement_force_field


class ForceFieldFitting:
    """
    fitting的主入口。根据不同的设置可以进行fitting,和validation
    输入：
        inchi_key_arr: List[inchi_key], 一系列的inchi_key。如果该项为空，后面的smiles才起作用
                        该项不为空，则从DB生成分子和获取QM数据
        smiles: List[smiles], 一系列的smiles。从smiles 生成分子，并从g09 log文件拿到QM数据
        atom_type_file: file, 自有力场的原子类型定义文件
        force_field_file: file, 自有力场的参数文件
        gaff_atom_type_file: file or None, gaff力场的原子类型定义文件
        gaff_cover_atom_type_file: file or None, gaff力场的原子类型转换文件
        gaff_force_field_file: file or None, gaff力场的参数文件
        empirical_atom_type_file: file or None, 经验力场的原子类型定义文件
        empirical_force_field_file: file or None, 经验力场的参数文件
        qm_data_path: None or str, qm log文件的路径，与smiles配合使用
        charge_method: None or str, 是否采用非力场中的电荷。None表示用力场中的电荷，其他还有esp, nn
        only_validation_flag: True or False, 是否只是做validation
        fitting_scan_para_flag: True or False, 是否只拟合有scan数据的参数
        fitting_terms: List[str], 需要拟合的参数类型，如bondterm, angleterm, dihedralterm,binc等
        fitting_fourier_n: None or List[int]。确定那些傅立叶展开形式的二面角参数需要拟合。如果为None则所有参数都拟合
        target_prop: List[str], 拟合的目标性质
        out_put_dir: path, 结果文件保存的路径
    输出：
        total_ff.ff: file, 拟合后的参数文件。本次拟合的参数与force_field_file的参数合并得到的文件
        fitting_result.json.gz: file, 拟合或验证的结果文件
    """
    
    def __init__(self,config=None) -> None:
        self.user_config = config
    
    @staticmethod
    def atom_type(inputs,atf=None,output_directory="."):
        if output_directory != ".":
            Path(output_directory).mkdir(exist_ok=True)
        molecules = MX.molecule_create(inputs)
        molecules = MX.molecule_structure(molecules)
        molecules = MX.atom_type(molecules,atf=atf)
        with open(f"{output_directory}/atom_type.txt",'w') as outf:
            for molecule in molecules:
                MX.molecule_show(molecule,attrs=["atom_type_name"],opath=output_directory)
                outf.write(f"{molecule.mole_name}\n")
                for atom in molecule.Atoms:
                    outf.write(f"{atom.ID} {atom.elem} {atom.atom_type_name}\n")
                outf.write("\n")

    @staticmethod
    def grasp_force_field(inputs,atf=None,fff=None,output_directory="."):
        if output_directory != ".":
            Path(output_directory).mkdir(exist_ok=True)
        molecules = MX.molecule_create(inputs)
        molecules = MX.molecule_structure(molecules)
        molecules = MX.atom_type(molecules,atf=atf)
        molecules = MX.grasp_force_field(molecules,atom_type_file=atf,force_field_file=fff)
        
        MX.format_convert(molecules,otype="mtx",opath=f"{output_directory}",extra_var="all")
        
    @staticmethod
    def am1bcc_charge(inputs,output_directory="."):
        if output_directory != ".":
            Path(output_directory).mkdir(exist_ok=True)
        molecules = MX.molecule_create(inputs)
        molecules = MX.molecule_structure(molecules)
        molecules = MX.atom_type(molecules, atf=None)
        molecules = MX._am1bcc_charge(molecules)
        with open(f"{output_directory}/am1bcc_charge.txt",'w') as outf:
            for molecule in molecules:
                MX.molecule_show(molecule,attrs=["ff_charge"],opath=output_directory)
                outf.write(f"{molecule.mole_name}\n")
                for atom in molecule.Atoms:
                    outf.write(f"{atom.ID} {atom.elem} {atom.ff_charge}\n")
                outf.write("\n")   
        
    def update_config(self):
        self.config = MX.update_configure(self.user_config)
        Path(self.config['EnvironmentSetting']['output_directory']).mkdir(exist_ok=True)
        self.molecule_paths = [self.config['path']['molecule']] + [f"{self.config['path']['molecule']}/{dd}" 
                              for dd in os.listdir(self.config["path"]["molecule"]) 
                              if os.path.isdir(f"{self.config['path']['molecule']}/{dd}" )]

        self.guess_fitting_para = self.config["FFFitting"]["using_force_field_file"] is None
        self.output_dir = self.config["EnvironmentSetting"]["output_directory"]
        
        
        if self.config["FFFitting"]["using_force_field_file"] is None:
            self.config["FFFitting"]["using_force_field_file"] = self.config["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"]

        
        if self.config["FFFitting"]["fitting_method"] == "st":
            self.torsion_constraint_step = None  # use default value
        elif self.config["FFFitting"]["fitting_method"] == "mst":
            self.torsion_constraint_step = [1500.0, 50.0, 4.0, 0.0]
        else:
            logger.error("fitting_method must be st or mst")
            raise Exception("Invalid fitting_method")
        self.fitting_terms = self.config["FFFitting"]["fitting_terms"]
        self.target_prop = self.config["FFFitting"]["target_prop"]
        self.optimizer = self.config["FFFitting"]["optimizer"]
        self.force_field_file = self.config["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"]
        self.validation_terms = self.config["FFFitting"]["validation_terms"]
        self.optimize_flag = self.config["FFFitting"]["optimize_flag"]
        self.used_vdw = self.config["ForceFieldSetting"]["use_scalevdw"]
        
        
        self.input_files = self.config["FFFitting"]["fitting_molecules"]

    def run(self):
        self.update_config()
        if self.input_files is None:
            return
        self.create_molecule()
        self.get_force_field()
        self.expend_molecule()
        if self.config["FFFitting"]["only_validation_flag"]:
            self.run_validation()
        else:
            self.get_fitting_parameter()
            if "binc" in self.fitting_terms:
                self.run_binc_fitting()
            self.run_intra_fitting()
        self.validation_figure()
        
    def create_molecule(self):
        if not isinstance(self.input_files[0],Molecule):
            extra_var ={"datasearch":{"data_type":"qmdata","compound_style":"molecule",}}
            self.training_set = MX.molecule_create(self.input_files,extra_var=extra_var)
        else:
            self.training_set = self.input_files

    def get_force_field(self):
        model_molecules = MX.find_stablest_conformer(self.training_set)
        
        model_molecules = MX.molecule_structure(model_molecules)
        
        self.model_molecules,total_ff = MX.grasp_force_field(model_molecules,
                           atom_type_file=self.config["ForceFieldSetting"]["DEFAULT_TYPING_FILE"],
                           force_field_file=self.config["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"],
                           reassign_atom_type=True,
                           charge_method=self.config["FFFitting"]["charge_method"],
                           ignore_existing=False,
                           empi_ff_flag = True,
                           use_scalevdw=True,
                           return_ff=True,
                           parallel=True,
                           )
        self.this_ff = MX.force_field_checkout(self.model_molecules,total_ff)
        
    def expend_molecule(self):
        if set(self.fitting_terms) - {"binc"}:
            self.total_molecules = []
            self.total_molecules = MX.conformer_expand(self.model_molecules,self.training_set,attrs=["coordinates",
                                                                                       "esp_charge",
                                                                                            "energy", 
                                                                                            "force", 
                                                                                            "hessian", 
                                                                                            "freq", 
                                                                                            "constrain", 
                                                                                            "scan_term", 
                                                                                            "conform_type",
                                                                                            "confID"],
                                                                                    )
            for molecule in self.total_molecules:
                if molecule.conform_type not in ["local minimum"]:
                    if hasattr(molecule,"force"):
                        delattr(molecule,"force")
        else:
            self.total_molecules = self.model_molecules[:]

        self.total_molecules = MX.update_structure_topol(self.total_molecules)
        
    def run_validation(self):
        MX._write_ff_file(self.this_ff,self.output_dir + "/validation_ff.ff")
        self.results = MX._mm_qm_analyze(
                                    self.total_molecules,
                                    results_path=self.output_dir,
                                    force_field=self.this_ff,
                                    optimizer=self.optimizer,
                                    done_fitting=[],
                                    init_this_ff= self.this_ff,
                                    validation_terms=self.validation_terms,
                                    optimize_flag = self.optimize_flag
                                    )
        
    def get_fitting_parameter(self):
        self.this_ff = MX.get_fitting_parameters(
                                        self.this_ff,
                                        fix_tag=["V", "Fit","amber99sb"],
                                        terms=self.config["FFFitting"]["fitting_terms"],
                                        preprocessing_fitting_parameter=self.guess_fitting_para,
                                        molecules=self.total_molecules,
                                        flag_fitting_scan_parameter=self.config["FFFitting"]["fitting_scan_para_flag"],
                                        atom_type_file=self.config["ForceFieldSetting"]["DEFAULT_TYPING_FILE"],
                                        )

        self.init_this_ff = deepcopy(self.this_ff)
        MX._write_ff_file(self.this_ff, self.output_dir + "/init_ff.ff")
    
    def run_binc_fitting(self):
        logger.info("Binc parameters Fitting ......")
        MX._binc_fitting(self.total_molecules,self.this_ff,target="esp")
        logger.info("Binc parameters Done")

    def run_intra_fitting(self):

        logger.info("Bonded parameters Fitting ......")

        self.fitting_results=MX._intra_fitting(
                                                self.this_ff,
                                                self.total_molecules,
                                                fitting_terms=self.fitting_terms,
                                                target_prop=self.target_prop,
                                                torsion_constraint_step=self.torsion_constraint_step,
                                                optimizer=self.optimizer
                                            )
        logger.info("Bonded parameters Done")
        MX._write_ff_file(self.this_ff, os.path.join(self.output_dir, "this_ff.ff"))
        MX._combine_ff_file(
                        self.force_field_file,
                        os.path.join(self.output_dir, "this_ff.ff"),
                        self.used_vdw,
                        os.path.join(self.output_dir, "total_ff.ff"),
                )
        new_total_molecules = MX.assign_force_field(self.total_molecules,self.this_ff)
        new_total_molecules = MX.update_structure_topol(new_total_molecules)
        self.results = MX._mm_qm_analyze(
                                    new_total_molecules,
                                    results_path=self.output_dir,
                                    force_field=self.this_ff,
                                    optimizer=self.optimizer,
                                    done_fitting=self.fitting_results,
                                    init_this_ff= self.init_this_ff,
                                    validation_terms=self.validation_terms,
                                    optimize_flag = self.optimize_flag
                                )
        
    def validation_figure(self):        
        if self.config["FFFitting"]["validation_figure_flag"]:
            MX._mm_qm_result_show(
                                    self.results["results"]["properties"],
                                    param_val_data=self.results["results"]["parameters"],
                                    save_path=self.output_dir
                                )

def run_auto_force_field_fitting(init_force_field_file,output_dir,jobs,validation_flag=False):
    input_force_field_file = [init_force_field_file]
    for ii,job in enumerate(jobs):
        ss = job.split("\/")
        input_force_field_file.append(f"{output_dir}/{ii}-{ss[-1]}/total_ff.ff")
    for ii,job in enumerate(jobs):
        ss = job.split("\/")
        config = {
                    "output_directory":f"{output_dir}/{ii}-{ss[-1]}",
                    "DEFAULT_FORCE_FIELD_FILE":init_force_field_file,
                    "fitting_molecules": job,
                    "only_validation_flag":validation_flag,
                    }
        FFF = ForceFieldFitting(config=config)
        FFF.run()

def meta_force_field_fitting(init_force_field_file,output_dir,validation_flag=False,nn=None):
    if nn is None:
        nn = MAX_TASK_NUM
    from craton.craton.chem.database.mongodb import MongoDB
    db = MongoDB()
    coll = db.ff_develop_coll
    ii = max(list(coll.distinct("order")))
    for i in range(ii):
        doc = coll.find_one({"order":i})
        if "fit" not in doc:
            break
        else:
            if doc["fit"] != "OK":
                break
    job_ids = [k for k in range(i,nn+1)]

    jobs = []
    input_ff_file = init_force_field_file

    for job_id in job_ids:
        doc = coll.find_one({"order":job_id})
        job_dir = f"{output_dir}/{job_id}.{doc['function_group'][:-4]}"
        os.system("mkdir %s" %job_dir)
        job_txt = f"{output_dir}/{job_id}.{doc['function_group'][:-4]}/{doc['function_group']}"
        with open(job_txt,'w') as outf:
            outf.write("inchi_key, smiles, frag_level, topol_label, function_group_label, frag_number, elements, \n")
            for rr in doc["molecules"]:
                for rrr in rr:
                    outf.write("%s,"%rrr)
                outf.write("\n")
        
        jobs.append({
                    "output_directory":job_dir,
                    "DEFAULT_FORCE_FIELD_FILE":input_ff_file,
                    "fitting_molecules": job_txt,
                    "only_validation_flag":validation_flag,
                    "job_id":job_id,
                    })
        input_ff_file = f"{job_dir}/total_ff.ff"

    for job in jobs:
        FFF = ForceFieldFitting(config=job)
        FFF.run()
        coll.update_one({"order":job["job_id"]},{"$set":{"fit":"OK"}})



if __name__ == "__main__":
    pass
