import textwrap
from pathlib import Path
import yaml

import click
import pandas as pd
#import simplejson

from ..application.force_field_fitting import ForceFieldFitting, meta_force_field_fitting
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="benchmark",context_settings=CONTEXT_SETTINGS)
def benchmark():
    """MD simulation benchmark

    \b
    * liquid, target property include density, HOV et al.
    * hfe, hydration free energy
    * rbfe, relative binding free energy for R group replace
    * rbfe-c, relative binding free energy for core hopping and charge change
    * abfe, abs binding free energy
    * pep-rbfe, relative binding free energy for peptide
    * mem-rbfe, relative binding free energy for membre protein, such as GPCR
    * cov-rbfe, relative binding free energy for covalence drug
    * mutation, the relative binding free energy of protein mutation
    * antibody-stble, the relative stability fo antibody
    * rna-rbfe, relative binding free energy for rna systems
    * rlogp, relative logp
    * alogp, abs logp
    * rlogs, relative solubility
    """


#try:
#    from craton.craton.fep.parse_gpickle import show_graph

#    fep_analyze.add_command(show_graph)
#except ModuleNotFoundError:
##    dash_found = False
#    logger.warn("dash cannot be found, the show_graph command cannot be used")


@benchmark.command("run")
@click.argument(
    "benchmark_type",
    type=click.Choice(
        ["liquid","hfe","rbfe","abfe","pep-rbfe","pep_rbfe","mem-rbfe","mem_rbfe",
         "mutation","rlogp",
         ]
         
    ),
)
@click.option("-i","--input_files",default=".",show_default=True,help="different job type base on different directory",)
@click.option("-o", "--output_directory",default=".", show_default=True,help="The directory of output",)
@click.option("-if", "--init_force_field_file",default="./0.ff", show_default=True,help="The init force field file",)
@click.option("-f", "--yaml_file", type=click.File("r"), help="if need more fined control, please specify a yaml file")
def fitting(fit_type,input_files, output_directory,init_force_field_file,yaml_file):
    """
    yaml file setting parameters：
        fitting_molecules: files or inchi keys or smiles, if smiles, get fitting data form gaussain log file
        atom_type_file: file, 自有力场的原子类型定义文件
        force_field_file: file, 自有力场的参数文件
        #gaff_atom_type_file: file or None, gaff力场的原子类型定义文件
        #gaff_cover_atom_type_file: file or None, gaff力场的原子类型转换文件
        #gaff_force_field_file: file or None, gaff力场的参数文件
        #empirical_atom_type_file: file or None, 经验力场的原子类型定义文件
        #empirical_force_field_file: file or None, 经验力场的参数文件
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
    if fit_type == "autofit":
        meta_force_field_fitting(init_force_field_file,output_directory)
        return

    elif fit_type == "yaml":
        if yaml_file is None:
            raise RuntimeError("When choose yaml, the yaml file must be supported followed by -f ")
        config = yaml.safe_load(yaml_file.read())
    else:
        config = {"fitting_molecules":input_files,"output_directory":output_directory}
        
    FFF = ForceFieldFitting(config=config)
    FFF.run()

