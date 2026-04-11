import textwrap
from pathlib import Path

import click
#import pandas as pd
#import simplejson

from craton import molxpert as MX
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(name="analyze",context_settings=CONTEXT_SETTINGS)
def analyze():
    """Analyzing the fep molecule_dynamics molecule_dynamics
    The analyzing sub-group containing following commands:

    \b
    * md, including bar, mbar, exchange, rmsd, torsion, interaction et al
    * rbfe, all pair directorys
    * all_extra, all pair directorys
    """


#try:
#    from craton.craton.fep.parse_gpickle import show_graph

#    fep_analyze.add_command(show_graph)
#except ModuleNotFoundError:
##    dash_found = False
#    logger.warn("dash cannot be found, the show_graph command cannot be used")


@analyze.command("gmx")
@click.argument(
    "analyze_type",
    type=click.Choice(
        ["bar", "ddg","rmsd", "energy","torsion","interaction","pair_interaction","pair-interaction","fep_exchange","fep-exchange",
         "block_ddg","block-ddg","accum_ddg","accum-ddg",
         "rbfe","abfe","complex","protein","normal",
         "all-rbfe","all-abfe","all_rbfe","all_abfe","dimer","all_dimer","all-dimer","total-dimer","total_dimer"]
    ),
)
@click.option("-i","--input_directory",default=".",show_default=True,help="different job type base on different directory",)
@click.option("-o", "--output_directory",default="md_result", show_default=True,help="The directory of output",)
@click.option("-exp", "--expt_file", default=None,show_default=True,help="Experiemnt dG file, this could be csv file or gpickle file")
@click.option("--two_stages/--no-two_stages", default=False, show_default=True, help="If two stages are used?")
@click.option("-pka", "--pka-file", default=None, show_default=True, help="pka file")
@click.option("-mol","--molecule_dir",default=None,show_default=True,help="the directory of molecule file")
@click.option("-attrs","--attributes",default=None,show_default=True,help="the attrs for ......")
@click.option("-dimer_type","--dimer_type",default=None,show_default=True,help="the attrs for ......")
def md_analyze(analyze_type,input_directory, output_directory, expt_file, two_stages, pka_file,molecule_dir,attributes,dimer_type):
    """
    gmx md result analyze:\n
        "bar": calculate the bar of free energy simulation\n
        "ddg": calculate deltal free energy from the bar of two states\n
        "block-ddg":\n
        "accum-ddg":\n
        "rmsd": rmsd calculate of the simulation include protein, default group = Protein,LIG\n
        "energy": energy analyze, default attributes = Potential, Temperature, Pressure\n
        "torsion": ligand torsion distribution analyze for the complex, rbfe, abfe et. al\n
        "interaction": interaction model analyze between ligand and protein\n
        "pair-interaction": compare the interaction model between A and B ligand for rbfe\n
        "fep-exchange": the exchange ratio for fep calculation\n
        "rbfe": data analyze for single pair rbfe calculation\n
        "abfe": data analyze for single abfe calculation\n
        "complex": data analyze for single complex calculation\n
        "protein": data analyze for single protein calculation\n
        "normal": data analyze for single normal calculation\n
        "all-rbfe": data analyze for a rabe task\n
        "all-abfe" data analyze for a abfe task\n
    """
    if attributes is not None:
        if attributes != "all":
            attributes = attributes.split(":")
    args = {
            "job_dir":Path(input_directory),
            "output_dir":Path(output_directory),
            "expt_file":expt_file,
            "pka_file":pka_file,
            "molecule_dir":molecule_dir,
            "attrs": attributes,
            "dimer_type":dimer_type
            }
    MX.analyzer_gmx(analyze_type,args)

@analyze.command("get-property")
@click.option("-i","--input_file",default=".",show_default=True,help="input_file")
@click.option("-o", "--output_directory", default=".", show_default=True,help="The directory of output")
@click.option("-t", "--temperature", default=298.15,show_default=True,help="the temperature to calculation property")
@click.option("-f", "--output_file", default="property", show_default=True, help="output file",)
def get_properties(input_file,output_directory, temperature, output_file):
    """
    不同温度（或压强）下的实验数据，如果数据与温度、压强无关，则不考虑温度等环境
    -i: experimental file\n
    -t: temperature\n
    -o: output idrectory\n
    -f: output file\n
    """
    MX.get_all_property(input_file,temperature=temperature,output_dir=output_directory,outfn=output_file)

@analyze.command("property-bin")
@click.option("-i","--input_file",default=".",show_default=True,help="input_file")
@click.option("-o", "--output_directory", default=".", show_default=True,help="The directory of output")
@click.option("-b", "--bin_number", default=10,show_default=True,help="bin_number")
def get_properties_value_anlayze(input_file,output_directory, bin_number):
    """
    分析实验数据的分布，用到get_properties生成的文件
    -i: experimental file\n
    -b: bin number\n
    -o: output idrectory\n
    """
    from craton.craton.property.expt_equ import property_value_analyze
    property_value_analyze(input_file,output_dir=output_directory,bin_num=bin_number)

@analyze.command("admet-bin")
@click.option("-i","--input_file",default=".",show_default=True,help="input_file")
@click.option("-o", "--output_directory", default=".", show_default=True,help="The directory of output")
@click.option("-b", "--bin_number", default=10,show_default=True,help="bin_number")
def get_properties_value_anlayze(input_file,output_directory, bin_number):
    """
    分析ADMET数据的分布
    -i: experimental file\n
    -b: bin number\n
    -o: output idrectory\n
    """
    from craton.craton.property.expt_equ import admet_value_analyze
    admet_value_analyze(input_file,output_dir=output_directory,bin_num=bin_number)

@analyze.command("property-info")
@click.option("-i","--input_file",show_default=True,help="chem info input_file")
@click.option("-p", "--property_file",show_default=True,help="property file")
@click.option("-o", "--output_file", show_default=True,help="output file")
def get_properties_value_anlayze(input_file,property_file, output_file):
    """
    在实验数据文件中，添加分子的信息，如重原子数、环大小等
    -i: chem info input_file\n
    -p: property file\n
    -o: output file\n
    """
    from craton.craton.property.expt_equ import get_expt_info_file
    get_expt_info_file(input_file,property_file,output_file)

@analyze.command("property-block")
@click.option("-i","--input_file",show_default=True,help="block property file")
@click.option("-p", "--property_file",show_default=True,help="property file")
@click.option("-o", "--output_file", show_default=True,help="output file")
def get_properties_value_anlayze(input_file,property_file, output_file):
    """
    使某些分子的实验数据失效
    -i: block property file\n
    -p: property file\n
    -o: output file\n
    """
    from craton.craton.property.expt_equ import unavil_record_of_property
    unavil_record_of_property(input_file,property_file,output_file)

@analyze.command("property-result")
@click.option("-i","--input_path",show_default=True,help="input path")
def get_properties_value_anlayze(input_path):
    """
    性质预测的结果生成
    -i: block property file\n
    -p: property file\n
    -o: output file\n
    """
    from craton.craton.property.expt_equ import get_results
    get_results(input_path)
    
@analyze.command("admet-result")
@click.option("-i","--input_path",show_default=True,help="input path")
@click.option("-color","--color_shift",default=0,show_default=True,help="color")
def get_properties_value_anlayze(input_path,color_shift):
    """
    ADMET性质预测的结果生成
    -i: input directory\n
    """
    from craton.craton.property.expt_equ import get_admet_results
    color = int(color_shift)
    get_admet_results(input_path,color_shift=color)    
    
@analyze.command("property-figure")
@click.option("-i","--input_path",show_default=True,help="input path")
def get_properties_value_anlayze(input_path):
    """
    性质预测的结果重新生成图
    -i: block property file\n
    """
    from craton.craton.property.expt_equ import get_figure
    get_figure(input_path)    
    
@analyze.command("property-script")
@click.option("-i","--input_path",show_default=True,help="input path")
def get_properties_value_anlayze(input_path):
    """
    生成微调的bash文件
    -i: block property file\n
    -p: property file\n
    -o: output file\n
    """
    from craton.craton.property.expt_equ import get_fine_tune_script
    get_fine_tune_script(input_path)

@analyze.command("chem-space")
@click.option("-i","--input_file",default=".",show_default=True,help="input_file")
@click.option("-o", "--output_directory", default="./", show_default=True,help="The directory of output")
def chem_space_analyze(input_file,output_directory):
    """
    分析分子的化学空间分布情况
    -i: experimental file\n
    -o: output idrectory\n
    """
    MX.analyze_chem_space(input_file,output_dir=output_directory)

@analyze.command("train-test")
@click.option("-i","--input_file",default=".",show_default=True,help="input_file")
@click.option("-p","--property",default="density_of_liquid",show_default=True,help="property")
@click.option("-s","--style",default="random",show_default=True,help="split style")
@click.option("-v","--svalue",default=None,show_default=True,help="value for some style")
@click.option("-o", "--output_directory", default="./", show_default=True,help="The directory of output")
@click.option("-test", "--test_flag", default="./", show_default=True,help="run test .......")
def train_test_split(input_file,property,style,svalue,output_directory,test_flag):
    """
    得到训练集和测试集
    -i: experimental file\n
    -o: output idrectory\n
    -s: split style\n
        random: random split testing set, none value\n
        arom: the molecule has aromatic ring as testing set, none value\n
        ring: the molecule has ring as testing set, none value\n
        halogen: the molecule has halogen element as testing set, none value\n
        element: the molecule has  special element as testing set, default value: S\n
        zelement: the molecule is in special zelement as testing set, default value: ZCONSP\n 
        ha: the molecule with larger heavy atoms as testing set, default value: 20\n
    -v: value\n
        list example: S:Cl:P\n
    -p: property\n
        list example: density_of_liquid:vapor_pressure:heat_capacity_of_liquid\n
        below is aviaiblie property\n
        density_of_liquid:\n
    -t: test flag\n
        run test .......
    """
    if property is not None:
        if property != "all":
            property = property.split(":")
    if test_flag == "no":
        test_flag = False
    else:
        test_flag = True
    MX.split_train_test(input_file,property,style,value=svalue,output_dir=output_directory,test_flag=test_flag)
    
@analyze.command("property-add")
@click.option("-i","--input_file",help="add property file")
@click.option("-p","--property_file",help="property file")
@click.option("-o", "--output_file", default="./", show_default=True,help="The directory of output")
def add_property_to_file(input_file,property_file,output_file):
    """
    增加新的实验数据
    -i: add property file\n
    -p: property file\n
    -o: output file\n
    """
    from craton.craton.property.expt_equ import add_expt_data
    add_expt_data(input_file,property_file,output_file) 


@analyze.command("cc")
@click.option("-i", "--input", help="input json file", default="analyze_basic.json")
@click.option("-o", "--output", help="output json file", default="new_analyze_basic.json")
def recalculate_cycle_closure(input, output):
    with open(input, "r") as f:
        cc_data = simplejson.load(f)
    ddg_df = pd.DataFrame.from_dict(cc_data["edge"], orient="index")
    dg_df = pd.DataFrame.from_dict(cc_data["node"], orient="index")
    ddg_df["name"] = ddg_df.index
    dg_df["name"] = dg_df.index
    ddg_df, dg_df = AllPairBasicAnalyze.run_cycle_closure(ddg_df, dg_df)

    recalculate_result = {}
    recalculate_result["edge"] = ddg_df.to_dict(orient="index")
    recalculate_result["node"] = dg_df.to_dict(orient="index")
    with open(output, "w") as f:
        simplejson.dump(recalculate_result, f)


@analyze.command("all_check")
@click.option(
    "-i",
    "--input",
    help="The directory of rbfe pair, the pair has the name '*_to_*' ",
)
@click.option(
    "-o",
    "--output",
    help="output directory",
    default="fep_result",
    show_default=True,
)
@click.option("-n", "--n_jobs", help="Number of cpus to analyze joblib", default=1, type=int, show_default=True)
@click.option(
    "--restart/--no-restart", help="For failed pairs, write the restart shell script", default=False, show_default=False
)
def all_pair_check(input, output, n_jobs, restart):
    """Check FEP Simulation status and print the failed pairs"""

    # for front-end
    check_result = {}
    output_dir = Path(output).resolve()
    if not output_dir.exists():
        output_dir.mkdir(exist_ok=True, parents=True)
    fep_check = CheckFEPRunningStatus(Path(input), n_jobs)
    fep_check.simulation_time.to_csv(str(output_dir / "simulation_time.csv"))

    if not len(fep_check.failed_pairs):
        check_result["status"] = "success"
        click.echo("Congratulation! All pairs are finished successfully!")
    else:
        if not len(fep_check.success_pairs):
            check_result["status"] = "failed"
        else:
            check_result["status"] = "partial success"

    for pair, error_stage in fep_check.error_messages.items():
        click.secho(pair, fg="green")
        for stage, error in error_stage.items():
            click.secho(f"  {stage}", fg="yellow")
            click.secho(f'{textwrap.indent(error, prefix="     ")}', fg="red")
        click.echo("\n")

    check_result["detail"] = fep_check.error_messages
    with open(output_dir / "check.json", "w") as f:
        simplejson.dump(check_result, f)

    if restart and len(fep_check.failed_pairs):
        restart_path = []
        for pair, stage_error in fep_check.error_messages.items():
            for stage in stage_error.keys():
                restart_path.append(str((fep_check._pair_directory / pair / stage).resolve()))
        write_restart_script(restart_path, output_dir / "restart.sh")

