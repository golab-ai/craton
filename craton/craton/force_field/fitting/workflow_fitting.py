import json
import sys
from copy import deepcopy

import yaml


class FittingWorkFlow:
    def __init__(self, settings, write_db_flag=False, ctx=None):
        #with open(setting_path) as f:
        #    self.setting = yaml.safe_load(f)
        self.setting = settings

        if ctx is None:
            self.ctx = MongoDB()
        self.write_db_flag = write_db_flag
        self.parser_config()

        # self.tasks = list(self.path_dict["tasks"].items())
        for ii, rr in enumerate(self.tasks):
            name, txt = rr
            if not Path(txt).exists():
                logger.error(f"Mole file not exist: {name}")
                sys.exit(1)
            # if Path(os.path.join(self.fitting_setting.output_dir, name)).exists():
            #    logger.error(f"Job dir already exists: {name}")
            #    sys.exit(1)

        self.last_output_path = None
        Path(self.config["output_dir"]).mkdir(parents=True, exist_ok=True)

        self.run()

    def parser_config(self):
        self.force_field_setting: ForceFieldSetting = ForceFieldSetting()
        self.force_field_setting.update(self.setting)
        self.force_field_setting = path_from_config(self.force_field_setting)

        self.fitting_setting: ForceFieldFitting = ForceFieldFitting()
        self.fitting_setting.update(self.setting)
        self.fitting_setting = path_from_config(self.fitting_setting)

        with open(self.fitting_setting.tasks_file) as inf:
            tasks_arr = [dd for dd in inf.readlines() if dd != "\n"]

        self.tasks = []
        for i, txt in enumerate(tasks_arr):
            path = os.path.join(self.fitting_setting.input_dir, txt.strip())
            filename = os.path.basename(path)
            name = "%02i.%s" % (i + 1, os.path.splitext(filename)[0])
            self.tasks.append([name, path])
        self.config = {kk: vv for kk, vv in self.force_field_setting.__dict__.items()}
        self.config.update({kk: vv for kk, vv in self.fitting_setting.__dict__.items()})

    def run(self):
        if self.config["ff_only"]:
            self.save_init_ff()
            return
        self.config["parent_output_dir"] = self.config["output_dir"]
        self.config["origin_force_field_file"] = self.config["force_field_file"]
        self.config["origin_using_force_field_file"] = self.config["using_force_field_file"]
        for i in range(len(self.tasks)):
            this_config = self.run_single_task(i)
            if self.write_db_flag:
                self.write_into_db(i, this_config)

    def run_fitting(self, configure, inchi_key_arr, ff_only=False):
        fitting = Fitting(configure=configure)
        fitting.load_training_set_from_db(inchi_key_arr, self.ctx)
        if len(fitting.training_set) == 0:
            return
        fitting.prepare()
        if not ff_only:
            fitting.run()

    def write_into_db(self, i_task, config):
        desc = deepcopy(config)
        name, txt = self.tasks[i_task]
        desc["job_name"] = name
        if desc["only_validation_flag"]:
            if desc["only_gaff"]:
                desc["job_type"] = "gaff validation"
            else:
                desc["job_type"] = "validation"
        else:
            if desc["charge_method"] is None:
                desc["job_type"] = "binc + intra"
            else:
                desc["job_type"] = "intra"

        result_file = os.path.join(desc["output_dir"], "fitting_results.json.gz")
        desc["result_file"] = result_file

        if desc["only_validation_flag"]:
            output_ff_file = None
        else:
            output_ff_file = os.path.join(desc["output_dir"], "total_ff.ff")
        desc["output_ff_file"] = output_ff_file
        logger.info("Writing fitting results to DB ...")
        with open(f"{name}.json", "w") as outf:
            outf.write(json.dumps(desc))

    def run_single_task(self, i_task):
        config = deepcopy(self.config)

        name, txt = self.tasks[i_task]
        logger.info(f"===== Fitting task {i_task + 1}/{len(self.tasks)}: {txt} =====")

        config["output_dir"] = os.path.join(config["parent_output_dir"], name)
        Path(config["output_dir"]).mkdir(parents=True, exist_ok=True)

        inchi_key_arr = [m["inchi_key"] for m in read_molecule_csv_file(txt)]

        if self.last_output_path and not config["only_validation_flag"]:
            config["force_field_file"] = os.path.join(self.last_output_path, "total_ff.ff")
            config["using_force_field_file"] = None

        self.last_output_path = config["output_dir"]

        self.run_fitting(config, inchi_key_arr)
        return config

    def save_init_ff(self):
        inchi_key_arr = []
        for rr in self.tasks:
            __, txt = rr
            inchi_key_arr.extend([m["inchi_key"] for m in read_molecule_csv_file(txt)])
        self.run_fitting(self.config, inchi_key_arr, ff_only=True)


if __name__ == "__main__":

    #try:
    #    with open(sys.argv[1]) as f:
    #        settings = yaml.safe_load(f)
    #except:
    #    with open("config.yml") as f:
    #        settings = yaml.safe_load(f)
    
    #fwk = FittingWorkFlow(settings)
    parser ={"input_setting_yaml":sys.argv[1],
             }
    if len(sys.argv) > 2:
        parser["job_list"] = sys.argv[2]
    else:
        parser["job_list"] = None
    with open(parser["input_setting_yaml"]) as f:
        settings = yaml.safe_load(f)
    if parser["job_list"] is None:
        fwk = FittingWorkFlow(settings)
    else:
        with open(parser["job_list"]) as f:
            joblist = [p.strip() for p in f.readlines()]
        for p in joblist:
            this_settings = {aa:bb for aa,bb in settings.items() if aa not in ["input_dir","tasks_file","output_dir"]}
            this_settings["input_dir"] = p
            this_settings["tasks_file"] = f"{p}/list.log"
            this_settings["output_dir"] = settings["output_dir"]
            for dd in  p.split("/"):
                if dd not in [".","input"]:
                    this_settings["output_dir"] += f"/{dd}"
            fwk = FittingWorkFlow(this_settings)
        