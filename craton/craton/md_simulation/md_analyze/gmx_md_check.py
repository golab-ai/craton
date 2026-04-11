import datetime
import os
#from collections import defaultdict
from datetime import timedelta
#from functools import cached_property
#from pathlib import Path
#from typing import Dict, List, Set, Tuple

#import pandas as pd
#from joblib import Parallel, delayed

from ...utils import logger
from ...utils.commons import parallel_run

class CheckGmxRun:
    def __init__(self,directory,batchfile=None,parallel=True):
        self.directory = directory
        self.batchfile = batchfile
        self.parallel = parallel

    @staticmethod
    def _get_info_from_log_file(filename):
        if not os.path.exists(filename[:-4] + ".tpr"):
            return False, False,False,False

        if not os.path.exists(filename):
            return True, False,False,False
        f_create_time = os.path.getctime(filename[:-4] + ".tpr")
        f_modify_time = os.path.getmtime(filename)
        with open(filename) as f:
            lines = f.readlines()
            run_done = True if lines[-2].startswith("Finished") else False
            done_err = True if lines[-1].startswith("------------------") else False
        now_step  = 1
        if done_err:
            for ii,line in enumerate(lines):
                if line.startswith("Fatal error"):
                    nn = ii + 1
                    break
            error_info = ""
            for line in lines[ii:]:
                if line == "":
                    break
                else:
                    line += "%s\n"%line 
            return True, True, False, False, error_info
        else:
            for kk,line in enumerate(lines):
                if line.startswith("Started"):
                    start_time = " ".join(line.split()[-5:])
                    start_datetime = datetime.datetime.strptime(start_time, "%a %b %d %H:%M:%S %Y")
                elif line.startswith("Command line"):
                    ss = lines[kk+1].split()
                    try:
                        ncore = int(ss[ss.index("-ntomp") + 1])
                    except:
                        ncore = 1
                elif line.startswith("   nsteps"):
                    nsteps = float(line.split("=")[1].strip())
                elif line.startswith("   dt"):
                    dt = float(line.split("=")[1].strip())
                elif line.startswith("           Step"):
                    now_step = int(lines[kk + 1].split()[0])
                elif line.startswith("Finished"):
                    end_time = " ".join(line.split()[-5:])
                    end_datetime = datetime.datetime.strptime(end_time, "%a %b %d %H:%M:%S %Y")
                    try:
                        performance = [float(lines[kk-1].split()[1]),float(lines[kk-1].split()[2])]
                    except:
                        performance = [None,None]
                    run_time = lines[kk-3].split()[0]
                    try:
                        cpu_time  = [float(lines[kk-4].split()[1]),float(lines[kk-4].split()[2])]                
                    except:
                        try:
                            cpu_time  = [float(lines[kk-3].split()[1]),float(lines[kk-3].split()[2])]
                        except:
                            cpu_time = [None,None]
        if run_done:
            return True, True, True, True,  nsteps*dt/1000 , nsteps*dt/1000, ncore, cpu_time,run_time, performance,end_datetime - start_datetime + timedelta(0.001)
        else:
            _time_ = f_modify_time - f_create_time
            per_time = (_time_ / now_step)
            res_step = nsteps - now_step
            res_time = per_time * res_step
            
            return True, True, True, False, nsteps*dt/1000, now_step*dt/1000, ncore,  nsteps*dt/1000 - now_step*dt/1000, res_time , per_time

    def _get_tasks(self,this_path):
        if os.path.exists(f"{this_path}/job.sh"):
            self.tasks.append(this_path)
        dds = [dd for dd in os.listdir(this_path) if os.path.isdir(f"{this_path}/{dd}")]
        for dd in dds:
            if dd.find("_to_") != -1:
                self.tasks.append(f"{this_path}/{dd}/rbfe")
                self.tasks.append(f"{this_path}/{dd}/rhfe")
            else:
                self._get_tasks(f"{this_path}/{dd}")


    def _get_info_job(self,filename):
        jobs = []
        with open(filename) as inf:
            for line in inf:
                if line.find("mdrun") != -1:
                    ss = line.split()
                    jobs.append(ss[ss.index("-deffnm") + 1])
        return jobs

    def _find_jobID_from_batch(self,task,batch_arr):
        for ii,bb in enumerate(batch_arr):
            if bb.find(task) != -1:
                return (ii)

    @staticmethod
    def _get_single_task_info(task,jobs=None,batch_arr=None):
        def get_single_job(task,jobs,jobID):
            done_job = []
            for job in jobs:
                status = CheckGmxRun._get_info_from_log_file(f"{task}/{job}.log")
                textt = ""
                if len(done_job) > 0:
                    textt = "  ."
                    textt += ",".join(done_job)
                    textt += "完成"
                
                if not status[0]:    
                    return "error_jobs", f"{jobID}{task}: 生成{job}.tpr错误{textt}",task
                if not status[1]:
                    return "error_jobs", f"{jobID}{task}: {job} 没有启动{textt}" ,task
                if not status[2]:
                    return "error_jobs", f"{jobID}{task}: {job} 错误{textt}\n 错误类型: {status[-1]}",task
                if not status[3]:
                    return "run_jobs",f"{jobID}{task}: {job} 计算中,共{status[4]} ns, 完成{status[5]} ns, 剩下{status[7]}, 预计还需要{status[8]}{textt}",task
                else:
                    done_job.append(job)
                    if len(done_job) == len(jobs):
                        return "done_jobs", f"{jobID}{task}: 计算完成，共{status[4]} ns, 运行时间{status[7][1]}, 费用(CFFF): {0.03*status[7][0]/3600} ¥。GMX效率 {status[9][0]} ns/day",task
        
        jobID = ""
        if batch_arr is not None:
            jobID = f"{job}{self._find_jobID_from_batch(task,batch_arr)}: "

        error_file = [f for f in os.listdir(task) if f[-4:] == ".err"]
        if len(error_file) == 0:
            return "wait_jobs", "{jobID}{task} 等待资源",task
        else:
            if task.find("_to_") != -1:
                dds = [dd for dd in os.listdir(f"{task}") if os.path.isdir(f"{task}/{dd}")]
                status = [get_single_job(f"{task}/{dd}",jobs,jobID) for dd in dds]
                this_infos = {"wait_jobs":[],"error_jobs":[],"run_jobs":[],"done_jobs":[]}
                for rr in status:
                    this_infos[rr[0]].append([rr[1],rr[2]])
                if len(this_infos["error_jobs"]) > 0:
                    return "error_jobs", this_infos["error_jobs"][0][0],task
                if len(this_infos["run_jobs"]) > 0:
                    return "run_jobs", this_infos["run_jobs"][0][0],task
                return "done_jobs",this_infos["done_jobs"][0][0],task
            else:
                return get_single_job(task,jobs,jobID)
                    
    def _get_info(self):
        self.parallel = False
        batch_arr = None
        if self.batchfile is not None:
            batch_arr = []
            with open(self.batchfile) as inf:
                for line in inf:
                    batch_arr.append(line.split()[-1].strip())
        self.tasks = []
        self._get_tasks(self.directory)
        jobs_arr = []
        for task in self.tasks:
            jobs_arr.append(self._get_info_job(f"{task}/job.sh"))
    
        infos = {"wait_jobs":[],"error_jobs":[],"run_jobs":[],"done_jobs":[]}

        if self.parallel:
            total_status = parallel_run(CheckGmxRun._get_single_task_info,
                                    self.tasks,
                                     kwds=[{"jobs":jobs_arr[ii],"batch_arr":batch_arr} for ii in range(len(self.tasks))],
                                     keep_order=False)

        else:
            total_status = []
            for task,jobs in zip(self.tasks,jobs_arr):
                total_status.append(self._get_single_task_info(task,jobs,batch_arr))

        for rr in total_status:
            infos[rr[0]].append([rr[1],rr[2]])
        return infos
