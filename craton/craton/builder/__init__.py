from .builder_manager import BuildManager
from pathlib import Path
from copy import deepcopy

#def builder_manager(molecules,config=None,style=None,parallel=True):
#    BM = BuilderManager(molecules,config=config,style=style,parallel=parallel)
#    sms = BM()
#    return sms

def builder_manager(molecules,config=None,style=None,parallel=True):
    return BuildManager(molecules,config,parallel=parallel)


def make_propety_system(sms):
        bulk_npt_properties = ["den","hov","cp","kt","ap","rdf","rg"]
        bulk_nvt_properties = ["dc","vis","td","er"]
        bulk_properties = bulk_npt_properties + bulk_nvt_properties
        slab_properties = ["st","vle","ct","cd","cb","svp","nbp"]
        tmp = []
        for sm in sms:
            this_property = sm.md_setting["property"]
            bulk_nvt_pp = list(set(this_property).intersection(set(bulk_nvt_properties)))
            slab_pp = list(set(this_property).intersection(set(slab_properties)))
            bulk_npt_pp = list(set(this_property).difference(set(bulk_nvt_pp+slab_pp)))
            _or_ = [len(bulk_npt_pp)>0,len(bulk_nvt_pp)>0,len(slab_pp)>0]
            use_this_sm = True
            change_this_sm_name = True
            if sum([1 for rr in _or_ if rr]) == 1:
                change_this_sm_name = False
            output_dir = sm.output_dir
            sm.info_dir = Path(sm.output_dir)
            for pp,dirname,flag in zip([bulk_npt_pp,bulk_nvt_pp,slab_pp],["_bulk_npt","_bulk_nvt","_slab"],_or_):
                if flag:
                    if use_this_sm:
                        this_sm = sm
                        use_this_sm = False
                        if change_this_sm_name:
                            this_sm.output_dir = f"{output_dir}/{dirname}"
                    else:
                        tmp.append(deepcopy(sm))
                        this_sm = tmp[-1]
                        this_sm.output_dir = f"{output_dir}/{dirname}"
                    this_sm.md_setting["property"] = pp
                    if dirname != "_bulk_npt":
                        if dirname == "_slab":
                            this_sm.lattics[-1] += 100.0
                        this_sm.md_setting["pressure"]["pressure_coupl"] = ["no" for __ in this_sm.md_setting["jobs"]]
                        _tmp_job = []
                        for job in this_sm.md_setting["jobs"]:
                            #if job == "eq_npt":
                            #    _tmp_job.append("eq_nvt")
                            # elif job == "prod_npt":
                            if job == "prod_npt":
                                _tmp_job.append("prod_nvt")
                            else:
                                _tmp_job.append(job)
                        this_sm.md_setting["jobs"] = _tmp_job
        sms.extend(tmp)
        return sms