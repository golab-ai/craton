from copy import deepcopy
from multiprocessing import Pool
import psutil

def parallel_run(func,args,kwds=None,objs=None,single_args_flag=True,keep_order=True,return_result=True):
    if objs is None:
        iterable_para = args
    else:
        if single_args_flag:
            iterable_para = objs
        else:
            iterable_para = args

    #ncore = min(psutil.cpu_count(logical=False) - 1, len(iterable_para))
    ncore = min(psutil.cpu_count(logical=True) - 2, len(iterable_para))
    pool = Pool(processes=ncore)
    jobs = []

    if objs is None:
        for ii,arg in enumerate(args):
            if kwds is not None:
                if isinstance(kwds,list):
                    kwd = kwds[ii]
                elif isinstance(kwds,dict):
                    kwd = deepcopy(kwds)
            else:
                kwd = None

            if keep_order:
                if kwd is None:
                    kwd = {"idx":ii}
                else:
                    kwd["idx"] = ii


                #jobs.append(pool.apply_async(func,args=[arg],kwds=kwd))
            if kwd is None:
                jobs.append(pool.apply_async(func,args=[arg]))
            else:
                jobs.append(pool.apply_async(func,args=[arg],kwds=kwd))
    else:
        for ii,__ in enumerate(iterable_para):
            if kwds is not None:
                if isinstance(kwds,list):
                    kwd = kwds[ii]
                elif isinstance(kwds,dict):
                    kwd = deepcopy(kwds)
            else:
                kwd = None

            if keep_order:
                if kwd is None:
                    kwd = {"idx":ii}
                else:
                    kwd["idx"] = ii

            if single_args_flag:
                arg = args
                this_func = getattr(objs[ii],func)
            else:
                arg = args[ii]
                this_func = getattr(objs,func)

            if kwd is None:
                jobs.append(pool.apply_async(this_func,args=[arg]))
            else:
                jobs.append(pool.apply_async(this_func,args=[arg],kwds=kwd))


    pool.close()
    pool.join()
    if return_result:
        if keep_order:
            results = [None] * len(args)
            for job in jobs:  
                data,idx = job.get()
                #results[idx] = (data,idx)
                results[idx] = data
        else:
            results = []
            for job in jobs:
                results.append(job.get()) 
        return results