import os
from pathlib import Path

from .qm_calculator import parse_qm_setting, QmCalculator

from ..utils.commons import parallel_run

class QMCalc:
    def __init__(self):
        pass

    @staticmethod
    def qm_input_file(
        molecules,
        qmpara=None,
        step=None,
        local_path="./",
        indexs=None,
        zmatrixs=None,
        fpath_pre="",
        parallel=True,
        ):
        """
        Generate for multi gjf file.
        """
        if not isinstance(molecules,list):
            molecules = [molecules]
        if indexs is not None:
            if not isinstance(indexs,list):
                indexs = [indexs]
        if zmatrixs is not None:
            if not isinstance(zmatrixs,list):
                zmatrixs = [zmatrixs]
        #if qm_setting_file is not None:
        #    qmpara = parse_qm_setting(qm_setting_file)


        this_path = os.path.join(local_path)
        Path(this_path).mkdir(parents=True, exist_ok=True)
        qmcal = QmCalculator(path=this_path)
        if parallel:
            kwds = [{"qmpara":qmpara,"step":step,"fpath_pre":fpath_pre} for _ in molecules]
            if indexs is not None:
                for ii,idx in enumerate(indexs):
                    kwds[ii]["index"] = idx
            if zmatrixs is not None:
                for ii,zmat in enumerate(zmatrixs):
                    kwds[ii]["zmatrix"] = zmat
            
            parallel_run('create_qm_input_files',molecules,kwds=kwds,objs=qmcal,single_args_flag=False,keep_order=False,return_result=False)
        else:
            for ii,molecule in enumerate(molecules):
                index = indexs[ii] if indexs is not None else None
                zmatrix = zmatrixs[ii] if zmatrixs is not None else None
                qmcal.create_qm_input_files(molecule, qmpara=qmpara, step=step, index=index, zmatrix=zmatrix,fpath_pre=fpath_pre)

