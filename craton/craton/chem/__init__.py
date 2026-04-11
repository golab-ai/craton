import os
from pathlib import Path
import tarfile
import zipfile


from .molecule import Molecule
from .atom import Atom
from .topology import Bond,Angle,Dihedral,Improper,Pair,Constrain
from ..utils.common.utils import parse_inchi_key

from .format.mol_parse import MolData, SdfData
from .format.mol2_parse import Mol2Data
from .format.db_parse import DBData
from .format.pdb_parse import PdbData
from .format.pdbqt_parse import PdbqtData
from .format.smiles_parse import SmilesData
from .format.xyz_parse import XyzData
from .format.csv_parse import CsvData
from .format.cif_parse import CifData
from .format.mtx_parse import MtxData
from ..software.gaussian import GauInputFile, GauOutputFile
from .format._rdkit import RdkitMol
from ..utils import logger
from ..utils.commons import parallel_run

from .image import MoleculeImage
from PIL import Image

def _show_molecule(molecule,attrs=None,extra=None,save_file=True,opath=None,TD_flag=False,remove_H_flag=True):
    molshow = MoleculeImage(molecule,target=attrs,extra=extra,save_file=save_file,fpath=opath,TD_flag=TD_flag,remove_H_flag=remove_H_flag)
    return molshow

def create_3D(molecule):
    molobj = MolData("normal")
    script = molobj._convert(molecule)
    rdkobj = RdkitMol()
    rdkobj._convert(script,normalization=False)
    rdkobj._get_3d()
    script_3d = rdkobj._get_script()
    datas =  molobj._parse(script_3d)
    for atom in molecule.Atoms:
        atom.coor = datas["coordinates"][atom.ID]

def pdbqt_file(molecule,extra_var=None):
    pdbqt = PdbqtData()
    return pdbqt.convert_info(molecule,extra_var=extra_var)

class FormatMolecule():
    __TYPE_FUNC = {
        ".mol":  MolData,
        ".mol2": Mol2Data,
        ".pdb":  PdbData,  # TODO, 很不完善
        ".pdbqt": PdbqtData,
        ".sdf":  SdfData,
        ".txt":  CsvData,
        ".csv":  CsvData,
        ".gjf":  GauInputFile,
        ".log":  GauOutputFile,
        ".xyz":  XyzData, # TODO, _parse不完善
        ".inchikey":DBData,
        ".mtx":  MtxData,
        "smiles": SmilesData,
        "inchi_key": DBData,
        ".fromdb":DBData,
        }
    def __init__(self):
        pass

    @staticmethod
    def _parse(input_files,extra_var=None,parallel=True):
        """
        extra_var:
            DBMol: conformer_selector,
            Smiles: create3d
            CsvFile: DBMol + CallRdkit
        """
        msg = ""
        if not isinstance(input_files,list):
            input_files = [input_files]

        molecules = []
        real_input_file = []
        for input_file in input_files:
            if input_file is None:
                return FormatMolecule._parse_single(None,"inchi_key",extra_var=extra_var)
            if isinstance(input_file,Molecule):
                molecules.append(input_file)
                continue
            if isinstance(input_file,dict):
                molecules.append(FormatMolecule._create_molecule(input_file)[0])
                continue
            

            if input_file.endswith(".zip"):
                this_scripts = FormatMolecule._read_tar_file(input_file,extra_var=extra_var,zip_file=True)
                real_input_file.extend(this_scripts)
            elif input_file.endswith(".tar.gz") or input_file.endswith(".tgz"):
                this_scripts = FormatMolecule._read_tar_file(input_file,extra_var=extra_var)
                real_input_file.extend(this_scripts)
            else:
                #datas = FormatMolecule._run_parse(input_file,extra_var=extra_var)
                real_input_file.append([input_file,{"extra_var":extra_var}])
        if len(real_input_file) > 0:
            if parallel:
                #total_datas = parallel_run(FormatMolecule._run_parse,real_input_file)
                total_datas = parallel_run(FormatMolecule._run_parse,[rr[0] for rr in real_input_file],kwds=[rr[1] for rr in real_input_file])
            else:
                total_datas = []
                for rr in real_input_file:
                    total_datas.append(FormatMolecule._run_parse(rr[0],**rr[1]))
            for datas in total_datas:
                if datas[1] is not None:
                    msg += f"{datas[1]}\n"
                else:
                    if isinstance(datas[0],list):
                        molecules.extend(datas[0])
                    else:
                        molecules.append(datas[0])

        if msg != "":
            logger.warning(msg)

        return molecules

    @staticmethod
    def _read_tar_file(input_file,extra_var=None,zip_file=False):
        def _get_molecule(script,extra_var,suffix):
            this_molecules = []
            datas = FormatMolecule._run_parse(script,extra_var=extra_var,suffix=suffix)
            if datas[1] != "":
                msg += f"{datas[1]}\n"
            else:
                if isinstance(datas[0],list):
                    this_molecules.extend(datas[0])
                else:
                    this_molecules.append(datas[0])

        #msg = ""
        #this_molecules = []

        tar_scripts = []
        if zip_file:
            zip = zipfile.ZipFile(input_file, "r")
            files = [name for name in zip.namelist() if name.endswith(".log")]
            for f in files:
                suffix = Path(f).suffix
                script = zip.open(f).read().decode()
                tar_scripts.append([script,{"extra_var":extra_var,"suffix":suffix,"script_filename":f,"parent_input_file":input_file}])
                #_get_molecule(script,extra_var,suffix)
        else:
            tar = tarfile.open(input_file, "r:gz")
            files = [m.name for m in tar.getmembers()]
            for f in files:
                suffix =  Path(f).suffix
                script = tar.extractfile(f).read().decode()
                tar_scripts.append([script,{"extra_var":extra_var,"suffix":suffix,"script_filename":f,"parent_input_file":input_file}])
                #_get_molecule(script,extra_var,suffix)
        return tar_scripts
        #return this_molecules, msg

    @staticmethod
    def _run_parse(input_file,extra_var=None,suffix=None,script_filename=None,parent_input_file=None,idx=None):
        datas = []
        msg = None
        __error_flag = False
        

        if os.path.isfile(input_file):
            if suffix is None:
                suffix = Path(input_file).suffix
            if suffix in FormatMolecule.__TYPE_FUNC:
                with open(input_file) as inf:
                    input_script = inf.read()
                #try:
                datas = FormatMolecule._parse_single(input_script,suffix,extra_var=extra_var)
                #except:
                #    __error_flag = True
            else:
                __error_flag = True
            if __error_flag:
                msg = f"create molelcule error:  {suffix if suffix is not None else 'Unknow file'} type -> {input_file}"
            if idx is not None:
                return [datas,msg],idx
            else:
                return [datas,msg]
            
        if suffix is not None:
            try:
                datas = FormatMolecule._parse_single(input_file,suffix,extra_var=extra_var)
            except:
                __error_flag = True 
        else:
            try:
                datas =  FormatMolecule._parse_single_try(input_file,extra_var=extra_var)
            except:
                __error_flag = True

        if __error_flag:
            msg = f"create molecule error:"
            msg += suffix if suffix is not None else 'Unknow file'
            msg += f" type -> {input_file}"
            if script_filename is not None:
                msg += f"in {script_filename}"
            if parent_input_file is not None:
                msg += f"in {parent_input_file}"
        if idx is not None:
            return [datas,msg], idx
        else:
            return [datas,msg]

    @staticmethod
    def _parse_single(input_script,suffix,extra_var=None):
        parser = FormatMolecule.__TYPE_FUNC[suffix]()
        datas = parser._parse(input_script,extra_var=extra_var)
        smiles_flag = True if extra_var is None else True if "smiles_flag" not in extra_var else extra_var["smiles_flag"]
        if suffix not in [".mtx"]:
            return FormatMolecule._create_molecule(datas,smiles_flag=smiles_flag)
        else:
            return FormatMolecule._set_molecule(datas)

    @staticmethod
    def _parse_single_try(input_file,extra_var):
        inchikey =  parse_inchi_key(input_file)
        if inchikey is not None:
            parser = DBData()
            data = parser._parse(input_file,extra_var=extra_var)
            return FormatMolecule._create_molecule(data)
        try:
            parser = SmilesData()
            data = parser._parse(input_file,extra_var=extra_var)
            return FormatMolecule._create_molecule(data)
        except:
            
            for __, _parser in FormatMolecule.__TYPE_FUNC.items():
                parser = _parser()
                try:
                    data = parser._parse(input_file,extra_var=extra_var)
                    return FormatMolecule._create_molecule(data)
                except:
                    pass
            raise KeyError("Cannot read input file %s" %input_file)

    @staticmethod
    def _parse_single_old(input_file,extra_var=None):
        suffix = Path(input_file).suffix
        if suffix in FormatMolecule.__TYPE_FUNC:
            parser = FormatMolecule.__TYPE_FUNC[suffix]()
            with open(input_file) as inf:
                input_script = inf.read()
            datas = parser._parse(input_script,extra_var=extra_var)
            return FormatMolecule._create_molecule(datas)
        
        inchikey =  parse_inchi_key(input_file)
        if inchikey is not None:
            parser = DBData()
            data = parser._parse(input_file,extra_var=extra_var)
            return FormatMolecule._create_molecule(data)
        
        try:
            parser = SmilesData()
            data = parser._parse(input_file,extra_var=extra_var)
            return FormatMolecule._create_molecule(data)
        except:
            try:
                parser = SdfData()
                data = parser._parse(input_file,extra_var=extra_var)
                return FormatMolecule._create_molecule(data)
            except:
                try:
                    parser = MolData()
                    data = parser._parse(input_file,extra_var=extra_var)
                    return FormatMolecule._create_molecule(data)
                except:
                    try:
                        parser = Mol2Data()
                        data = parser._parse(input_file,extra_var=extra_var)
                        return FormatMolecule._create_molecule(data)
                    except:
                        raise KeyError("Cannot read input file %s" %input_file)

    @staticmethod
    def _set_molecule(datas):
        __label = {"Atoms":Atom,"Bonds":Bond,"Angles":Angle,"Dihedrals":Dihedral,
                   "Impropers":Improper,"Pair12":Pair,"Pair13":Pair,"Pair14":Pair,"Pair1n":Pair,
                   "constrain":Constrain}
        __atoms = ["a1","a2","a3","a4"]
        if not isinstance(datas,list):
            datas = [datas]
        molecules = []
        for data in datas:
            molecule=Molecule("molecule")
            if "Molecule" in data:
                for aa,bb in data["Molecule"].items():
                    if aa not in ["formula","torsion_number","net_charge","heavy_atoms"]:
                        setattr(molecule,aa,bb)
            for attr in __label:
                if attr in data:
                    _tmp = []
                    if attr == "Atoms":
                        for rr in data[attr]:
                            _tmp.append(Atom("atom"))
                            for aa,bb in rr.items():
                                setattr(_tmp[-1],aa,bb)
                            if hasattr(_tmp[-1],"parameter"):
                                _ff_parameter = {tt:getattr(_tmp[-1],tt,None) for tt in ["pstyle","mass","fix_parameter","parameter","ptag","pscore","pcount"]}
                                #####_ff_parameter["parameter"] = getattr(_tmp[-1],"vdw_parameter",None)
                                _tmp[-1]._ff_parameter = _ff_parameter
                    elif attr == "constrain":
                        for rr in data[attr]:
                            atoms = [rr[an] for an in __atoms if an in rr]
                            _tmp.append(__label[attr](atoms,rr["fix_value"]))
                            for aa,bb in rr.items():
                                
                                if aa not in __atoms + ["fix_value"]:
                                    setattr(_tmp[-1],aa,bb)
                    else:
                        for rr in data[attr]:
                            atoms = [rr[an] for an in __atoms if an in rr]
                            _tmp.append(__label[attr]("normal",*atoms))
                            for aa,bb in rr.items():
                                if aa not in __atoms:
                                    setattr(_tmp[-1],aa,bb)
                            if hasattr(_tmp[-1],"parameter"):
                                _ff_parameter = {tt:getattr(_tmp[-1],tt,None) for tt in ["pstyle","fix_parameter","parameter","ptag","pscore","pcount"]}
                                _tmp[-1]._ff_parameter = _ff_parameter
                    setattr(molecule,attr,_tmp)
            molecules.append(molecule)
        return molecules

    @staticmethod
    def _create_molecule(datas,smiles_flag=True):
        if not isinstance(datas,list):
            datas = [datas]
        molecules  = []
        for data in datas:
            molecules.append(Molecule("molecule"))
            molecules[-1].get_mole_info(data)
            molecules[-1].get_atoms_info(data)
            if not hasattr(molecules[-1].Atoms[0],"connectivity"):
                pass
                #logger.warning("the molecule has not connectivity")
            else:
                molecules[-1].create_topols(smiles_flag=smiles_flag)
                if molecules[-1].style not in ["molecule"]:
                    molecules[-1].create_intra_nonbond_macromole()
                else:
                    molecules[-1].create_intra_nonbond()
                molecules[-1].steps.append("topol")
        return molecules

    @staticmethod
    def _check_molecule(molecules):
        for molecule in molecules:
            if hasattr(molecules[-1].Atoms[0],"bond_type"):
                molecules[-1].append("bond_type")
            
    @staticmethod
    def _write_file(texts,molecule_names,otype,ofilename,opath):
        if opath is None:
            opath = "./"
        if not isinstance(texts,list) or len(texts) == 1:
            if ofilename is None:
                ofilename = "compound"
            fn = f"{opath}/{ofilename}{otype}"
            if len(texts) == 1:
                script = texts[0]
            else:
                script = texts
            with open(fn,'w') as outf:
                outf.write(script)
        else:
            for ii,text in enumerate(texts):
                fname = molecule_names[ii] if ofilename is None else ofilename
                fn = f"{opath}/{fname}-{str(ii)}{otype}"
                with open(fn,'w') as outf:
                    outf.write(text)

    @staticmethod
    def _convert(molecules,otype=None,ofilename=None,opath=None,extra_var=None,parallel=True):
        """
        csv: atrributions,list
        mol: has3d,yes or no; vsflag,bool
        mol2: has3d,yes or no; vsflag,bool
        sdf: elem_atomtype,bool; real_bond_type, bool; vsflag,bool
        """
        if not isinstance(molecules,list):
            molecules = [molecules]
        if otype in ["png", "svg", "figure"]:
            from .image import MoleculeImage
            for molecule in molecules:
                MoleculeImage(molecule,save_file=True,fpath=opath)
            return
        
        if otype in ["db"]:
            converter = DBData()
            converter._convert(molecules,extra_var=extra_var)
            return

        otype = ".csv" if otype is None or otype in ["inchi_key","smiles"] else f".{otype}"
        converter = FormatMolecule.__TYPE_FUNC[otype]()
        if otype in [".csv",".db",".sdf",]:
            texts = converter._convert(molecules,extra_var=extra_var)
        else:
            texts = []
            for molecule in molecules:
                texts.append(converter._convert(molecule,extra_var=extra_var))
        molecule_names = []
        for molecule in molecules:
            fname = molecule.mole_name
            if hasattr(molecule,"constrain"):
                for cons in molecule.constrain:
                    fname += f"_{cons.name}-{round(cons.fix_value,2)}"
            molecule_names.append(fname)
        FormatMolecule._write_file(texts,molecule_names,otype,ofilename,opath)
        return texts
