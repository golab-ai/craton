import os
import tarfile
import zipfile
import copy
import math
import datetime
import numpy as np
from pathlib import Path
from copy import deepcopy

CRATON_DIR = Path(__file__).parent
ROOT_DIR = CRATON_DIR.parent
DATABASE_EXISTS = Path(f"{ROOT_DIR}/database").exists()
if DATABASE_EXISTS:
    from ..database.mongodb import MongoDB

    qmcalc_coll = MongoDB().compound_qmcalc_coll

TASK_MAX_MOLs = 5000

class2=[
        ["first_train_set","Q0"],
        ["first_train_set","Q1"],
        ["second_train_set","Q0"],
        ["second_train_set","Q1"],
        ["first_train_set","Q4"],
        ["first_test_set","Q0"],
        ["first_test_set","Q1"],
        ['second_test_set',"Q0"],
        ['second_test_set',"Q1"],
        ["second_train_set","Q4"],
        ['third_test_set',"Q0"],
        ['third_test_set',"Q1"],
        ['forth_test_set',"Q0"],
        ['forth_test_set',"Q1"],
        ]
class3 = [
            'ZC','ZCO_chain','ZCO_ring','ZCN_chain', 'ZCN_ring','ZCP', 'ZCS', 
            'ZCON_chain', 'ZCON_ring', 'ZCOS_chain', 'ZCOS_ring', 'ZCOS',
            'ZCOP', 'ZCOP_chain', 'ZCOP_ring','ZCNS_chain', 'ZCNS_ring',
            'ZCNP', 'ZCNS','ZCSP',
            'ZCONS_chain', 'ZCONS_ring', 'ZCONS','ZCONP', 'ZCONP_chain', 'ZCONP_ring', 'ZCOSP', 'ZCNSP', 
            'ZCONSP', 
            'F', 'Cl', 'Br', 'halogen','mix_X',
            'CONSP_fr1', 'F_fr1', 'Cl_fr1', 'Br_fr1','mix_X_fr1', 
            'CONSP_fr2', 'X_fr2','fr2', 
            'CONSP_mfr', 'mfr', 'X_mfr',
            'tf','scaffold', 
            ]


from .. import molxpert as MX

class QMConformationSearch:
    """
    Q0: initial structure opt, freq, charge calculation
    Q1: torsion scan
    Q2: bond stretch and angle bending
    Q3: torsion 2D scan
    Q4: local minimum search from torsion scan
    Q5: random conformation search
    Q6: oligomer (intermolecular)
    Q7: transition state
    Q10: SP calculation
    Q100: CCSD(T)/CBS calculation
    """

    def __init__(self,config) -> None:
        self.config = config
        self.config = MX.update_configure(self.config)
        self.stage = self.config["QMSetting"]["stage"]
        self.opath = self.config["EnvironmentSetting"]["output_directory"]
        self.qmpara = self.config[f"{self.stage}QMSetting"]
        self.ignore_alkane = self.config["QMSetting"]["ignore_alkane"]
        self.datasearch = self.config["datasearch"]
        self.extra_var = {"datasearch":self.datasearch}
        if self.stage == "Q6":
            self._run_initial()
            self.molecules = MX.get_from_db(self.config)
            self._run()
        else:
            if not isinstance(self.config["MoleculeFileSetting"]["molecules"],list):
                self.input_list = [self.config["MoleculeFileSetting"]["molecules"]]
            else:
                self.input_list = self.config["MoleculeFileSetting"]["molecules"]
            #self.stage = self.config["QMSetting"]["stage"]
            
            self._run_initial()
            #self._get_extra_var()

            self._generate_molecules()
            self._run()


    def _run_initial(self):
        self.max_heavy = None
        #self.conform_selector_extra = None
        self.conform_flag = False
        
        self.molecules = []

        self.conform_flag = True
        self.qm_method = "RB3LYP"
        self.qm_basis_set = "def2SVP"
        self.conform_selector_extra = None
        self.all_conform_flag = False
        self.get_attrs = ["energy", "force", "hessian", "freq", "esp_charge", "conform_type"]
        self.return_origin_conformation=True

    def _get_extra_var(self):
        if self.stage == "Q0":
            self.extra_var = None

        elif self.stage == "Q10":
            self.extra_var = {"conformer_setting":
                              {"conformer_selector":
                               {"confID":{"$exists":True},"conform_type":{"$ne":"optimizing"}}
                               }
                               }
        elif self.stage == "Q1":
            self.extra_var = {"conformer_setting":
                              {"conformer_selector":
                               {"confID":{"$exists":True},"conform_type":{"$ne":"optimizing"}}
                               },
                               "select_setting":"stablest_conformer"
                               }
        elif self.stage == "Q2":
            self.extra_var = {"conformer_setting":
                              {"conformer_selector":
                               {"confID":{"$exists":True},"conform_type":{"$ne":"optimizing"}}
                               },
                               "select_setting":"stablest_conformer"
                               }
        elif self.stage == "Q4":
            #self.extra_var = {"conformer_setting":
            #                  {"conformer_selector":
            #                   {"confID":{"$exists":True},"conform_type":{"$ne":"optimizing"}}
            #                   },
            #                   }
            self.extra_var = {"conformer_setting":
                              {"conformer_selector":
                               {"conform_type":{"$ne":"optimizing"}}
                               },
                               }
        else:
            self.conform_flag = False
            self.extra_var = {}

    def _generate_molecules(self):
        from ..craton.chem.molecule import Molecule
        for ii, _input in enumerate(self.input_list):
            if isinstance(_input,Molecule):
                self.molecules.append(_input)
            else:
                molecules = MX.molecule_create(_input,extra_var=self.extra_var,show_figure=False,parallel=False)
                self.molecules.extend(molecules)
        
    def _run(self):
        if self.stage == "Q0":
            self._run_Q0()
        elif self.stage == "Q1":
            self._run_Q1()
        elif self.stage == "Q4": # Q4 r2
            self._run_Q4()
        elif self.stage == "Q2": # Q2 r4
            self._run_Q2()
        elif self.stage == "Q6": # Q2 r4
            self._run_Q6()
        elif self.stage == "Q10":
            self._run_Q10()
        elif self.stage == "Q8":
            self._run_Q8()

    def _run_Q10(self):
        MX.qm_input_file(self.molecules,self.qmpara,"Q10",local_path=self.opath,indexs=[molecule.confID for molecule in self.molecules])

        #for molecule in self.molecules:
        #    MX._qm_input_file(molecule,self.setting_file,"Q10",self.opath,index=molecule.confID,)

    def _run_Q0(self):
        MX.qm_input_file(self.molecules,self.qmpara,"Q0",local_path=self.opath)
        #for molecule in self.molecules:
        #    MX._qm_input_file(molecule,self.setting_file,"Q0",self.opath,)

    def _run_Q1(self):
        self.molecules = MX.molecule_structure(self.molecules)
        self.molecules = MX.molecule_torsion(self.molecules)
        if self.ignore_alkane:
            self.molecules = MX.atom_type(self.molecules)
            self.molecules = MX._ignore_alkane_torsion_(self.molecules)
        
        MX.qm_input_file(self.molecules,self.qmpara,"Q1",local_path=self.opath)
        #for molecule in self.molecules:
        #    MX._qm_input_file(molecule,self.setting_file,"Q1",self.opath)

    def _run_Q2(self):
        self.molecules = MX.molecule_structure(self.molecules)
        self.molecules = MX.atom_type(self.molecules)
        molecule_dict = MX.Q2_bond_angle_conformer_(self.molecules,ignore_alkane=self.ignore_alkane)
        for molecule in self.molecules:
            indexs = []
            zmatrixs = []
            for name, zmat in molecule_dict[molecule.inchi_key].items():
                indexs.append(name)
                zmatrixs.append(zmat)
            MX.qm_input_file([molecule for _ in molecule_dict[molecule.inchi_key]],
                              self.qmpara,
                              "Q2",
                              local_path=self.opath,
                              indexs=indexs,
                              zmatrixs=zmatrixs
                              )
            #for name,zmat in molecule_dict[molecule.inchi_key].items():
            #    MX._qm_input_file(molecule,self.setting_file,"Q2",self.opath,index=name,zmatrix=zmat)

    def _run_Q8(self):
        self.molecules = MX.update_structure_topol(self.molecules)
        self.molecules = MX.molecule_structure(self.molecules)
        if self.qmpara["exists_type"] is not None:
            self.molecules = MX.atom_type(self.molecules)
        for molecule in self.molecules:

            scan_term = MX.Q8_bond_angle_scan(molecule,
                                              inter_val=self.qmpara["inter_val"],
                                              ignore_ring=self.qmpara["ignore_ring"],
                                              #exists_type=self.qmpara["exists_type"]
                                              exists_type=[],
                                              )
            for tt in scan_term:
                this_qmpara = deepcopy(self.qmpara)
                this_qmpara["scan_term"] = tt[0]
                this_qmpara["scan_setting"] = tt[1]
                if tt[1][1] > 0:
                    _tmp = "-f"
                else:
                    _tmp = "-r"
                MX.qm_input_file(molecule,this_qmpara,"Q8",local_path=self.opath,fpath_pre=_tmp,parallel=False)

    def _run_Q6(self,):
        molecules = []
        for chemsystem in self.molecules:
            tmp = MX.molecule_create([chemsystem["molecule1"],chemsystem["molecule2"]])
            molecule = tmp[0]
            molecule.mole_name = chemsystem["name"]
            molecule.smiles = chemsystem["name"]
            molecule.inchi_key = chemsystem["name"]
            molecule.Atoms.extend(tmp[1].Atoms)
            molecules.append(molecule)
        
        MX.qm_input_file(molecules,self.qmpara,"Q6",local_path=self.opath)

    def old_run_Q6(self,optimize_structure=True):
        parent_path = self.config["MoleculeFileSetting"]["molecules"]
        parent_dds = [dd for dd in os.listdir(parent_path) if os.path.isdir(f"{parent_path}/{dd}")]
        for dd in parent_dds:
            print(dd)
            this_path = f"{parent_path}/{dd}"
            this_output_path = f"{self.opath}/{dd}"
            Path(this_output_path).mkdir(exist_ok=True)
            this_dds = [ddd for ddd in os.listdir(this_path)]
            molecules = []
            for ddd in this_dds:
                print(ddd)
                this_this_path = f"{this_path}/{ddd}"
                residue1 = MX.molecule_create(f"{this_this_path}/{ddd}.mtx")[0]
                mtxs = [f"{this_this_path}/dimer/{fff}" for fff in os.listdir(f"{this_this_path}/dimer") if fff.find(".mtx") != -1]
                residue2s = MX.molecule_create(mtxs)
                for res2 in residue2s:
                    _tmp = deepcopy(residue1)
                    _tmp.Atoms.extend(res2.Atoms)
                    _tmp.mole_name = f"{residue1.mole_name}_{res2.mole_name}"
                    _tmp.smiles = _tmp.mole_name
                    _tmp.inchi_key = _tmp.mole_name
                    molecules.append(_tmp)
            MX.qm_input_file(molecules,self.qmpara,"Q6",local_path=this_output_path,parallel=False)


    def _run_Q4(self,optimize_structure=True):
        self.molecules = MX.molecule_structure(self.molecules,parallel=False)
        scan_curve = MX.scan_curve(self.molecules)
        stable_molecules = MX.find_stablest_conformer(self.molecules)
        stable_molecules = MX.atom_type(stable_molecules)

        stable_molecules = MX.grasp_force_field(stable_molecules)



        rlm_dict = MX.pes_local_minimum(scan_curve)
        rlm_dict = [rlm_dict[molecule.inchi_key] for molecule in stable_molecules]
        total_lm_molecules = MX.lm_by_combine_scan_curve(stable_molecules,rlm_dict,n=64,create_constrain=False,parallel=False)

        for ii,molecule in enumerate(stable_molecules):
            lm_molecules = MX.conformer_expand(molecule,total_lm_molecules[ii])
            lm_molecules = MX.update_structure_topol(lm_molecules)
            #lm_molecules = total_lm_molecules[ii]
            if optimize_structure:
                lm_molecules = MX._optimize(lm_molecules,optimizer="openmm", all_torsion_constraint = 0.0, write_mol=None)

            after_molecules = MX.remove_similar_conformer(lm_molecules,target_molecule=molecule)
            
            MX.qm_input_file(after_molecules,self.qmpara,"Q4",local_path=self.opath,indexs=[ii for ii in range(len(after_molecules))])
            #for ii, mm in enumerate(after_molecules):
            #    MX._qm_input_file(mm,self.setting_file,"Q4",self.opath,index=ii)

def get_calculated_molecules(nn=None):
    if nn is None:
        nn = TASK_MAX_MOLs
    flag = False
    count = 0
    select_moles = {"Q0":[],"Q1":[],"Q2":[],"Q3":[],"Q4":[]}
    for items in class2:
        selector = {"class2":items[0], "steps":{"$ne":items[1]}}
        for fg in class3:
            selector["class3"] = fg
            for doc in qmcalc_coll.find(selector,{"name":1}):
                if count <= nn:
                    select_moles[items[1]].append(doc["name"])
                    qmcalc_coll.update_one({"name":doc["name"]},{"$push":{"steps":items[1]}})
                    count += 1
                else:
                    flag = True
                    break
            if flag:
                break
        if flag:
            break
    return select_moles

def create_qm_jobs(select_moles,fpath=None):
    if fpath is None:
        fpath = "."
    for kk,vv in select_moles.items():
        if len(vv) > 0:
            output_dir = f"{fpath}/{kk}"
            config = {"output_directory":output_dir,"stage":kk,"molecules":vv}
            qmc = QMConformationSearch(config)

def meta_qm_calc(fpath=None,nn=None):
    if fpath is None:
        fpath = f"./{datetime.date.today()}"
        if os.path.exists(fpath):
            for ii in range(1,10000):
                if not os.path.exists(f"{fpath}-{ii}"):
                    fpath = f"{fpath}-{ii}"
                    break
        os.mkdir(fpath)

    select_moles = get_calculated_molecules(nn=nn)
    create_qm_jobs(select_moles,fpath=fpath)