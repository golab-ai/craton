
from .. import CRATON_CONFIGURE

import json
from ..utils.commons import parallel_run

from .atom_type import TypingEngine, assign_atom_type_to_term, TypingDefine
from .forcefield_manager import ForceFieldManager as FFM
from .fitting import ForceFieldFitting
from .analyze import FittingAnlayze

def get_atoms_arr_and_this_terms(molecules,this_terms):
    _this_terms = this_terms
    _atoms_arr = []
    for jj,terms in enumerate(this_terms):
        if terms is None:
            _atoms_arr.append(None)
        else:
            tmp = []
            for kk,vv in terms.items():
                if kk == "binc":
                    for kkk,vvv in vv.items():
                        tmp.append(kkk)
                        tmp.extend([molecules[jj].Atoms[kkk].connect[vvvv] for vvvv in vvv])
                elif kk == "pair_charge":
                    for kkk,vvv in vv.items():
                        items = getattr(molecules[jj],kkk,[])
                        for jjj in vvv:
                            item = items[jjj]
                            tmp.append(item.a1)
                            tmp.append(item.a2)
                else:
                    items = getattr(molecules[jj],kk,[])
                    for jjj in vv:
                        item = items[jjj]
                        for an in range(1,5):
                            if hasattr(item,f"a{an}"):
                                tmp.append(getattr(item,f"a{an}"))
            _atoms_arr.append(list(set(tmp)))
    return _this_terms, _atoms_arr

class AtomType:
    def __init__(self,atf=None,check_at_types=False) -> None:
        self.at = TypingEngine(atversion=atf,check_at_types=check_at_types)

    def _atom_type(self,molecules,
                   atoms_arr=None,
                   equtable_flag=True,
                   create_improper=True,
                   raise_invalid_charge=False,
                   raise_unsupported_coordination=False,
                   assign_atom_type_flag=True,
                   ignore_existing=False,
                   ignore_ff_existing=False,
                   this_terms=None,
                   parallel=True,
                   ):
        if not isinstance(molecules,list):
            molecules = [molecules]
            if atoms_arr is not None:
                atoms_arr = [atoms_arr]
            if this_terms is not None:
                this_terms = [this_terms]
        _atoms_arr = [None for ii in range(len(molecules))]
        _this_terms = [None for ii in range(len(molecules))]
        
        if atoms_arr is not None:
            _atoms_arr = atoms_arr
        
        if this_terms is not None:
            _this_terms,_atoms_arr = get_atoms_arr_and_this_terms(molecules,this_terms)              
        if parallel:
            molecules = parallel_run(self.at.assign_mole_at,molecules,kwds=[{"atoms_arr":rr,"ignore_existing":ignore_existing} for rr in _atoms_arr],keep_order=False)
            if assign_atom_type_flag:
                #args = {"atom_type_rule":self.at.at_types["root"],
                #        "equtable_flag":equtable_flag,
                #        "create_improper":create_improper,
                #        "raise_invalid_charge":raise_invalid_charge,
                #        "raise_unsupported_coordination":raise_unsupported_coordination,
                #        "ignore_existing":ignore_existing,
                #        "this_terms": _this_terms,
                #        "atom_arr": _atoms_arr
                #        }
                args = [{"atom_type_rule":self.at.at_types["root"],
                        "equtable_flag":equtable_flag,
                        "create_improper":create_improper,
                        "raise_invalid_charge":raise_invalid_charge,
                        "raise_unsupported_coordination":raise_unsupported_coordination,
                        "ignore_existing":ignore_existing,
                        "this_terms": _this_terms[ii],
                        "atom_arr": _atoms_arr[ii]} for ii in range(len(_atoms_arr))]
                molecules = parallel_run(assign_atom_type_to_term,molecules,kwds=args,keep_order=False)
        else:
            for ii,molecule in enumerate(molecules):
                self.at.assign_mole_at(molecule,atoms_arr=_atoms_arr[ii],ignore_existing=ignore_existing,ignore_ff_existing=ignore_ff_existing)
                if assign_atom_type_flag:
                    assign_atom_type_to_term(molecule,atom_type_rule=self.at.at_types["root"],
                                     equtable_flag=equtable_flag,
                                     create_improper=create_improper,
                                     raise_invalid_charge=raise_invalid_charge,
                                     raise_unsupported_coordination=raise_unsupported_coordination,
                                     ignore_existing=ignore_existing,
                                     this_terms=_this_terms[ii],
                                     atom_arr=_atoms_arr[ii]
                                     )
        for molecule in molecules:
            if "atom type" not in molecule.steps:
                molecule.steps.append("atom type")
        return molecules


    @staticmethod
    def _set_ats_term(molecule):
        assign_atom_type_to_term(molecule,create_improper=False)


    @staticmethod
    def _convert_atom_type_file_to_json(atf):
        atdefine = TypingDefine(atversion=atf,check_at_types=False,convert_flag=True)
        new_file = atf[:-4] + ".json"
        with open(new_file,'w') as outf:
            outf.write(json.dumps(atdefine.at_types))

class MolForceField:
    def __init__(self) -> None:
        pass

    @staticmethod
    def get_force_field(molecules,
                           atom_type_file,
                           force_field_file,
                           this_terms=None,
                           reassign_atom_type=False,
                           compensating_at_file=None,
                           compensating_ff_file=None,
                           charge_method=None,
                           charge_ff = None,
                           ignore_existing=False,
                           return_ff=False,
                           empi_ff_flag = True,
                           use_scalevdw=True,
                           parallel=True,):
        if not isinstance(molecules,list):
            molecules = [molecules]
        if atom_type_file is None:
            atom_type_file=CRATON_CONFIGURE["ForceFieldSetting"]["DEFAULT_TYPING_FILE"]
        if force_field_file is None:
            force_field_file=CRATON_CONFIGURE["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"]
        if empi_ff_flag:
            empi_ff = CRATON_CONFIGURE["ForceFieldSetting"]["EMPIRICAL_FORCE_FIELD_FILE"]
        else:
            empi_ff = None
        molecules = FFM.get_force_field(molecules,
                           atom_type_file,
                           force_field_file,
                           this_terms=this_terms,
                           reassign_atom_type=reassign_atom_type,
                           compensating_at_file=compensating_at_file,
                           compensating_ff_file=compensating_ff_file,
                           charge_method=charge_method,
                           charge_ff = charge_ff,
                           ignore_existing=ignore_existing,
                           return_ff=return_ff,
                           empi_ff = empi_ff,
                           use_scalevdw=use_scalevdw,
                           parallel=parallel,
                               )
        return molecules
        
    @staticmethod
    def grasp_force_field(molecules,
                          atom_type_file,
                          force_field_file,
                          reassign_atom_type=False,
                          charge_method=None,
                          charge_ff=None,
                          ignore_existing=False,
                          empi_ff_flag = True,
                          use_scalevdw=True,
                          return_ff=False,
                          parallel=True,
                          ):
        if not isinstance(molecules,list):
            molecules = [molecules]
        if atom_type_file is None:
            atom_type_file=CRATON_CONFIGURE["ForceFieldSetting"]["DEFAULT_TYPING_FILE"]
        if force_field_file is None:
            force_field_file=CRATON_CONFIGURE["ForceFieldSetting"]["DEFAULT_FORCE_FIELD_FILE"]
        if empi_ff_flag:
            empi_ff = CRATON_CONFIGURE["ForceFieldSetting"]["EMPIRICAL_FORCE_FIELD_FILE"]
        else:
            empi_ff = None
        molecules = FFM.assign_force_field(molecules,
                               atom_type_file,
                               force_field_file,
                               reassign_atom_type=reassign_atom_type,
                               charge_method=charge_method,
                               charge_ff=charge_ff,
                               ignore_existing=ignore_existing,
                               empi_ff = empi_ff,
                               use_scalevdw=use_scalevdw,
                               return_ff=return_ff,
                               parallel=parallel,
                               )
        return molecules

    @staticmethod
    def checkout_force_field(molecules,ffjson):
        return FFM.checkout_force_field(molecules,ffjson)

    @staticmethod
    def get_fitting_parameters(this_ff,
            fix_tag=["V", "Fit"],
            terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
            preprocessing_fitting_parameter=True,
            molecules=None,
            flag_fitting_scan_parameter=False,
            atom_type_file=None,
            ):

        return FFM.get_fitting_parameters(
            this_ff,
            fix_tag=fix_tag,
            terms = terms,
            preprocessing_fitting_parameter=preprocessing_fitting_parameter,
            molecules=molecules,
            flag_fitting_scan_parameter=flag_fitting_scan_parameter,
            atom_type_file=atom_type_file,
        )

    @staticmethod
    def assign_force_field_parameter(molecules,this_ff,this_terms=None,parallel=True):
        return FFM.assign_force_field_parameter(
            molecules,
            this_ff,
            this_terms=this_terms,
            parallel=parallel
        )

    @staticmethod
    def combine_ff(ff_ref, ff_new, prior_tags=None):
        return FFM.combine_ff(ff_ref,ff_new,prior_tags=prior_tags)

    @staticmethod
    def combine_ff_file(ff_file_ref, ff_file_new, used_vdw,out_path):
        FFM.combine_ff_file(ff_file_ref, ff_file_new, used_vdw, out_path)

    @staticmethod
    def write_ff_file(this_ff,opath):
        FFM.transfer(this_ff,opath)

    @staticmethod
    def _special_dihedrals_infos(new_molecules,this_ff):
        FFM.get_special_infos(new_molecules,this_ff)

    @staticmethod
    def _convert_force_field_file_to_json(fff):
        import json
        ffjson = FFM._convert_force_field_format_to_json(fff)
        new_file = fff + ".json"
        with open(new_file,'w') as outf:
            outf.write(json.dumps(ffjson))

    @staticmethod
    def create_empirical_force_field(atf,fff1,fff2):
        FFM._get_empirical_force_field(atf,fff1,fff2)

    @staticmethod
    def read_force_field(fff,use_scalevdw):
        return FFM.read_force_field_file(fff,use_scalevdw)
    
    