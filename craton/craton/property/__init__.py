from .property import calc_inertia
from .property import calc_dipole
from .property import calc_multipolar
#from .expt_equ import require_expt_data as RED
#from .expt_equ import get_all_json_property 
#from .expt_equ import SplitDataSet

class MolProperty:
    def __init__(self) -> None:
        pass

    @staticmethod
    def _inertia_calculate(molecule,ignore_hydrogen=False):
        return calc_inertia(molecule,ignore_hydrogen=False)
    
    @staticmethod
    def _dipole_calculate(molecule):
        return calc_dipole(molecule)
    
    @staticmethod
    def thermodyna_property(inputs,props,molecule_type=None,sources=None,temperatures=None,pressures=None,condinations=None):
        return RED(inputs,props,molecule_type=molecule_type,sources=sources,temperatures=temperatures,pressures=pressures,condinations=condinations)
    
    @staticmethod
    def get_all_property_at_temperature(inf,temperature,output_dir=".",outfn = "property"):
        results,error = get_all_json_property(inf,temperature,fn=f"{output_dir}/{outfn}_{temperature}K.json")
        return results, error
        
    @staticmethod
    def split_train_set(inf,prop,style,value=None,output_dir=".",test_flag=True):
        SDS = SplitDataSet(inf,prop,style,value=value,output_dir=output_dir,test_flag=test_flag)
        SDS.run()
    
    
