import statistics
from ..utils import logger

def get_am1bcc_charge(molecules):
    from ..mm_calculator import Calculator as Calc
    charges = {}
    for molecule in molecules:
        #try:
        Calc.am1bcc_charge(molecule)
        charges[molecule.mole_name] = molecule.ff_charge
        #except:
        #    logger.error("Mol2 with am1bcc information not found")
        #    raise Exception("Calc am1bcc charge failed")
    return charges

def get_QM_charge(molecules, charge_type="esp"):
    """
    得到qm计算的电荷
    """

    charges = {}
    for molecule in molecules:
        if getattr(molecule,charge_type) is not None:
            if molecule.mole_name not in charges:
                charges[molecule.mole_name] = []
            charges[molecule.mole_name].append(getattr(molecule,charge_type))
    for mname,rr in charges.items():
        ci = len(rr[0])
        charges[mname] = [statistics.mean([rrr[ii] for rrr in rr]) for ii in range(ci)]
    return charges
    
def get_nn_charge(molecules):
    return {}

def get_nonbinc_charge(molecules,charge_method,ignore_existing=False):
    if charge_method not in ["esp", "mulliken","am1bcc", "nn","qm"]:
        raise Exception("Invalid charge_method. Should be esp, am1bcc, nn, mulliken")

    this_molecules = []
    for molecule in molecules:
        if not ignore_existing and "force field" in molecule.steps:
            pass
        else:
            this_molecules.append(molecule)

    if charge_method == "qm":
        charge_method = "esp"
            
    if charge_method == "am1bcc":
        charges= get_am1bcc_charge(this_molecules)
    elif charge_method in ["nn","ai"]:
        charges = get_nn_charge(this_molecules)
    else:
        charges = get_QM_charge(this_molecules,charge_type=charge_method)

    return charges