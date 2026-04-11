from .fitting import binc_fitting, fitting
from .validation import validation

class ForceFieldFitting:
    
    @staticmethod
    def _validation(
                    molecules,
                    this_ff,
                    output_dir="./",
                    optimizer="openmm",
                    hessian_flag=False,
                    fitting_info=None,
                    init_this_ff=None,
                    ):
        return validation(
                    molecules,
                    this_ff,
                    output_dir,
                    optimizer=optimizer,
                    hessian_flag=hessian_flag,
                    fitting_info=fitting_info,
                    init_this_ff=init_this_ff
        )
    
    @staticmethod
    def _binc_fitting(molecules, this_ff, target= "esp"):
        binc_fitting(molecules,this_ff,target=target)

    @staticmethod
    def _intra_fitting(this_ff,
                        molecules,
                        fitting_terms=["bondterm", "angleterm", "dihedralterm", "improperterm", "binc"],
                        target_prop = ["energy", "force", "hessian", "penalty_torsion"],
                        torsion_constraint_step=None,
                        optimizer="openmm",):
        
        return fitting(this_ff,
                        molecules,
                        fitting_terms=fitting_terms,
                        target_prop = target_prop,
                        torsion_constraint_step=torsion_constraint_step,
                        optimizer=optimizer,
        )
