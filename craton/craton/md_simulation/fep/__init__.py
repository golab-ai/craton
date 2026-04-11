from .lambda_schedule import LambdaSchedule
from .abfe import get_intermolecular_restrain

def get_lambda_schedule(fep_setting,fep_type="r_group",mixed_lambda=False,is_relative=False):
    if fep_type == "r_group":
        return LambdaSchedule(fep_setting=fep_setting,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()
    elif fep_type == "charge_hopping":
        return LambdaSchedule(fep_setting=fep_setting, is_charge_hopping=True,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()
    elif fep_type == "core_hopping":
        return LambdaSchedule(fep_setting=fep_setting, is_core_hopping=True,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()
    elif fep_type == "couple":
        return LambdaSchedule.hfe_lambda()
    else:
        return LambdaSchedule(fep_setting=fep_setting,mixed_lambda=mixed_lambda,is_relative=is_relative).generate_lambdas()

def abfe_intermolecule(system):
    return get_intermolecular_restrain(system)
