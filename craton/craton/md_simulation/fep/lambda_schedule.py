from typing import Union

import numpy as np
import pandas as pd

from ...utils import logger


def smooth_function(function_type):
    if function_type == 0:
        return lambda x: x
    if function_type == 1:
        return lambda x: -2 * x**3 + 3 * x**2
    if function_type == 2:
        return lambda x: 6 * x**5 - 15 * x**4 + 10 * x**3
    if function_type == 3:
        return lambda x: -20 * x**7 + 70 * x**6 - 84 * x**5 + 35 * x**4
    if function_type == 4:
        return lambda x: 70 * x**9 - 315 * x**8 + 540 * x**7 - 420 * x**6 + 126 * x**5


class LambdaSchedule:
    @staticmethod
    def hfe_lambda():
        half_window = 8
        half_vdw = smooth_function(0)(np.linspace(0, 1, half_window))
        half_coul = smooth_function(2)(np.linspace(0, 1, half_window))
        bonded = smooth_function(2)(np.linspace(0, 1, half_window * 2))
        coul = np.concatenate((half_coul, np.ones(half_window)))
        vdw = np.concatenate((np.zeros(half_window), half_vdw))
        df = pd.DataFrame.from_records({"bonded": bonded, "vdw": vdw, "coul": coul})
        df = df.round(3)
        return df

    def __init__(self, fep_setting=None, is_charge_hopping=False, is_core_hopping=False,is_abfe=False,mixed_lambda=False,is_relative=False):
        self.fep_setting = fep_setting
        self.mixed_lambda = mixed_lambda
        #self.is_relative = isinstance(self.fep_setting, RBFEFEPSetting)
        self.is_relative = is_relative
        self.is_charge_hopping = is_charge_hopping
        self.is_core_hopping = is_core_hopping

    def generate_lambdas(self):
        if self.fep_setting is None:
            coul_lambdas = None
            vdw_lambdas = None
            bond_lambdas = None
        else:
            coul_lambdas, vdw_lambdas, bond_lambdas = (
                self.fep_setting["coul_lambdas"],
                self.fep_setting["vdw_lambdas"],
                self.fep_setting["bonded_lambdas"],
            )
        if coul_lambdas is None or vdw_lambdas is None:
            if self.is_charge_hopping:
                return self.default_charge_rule()
            elif self.is_core_hopping:
                return self.default_core_hopping_rule()
            else:
                return self.default_rgroup_rule()
        else:
            #logger.debug("read lambda from user input file")
            return self.read_lambda(coul_lambdas, vdw_lambdas, bond_lambdas)

    def default_rgroup_rule(self):
        bond_smooth_type = 0
        vdw_smooth_type = 2
        coul_smooth_type = 2

        if self.is_relative:
            half_window = 8

            half_vdw = smooth_function(vdw_smooth_type)(np.linspace(0, 1, half_window))
            half_coul = smooth_function(coul_smooth_type)(np.linspace(0, 1, half_window))
            bonded = smooth_function(bond_smooth_type)(np.linspace(0, 1, half_window * 2))

            if self.mixed_lambda:  # change the vdw and coul same
                perturbB_vdw = np.concatenate((np.ones(half_window), np.zeros(half_window)))
                perturbA_vdw = np.concatenate((np.zeros(half_window), np.ones(half_window)))
                vdw = np.concatenate((half_vdw, half_vdw))
                coul = np.concatenate((half_coul, half_coul))
            else:
                perturbB_vdw = np.concatenate((np.ones(16), np.zeros(8)))
                perturbA_vdw = np.concatenate((np.zeros(16), np.ones(8)))
                vdw = np.concatenate((half_vdw, np.ones(8), half_vdw))
                coul = np.concatenate((np.zeros(8), half_coul, np.ones(8)))
                bonded = np.concatenate((bonded, np.ones(half_window)))
            df = pd.DataFrame.from_records(
                {"bonded": bonded, "vdw": vdw, "coul": coul, "perturbB_vdw": perturbB_vdw, "perturbA_vdw": perturbA_vdw}
            )

        else:
            half_window = 24
            half_vdw = smooth_function(vdw_smooth_type)(np.linspace(0, 1, half_window))
            half_coul = smooth_function(coul_smooth_type)(np.linspace(0, 1, half_window))
            bonded = smooth_function(bond_smooth_type)(np.linspace(0, 1, half_window * 2))

            vdw = np.concatenate((np.zeros(half_window), half_vdw))
            coul = np.concatenate((half_coul, np.ones(half_window)))
            df = pd.DataFrame.from_records({"bonded": bonded, "vdw": vdw, "coul": coul})
        df = df.round(3)
        return df

    def read_lambda(self, coul_lambdas, vdw_lambdas, bond_lambdas) -> pd.DataFrame:
        
        n = len(vdw_lambdas)
        if self.is_relative:
            perturb_idx = 0
            if self.mixed_lambda:
                perturb_idx = np.argmax(vdw_lambdas) + 1
                if not np.all(np.diff(vdw_lambdas[:perturb_idx]) > 0):
                    raise RuntimeError("vdw lambda error!")
                if not np.all(np.diff(vdw_lambdas[perturb_idx:]) > 0):
                    raise RuntimeError("vdw lambda error!")
                if not bond_lambdas:
                    bond_lambdas = smooth_function(1)(np.linspace(0, 1, n))
            else:
                for i in range(len(vdw_lambdas)):
                    if vdw_lambdas[i] == 1 and vdw_lambdas[i + 1] == 0:  # two stages
                        perturb_idx = i + 1
                        break
                if not bond_lambdas:
                    bond_lambdas = np.concatenate(
                        (smooth_function(1)(np.linspace(0, 1, perturb_idx)), np.ones(n - perturb_idx))
                    )
            if not perturb_idx:
                raise RuntimeError("vdw lambda error!")
            perturbB_vdw = np.concatenate((np.ones(perturb_idx), np.zeros(n - perturb_idx)))
            perturbA_vdw = np.concatenate((np.zeros(perturb_idx), np.ones(n - perturb_idx)))
            df = pd.DataFrame.from_records(
                {
                    "bonded": bond_lambdas,
                    "vdw": vdw_lambdas,
                    "coul": coul_lambdas,
                    "perturbB_vdw": perturbB_vdw,
                    "perturbA_vdw": perturbA_vdw,
                }
            )
        else:
            if bond_lambdas is None:
                #print("bond_lambdas:",bond_lambdas,"abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*abfe*")
                #bond_lambdas = smooth_function(1)(np.linspace(0, 1, len(vdw_lambdas)))
                df = pd.DataFrame.from_records({"vdw": vdw_lambdas, "coul": coul_lambdas})
            else:
                df = pd.DataFrame.from_records({"bonded": bond_lambdas, "vdw": vdw_lambdas, "coul": coul_lambdas})
            # do not allow coul and vdw change at same time
        if not self.is_relative or not self.mixed_lambda:
            for i in range(1, n):
                vdw_perturb = vdw_lambdas[i] - vdw_lambdas[i - 1] > 0
                coul_perturb = coul_lambdas[i] - coul_lambdas[i - 1] > 0
                if vdw_perturb and coul_perturb:
                    raise RuntimeError("Vdw and coul lambdas are not allowed change at the same time!")
        return df

    def default_charge_rule(self):
        half_window = 12
        vdw_smooth_type = 2
        coul_smooth_type = 2
        bond_smooth_type = 0

        half_vdw = smooth_function(vdw_smooth_type)(np.linspace(0, 1, half_window))
        half_coul = smooth_function(coul_smooth_type)(np.linspace(0, 1, half_window))

        perturbB_vdw = np.concatenate((np.ones(half_window), np.zeros(half_window)))
        perturbA_vdw = np.concatenate((np.zeros(half_window), np.ones(half_window)))
        vdw = np.concatenate((half_vdw, half_vdw))
        coul = np.concatenate((half_coul, half_coul))
        bonded = smooth_function(bond_smooth_type)(np.linspace(0, 1, half_window * 2))
        df = pd.DataFrame.from_records(
            {"bonded": bonded, "vdw": vdw, "coul": coul, "perturbB_vdw": perturbB_vdw, "perturbA_vdw": perturbA_vdw}
        )
        df = df.round(3)
        return df

    def abfe_lambda(self):
        if "vdw_lambdas" not in self.fep_setting:
            return LambdaSchedule.hfe_lambda()
        else:
            return pd.DataFrame({"vdw":self.fep_setting["vdw_lambdas"],"coul":self.fep_setting["coul_lambdas"],})

    def default_core_hopping_rule(self):
        return self.default_rgroup_rule()


if __name__ == "__main__":
    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.use("macosx")

    def plot_smooth():
        plt.figure()
        x = np.linspace(0, 1, 100)
        for i in range(5):
            y = smooth_function(i)(x)
            plt.plot(x, y, label=f"type {i}")
        plt.legend()
        plt.show()

    def plot_mixed_lambda(mixed_lambda):
        plt.figure()
        df = LambdaSchedule(mixed_lambda=mixed_lambda).default_rgroup_rule()
        df.plot(kind="line")
        plt.show()

    plot_mixed_lambda(True)
    # print(df)
