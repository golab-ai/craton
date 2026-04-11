import os
from typing import Any

AVOGADRO = 6.02 #6.02e23
CAL_TO_J = 4.184
HA_TO_KCAL_MOL = 627.51
EV_TO_KCAL_MOL = 23.0605
BOLTZMANN = 1.38e-23 #J/K
PLANCK = 6.626e-34 #J.s
REDUCED_PLANCK = 1.05457e-34
IDEA_GAS_CONSTANT = 8.314

class Unit:
    def __init__(self,value,reference,target,extra=None):
        self.value = value
        self.reference = reference
        self.target = target
        self.extra = extra
        self.style = "concentration"

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        func = self.__FUNC[self.style]
        return func(self)

    __units = {
        "concentration":["n/A^3","n/nm^3","mol/ml","mol/l","mmol/l","g/ml","mg/ml","mg/l","g/l","kg/l","kg/m^3","%"]
    }

    def _concentration(self):
        func = {
            "n/A^3":[lambda x : x * 1, lambda y : y / 1.0],
            "n/nm^3":[lambda x : x / 0.001, lambda y : y * 0.001],
            "mol/ml":[lambda x : x * AVOGADRO * 10.0, lambda y : y / 10.0 / AVOGADRO],
            "mol/l":[lambda x : x * 0.0001 * AVOGADRO,lambda y : y * 10000.0 / AVOGADRO],
            "mmol/l":[lambda x : x * 0.000001 * AVOGADRO, lambda y : y / 0.000001 / AVOGADRO],
            "g/ml" :[lambda x : x * AVOGADRO * 0.1 / self.extra["mass"], lambda y : y * self.extra["mass"] / AVOGADRO / 0.1],
            "mg/ml":[lambda x : x * AVOGADRO * 0.0001 / self.extra["mass"], lambda y : y * self.extra["mass"] / AVOGADRO / 0.0001],
            "mg/l":[lambda x : x * AVOGADRO * 0.000001 / self.extra["mass"], lambda y : y * self.extra["mass"] / AVOGADRO / 0.000001],
            "g/l":[lambda x : x * AVOGADRO * 0.0001 / self.extra["mass"], lambda y : y * self.extra["mass"] / AVOGADRO / 0.0001],
            "kg/l":[lambda x : x * AVOGADRO * 0.1 / self.extra["mass"], lambda y : y * self.extra["mass"] / AVOGADRO / 0.1],
            "kg/m^3":[lambda x : x * AVOGADRO * 0.0001 / self.extra["mass"], lambda y : y * self.extra["mass"] / AVOGADRO / 0.0001],
            "%1":[
                  lambda x : x * self.extra["density"] * AVOGADRO * 0.1 / self.extra["mass"],
                  lambda y : y * self.extra["mass"] / self.extra["density"] / AVOGADRO / 0.1
                  ],
            "%2":[
                  lambda x : x * self.extra["solvent_density"] * AVOGADRO * 0.1 / self.extra["mass"] / (1 - x),
                  lambda y : y * self.extra["mass"] * (1 - y) / AVOGADRO / 0.1 / self.extra["solvent_density"]
            ],
        }
        if self.reference == "%":
            if "density" in self.extra:
                self.reference = "%1"
            else:
                if "solvent_density" in self.extra:
                    self.reference = "%2"

        if self.target == "%":
            if "density" in self.extra:
                self.target = "%1"
            else:
                if "solvent_density" in self.extra:
                    self.target = "%2"
        v0 = func[self.reference][0](self.value)
        return func[self.target][1](v0)
    __FUNC = {
        "concentration": _concentration,
        }

