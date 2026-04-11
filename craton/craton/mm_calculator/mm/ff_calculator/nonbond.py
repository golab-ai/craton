import math

class nonbond_calculator:
    def __init__(
        self, doing, funcid, x, para, style="energy", combination_rule="pair", scale_factor=1.0
    ):
        """
        calculate non-bonded energy
        include: coul and vdw
        """

        self.funcid = funcid
        self.style = style
        self.x = x
        self.para = para
        self.doing = doing
        self.combination_rule = combination_rule
        self.scale_factor = scale_factor

    def coul_vdw(self):
        """
        coul and vdw
        """
        wow = charge_calculator(self.funcid[0],self.x,self.para[0],style=self.style)
        wow()
        self.charge_value = wow.value
        wow = vdw_calculator(self.funcid[1],self.x,self.para[1],style=self.style,combination_rule=self.combination_rule,scale_factor=self.scale_factor)
        wow()
        self.vdw_value = wow.value
        ####return charge_value, vdw_value

    def coul(self):
        """
        only coul
        """
        wow = charge_calculator(self.funcid[0], self.x,self.para[0],style=self.style)
        charge_value = wow()
        vdw_value = 0
        return charge_value, vdw_value

    def vdw(self):
        """
        only vdw
        """
        charge_value = 0
        wow = vdw_calculator(self.funcid[1],self.x,self.para[1],style=self.style,combination_rule=self.combination_rule,scale_factor=self.scale_factor,)
        vdw_value = wow()
        return charge_value, vdw_value

    __Func_ID = {
        "coul": coul,
        "vdw": vdw,
        "coul_vdw": coul_vdw,
    }

    def __call__(self):
        func = self.__Func_ID[self.doing]
        func(self)


class vdw_calculator:
    """
    """

    def __init__(self, funcid, x, para, style="energy", combination_rule="pair", scale_factor=1.0):
        """
        calculate vdw
        """

        self.funcid = funcid
        self.x = x
        self.style = style
        self.para = para
        self.combination_rule = combination_rule
        self.scale_factor = scale_factor

    def LJ12_6(self):
        """
        LJ12_6 function: v = 4 * espl * [(sigma/r)**12 - (sigma/r)**6]
        """
        r1 = 1.0 / self.x
        r3 = r1**3
        r6 = r3 * r3
        q3 = self.para[0] ** 3
        q6 = q3**2
        v = q6 * r6
        if self.style == "energy":
            self.value =  4 * self.para[1] * (v**2 - v)

    __Func_ID = {
        "LJ12_6": LJ12_6,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)



class charge_calculator:
    """
    coul factor: 332.063714 kcal mol(-1) A e(-2)
    """

    def __init__(self, funcid, x, para, style="energy", charge_type="point", n=0, m=0):
        self.funcid = funcid
        self.x = x
        self.para = para
        self.style = style
        self.charge_type = charge_type
        self.n = n
        self.m = m

    def coul(self):
        """
        coul interaction: v = 332.063714 * q1 * q2 / r
        """
        q = self.para[0] * self.para[1]
        r = 1.0 / self.x
        if self.style == "energy":
            self.value = 332.063714 * q * r

    __Func_ID = {
        "coul": coul,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)