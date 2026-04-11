import math

class bond_calculator:
    """
    """

    def __init__(self, funcid, x, para, style="energy"):
        """
        """

        self.funcid = funcid
        self.x = x
        self.para = para
        self.style = style

    def harmonic(self):
        """
        harmonic形式的计算: v = k * (r - r0)**2
        """
        if self.style == "energy":
            self.value = self.para[1] * (self.x - self.para[0])**2
            #return self.para[1] * (self.x - self.para[0])**2

    def _harmonic_(self):
        """
        harmonic形式的计算: v = k * (r - r0)**2
        """
        if self.style == "energy":
            if self.outtype == "value":
                self.value = self.para[1] * (self.x - self.para[0]) ** 2
            elif self.outtype == "operator":
                self.value = [[self.x**2, 1, -2 * self.x], [[1], [0, 0, 1], [0, 1]]]

    __Func_ID = {
        "harmonic": harmonic,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)


class angle_calculator:
    """
    """

    def __init__(self, funcid, x, para, style="energy"):
        """
        """

        self.funcid = funcid
        self.x = x
        self.para = para
        self.style = style

    def harmonic(self):
        """
        harmonic function: v = k * (radian - radian0)**2
        """
        if self.style == "energy":
            self.value = self.para[1] * (math.radians(self.x) - math.radians(self.para[0])) ** 2

    def _harmonic_(self):
        """
        harmonic形式的计算
        """

        if self.style == "energy":
            if self.outtype == "value":
                self.value = self.para[1] * (self.x - math.radians(self.para[0])) ** 2
            elif self.outtype == "operator":
                self.value = [[self.x**2, 1, -2 * self.x], [[1], [0, 0, 1], [0, 1]]]

    __Func_ID = {
        "harmonic": harmonic,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)


class dihedral_calculator:
    """
    """

    def __init__(self, funcid, x, para, style="energy"):
        """
        """

        self.funcid = funcid
        self.x = x
        self.para = para
        self.style = style

    def fourier_expansion(self):
        """
        fourier_expansion funciton (funcid = fourier, opls,amber)
        v = k1*(1+cos(angle-phase1)) + k2*(1+cos(2*angle-phase2)) + k3*(1+cos(3*angle-phase3)) + k4*(1+cos(4*angle-phase4))
        """
        if self.style == 'energy':
            value = 0
            for i in range(0,8,2):
                if self.para[i] is not None:
                    value += self.para[i]*(1 + math.cos((i / 2 + 1)*math.radians(self.x) - math.radians(self.para[i + 1])))
            self.value = value


    def fourier_expansion_new(self):
        """
        fourier_expansion形式的计算。包括 funcid = fourier, opls, amber
        v = k1*(1+cos(angle-phase1)) + k2*(1+cos(2*angle-phase2)) + k3*(1+cos(3*angle-phase3)) + k4*(1+cos(4*angle-phase4))
        """
        if self.style == "energy":
            if self.outtype == "value":
                value = 0
                for i in range(0,8,2):
                    if self.para[i] is not None:
                        value += self.para[i]*(1 + math.cos((i / 2 + 1)*self.x - self.para[i + 1]))
                        #value += self.para[i] * (1 + ((-1) ** i) * math.cos((i + 1) * self.x))
                self.value = value
            elif self.outtype == "operator":
                tmp = []
                for i in range(4):
                    if self.para[i] is not None:
                        tmp.append(1 + ((-1) ** i) * math.cos((i + 1) * self.x))
                    else:
                        tmp.append(0)
                self.value = [tmp, [[0], [1], [2], [3]]]

    def _fourier_expansion_(self):
        """
        fourier_expansion形式的计算。包括 funcid = fourier, opls, amber
        v = k1*(1+cos(angle-0.0)) + k2*(1+cos(2*angle-180.0)) + k3*(1+cos(3*angle-0.0)) + k4*(1+cos(4*angle-180.0))
        """
        if self.style == "energy":
            if self.outtype == "value":
                value = 0
                for i in range(4):
                    if self.para[i] is not None:
                        value += self.para[i] * (1 + ((-1) ** i) * math.cos((i + 1) * self.x))
                self.value = value
            elif self.outtype == "operator":
                tmp = []
                for i in range(4):
                    if self.para[i] is not None:
                        tmp.append(1 + ((-1) ** i) * math.cos((i + 1) * self.x))
                    else:
                        tmp.append(0)
                self.value = [tmp, [[0], [1], [2], [3]]]

    __Func_ID = {
        "fourier": fourier_expansion,
        "opls": fourier_expansion,
        "amber": fourier_expansion,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)


class improper_calculator:
    """
    """

    def __init__(self, funcid, x, para, style="energy"):
        """
        """

        self.funcid = funcid
        self.x = x
        self.para = para
        self.style = style

    def fourier_2n(self):
        """
        fourier_2n funciton(funcid = fourier, opls, amber,cfvv):
            v = k * (1 - cos(2 * angle))
        """
        if self.style == "energy":
            self.value = self.para[0] * (1 - math.cos(2 * math.radians(self.x)))

    def harmonic(self):
        """
        harmonic function (funcid = harmoinc,charmm):
            v = k * (angle)**2
        """
        if self.style == "energy":
            self.value = self.para[0] * (math.radians(self.x))**2   

    def _fourier_2n_(self):
        """
        fourier_2n形式的计算。包括 funcid = fourier, opls, amber,cfvv
        """
        if self.style == "energy":
            if self.outtype == "value":
                self.value = self.para[0] * (1 - math.cos(2 * self.x))
            elif self.outtype == "operator":
                self.value = [[1 - math.cos(2 * self.x)], [[0]]]

    def _harmonic_(self):
        """
        harmonic形式的计算。包括 funcid = harmoinc,charmm
        """
        if self.style == "energy":
            if self.outtype == "value":
                self.value = self.para[0] * (self.x) ** 2
            elif self.outtype == "operator":
                self.value = [[self.x**2], [[0]]]

    __Func_ID = {
        "fourier": fourier_2n,
        "opls": fourier_2n,
        "amber": fourier_2n,
        "harmonic": harmonic,
        "cvff": fourier_2n,
        "charmm": harmonic,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)


class constrain_calculator:
    """
    constrain计算器.优化结构时起效，能量性质计算时不起作用
    """

    def __init__(self, funcid, x, para, k_para="default"):
        """
            funcid: constrain的类型，如bond,angle,dihedral等
            x: float, 作用项的数值.根据funcid不同而不同
            para: List[float], 该作用项的参数
            k_para: constrain k值的大小
        """

        self.funcid = funcid
        self.x = x
        self.para = para

    def calc_bond(self):
        """
        constrain一个键
        """
        self.value = 10000 * (self.x - self.para) ** 2

    def calc_angle(self):
        """
        constrain一个角
        """
        xx = self.para * math.pi / 180.0
        self.value = 5000 * (math.radians(self.x) - xx) ** 2

    def calc_dihedral(self):
        """
        constrain一个二面角
        """
        xx = self.para * math.pi / 180.0
        self.value = 1500 * (1 - math.cos(math.radians(self.x) - xx))

    __Func_ID = {
        "bond": calc_bond,
        "angle": calc_angle,
        "dihedral": calc_dihedral,
    }

    def __call__(self):
        func = self.__Func_ID[self.funcid]
        func(self)
