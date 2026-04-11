#from ..chemkit.structure import determine_hydrogen_bond
from ..utils import numerical_algorithm


class Ensemble:
    __default_attributs = [
        "system",  # the System class
        "trj",  # 轨迹文件 “coor", "volicity", "force", "lattics"
        "property",
        "thermo",
    ]  # "temperature", "pressure","potential",

    def __init__(self, style):
        self.s = style

    @property
    def density(self):
        if "Density" not in self.thermo:
            self.calculate_density()
        return sum(self.thermo["Density"]) / len(self.thermo["Density"])

    def calculate_density(self):
        self.thermo["Density"] = []
        mass = 0.0
        for aa in self.system.Atoms:
            mass += aa.mass
        for ll in self.trj["lattics"]:
            self.thermo["Density"].append(mass * 10 / 6.02 / ll[0] / ll[1] / ll[2])

    @property
    def volume(self):
        if "volume" not in self.thermo:
            self.calculate_volume()
        return sum(self.thermo["volume"]) / len(self.thermo["volume"])

    def calculate_volume(self):
        self.thermo["volume"] = []
        for ll in self.trj["lattics"]:
            self.thermo["volume"].append(ll[0] * ll[1] * ll[2])

    @property
    def potential(self):
        return sum(self.thermo["Potential"]) / len(self.thermo["Ppotential"])

    @property
    def temperature(self):
        return sum(self.thermo["Temperature"]) / len(self.thermo["Temperature"])

    @property
    def pressure(self):
        return sum(self.thermo["Pressure"]) / len(self.thermo["Ppressure"])

    def get_property(self, pp):
        if pp in self.thermo:
            return sum(self.thermo[pp]) / len(self.thermo[pp])
        elif pp in self.property:
            return sum(self.thermo[pp]) / len(self.thermo[pp])
        else:
            return None

    @property
    def density_profile(self, interval=1):
        den_profile_tmp = {
            "x": {},
            "y": {},
            "z": {},
        }
        for ii in range(len(self.trj["lattics"])):
            lattic = self.trj["lattics"][ii]
            area = {
                "x": interval * lattic[1] * lattic[2],
                "y": lattic[0] * interval * lattic[2],
                "z": lattic[0] * lattic[1] * interval,
            }
            coor_arr = self.trj["coor"][ii]
            den = {"x": {}, "y": {}, "z": {}}
            for jj in range(len(coor_arr)):
                coor = self.trj["coor"][ii][jj]
                mass = self.system.Atoms[jj].mass
                nx = int(coor[0] / interval)
                ny = int(coor[1] / interval)
                nz = int(coor[2] / interval)
                if nx not in den["x"]:
                    den["x"][nx] = 0.0
                if ny not in den["y"]:
                    den["y"][ny] = 0.0
                if nz not in den["z"]:
                    den["z"][nz] = 0.0
                den["x"][nx] += mass
                den["y"][ny] += mass
                den["z"][nz] += mass
            for term in ["x", "y", "z"]:
                for aa, bb in den[term].items():
                    if aa not in den_profile_tmp[term]:
                        den_profile_tmp[term][aa] = [0.0, 0]
                    den_profile_tmp[term][aa][0] += bb * 10 / area[term]
                    den_profile_tmp[term][aa][1] += 1
        den_profile = {
            "x": [[], []],
            "y": [[], []],
            "z": [[], []],
        }
        for term in ["x", "y", "z"]:
            for aa, bb in den_profile_tmp[term].items():
                den_profile[term][0].append(aa)
                den_profile[term][1].append(bb[0] / bb[1])
        for term in ["x", "y", "z"]:
            for aa, bb in den_profile.items():
                bb = sorted(bb, key=lambda x: x[0])
        return den_profile

    def check_equilibrium(self, pp=["temperature", "pressure", "Potential"]):
        equilibrium = {}
        for p in pp:
            equilibrium[p] = [[], []]
            p_arr = self.thermo[p]
            for n, vv in enumerate(p_arr):
                equilibrium[p][0].append(n)
                equilibrium[p][1].append(vv)
                equilibrium[p].append(numerical_algorithm.linear_fitting(equilibrium[p][0], equilibrium[p][1]))
        return equilibrium

    def find_hydrogen_bond(self, water_flag=False):
        hbond = {}
        # ligand_donor = [[],[]]
        # ligand_acceptor = []
        # other_donor = [[],[]]
        # other_acceptor =[]
        donors = []
        acceptors = []
        nn = len(self.system.mole[0].Atoms)
        for ii in range(nn):
            if self.system.Atoms[ii].elem == "H":
                if not water_flag:
                    if self.system.Atoms[ii].mole_type == self.system.mole[0].mole_name:
                        ca = self.system.Atoms[ii].connect[0]
                        if self.system.Atoms[ca].elem in ["O", "N", "S", "P"]:
                            donors.append([ii, ca])
                else:
                    ca = self.system.Atoms[ii].connect[0]
                    if self.system.Atoms[ca].elem in ["O", "N", "S", "P"]:
                        donors.append([ii, ca])
            if self.system.Atoms[ii].elem in ["O", "N", "S", "P"]:
                acceptors.append(ii)
        for fram in self.trj["coor"]:
            for donor in donors:
                for acceptor in acceptors:
                    if acceptor != donor[1]:
                        if determine_hydrogen_bond(fram[acceptor], fram[donor[0]], fram[donor[1]]) == "yes":
                            hbond_key = "%s-%s-%s-%s-%s-%s$%s-%s-%s-%s-%s-%s" % (
                                self.system.Atoms[donor[0]].mole_id,
                                self.system.Atoms[donor[0]].mole_type,
                                self.system.Atoms[donor[0]].residu_number,
                                self.system.Atoms[donor[0]].residu,
                                self.system.Atoms[donor[0]].No_offset,
                                self.system.Atoms[donor[1]].elem,
                                self.system.Atoms[acceptor].mole_id,
                                self.system.Atoms[acceptor].mole_type,
                                self.system.Atoms[acceptor].residu_number,
                                self.system.Atoms[acceptor].residu,
                                self.system.Atoms[acceptor].No_offset,
                                self.system.Atoms[acceptor].elem,
                            )
                            if hbond_key not in hbond:
                                hbond[hbond_key] = 0
                            hbond[hbond_key] += 1
        return hbond


class GEnsemble:
    def __init__(self, style):
        self.style = style
