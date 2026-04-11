class XyzData:
    def __init__(self):
        self.elem = []
        self.coor = []

    def _parse(self,input_script,extra_var=None):
        data = {"molecule_name"  : "",
            "atom_count"     : 0,
            "elements"       : [],
            "coordinates"    : [],
            }
        for line in input_script.splitlines()[2:]:
            ss = line.strip().split()
            data["elements"].append(ss[0])
            data["coordinates"].append([float(ss[1]), float(ss[2]), float(ss[3])])
        data["atom_count"] = len(data["elements"])
        return data

    def _convert(self,molecule,extra_var=None):
        text = "%d\n" % len(molecule.Atoms)
        text += "craton create xyz file\n"
        for atom in molecule.Atoms:
            text += "%s %.3f %.3f %.3f \n" %(atom.element,atom.coordinates[0],atom.coordinates[1],atom.coordinates[2])
        return text

