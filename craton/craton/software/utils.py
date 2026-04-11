import string

def set_gromacs_atom_info(self, resname="UNK"):
    """
    生成gromacs所需要的特殊的性质
    ###TODO: 该方法需要UT###
    输入：
    输出：
    """
    self.residu_n = 1
    for i, atom in enumerate(self.Atoms):
        atom.residu = resname
        atom.residue_ID = 1
        atom.charge_group = i + 1


