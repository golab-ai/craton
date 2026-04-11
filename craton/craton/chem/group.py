
class Group:
    def __init__(self, style="residue",group_str=None,group_name=None,group_idx=None,atoms=None,group_chain_name=None,net_charge=0):
        self.style = style
        self.group_str = group_str
        self.group_name = group_name
        self.group_idx = group_idx
        self.atoms = atoms
        self.group_chain_name = group_chain_name
        self.net_charge = net_charge

