from typing import Any
from .elements import Element,get_elem_property

class Atom:
    def __init__(self, style="atom", aid=None):
        self.s = style
        self.formal_charge = 0
        if aid is not None:
            self.ID = aid

    @property
    def vdw_radius(self):
        return Element.get(self._symbol).vdw_radius
    
    @property
    def valent_radius(self):
        return Element.get(self._symbol).valent_radius
    
    @property
    def doubel_bond_radius(self):
        return Element.get(self._symbol).double_bond_radius
    
    @property
    def triple_bond_radius(self):
        return Element.get(self._symbol).triple_bond_radius

    @property
    def elem(self):
        return self._symbol

    @elem.setter
    def elem(self,x):
        self._symbol = x
        self._mass = Element.get(x).mass
        self._number = Element.get(x).number

    @property
    def element(self):
        return self._symbol

    @element.setter
    def element(self,x):
        self._symbol = x
        self._mass = Element.get(x).mass
        self._number = Element.get(x).number

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self,x):
        self._symbol = x
        self._mass = Element.get(x).mass
        self._number = Element.get(x).number

    @property
    def elements(self):
        return self._symbol
    
    @elements.setter
    def elements(self,x):
        self._symbol = x
        self._mass = Element.get(x).mass
        self._number = Element.get(x).number

    
    @property
    def number(self):
        return self._number

    @number.setter
    def number(self,x):
        self._symbol = Element.get(x).symbol
        self._mass = Element.get(x).mass
        self._number = x

    @property
    def atom_number(self):
        return self._number

    @atom_number.setter
    def atom_number(self,x):
        self._symbol = Element.get(x).symbol
        self._mass = Element.get(x).mass
        self._number = x

    @property
    def atom_mass(self):
        return self._mass

    @atom_mass.setter
    def atom_mass(self,x):
        self._symbol = Element.get_by_mass(x).symbol
        self._mass = x
        self._number = Element.get_by_mass(x).number

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self,x):
        #self._symbol = Element.get_by_mass(x).symbol
        self._mass = x
        #self._number = Element.get_by_mass(x).number

    #def __setattr__(self, __name: str, __value: Any) -> None:
    #    if __name in ["elem","element","symbol","elements"]:
    #        self._symbol = __value
    #        self._mass = Element.get(__value).mass
    #        self._number = Element.get(__value).number
    #
    #    elif __name in ["number","atom_number"]:
    #        self._symbol = Element.get(__value).symbol
    #        self._mass = Element.get(__value).mass
    #        self._number = __value
    #
    #    elif __name in ["mass","atom_mass"]:
    #        self._symbol = Element.get_by_mass(__value).symbol
    #        self._mass = __value
    #        self._number = Element.get_by_mass(__value).number
    #    else:
    #        object.__setattr__(self, __name, __value)

    #def __getattribute__(self, __name: str) -> Any:
    #    if __name in ["elem","element","symbol","elements"]:
    #        return self._symbol
    #    elif __name in ["atom_number","number"]:
    #        return self._number
    #    elif __name in ["mass","atom_mass"]:
    #        return self._mass
    #    else:
    #        return object.__getattribute__(self, __name)


    #######暂时用到#####################
    @property
    def No(self):
        return self.ID
    
    @No.setter
    def No(self,i):
        self.ID = i

    @property
    def primitive_formal_charge(self):
        return self.partial_formal_charge
    
    @primitive_formal_charge.setter
    def primitive_formal_charge(self,f):
        self.partial_formal_charge  = f

    @property
    def mole_ID(self):
        return self.molecule_ID

    @mole_ID.setter
    def mole_ID(self,s):
        self.molecule_ID = s

    @property
    def mole_type(self):
        return self.molecule_type

    @mole_type.setter
    def mole_type(self,s):
        self.molecule_type = s

    @property
    def residu(self):
        return self.residue

    @residu.setter
    def residu(self,s):
        self.residue = s

    @property
    def residu_number(self):
        return self.residue_ID

    @residu_number.setter
    def residu_number(self,i):
        self.residue_ID = i 

    @property
    def ff_charge(self):
        return self.point_charge

    @ff_charge.setter
    def ff_charge(self,f):
        self.point_charge = f

    @property
    def ff_charge_base(self):
        return self.point_charge_base

    @ff_charge_base.setter
    def ff_charge_base(self,f):
        self.point_charge_base = f

    @property
    def binc_tag(self):
        return self.binc_str

    @binc_tag.setter
    def binc_tag(self,s):
        self.binc_str = s

    @property
    def binc_para(self):
        return self.binc_parameters

    @binc_para.setter
    def binc_para(self,a):
        self.binc_parameters = a

    @property
    def style(self):
        return self.pstyle

    @style.setter
    def style(self,s):
        self.pstyle = s

    @property
    def para(self):
        return self.parameter

    @para.setter
    def para(self,s):
        self.parameter = s

    @property
    def tag(self):
        return self.vdw_str

    @tag.setter
    def tag(self,s):
        self.vdw_str = s

    @property
    def pscore(self):
        return self.vdw_score

    @pscore.setter
    def pscore(self,f):
        self.vdw_score = f

    @property
    def inring(self):
        return self.has_ring

    @inring.setter
    def inring(self,a):
        self.has_ring = a

    @property
    def ring_size(self):
        return self.has_ring_size

    @ring_size.setter
    def ring_size(self,a):
        self.has_ring_size = a

    @property
    def ring_prop(self):
        return self.has_ring_property

    @ring_prop.setter
    def ring_prop(self,a):
        self.has_ring_property = a

    @property
    def coor(self):
        return self.coordinates

    @coor.setter
    def coor(self,a):
        self.coordinates = a

    @property
    def connect(self):
        return self.connectivity

    @connect.setter
    def connect(self,a):
        self.connectivity = a

    @property
    def bond_type_detail(self):
        return self.connectivity_type

    @bond_type_detail.setter
    def bond_type_detail(self,a):
        self.connectivity = a

    @property
    def atom_name(self):
        if hasattr(self,"_name"):
            return self._name
        else:
            self._name = f"{self.element}{self.ID}"
    
    @atom_name.setter
    def atom_name(self,a):
        self._name = a


    @property
    def name(self):
        if hasattr(self,"_name"):
            return self._name
        else:
            self._name = f"{self.element}{self.ID}"
    
    @name.setter
    def name(self,a):
        self._name = a

    @property
    def atom_type(self):
        return [self.atom_type_name]
    
    @property
    def atom_type_used(self):
        return [self.nonb_atom_type]

    @property
    def atom_type_names(self):
        return [self.atom_type_name]
    
    @property
    def atom_type_used_names(self):
        return [self.nonb_atom_type]

    @property
    def atom_type_used_name(self):
        return self.nonb_atom_type
            
    @atom_type_used_name.setter
    def atom_type_used_name(self,ss):
        pass


    def to_dict(self):
        return {
            "residue_index": self.residue_ID,
            "residue_name": self.residue,
            "name": self.atom_name,
            "atom_index": self.ID,
        }
        
    #bond_type, bond_type_old

    ###########################
    #def __repr__(self):
    #    return f"<Atom: {self.elem} {self.No}>"


class VS:
    def __init__(self, style, center_atom, connect_atoms, psets):
        self.style = style
        self.center_atom = center_atom
        self.connect_atoms = connect_atoms
        self.psets = psets
        self.coor = []

    def calc_value(self, this_coor):
        pass

class CG:
    def __init__(self):
        pass