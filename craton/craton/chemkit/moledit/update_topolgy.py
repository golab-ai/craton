from copy import deepcopy

def reset_atom(molecule1,molecule2,attrs):
    from ...chem.molecule import Molecule
    from ...chem.constants import ATOM_DEFAULT_ATTRIBUTES

    if isinstance(molecule2,Molecule):
        for attr in attrs:
            if attr in ATOM_DEFAULT_ATTRIBUTES:
                if hasattr(molecule2.Atoms[0],attr):
                    for ii,atom in enumerate(molecule1.Atoms):
                        setattr(atom,attr,getattr(molecule2.Atoms[ii],attr))
                else:
                    if hasattr(molecule1.Atoms[0],attr):
                        for ii,atom in enumerate(molecule1.Atoms):
                            delattr(atom,attr)
            else:
                if hasattr(molecule2,attr):
                    setattr(molecule1,attr,getattr(molecule2,attr))
                else:
                    if hasattr(molecule1,attr):
                        delattr(molecule1,attr)
    elif isinstance(molecule2,dict):
        for ii,atom in enumerate(molecule1.Atoms):
            for attr in attrs:
                setattr(atom,attr,molecule2[attr][ii])
    elif isinstance(molecule2,list):
        for ii,atom in enumerate(molecule1.Atoms):
            for jj,attr in enumerate(attrs):
                setattr(atom,attr,molecule2[jj][ii])

def template_molecule_topolgy(molecules1,molecules2,match_key=None):
    from ...chem.molecule import Molecule
    """
    通过molecules2补全或更新molecules1中缺失的拓扑结构.
    通常情况下，molecules1缺失部分topolgy信息，但我们想保留molecules1中的3D信息（原子坐标）
    如果molecules2不是一个数组，而是一个Molecule,所有molecules1中分子以相同的molecules2为模板
    match_key指定根据什么方式匹配，
        None: one by one
        string: 必须是inchi key, canonic smiles, inchi, molecule name等能分子唯一标志性质 
    """
    def _update_atom(molecule1,molecule2):
        reset_atom(molecule1,molecule2,["connectivity","bond_type","formal_charge"])
        molecule.create_topols()
        molecule.steps = ["topol"]

    if not isinstance(molecules1,list):
        molecules1 = [molecules1]

    if isinstance(molecules2,Molecule):
        for molecule in molecules1:
            _update_atom(molecule,molecules2)
        return
    if match_key is None:
        for ii,molecule in enumerate(molecules1):
            _update_atom(molecule,molecules2[ii])
    else:
        tmp_dict = {getattr(mm,match_key):mm for mm in molecules2}
        for molecule in molecules1:
            _update_atom(molecule,tmp_dict[getattr(molecule,match_key)])


def expand_conformers(molecules1,molecules2,attrs=["coordinates"],match_key=None):
    from ...chem.molecule import Molecule
    """
    以molecules1为模板，扩展更多的molecules2中的构象
    通常情况下，molecules2中有很多不同的分子构象(3D坐标)，为相同的分子提供更多的构象
    molecules可以不是Molecule的数组，而是坐标的集合，如：
    [[coordinates0],[coordinates1],[coordinates2],......]
    {
        "molecule_name0":[[coordinates0],[coordinates1],[coordinates2],......],
        "molecule_name1":[[coordinates0],[coordinates1],[coordinates2],......],
    }
    coordinates = [[x0,y0,z0],[x1,y1,z1],[x2,y2,z2],......]
    如果molecules1不是一个数组，而是一个Molecule,所有molecules2中分子以相同的molecules1为模板
    match_key指定根据什么方式匹配，
        None: one by one
        string: 必须是inchi key, canonic smiles, inchi, molecule name等能分子唯一标志性质 
    """
    def get_molecules_type(molecules2):
        if isinstance(molecules2,dict):
            return "DICT"
        elif isinstance(molecules2,list):
            if isinstance(molecules2[0],Molecule):
                return "MOL_LIST"
            else:
                return "ARR_LIST"
    
    def _gen_new_conformer(molecule,mm2):
        this_molecules = []
        for m in mm2:
            this_molecules.append(deepcopy(molecule))
            if not isinstance(m,Molecule):
                m = [m]
            reset_atom(this_molecules[-1],m,attrs)
        return this_molecules

    if isinstance(molecules1,Molecule):
        molecules1 = [molecules1]

    if isinstance(molecules2,Molecule):
        molecules2 = [molecules2]

    _m2_type = get_molecules_type(molecules2)

    new_molecules = []
    if _m2_type == "DICT":
        for molecule in molecules1:
            new_molecules.extend(_gen_new_conformer(molecule,molecules2[getattr(molecule,match_key)]))
    elif _m2_type == "MOL_LIST":
        for molecule in molecules1:
            new_molecules.extend(_gen_new_conformer(molecule,[m for m in molecules2 if getattr(m,match_key) == getattr(molecule,match_key)]))
    elif _m2_type == "ARR_LIST":
        for molecule in molecules1:
            new_molecules.extend(_gen_new_conformer(molecule,molecules2))
    
    return new_molecules
  

