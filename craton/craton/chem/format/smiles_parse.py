"""
测试：
    1.含有手性的smiles
    2. 使用3D结构 和不用3D结构
    3. normalization Tru or False
"""
from ._rdkit import RdkitMol
from .mol_parse import MolData
from ...utils import logger

class SmilesData:
    def __init__(self, style=""):
        self.style = style

    def _parse(self,smiles,extra_var=None):
        structure_3d = extra_var["structure_3d"] if extra_var is not None and "structure_3d" in extra_var \
            else True
        return self._parse_script(smiles,structure_3d)

    def _parse_script(self,smiles,structure_3d):

        self.rdkmol = RdkitMol()
        self.rdkmol._parse_smiles(smiles)
        if structure_3d:
            self.rdkmol._get_3d()
        datas = self._parse_molecule_datas(self.rdkmol._get_script())
        return datas

    def _parse_molecule_datas(self, script):
        mol_data = MolData()
        datas =  mol_data._parse(script,extra_var=None)
        datas["molecule_name"] = self.rdkmol.inchi_key
        datas["inchi_key"] = self.rdkmol.inchi_key
        datas["inchi"] = self.rdkmol.inchi
        datas["smiles"] = self.rdkmol.smiles
        return datas

    def _convert(self,molecule,extra_var=None):
        structure_3d = extra_var["structure_3d"] if extra_var is not None and "structure_3d" in extra_var \
            else False
        normalization = extra_var["normalization"] if extra_var is not None and "normalization" in extra_var \
            else True
        
        return self._convert_molecule(molecule,structure_3d,normalization)

    def _convert_molecule(self, molecule, structure_3d, normalization):
        """
        if structure_3d, the smiles includ "chirality","cis_trans". Usually this is False
        usually, normalization is True. False for figure
        """

        mol_script = MolData("normal")
        script = mol_script._convert(molecule,extra_var = {"structure_3d": structure_3d})

        self.rdkmol = RdkitMol()
        self.rdkmol._convert(script)
        self.smiles = self.rdkmol.smiles
        self.inchi_key = self.rdkmol.inchi_key
        self.inchi = self.rdkmol.inchi

