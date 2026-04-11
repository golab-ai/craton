from rdkit import Chem

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AllChem import AlignMol

from ...utils import logger

class RdkitMol:

    def __init__(self,style="") -> None:
        self.style = style

    def _parse_smiles(self, string):
        self.smiles = AllChem.CanonSmiles(string)
        self.rdkm = Chem.MolFromSmiles(self.smiles)
        self.rdkmh = Chem.AddHs(self.rdkm)
        self._parse_inchi()

    def _parse_inchi(self):
        #self.smiles = Chem.MolToSmiles(self.rdkm)
        self.inchi = Chem.inchi.MolToInchi(self.rdkm)
        self.inchi_key = Chem.inchi.InchiToInchiKey(self.inchi)

    def _get_3d(self):
        AllChem.EmbedMolecule(self.rdkmh)

    def _get_script(self):
        return Chem.MolToMolBlock(self.rdkmh)

    def _convert(self, script,normalization=True,draw_image=False,TD_flag=False):
        _rdkm = Chem.MolFromMolBlock(script)
        Chem.rdmolops.AssignAtomChiralTagsFromStructure(_rdkm)
        _smiles = Chem.MolToSmiles(_rdkm)
        smiles = AllChem.CanonSmiles(_smiles)

        if normalization:
            try:
                self.rdkm = Chem.MolFromSmiles(smiles)
                self.rdkmh = Chem.AddHs(self.rdkm)
                self.smiles = smiles
            except:  # noqa
                logger.error("rdkit can not explanate the caonsmiles %s which create by rdkit" % self.smiles)
                logger.error("we use the new smiles %s" % _smiles)
                self.rdkm = Chem.MolFromSmiles(_smiles)
                self.rdkmh = Chem.AddHs(self.rdkm)
                self.smiles = _smiles 
        else:
            self.rdkmh = Chem.MolFromMolBlock(script, removeHs=False)
            self.rdkm = _rdkm
            self.siles = smiles
        
        self._parse_inchi()
        if draw_image:
            if not TD_flag:
                AllChem.Compute2DCoords(self.rdkm)

    def _get_RMSD(self,script1, script2):
        
        rdkm1 = Chem.MolFromMolBlock(script1)
        rdkm2 = Chem.MolFromMolBlock(script2)

        rmsd = AlignMol(rdkm1, rdkm2, -1, -1)
        return rmsd

    def draw_figure(self, fname=None, mole_name=None,highlighttype=None, highlightarrs=None, removeH_flag=True):
        if mole_name is None:
            mole_name = self.inchi_key
        
        if highlighttype is None:
            imp = Chem.Draw.MolToImage(self.rdkm, legend=mole_name)
        elif highlighttype == "atoms":
            if removeH_flag:
                imp = Chem.Draw.MolToImage(self.rdkm, highlightAtoms=highlightarrs, legend=mole_name)
            else:
                imp = Chem.Draw.MolToImage(self.rdkmh, highlightAtoms=highlightarrs, legend=mole_name)
        elif highlighttype == "bonds":
            imp = Chem.Draw.MolToImage(self.rdkm, highlightBonds=highlightarrs, legend=mole_name)
        elif highlighttype == "label":
            drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
            AllChem.Compute2DCoords(self.rdkm, canonOrient=True)
            for i in range(len(self.rdkm.GetAtoms())):
                self.rdkm.GetAtomWithIdx(i).SetProp("molAtomMapNumber", highlightarrs[i])
            drawer.DrawMolecule(self.rdkm,legend=mole_name)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
        
        if fname is not None:
            if highlighttype == "label":
                with open(fname, "w") as outf:
                    outf.write(svg)
            else:
                imp.save(fname)
        return [svg,"svg"] if highlighttype == "label" else [imp,"img"]
    
