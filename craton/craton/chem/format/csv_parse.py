import csv
from .db_parse import DBData
from .smiles_parse import SmilesData
from ...utils.common.utils import parse_string
from ...utils import logger
from io import StringIO


DEFAULT_MOLECULE_ATTR_TO_CSV = [
    "smiles","inchi_key","inchi","name","mole_name",
    "iupac_name","nick_name","drug_name","internal_name",
    "mass","formula","heavy_atoms","net_charge","multiple","element_set",
    "function_group","function_group_label",
    "torsions","torsion_number","constrain_term","constrain_value",
    "rings","ring_number","ring_size","ring_property","ring_blocks",
    "energy","dipole","inertia","density",
    "source","count","topol_label","frag_type",
]

class CsvData:
    def __init__(self,style=".csv"):
        self.style = style
        pass

    def _parse(self,input_script,extra_var=None):
        script = input_script.splitlines()
        return self.read_file(script,extra_var=extra_var)

    def _convert(self,molecules,extra_var=None):
        if not isinstance(molecules,list):
            molecules = [molecules]
        if extra_var is None:
            extra_var = []
        attributions = ["mole_name","inchi_key","smiles"] if "mole_name" in extra_var \
                else ["inchi_key","smiles"]
        for attr in extra_var:
            if attr in DEFAULT_MOLECULE_ATTR_TO_CSV and attr not in ["inchi_key","smiles","mole_name"]:
                attributions.append(attr)
        text = " ".join(attributions)
        text += "\n"
        for molecule in molecules:
            for attr in attributions:
                if isinstance(getattr(molecule,attr,None),list) or isinstance(getattr(molecule,attr,None),set):
                    if len(getattr(molecule,attr,None)) == 0:
                        text += "None "
                    else:
                        text += "%s " %":".join([str(sss) for sss in getattr(molecule,attr)])
                else:
                    text += "%s " % getattr(molecule,attr,None)
            text += "\n"
        return text

    def read_file(self,input_script,extra_var=None):
        smiles_flag = False
        inchikey_flag = False
        name_flag = False
        input_script = [line for line in input_script if line[0] != "#"]
        if "," in input_script[0]:
            titles = input_script[0].split(",")
        else:
            titles = input_script[0].split()
        titles = [tt.strip(",") for tt in titles]

        if "smiles" in titles or "SMILES" in titles:
            smiles_flag = True
        if "smiles_base" not in titles:
            if "inchi_key" in titles:
                inchikey_flag = True
            if "name" in titles:
                name_flag = True

        dicts = []
        if "," in input_script[0]:
            scripts = "\n".join(input_script)
            f = StringIO(scripts)
            reader = csv.DictReader(f)
            for row in reader:
                dicts.append(row)
        else:
            for line in input_script[1:]:
                ss = line.split()
                tmp = {}
                for ii,key in enumerate(titles):
                    if key != "smiles_base":
                        tmp[key] = ss[ii].strip(",")
                dicts.append(tmp)
        if name_flag:
            names = [line["name"] for line in dicts]
            doc = DBData()
            return doc._parse(names,extra_var=extra_var)
            
        if inchikey_flag:
            inchikeys = [line["inchi_key"] for line in dicts]
            doc = DBData()
            return doc._parse(inchikeys,extra_var=extra_var)

        docs = []
        smilesobj = SmilesData()
        
        for ii,row in enumerate(dicts):
            smiles = row["smiles"].strip()
            #try:
                    #datas = smilesobj._parse(smiles,extra_var=extra_var)
                
            docs.append(smilesobj._parse(smiles,extra_var=extra_var))
            for kk,vv in row.items():
                if kk not in ["smiles","smiles_base"]:
                    values  = vv.split(":")
                    if set(vv) != set(["None"]) and vv != "":
                        docs[-1][kk] = parse_string(values[0]) if len(values) == 1 else parse_string(values)
                   
            #except:
            #    logger.warning(f"{smiles} read smiles error: {ii+1} line in ")
        return docs

    def old_read_file(self,input_script,extra_var=None):
        if "," in input_script[0]:
            titles = input_script[0].split(",")
        else:
            titles = input_script[0].split()
        titles = [tt.strip(",") for tt in titles]
        
        inchikey_flag = True if "inchi_key" in titles else False
        name_flag = True if "name" in titles else False
        inchikey_index = None if not inchikey_flag else titles.index("inchi_key")
        name_index = None if not name_flag else titles.index("name")
        smiles_index = titles.index("smiles") if "smiles" in titles else None
        SMILES_index = titles.index("SMILES") if "SMILES" in titles else None
        smiles_index = smiles_index if smiles_index is not None else SMILES_index
        
        #mole_name_index = None
        if "smiles_base" in titles:
            inchikey_flag = False
            name_flag = False 
        if name_flag:
            names = [line.split()[name_index].strip(",") for line in input_script]
            doc = DBData()
            return doc._parse(names,extra_var=extra_var)
            
        if inchikey_flag:
            inchikeys = [line.split()[inchikey_index].strip(",") for line in input_script[1:]]
            doc = DBData()
            return doc._parse(inchikeys,extra_var=extra_var)

        if smiles_index is None:
            return []

        docs = []
        smilesobj = SmilesData()
        
        for ii,line in enumerate(input_script[1:]):
            
            if line[0] != "#":
                ss = line.split()
                smiles = ss[smiles_index].strip(",")
                try:
                    #datas = smilesobj._parse(smiles,extra_var=extra_var)
                
                    docs.append(smilesobj._parse(smiles,extra_var=extra_var))
                    for jj,title in enumerate(titles):
                        if title not in ["smiles_base"]:
                            values = ss[jj].strip(",").split(":")
                            if set(values) != set(["None"]):
                                docs[-1][title] = parse_string(values[0]) if len(values) == 1 else parse_string(values)
                except:
                    logger.warning(f"read smiles error: {ii+1} line in ")
        return docs
    

