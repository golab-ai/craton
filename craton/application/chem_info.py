from craton import molxpert as MX
import csv,os
import json
from pathlib import Path
from copy import deepcopy    

class ChemInfo:
    def __init__(self,inputs,output_directory="."):
        self.molecules = MX.molecule_create(inputs)
        self.molecules = MX.molecule_structure(self.molecules)
        self.output_directory = output_directory

    def ring(self):
        for molecule in self.molecules:
            print(molecule.ring_dict)
            MX.molecule_show(molecule,attrs=["ring"],opath=self.output_directory)

    def chiral_atom(self):
        molecules = MX.molecule_chiral(self.molecules)
        text = ""
        chiral_atoms = [[atom for atom in molecule.Atoms if atom.chirality_flag] for molecule in molecules]
        for ii,atoms in enumerate(chiral_atoms):
            text += f"{molecules[ii].mole_name} {molecules[ii].smiles}\n"
            if len(atoms) > 0:
                for atom in atoms:
                    text += f"{atom.ID} {atom.elem} {atom.has_ring}"
            text += f"\n"
        print(text)
        with open(f"{self.output_directory}/chiral_atom.txt",'w') as outf:
            outf.write(text)

        for molecule in molecules:
            MX.molecule_show(molecule,attrs=["chirality_flag"],opath=f"{self.output_directory}")

    def torsion(self):
        molecules = MX.molecule_torsion(self.molecules)
        with open(f"{self.output_directory}/torsion.txt",'w') as outf:
            for molecule in molecules:
                outf.write(f"{molecule.mole_name} {molecule.smiles}\n")
                for torsion in molecule.torsions:
                    outf.write(f"{torsion[0]} {torsion[1]} {torsion[2]} {torsion[3]}\n")
                outf.write("\n\n")
                MX.molecule_show(molecule,attrs=["torsion"],opath=f"{self.output_directory}")
                MX.molecule_show(molecule,attrs=["ID"],opath=f"{self.output_directory}")
    
    def hybrid(self):
        molecules = MX.molecule_hybrid(self.molecules)
        for molecule in molecules:
            MX.molecule_show(molecule,attrs=["hybrid"],opath=f"{self.output_directory}")

    def interaction_site(self):
        molecules = MX.molecule_model(self.molecules)
        with open(f"{self.output_directory}/interaction_model.txt",'w') as outf:
            for molecule in molecules:
                MX.molecule_show(molecule,attrs=["ID"],opath=f"{self.output_directory}")
                outf.write(f"{molecule.mole_name} {molecule.smiles}\n")
                for kk,vv in molecule._interaction_model.items():
                    if kk not in ["hydrophobic"]:
                        outf.write(f"{kk}:\n")
                        for vvv in vv:
                            outf.write(f"{vvv}\n")
                outf.write("\n\n")

    def fragmentation(self):
        pass

    def atom_cluster(self):
        pass

    def function_group(self):
        pass

    def molecule_image(self):
        pass

class MoleEdit:
    def __init__(self,inputs,output_directory="."):
        molecules = MX.molecule_create(inputs)[0]
        self.molecule = MX.molecule_structure(molecules)[0]
        self.output_directory = output_directory

    def calculate_structure_parameter(self,atoms):
        value = MX.structure_calculate(self.molecule,atoms)
        if len(atoms) == 2:
            print(f"the distance between {atoms[0]} and {atoms[1]} is : {value}Å")
        if len(atoms) == 3:
            print(f"the angle between {atoms[0]} {atoms[1]} and {atoms[2]} is : {value} deg")
        if len(atoms) == 4:
            print(f"the torsion between {atoms[0]} {atoms[1]} {atoms[2]} and {atoms[3]} is : {value} deg")

    def change_structure_parameter(self,atoms,value,del_value=False):
        #self.molecule = MX.update_structure_topol(self.molecule)[0]
        value_tmp = MX.structure_calculate(self.molecule,atoms)
        molecule = deepcopy(self.molecule)
        molecule = MX.structure_change(molecule,atoms,value,del_value=del_value)
        molecule = MX.update_structure_topol(molecule)[0]
        value_ch = MX.structure_calculate(molecule,atoms)
        molecule.mole_name = f"{molecule.mole_name}_changed"
        MX.format_convert([molecule,self.molecule],otype="sdf",opath=self.output_directory)
        print(f"{'-'.join([str(an) for an in atoms])}: {value_ch} ({value_tmp})")

    def get_RMSD(self,input2):
        molecule2 = MX.molecule_create(input2)[0]
        print(MX.conformer_RMSD(self.molecule,molecule2))

def molecule_descript(input,fpath):
    molecules = MX.molecule_create(input)
    molecules = MX.molecule_structure(molecules)
    molecules = MX.molecule_torsion(molecules)
    molecules = MX.molecule_hybrid(molecules)
    molecules = MX.molecule_function_group(molecules)
    molecules = MX.function_group(molecules)
    molecules = MX.molecule_model(molecules)
    MX.save_molecule_info(molecules,fpath=fpath)
    os.mkdir(f"{fpath}/PNG")
    MX.format_convert(molecules,otype="png",opath=f"{fpath}/PNG")
    os.mkdir(f"{fpath}/PNG_ID")
    for molecule in molecules:
        MX.molecule_show(molecule,attrs=["ID"],save_file=True,opath=f"{fpath}/PNG_ID",)

def register_molecule(input,logf=None,configure=None):
    if configure is None:
        configure = {}
    config = MX.update_configure(configure)
    #else:
    #    config = MX.
    molecules = MX.molecule_create(input)
    if logf is not None:
        charge_molecules = MX.molecule_create(logf)
        _dicts_ = {molecule.name:molecule for molecule in charge_molecules}
    
        for molecule in molecules:
            charge_molecule = _dicts_[molecule.name]
            for ii,atom in enumerate(molecule.Atoms):
                atom.ff_charge = charge_molecule.Atoms[ii].esp_charge
                
    
    molecules = MX.molecule_structure(molecules)
    atf = config["ForceFieldSetting"]["PROTEIN_TYPING_FILE"]
    fff = config["ForceFieldSetting"]["PROTEIN_FORCE_FIELD_FILE"]
    molecules = MX.atom_type(molecules,atf=atf)
    MX.AA_registered(molecules)

def chem_info(ref,attr_style=0,f_style="json",index_key="inchi_key"):
    total_attrs = {
                0:["IUPAC_Name","CAS_number","inchi_key","formula",
                   "smiles","mass","heavy_atoms","net_charge","torsion_number",
                   "ring_number","ring_size","ring_property","element_count",
                   "zelement","function_group","function_group_label"],
                1:[],
                2:[],
            }
    suffix = Path(ref).suffix    
    attrs = total_attrs[attr_style]
    
    if suffix in ["csv","txt"]:
        tmp_datas = list(csv.reader(open(ref)))
        labels = tmp_datas[0]
    
        datas = [{s:data[ii] for ii,s in enumerate(labels)} for data in tmp_datas[1:]]
    
        smiles_arr = [data["smiles"].strip() for data in datas]
    else:
        smiles_arr = ref
    molecules = MX.molecule_create(smiles_arr,show_figure=False)
    molecules = MX.molecule_structure(molecules)
    molecules = MX.molecule_torsion(molecules)
    molecules = MX.molecule_hybrid(molecules)
    molecules = MX.molecule_function_group(molecules)
    molecules = MX.molecule_model(molecules)
    text = write_txt(molecules[0])
    # print()
    # print(text)
    # print(molecules[0].__dict__.keys())
    # print(molecules[0].function_group)
    # print(molecules[0].torsions)
    if suffix in ["csv","txt"]:
        for ii,m in enumerate(molecules):
            data = datas[ii]
            for aa,bb in data.items():
                if aa not in ["inchi_key","smiles"]:
                    setattr(m,aa.strip(),bb.strip())
    #for molecule in molecules:
    #    for atom in molecule.Atoms:
    #        print(atom.__dict__)
    if f_style == "json":
        docs = {}
        for m in molecules:
            kk = getattr(m,index_key)
            docs[kk] ={}
            for attr in attrs:
                docs[kk][attr] = getattr(m,attr,"-")
        with open("chem_info.json",'w') as outf:
            outf.write(json.dumps(docs))

