import os,sys
from pathlib import Path
import csv

from craton import molxpert as MX

file_path = Path(__file__).resolve()
base_dir = file_path.parent.parent
py_file = f"{base_dir}/craton/dock/pocket_analyzer.py"
##tmp_dir = f"{base_dir}/tmp"

class MolDock:
    def __init__(self,protein_file,ligand_file,center=None,box_size=None,charge_method="binc", dock_atf=None,output_directory=".",configure=None):
        """
        __init__ 的 Docstring
        
        :param protein_file: 蛋白的pdb文件
        :param ligand_file: 配体的sdf文件
        :param charge_method: 配体的电荷方法
        :param dock_atf: 配体的原子类型
        :param center: 对接的中心
        :param box_size: 对接的box size
        """
        self.user_config = configure
        if self.user_config is not None:
            self.config = MX.update_configure(self.user_config)
        else:
            self.config = MX.update_configure({})
        
        self.protein_file = protein_file
        self.ligand_file = ligand_file
        if center is None:
            self.center = None
            self.box_size = None
        else:
            self.center = center
            if box_size is None:
                self.box_size = [30.0, 30.0, 30.0]
            else:
                self.box_size = box_size

        self.charge_method = charge_method
        self.dock_atf = self.config["ForceFieldSetting"]["DOCK_ATOM_TYPING_FILE"]
        self.output_directory = output_directory

        Path(output_directory).mkdir(exist_ok=True)
        

    @staticmethod
    def pocket(protein_file,output_directory="."):
        Path(output_directory).mkdir(exist_ok=True)
        protein_file_name = protein_file.split("/")[-1]
        pre_file_name = protein_file_name[:-4]
        csv_file = f"{output_directory}/{pre_file_name}_pocket.csv"
        cavity_file = f"{output_directory}/{pre_file_name}_pocket_cavity.pdb"
        os.system(f"python {py_file} {protein_file} -o {csv_file} --cavity-pdb {cavity_file}")
        return csv_file
        
        #return MX.pocket_analyze(protein_file,csv_file_flag=True,cavity_file_flag=True)
    
    def assign_ring_for_protein(self,protein):
        for atom in protein.Atoms:
            atom.has_ring = []
            atom.has_ring_property=[]
            atom.has_ring_size =[]
        for kk,vv in protein.rings.items():
            size=len(vv)-1
            ring_p = vv[-1]
            for an in vv[:-1]:
                protein.Atoms[an].has_ring.append(kk)
                protein.Atoms[an].has_ring_property.append(ring_p)
                protein.Atoms[an].has_ring_size.append(size)

    def reassign_torsion(self,ligand):
        torsions = []
        for term in ligand.torsions:
            elems = [ligand.Atoms[an].elem for an in term[1:3]]
            if set(elems) == set(["N","C"]):
                ii = elems.index("C") + 1
                if "2" in ligand.Atoms[term[ii]].bond_type:
                    jj = ligand.Atoms[term[ii]].bond_type.index("2")
                    ann = ligand.Atoms[term[ii]].connect[jj]
                    if ligand.Atoms[ann].elem in ["O"]:
                        continue
            torsions.append(term)
        ligand.torsions = torsions

    def read_coor_from_pdbqt(self,pdbqt_file):
        coors = []
        with open(pdbqt_file) as outf:
            lines = outf.readlines()
        for line in lines:
            if line[:4] == "ATOM":
                coors.append([float(line[30:38].strip()),float(line[38:46].strip()),float(line[46:54].strip())])
            elif line[:6] == "ENDMDL":
                break
        return coors

    def vina_prepare(self):
        
        protein = MX.molecule_create(self.protein_file)[0]
        protein = MX.molecule_structure(protein)[0]
        protein = MX.protein_structure(protein)
        self.assign_ring_for_protein(protein)
        protein.scan_term = []
        ligands = MX.molecule_create(self.ligand_file)
        ligands = MX.molecule_structure(ligands)
        ligands = MX.molecule_torsion(ligands)
        for ligand in ligands:
            self.reassign_torsion(ligand)
        if self.charge_method == "binc":
            ligands = MX.atom_type(ligands)
            ligands = MX.grasp_force_field(ligands)
        elif self.charge_method == "am1bcc":
            ligands = MX._am1bcc_charge(ligands)
        

        ligands = MX.atom_type(ligands,atf=self.dock_atf,ignore_existing=True)
        protein = MX.atom_type(protein,atf=self.dock_atf,
                               assign_atom_type_flag=False,
                               ignore_existing=True,
                               ignore_ff_existing=True,
                               parallel=False)[0]

        if self.center is None:
            if not hasattr(protein,"pocket_center"):
                csv_file = MolDock.pocket(self.protein_file,output_directory=self.output_directory)
                with open(csv_file) as inf:
                    reader = csv.DictReader(inf)
                    for row in reader:
                        center = [float(row["center_x"]),float(row["center_y"]),float(row["center_z"])]
                        box_size = [float(row["size_x"]),float(row["size_y"]),float(row["size_z"])]
                        break
                protein.pocket_center = center
                protein.pocket_box = box_size
        else:
            protein.pocket_center = self.center
            protein.pocket_box = self.box_size
        
        Path(f"{self.output_directory}/input").mkdir(exist_ok=True)
        Path(f"{self.output_directory}/output").mkdir(exist_ok=True)
        protein_pdbqt_file = f"{protein.mole_name}".split("/")[-1]

        ligands_mappings = []
        self.ligands_files = []
        for ii, ligand in enumerate(ligands):
            ligand_infos = MX.pdbqt_file(ligand,extra_var={"ignore_connect":True})
            ligand_pdbqt_file = f"{self.output_directory}/input/{ligand.mole_name}_{ii}.pdbqt"
            ligand_pdbqt_file_output_head = f"{self.output_directory}/output/{ligand.mole_name}_{ii}"
            self.ligands_files.append([ligand_pdbqt_file_output_head, ligand_pdbqt_file])
            with open(ligand_pdbqt_file,'w') as outf:
                outf.write(ligand_infos[0])
            ligands_mappings.append(ligand_infos[1])

        self.protein_script = MX.format_convert(protein,otype="pdbqt",
                                                ofilename=protein_pdbqt_file,
                                                opath=f"{self.output_directory}/input",extra_var={"ignore_connect":True})
        
        self.protein_files = f"{self.output_directory}/input/{protein_pdbqt_file}.pdbqt"
        
        self.protein = protein
        self.ligands = ligands
        self.ligands_mappings = ligands_mappings

    def run_vina(self):
        
        self.vina_prepare()
        all_results = MX.molecule_docking(self.protein_files,
                                                  self.ligands_files,
                                                  self.protein.pocket_center,
                                                  self.protein.pocket_box,
                                                  output_directory=self.output_directory,
                                            )
        mole_name = self.protein.mole_name.split("/")[-1].replace(".pdb", "")
        with open(f"{self.output_directory}/{mole_name}_score.csv",'w') as outf:
            outf.write("name, smiles, score\n")
            for ii,results in enumerate(all_results):
                score = results
                ligand = self.ligands[ii]
                atom_mappings = self.ligands_mappings[ii]
                coors = self.read_coor_from_pdbqt(f"{self.ligands_files[ii][0]}_docked.pdbqt")
                ligand.associated_data = {"smiles":ligand.smiles,"score":score}
                for jj,atom in enumerate(ligand.Atoms):
                    
                    atom.coor = coors[atom_mappings[jj]-1]
                outf.write(f"{ligand.mole_name}, {ligand.smiles}, {score}\n")
        MX.format_convert(self.ligands,otype="sdf",ofilename=f"{self.protein.mole_name}_docked",
                          opath=self.output_directory)