import csv
from vina import Vina

import argparse,os,datetime

class VinaRun:
    def __init__(self,protein,ligand,center=None,box_size=None,output_pdbqt=None,file_input=False,):
        self.protein = protein
        self.ligand = ligand
        self.center = center
        self.box_size = box_size
        #self.file_input = file_input
        self.output_pdbqt = output_pdbqt

    def run_vina(self):
    
        # 1. 初始化 Vina 引擎（指定 CPU 核心数和搜索详尽度）
        v = Vina(sf_name='vina', cpu=4)

        # 2. 设置受体和配体
        #if self.file_input:
        v.set_receptor(self.protein)
        v.set_ligand_from_file(self.ligand)

        #else:
        #v.set_receptor_from_string(self.protein)
        #v.set_ligand_from_string(self.ligand)

        # 3. 定义对接盒子 (Center & Size)
        v.compute_vina_maps(center=self.center, box_size=self.box_size)

        # 4. 执行对接
        v.dock(exhaustiveness=8, n_poses=9)

        # 5. 输出结果
        if self.output_pdbqt is not None:
            v.write_poses(self.output_pdbqt, n_poses=5, overwrite=True)

        # 6. 获取得分
        energies = v.energies()
        #coors = v.poses(n_poses=1,coordinates_only=True)
        #for ii,coor in enumerate(coors):
        #    self.molecule.Atoms[ii].coor = coor
        #print(f"最高结合亲和力: {energies[0][0]:.2f} kcal/mol")
        return f"{energies[0][0]:.2f}"








def run():
    parser = argparse.ArgumentParser(description="run docking")
    parser.add_argument("--protein_file", type=str, default="protein.pdb", help="pdb file")
    parser.add_argument("--ligand_file", type=str, default="ligand.sdf", help="ligand file")
    parser.add_argument("--configure_file", type=str, default="pocket.csv", help="configure file")
    parser.add_argument("--results_dir", type=str, default="./Docking", help="output directory")
    args = parser.parse_args()
    protein_file = args.protein_file
    ligand_file = args.ligand_file
    configure_file = args.configure_file
    results_dir = args.results_dir
    os.makedirs(f"{results_dir}/input",exist_ok=True)
    os.makedirs(f"{results_dir}/pdbqt_output",exist_ok=True)
    os.makedirs(f"{results_dir}/output",exist_ok=True)
    rec_pdbqt = f"{results_dir}/receptor.pdbqt"
    pybel_convert_to_pdbqt(protein_file,rec_pdbqt,is_receptor=True)
    center,box_size = read_config(configure_file)

    arrs = read_sdf_file(ligand_file)
    infos = []
    for arr in arrs:
        ligand_pdbqt = f"{results_dir}/input/{arr[0]}.pdbqt"
        output_pdbqt = f"{results_dir}/pdbqt_output/{arr[0]}_docking.pdbqt"
        meeko_convert_sdf_to_pdbqt(arr[1], ligand_pdbqt)
        score = run_vina(rec_pdbqt,ligand_pdbqt,output_pdbqt,center,box_size)
        output_sdf = f"{results_dir}/output/{arr[0]}_docking.sdf"
        convert_pdbqt_to_sdf(output_pdbqt,output_sdf)
        with open(output_sdf) as inf:
            lines = inf.readlines()
        lines[0] = f"{arr[0]}\n"
        lines.append(f"><score>")
        lines.append(f"{score}")
        scripts = "".join(lines)
        infos.append([score,scripts])
    infos = sorted(infos,key = lambda x:x[0], reverse=True)
    total_scripts = "$$$$\n".join([s[1] for s in infos])
    ligand_name = os.path.basename(ligand_file)
    ligand_base = os.path.splitext(ligand_name)[0]
    total_output_sdf = f"{results_dir}/output/{ligand_base}_docking.sdf"
    with open(total_output_sdf, 'w') as outf:
        outf.write(total_scripts)

if __name__ == "__main__":
    run()






