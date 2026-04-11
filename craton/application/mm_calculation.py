from craton import molxpert as MX
from copy import deepcopy
from pathlib import Path

class MMCalculation:
    
    def __init__(self,inputs,prop="energy",output_directory="."):
        self.molecules = MX.molecule_create(inputs)
        self.molecules = MX.molecule_structure(self.molecules)
        self.molecules = MX.atom_type(self.molecules)
        self.molecules = MX.grasp_force_field(self.molecules)
        self.molecules = MX.update_structure_topol(self.molecules)
        self.prop = prop
        self.output_directory = output_directory
        Path(self.output_directory).mkdir(exist_ok=True)

    def energy_calculate(self):
        molecules =  MX.energy(self.molecules,prop=self.prop)
        MX.format_convert(molecules,otype="mtx",opath=f"{self.output_directory}/energy",extra_var="all")

    
    def molecule_optimize(self):
        MX.format_convert(self.molecules,otype="sdf",opath=f"{self.output_directory}/origin_structure")
        molecules = MX._optimize(deepcopy(self.molecules))
        MX.format_convert(molecules,otype="sdf",opath=f"{self.output_directory}/optimized_structure")
        with open(f"{self.output_directory}/rmds.csv",'w') as outf:
            outf.write("name,smiles,rmsd\n")
            for ii,molecule in enumerate(molecules):
                outf.write(f"{molecule.mole_name},{molecule.smiles},{MX.conformer_RMSD(self.molecules[ii],molecule)}\n")

    def molecule_optimize_with_fix(self,atoms,value=None):
        pass    


    def torsion_scan(self,show_figure=True,vs_mm=True,vs_qm=False,base_qm=False):
        molecules = MX.molecule_torsion(self.molecules)
        qm_datas = {}
        if vs_qm:
            inchi_keys = [m.inchi_key for m in molecules]
            qm_molecules = deepcopy(self.molecules)
            if len(qm_molecules) > 0:
                qm_scan_curve = MX.scan_curve(qm_molecules)
                if qm_scan_curve:
                    qm_datas = MX.scan_curve_data(qm_scan_curve)
    
        mm_datas = {}
        if vs_mm:
            if base_qm:
                if qm_datas:
                    mm_molecules = [deepcopy(molecule) for molecule in qm_molecules]
                    mm_molecules = MX._optimize(mm_molecules)
                    mm_molecules = MX.energy(mm_molecules)
                    mm_scan_curve = MX.scan_curve(mm_molecules)
                    mm_datas = MX.scan_curve_data(mm_scan_curve)
            else:
                mm_molecules = deepcopy(molecules)
                mm_scan_curve = MX._torsion_scan(mm_molecules)
                mm_datas = MX.scan_curve_data(mm_scan_curve)

        qm_mm_data = {name:[qm_datas[name],mm_datas[name]] for name in mm_datas if name in qm_datas}
        only_mm_data = {name:data for name,data in mm_datas.items() if name not in qm_datas}
        only_qm_data = {name:data for name,data in qm_datas.items() if name not in mm_datas}

        qm_mm_args = {"labels":["QM","MM"],"rmse":True,"fitting_curve":True,"xylabels":["Dihedral","Energy"],"save_path":self.output_directory}
        only_mm_args = {"labels":["MM"],"rmse":True,"fitting_curve":True,"xylabels":["Dihedral","Energy"],"save_path":self.output_directory}
        only_qm_args = {"labels":["QM"],"rmse":True,"fitting_curve":True,"xylabels":["Dihedral","Energy"],"save_path":self.output_directory}
        datas = [[name,qm_mm_args,data] for name,data in qm_mm_data.items()]
        datas.extend([[name,only_mm_args,data] for name,data in only_mm_data.items()])
        datas.extend([[name,only_qm_args,data] for name,data in only_qm_data.items()])
        qm_datas_vs = {}
        for kk,vv in qm_datas.items():
            vs = [[angle,vv[1][ii]] for ii,angle in enumerate(vv[0])]
            vs = sorted(vs,key=lambda x:x[0])
            qm_datas_vs[kk] = vs

        mm_datas_vs = {}
        for kk,vv in mm_datas.items():
            vs = [[angle,vv[1][ii]] for ii,angle in enumerate(vv[0])]
            vs = sorted(vs,key=lambda x:x[0])
            mm_datas_vs[kk] = vs

        for rr in datas:
            rr[1]["name"] = rr[0]
            MX.figure_show(rr[2],"pes",rr[1])
        if vs_mm:
            if vs_qm:
                with open(f"{self.output_directory}/torsion_scan.txt",'w') as outf:
                    for kk,vv in qm_datas_vs.items():
                        ss = kk.split("_")
                        outf.write(f"{ss[0]} {ss[1]}\n")
                        outf.write("dihedral qm mm\n")
                        mm_vv = mm_datas_vs[kk]
                        for ii,kv in enumerate(vv):
                            outf.write(f"{kv[0]} {kv[1]} {mm_vv[ii][1]}\n")
                        outf.write("\n\n")
                return mm_datas, qm_datas
            else:
                with open(f"{self.output_directory}/torsion_scan.txt",'w') as outf:
                    for kk,vv in mm_datas_vs.items():
                        ss = kk.split("_")
                        outf.write(f"{ss[0]} {ss[1]}\n")
                        outf.write("dihedral mm\n")
                        for kv in vv:
                            outf.write(f"{kv[0]} {kv[1]}\n")
                        outf.write("\n\n")
                return mm_datas
        else:
            with open(f"{self.output_directory}/torsion_scan.txt",'w') as outf:
                    for kk,vv in qm_datas_vs.items():
                        ss = kk.split("_")
                        outf.write(f"{ss[0]} {ss[1]}\n")
                        outf.write("dihedral qm\n")
                        for ii,kv in enumerate(vv):
                            outf.write(f"{kv[0]} {kv[1]}\n")
                        outf.write("\n\n")
            return qm_datas
        
    def volume_surface(self):
        molecules = MX.molecule_volume_surface(self.molecules)
        with open(f"{self.output_directory}/volume_surface.csv",'w') as outf:
            outf.write("name, smiles, volume(Å³), surface(Å²)\n")
            for molecule in molecules:
                outf.write(f"{molecule.mole_name}, {molecule.smiles}, {molecule.volume}, {molecule.surface}\n")
    
    def multipole_moment(self):
        if self.prop not in ["quadrupole","dipole","octupole"]:
            self.prop = "multipole"
        molecules = MX.molecule_multipole_moment(self.molecules)
        if self.prop == "multipole":
            props = ["dipole","quadrupole","octupole"]
        else:
            props = [self.prop]
        
        with open(f"{self.output_directory}/{self.prop}.txt",'w') as outf:
            for molecule in molecules:
                outf.write(f"{molecule.mole_name} {molecule.smiles}\n")
                if "dipole" in props:
                    outf.write("\n偶极矩 (单位: Debye):\n")
                    outf.write(f"  |μ| = {molecule.dipole:.6f}\n")
                    outf.write(f"  μx = {molecule.dipole_moment[0]:.6f}\n")
                    outf.write(f"  μy = {molecule.dipole_moment[1]:.6f}\n")
                    outf.write(f"  μz = {molecule.dipole_moment[2]:.6f}\n")
                if "quadrupole" in props:
                    outf.write("\n四极矩特征值 (单位: Buckingham):\n")
                    outf.write(f"  λ1 = {molecule.quadrupole[0]:.6f}\n")
                    outf.write(f"  λ2 = {molecule.quadrupole[1]:.6f}\n")
                    outf.write(f"  λ3 = {molecule.quadrupole[2]:.6f}\n")
                    outf.write("\n四极矩张量 (单位: Buckingham):\n")
                    outf.write("  Q_ij =\n")
                    for i in range(3):
                        outf.write("    [" + "  ".join([f"{molecule.quadrupole_moment[i, j]:10.4f}" for j in range(3)]) + "]\n")
                if "octupole" in props:
                    outf.write("\n八极矩张量 (单位: e * Å^3):\n")
                    outf.write("  O_ijk =\n")
                    for i in range(3):
                        outf.write(f"  i={i}:\n")
                        for j in range(3):
                            outf.write("    [" + "  ".join([f"{molecule.octupole_moment[i, j, k]:10.4f}" for k in range(3)]) + "]\n")
                outf.write("\n\n")   
        
    def inertia_moment(self):
        molecules = MX.molecule_inertia(self.molecules)
        with open(f"{self.output_directory}/inertia.csv",'w') as outf:
            outf.write("name, smiles, A(cm⁻¹), B(cm⁻¹), C(cm⁻¹), I_a(amu·Å²), I_b(amu·Å²), I_c(amu·Å²),\n")
            for molecule in molecules:
                outf.write(f"{molecule.mole_name}, {molecule.smiles}, {molecule.inertia_cm[0]: .6f}, {molecule.inertia_cm[1]: .6f}, {molecule.inertia_cm[2]: .6f}, {molecule.inertia[0]: .6f}, {molecule.inertia[1]: .6f}, {molecule.inertia[2]: .6f}\n")
    
    def center(self):
        if self.prop not in ["cog","com","cob","size"]:
            self.prop = "center"
        molecules = MX.molecule_center(self.molecules)
        if self.prop == "center":
            props = ["cog","com","cob","size"]
        else:
            props = [self.prop]

        with open(f"{self.output_directory}/center.txt",'w') as outf:
            for molecule in molecules:
                outf.write(f"{molecule.mole_name} {molecule.smiles}\n")
                if "cog" in props:
                    outf.write(f"几何中心：{molecule.cog[0]: .4f} {molecule.cog[1]: .4f} {molecule.cog[2]: .4f}\n")
                if "com" in props:
                    outf.write(f"质心：{molecule.com[0]: .4f} {molecule.com[1]: .4f} {molecule.com[2]: .4f}\n")
                if "cob" in props:
                    outf.write(f"盒子中心：{molecule.cob[0]: .4f} {molecule.cob[1]: .4f} {molecule.cob[2]: .4f}\n")
                if "size" in props:
                    outf.write(f"分子大小：{molecule.size[0][1]: .4f}-{molecule.size[0][0]: .4f} {molecule.size[1][1]: .4f}-{molecule.size[1][0]: .4f} {molecule.size[2][1]: .4f}-{molecule.size[2][0]: .4f}\n")
                outf.write("\n\n")
