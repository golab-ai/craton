#!/usr/bin/env python
"""

"""
import inspect
import os
from pathlib import Path

#from glx_cadd.commons import chem_const, ff_logger, mole_analy
#from glx_cadd.commons.app_context import AppContext  # noqa
#from glx_cadd.commons.mole_system import Atom, Molecule, System
#from glx_cadd.commons.third_party.read_file_utils import extra_info, search_line

from ..chem.molecule import Molecule
from ..chem.atom import Atom

current_dir = os.path.dirname(os.path.abspath(os.path.realpath(inspect.getfile(inspect.currentframe()))))
parent_dir = os.path.dirname(current_dir)

INTERMEDIATE_DIR = f"{parent_dir}/intermediate/"
Path(INTERMEDIATE_DIR).mkdir(parents=True, exist_ok=True)


def LmpFF(style,para):
    def bond_harmonic(para):
        return [para[1],para[0]]

    def angle_harmonic(para):
        return [para[1],para[0]]
    
    def dihedral_opls(para):
        return [para[0]*2.0,para[2]*2.0,para[4]*2.0,para[6]*2.0]
    
    def improper_cvff(para):
        return [para[0],1,2]
    
    def pair_lj12_6(para):
        return [para[1],para[0]]
    
    _Func = {
        "atomtype_LJ12_6": pair_lj12_6,
        "bondterm_harmonic":bond_harmonic,
        "angleterm_harmonic":angle_harmonic,
        "dihedralterm_amber":dihedral_opls,
        "dihedralterm_opls":dihedral_opls,
        "improperterm_amber": improper_cvff,
    }

    return _Func[style](para)

class LmpInputFile:
    # default properties: Head_arr, Body_dict

    _name = "Lammps"

    def __init__(self, sm, style="normal"):
        self.style = style
        self.sm = sm
        self.init_system()

    def init_system(self):
        #self.atom_style = self.sm.md_setting["atom_style"]
        self.atom_style = "full"
        self.output_dir = self.sm.output_dir
        self.head_arr = [0 for __ in self.__Head_Keywords]
        self.head_arr[0] = len(self.sm.coordinates)

        self.term_style = {"atomtype":None,"bondterm":None,"angleterm":None,"dihedralterm":None,"improperterm":None}

        for ii,term in enumerate(["Bonds","Angles","Dihedrals","Impropers"]):
            self.head_arr[ii + 1] = sum([len(getattr(molecule,term,[])) * nn for molecule,nn in zip(self.sm.molecules,self.sm.molecule_number)])

        for ii,term in enumerate(["atomtype", "bondterm", "angleterm", "dihedralterm", "improperterm"]):
            if term in self.sm.ff:
                self.head_arr[ii+9] = len(self.sm.ff[term])
                nn = 0
                for item,vv in self.sm.ff[term].items():
                    nn += 1
                    vv["type_id"] = nn
                self.term_style[term] = list(set([vv["pstyle"] for item,vv in self.sm.ff[term].items()]))
        if len(self.sm.lattics) == 9:
            self.head_arr[19:28] = self.sm.lattics 
        elif len(self.sm.lattics) == 6:
            self.head_arr[19:25] = self.sm.lattics
        elif len(self.sm.lattics) == 3:
            self.head_arr[19] = 0.00000
            self.head_arr[20] = self.sm.lattics[0]
            self.head_arr[21] = 0.00000
            self.head_arr[22] = self.sm.lattics[1]
            self.head_arr[23] = 0.00000
            self.head_arr[24] = self.sm.lattics[2]     

    __Head_Keywords = [
        "atoms",
        "bonds",
        "angles",
        "dihedrals",
        "impropers",

        "ellipsoids",
        "lines",
        "triangles",
        "bodies",

        "atom types",
        "bond types",
        "angle types",
        "dihedral types",
        "improper types",

        "extra bond per atom",
        "extra angle per atom",
        "extra_dihedral per atom",
        "extra improper per atom",
        "extra special per atom",
        
        "xlo",
        "xhi",
        "ylo",
        "yhi",
        "zlo",
        "zhi",
        "xy",
        "xz",
        "yz",
    ]
    __Body_Keywords = [
        "Atoms",
        "Bonds",
        "Angles",
        "Dihedrals",
        "Impropers",
        
        "Ellipsoids",
        "Lines",
        "Triangles",
        "Bodies",
        
        "Velocities",

        "Masses",
        "Pair Coeffs",
        "Bond Coeffs",
        "Angle Coeffs",
        "Dihedral Coeffs",
        "Improper Coeffs",
        "PairIJ Coeffs",
        "BondBond Coeffs",
        "BondAngle Coeffs",
        "MiddleBondTorsion Coeffs",
        "EndBondTorsion Coeffs",
        "AngleTorsion Coeffs",
        "AngleAngleTorsion Coeffs",
        "BondBond13 Coeffs",
        "AngleAngle Coeffs",
    ]
    __Atom_Style = {
        "angle": ["atom-ID", "molecule-ID", "atom-type", "x", "y", "z"],
        "atomic": ["atom-ID", "atom-type", "x", "y", "z"],
        "body": ["atom-ID", "atom-type", "bodyflag", "mass", "x", "y", "z"],
        "bond": ["atom-ID", "molecule-ID", "atom-type", "x", "y", "z"],
        "charge": ["atom-ID", "atom-type", "q", "x", "y", "z"],
        "dipole": ["atom-ID", "atom-type", "q", "x", "y", "z", "mux", "muy", "muz"],
        "electron": ["atom-ID", "atom-type", "q", "spin", "eradius", "x", "y", "z"],
        "ellipsoid": ["atom-ID", "atom-type", "ellipsoidflag", "density", "x", "y", "z"],
        "full": ["atom-ID", "molecule-ID", "atom-type", "q", "x", "y", "z"],
        "line": ["atom-ID", "molecule-ID", "atom-type", "lineflag", "denstiy", "x", "y", "z"],
        "meso": ["atom-ID", "atom-type", "rho", "e", "cv", "x", "y", "z"],
        "molecular": ["atom-ID", "molecule-ID", "atom-type", "x", "y", "z"],
        "peri": ["atom-ID", "atom-type", "volume", "density", "x", "y", "z"],
        "sphere": ["atom-ID", "atom-type", "diameter", "density", "x", "y", "z"],
        "template": ["atom-ID", "molecule-ID", "template-index", "template-atom", "atom-type", "x", "y", "z"],
        "tri": ["atom-ID", "molecule-ID", "atom-type", "triangleflag", "density", "x", "y", "z"],
        "wavepacket": ["atom-ID", "atom-type", "charge", "spin", "eradius", "etag", "cs_re", "cs_im", "x", "y", "z"],
        "hybrid": ["atom-ID", "atom-type", "x", "y", "z"],
    }
    # __Body_Keywords_index=[5,0,0,15,16,17,18,1,2,3,4,5,5,6,7,8,9,7,7,8,8,8,8,6,7]

    def write_data_head(self):
        text = ""
        for i in range(0, 19):
            if self.head_arr[i] > 0:
                text += "%8d %s\n" % (self.head_arr[i], self.__Head_Keywords[i])
        text += "\n"
        text += "%12.4f%12.4f   %s %s\n"% (self.head_arr[19], self.head_arr[20], self.__Head_Keywords[19], self.__Head_Keywords[20])
        text += "%12.4f%12.4f   %s %s\n"% (self.head_arr[21], self.head_arr[22], self.__Head_Keywords[21], self.__Head_Keywords[22])
        text += "%12.4f%12.4f   %s %s\n"% (self.head_arr[23], self.head_arr[24], self.__Head_Keywords[23], self.__Head_Keywords[24])

        if self.head_arr[25] != 0.0 or self.head_arr[26] != 0.0 or self.head_arr[27] != 0.0:
            text += "%12.4f%12.4f%12.4f   %s %s %s\n" % (self.head_arr[25],self.head_arr[26],self.head_arr[27],
                                                        self.__Head_Keywords[25],self.__Head_Keywords[26],self.__Head_Keywords[27])

        text += "\n"
        return text

    def write_data_structure(self):
        __label = {"Atoms":[0,"atomtype"],
                   "Bonds":["a1","a2","bondterm"],
                   "Angles":["a1","a2","a3","angleterm"],
                   "Dihedrals":["a1","a2","a3","a4","dihedralterm"],
                   "Impropers":["a1","a2","a3","a4","improperterm"]
                   }
        text  = ["","","","","","","","",""]
        IDs = [0,0,0,0,0,0,0,0,0,0]
        coords = self.sm.coordinates
        start_id = 1
        for nn,molecule in enumerate(self.sm.molecules):
            for ii in range(self.sm.molecule_number[nn]):
                IDs[-1] += 1
                for jj,term in enumerate(["Atoms","Bonds","Angles","Dihedrals","Impropers"]):
                    if term == "Atoms":
                        for atom in getattr(molecule,term,[]):
                            IDs[jj] += 1
                            if self.atom_style == "full":
                                text[jj] += "%10d %10d %10d %10.3f %15.3f %15.3f %15.3f\n" %(IDs[jj],IDs[-1], 
                                                                                             self.sm.ff[__label[term][-1]][atom.atom_type_used_name]["type_id"],
                                                                                             atom.point_charge,*coords[IDs[jj]-1]
                                                                                             )

                    else:
                        for item in getattr(molecule,term,[]):
                            IDs[jj] += 1
                            text[jj] += "%10d %10d "%(IDs[jj],self.sm.ff[__label[term][-1]][item.atom_type_used_name]["type_id"])
                            for ans in __label[term][:-1]:
                                text[jj] += "%10d " %(getattr(item,ans) + start_id)
                            text[jj] += "\n"
                start_id += len(molecule.Atoms)
        return text

    def write_data_parameter(self):
        
        text = "Masses\n\n"
        for at,vv in self.sm.ff["atomtype"].items():
            text += "%5d %8.3f # %s\n" % (vv["type_id"],vv["mass"],at)
        text += "\n"
        for ii,term in enumerate(["atomtype","bondterm","angleterm","dihedralterm","improperterm"]):
            if term in self.sm.ff.keys():
                text += "%s\n\n" %self.__Body_Keywords[ii+11]
                for at,vv in self.sm.ff[term].items():
                    text += "%5d "%vv["type_id"]
                    para = LmpFF(f"{term}_{vv['pstyle']}",vv["parameter"])
                    for pp in para:
                        if isinstance(pp,int):
                            text += "%15d" %pp
                        else:
                            text += "%15.3f "%pp
                    text += "\n"
                    
                text += "\n"
        return text 

    def write_data(self,file_name="lmp.data"):
        
        outf = open(f"{self.output_dir}/{file_name}", "w")
        outf.write("# cpy generated LAMMPS data file\n\n")

        ## write head section
        outf.write(self.write_data_head())
        
        ####write structure section
        
        texts = self.write_data_structure()
        for ii,text in enumerate(texts):
            if text != "":
                outf.write("%s\n\n"%self.__Body_Keywords[ii])
                outf.write(text)
                outf.write("\n")
       
        ###write parameters section
        outf.write(self.write_data_parameter())
        
        outf.close()

    def write_in_para_style(self):
        text = ""
        __label = {
                   "atomtype":{
                        "style":"pair_style",
                        "LJ12_6":"lj/cut/coul/long",
                        },
                   "bondterm":{
                       "style":"bond_style",
                       "harmonic":"harmonic",
                       "morse":"morse",
                       "class2":"class2",
                       },
                   "angleterm":{
                       "style":"angle_style",
                       "harmonic":"harmonic",
                       "class2":"class2",
                       },
                   "dihedralterm":{
                       "style":"dihedral_style",
                       "amber":"opls",
                       "fourier":"fourier",
                       "opls":"opls",
                       },
                   "improperterm":{
                       "style":"improper_style",
                       "amber":"cvff",
                       },
                   }
        for term,typ in self.term_style.items():
            if typ is None:
                text += "%s none\n" %__label[term]["style"]
            else:
                if len(typ) == 1:
                    text += "%s %s" %(__label[term]["style"],__label[term][typ[0]])
                    if term == "atomtype":
                        text += "    12.0"
                    text += "\n"
        text += "pair_modify        mix arithmetic\n"
        text += "pair_modify   tail yes\n"
        text += "kspace_style  pppm 1.0e-4\n"
        text += "special_bonds  lj 0.0 0.0 0.5 coul 0.0 0.0 0.8333\n"
        text += "dielectric     1.0\n\n"
        return text

    def write_in_jobs(self):
        text = ""
        velocity_flag = True
        for ii,job in enumerate(self.sm.md_setting["jobs"]):
            
            if job in ["mini","opt"]:
                text += "min_style       cg\n"
                text += f'minimize 1.0e-4 1.0e-6 100 {self.sm.md_setting["nsteps"][ii]}\n'
            elif job in ["eq_nvt"]:
                text += f'fix {job} all nvt temp {self.sm.md_setting["temperature"]["temperature"]} {self.sm.md_setting["temperature"]["temperature"]} 100.0\n'
                text += f"timestep        1.0\n"
                if velocity_flag:
                    text += f'velocity       all create {self.sm.md_setting["temperature"]["temperature"]} 1166140691 mom yes rot yes dist gaussian\n'
                    velocity_flag = False
                text += f'run      {self.sm.md_setting["nsteps"][ii]}\n'
                text += f"unfix  {job}\n"
            elif job in ["eq_npt"]:
                text += f'fix {job} all npt temp {self.sm.md_setting["temperature"]["temperature"]} {self.sm.md_setting["temperature"]["temperature"]} 100.0 '
                text += f'iso {self.sm.md_setting["pressure"]["pressure"]} {self.sm.md_setting["pressure"]["pressure"]} 1000.0\n'
                text += "timestep        1.0\n"
                text += f'run      {self.sm.md_setting["nsteps"][ii]}\n'
                text += f"unfix {job}\n"
            elif job in ["prod_npt"]:
                text += f'fix {job} all npt temp {self.sm.md_setting["temperature"]["temperature"]} {self.sm.md_setting["temperature"]["temperature"]} 100.0 '
                text += f'iso {self.sm.md_setting["pressure"]["pressure"]} {self.sm.md_setting["pressure"]["pressure"]} 1000.0\n'
                text += "timestep        1.0\n"
            elif job in ["prod_nvt"]:
                text += f'fix {job} all nvt temp {self.sm.md_setting["temperature"]["temperature"]} {self.sm.md_setting["temperature"]["temperature"]} 100.0\n'
                text += "timestep        1.0\n"
                text += f'run      {self.sm.md_setting["nsteps"][ii]}\n'
        text += "\n"
        return text
        
    def write_in_property(self):
        _nn = int(self.sm.md_setting["nsteps"][-1] / 10)
        block_size = _nn if _nn <= 1000000 else 1000000
        n_repeat = 100
        n_every = int(block_size / n_repeat)

        text = "variable  T equal temp\n"
        text += "variable P equal press\n"
        text += "variable V equal vol\n"
        text += "variable nmole equal %d\n" %sum(self.sm.molecule_number)
        text += "variable Tref equal %.3f\n" %self.sm.md_setting["temperature"]["temperature"]
        text += "variable Pref equal %.3f\n" %self.sm.md_setting["pressure"]["pressure"]
        text += "variable lz equal lz\n"
        # text += "group mole molecule <> 1 %d\n" %sum(self.sm.molecule_number)
        text += "group mole id 1:%d\n" %sum(self.sm.molecule_number)
        #text += "compute inter all inter \n"

        flag_0 = True
        flag_1 = True

        if len(set(["hov", "cp","kt","ap","heat capacity","compressibility","thermal expansion"]) & set(self.sm.md_setting["property"])) > 0:
            text += "compute pair_en all pe pair\n"
        for pp in self.sm.md_setting["property"]:
            if pp == "den":
                text += f"compute   atom_mass all property/atom mass\n"
                text += f"compute   total_mass all reduce sum c_atom_mass   # 单位 g/mol，数值上等于总分子量\n"
                text += f"variable  density equal (c_total_mass / 6.022e23) / (vol * 1e-24)   # g/cm^3\n"
                text += f"fix denout all ave/time {n_every} {n_repeat} {block_size} v_density file denout.log\n"
            elif pp in ["hov"]:
                text += f"variable hov equal -c_pair_en/v_nmole+8.314*v_Tref/4184\n"
                text += f"fix hovout all ave/time {n_every} {n_repeat} {block_size} v_hov file hov.log"
            elif pp in ["cp","kt","ap","heat capacity","compressibility","thermal expansion"]:
                if flag_0:
                    text += "variable Hconf equal c_pair_en+1.01325*0.0014388*v_Pref*v_V\n" 
                    text += "variable VHconf equal v_V*v_Hconf\n"
                    text += "variable UHconf equal c_pair_en*v_Hconf\n"

                    ####cp
                    text += f"variable cpfactor1 equal 2105587.7/v_Tref/v_Tref\n" 
                    text += f"variable cpfactor2 equal v_Pref*30.3188*1.01325/v_Tref/v_Tref\n" 
                    text += f"variable cpfactor3 equal v_nmole\n" 
                    text += f"fix cpout all ave/time 5 1 5 c_pair_en v_Hconf v_V v_cpfactor1 v_cpfactor2 v_cpfactor3 file cp.log\n"
                
                    #### kt
                    text += f"variable          VV equal v_V*v_V\n" 
                    text += f"variable          ktfactor equal 0.072464/v_Tref\n" 
                    text += f"fix               ktout all ave/time 5 1 5 v_V v_ktfactor file kt.log\n"

                    #### ap
                    text += f"variable          apfactor equal 503.2475/v_Tref/v_Tref\n" 
                    text += f"fix               apout all ave/time 5 1 5 v_V v_Hconf v_apfactor file ap.log\n"
                    flag_0 = False

            elif pp == "rdf":
                text += "compute           11rdf  all rdf 100 1 1\n" 
                text += f"fix               rdfout all ave/time {n_every} {n_repeat} {block_size} c_11rdf file rdf_$NAME$_$TEMP$_$PRESS$.log mode vector\n" 
            elif pp == "rg":
                text += "compute           chunkrg all chunk/atom molecule nchunk once ids once\n" 
                text += "compute           rg all gyration/chunk chunkrg\n"
                text += f"fix               rgout all ave/time {n_every} {n_repeat} {block_size} c_rg file rg_$NAME$_$TEMP$_$PRESS$.log mode vector\n"

            elif pp in ["dc"]:
                text += f"compute           chunkmsd all chunk/atom molecule nchunk once ids once\n" 
                text += f"compute           msd all msd/chunk chunkmsd\n"
                text += f"fix               msdout all ave/time 100 1 100 c_msd[4] file dc.log mode vector\n"
                #text += "#compute           vacf all vacf\n"
                #text += "#fix               vacfout all ave/time 100 1 100 c_vacf[4] file vacf.log\n"
            elif pp in ["viscosity","vis"]:
                text += f"variable          visfactor equal 0.000000001*0.74397*{n_every}*v_V/v_Tref\n"
                text += "variable          pxy equal pxy\n"
                text += "variable          pyz equal pyz\n"
                text += "variable          pxz equal pxz\n"
                text += "fix               nonPout all ave/time 5 1 5 v_visfactor  v_pxy v_pyz v_pxz file nonp.log\n"
                text += f"fix               visacfout1 all ave/correlate {n_every} {n_repeat} {block_size} v_pxy v_pxy type upper file visacf1.log\n" 
                text += f"fix               visacfout2 all ave/correlate {n_every} {n_repeat} {block_size} v_pyz v_pyz type upper file visacf2.log\n" 
                text += f"fix               visacfout3 all ave/correlate {n_every} {n_repeat} {block_size} v_pxz v_pxz type upper file visacf3.log\n" 
            elif pp in ["td"]:
                text += "variable          tdfactor equal 10000000*3.50035*100*v_V/v_Tref/v_Tref\n"
                text += "compute           ke all ke/atom\n" 
                text += "compute           pe all pe/atom\n" 
                text += "compute           stress all stress/atom NULL virial\n" 
                text += "compute           hflux all heat/flux ke pe stress\n" 

                text += "fix               hfluxout all ave/time 5 1 5 v_tdfactor c_hflux[1] c_hflux[2] c_hflux[3] file hflux.log\n"
                text += f"fix               tdacfout1 all ave/correlate {n_every} {n_repeat} {block_size} c_hflux[1] c_hflux[1] type upper file tdacf1.log\n" 
                text += f"fix               tdacfout2 all ave/correlate {n_every} {n_repeat} {block_size} c_hflux[2] c_hflux[2] type upper file tdacf2.log\n" 
                text += f"fix               tdacfout3 all ave/correlate {n_every} {n_repeat} {block_size} c_hflux[3] c_hflux[3] type upper file tdacf3.log\n" 
            elif pp in ["surface tension","st"]:
                text += f"variable          st equal 0.01*(lz/2.0)*(pzz-(pxx+pyy)/2.0)*1.01325\n"
                text += f"fix               stout all ave/time {n_every} {n_repeat} {block_size} v_st file st.log\n"
            elif pp in ["vle","ct","cd","cb","svp","nbp"]:
                if flag_1:
                    text += "variable          svp equal pzz\n"
                    text += "variable          st equal 0.01*(lz/2.0)*(pzz-(pxx+pyy)/2.0)*1.01325\n"
                    text += f"fix               stout all ave/time {n_every} {n_repeat} {block_size} v_st file st_$NAME$_$TEMP$.log\n"
                    text += f"fix               svpout all ave/time {n_every} {n_repeat} {block_size} v_svp file svp_$NAME$_$TEMP$.log\n"
                    text += f"fix               vleout all ave/chunk {n_every} {n_repeat} {block_size} z lower 1.0 density/mass ave running file vle_$NAME$_$TEMP$.log\n"
                    falg_1 = False

        text += "\n"
        return text

    def write_in_output(self):
        self.get_atom_type_element()
        text = ""
        text += f'thermo          {self.sm.md_setting["output"]["nstenergy"]}\n'
        text += f"thermo_style    custom step press etotal ke pe emol ebond eangle edihed eimp epair ecoul evdwl\n"
        text += f'dump            1 all custom {self.sm.md_setting["output"]["nstxout-compressed"]} dump.lmptrajectory id type x y z element\n'
        text += f"dump_modify     1 sort id \n"   
        text += f'dump_modify     1 element {" ".join(self.atom_type_element)}\n'
        text += f'run             {self.sm.md_setting["nsteps"][-1]}\n'
        text += "\n"
        return text
    
    def get_atom_type_element(self):
        self.atom_type_element = []
        atomtypes_sorted = sorted(self.sm.ff["atomtype"].items(), key=lambda x: x[1]["type_id"])
        for at_name, _ in atomtypes_sorted:
            # 从名称中提取元素符号：
            # 简单方法：取第一个字符并大写（如 o_2n -> O, n_3o -> N, h_1s -> H）
            # 特殊处理：cg -> CG (或者保持 C，但根据实际)
            if at_name.startswith("cg"):
                element = "C"
            else:
                element = at_name[0].upper()
            self.atom_type_element.append(element)

    def write_in(self,file_name="lmp.in"):
        outf = open(f"{self.output_dir}/{file_name}", "w")
        outf.write("# cpy generated LAMMPS data file\n\n")
        outf.write("units         real\n")
        outf.write("boundary      p p p\n")
        outf.write("atom_style   full\n")
        outf.write("\n")
        outf.write(self.write_in_para_style())
        
        outf.write("read_data       lmp.data\n\n")
        outf.write(self.write_in_jobs())
        outf.write(self.write_in_property())
        outf.write(self.write_in_output())
        outf.close()

class LmpDumpFile:
    def __init__(self, style):
        self.s = style

    def read_script(self, script):
        self.script = script

    def read_info(self):
        self.lattic = []
        self.timestep = []
        self.coor = []
        self.force = []
        self.vel = []
        for i in range(len(self.script)):
            if self.script[i].strip() == "ITEM: NUMBER OF ATOMS":
                break
        atom_n = int(self.script[i + 1].strip())
        record_key = {}
        tmp = self.script[i + 6].strip().split()[2:]
        for i in range(len(tmp)):
            record_key[tmp[i]] = i
        for i in range(0, len(self.script), atom_n + 9):
            if i % (atom_n + 9) == 1:
                self.timestep.append(int(self.script[i]))
                if "x" in record_key.keys():
                    self.coor.append([None for i in range(atom_n)])
                if "fx" in record_key.keys():
                    self.force.append([None for i in range(atom_n)])
                if "vx" in record_key.keys():
                    self.vel.append([None for i in range(atom_n)])
            elif i % (atom_n + 9) == 5:
                for j in range(3):
                    self.lattic.append(self.script[i + j].strip().split())
            elif i % (atom_n + 9) >= 9:
                ss = self.script[i].strip().split()
                n = int(ss[record_key["id"]])
                if "x" in record_key.keys():
                    self.coor[-1][n - 1] = [
                        float(ss[record_key["x"]]),
                        float(ss[record_key["y"]]),
                        float(ss[record_key["z"]]),
                    ]
                if "fx" in record_key.keys():
                    self.force[-1][n - 1] = [
                        float(ss[record_key["fx"]]),
                        float(ss[record_key["fy"]]),
                        float(ss[record_key["fz"]]),
                    ]
                if "vx" in record_key.keys():
                    self.force[-1][n - 1] = [
                        float(ss[record_key["vx"]]),
                        float(ss[record_key["vy"]]),
                        float(ss[record_key["vz"]]),
                    ]

class LmpOutputFile:
    def __init__(self, style):
        self.s = style

    def read_script(self, script):
        self.script = script
        self.read_info()

    def read_info(self):
        lp = []
        kws = []
        self.lmp_thermo = []
        for i in range(0, len(self.script)):
            if self.script[i][:4] == "Step":
                string = self.script[i].split()
                lp.append(i)
                kws.append(string)
            if self.script[i][:9] == "Loop time":
                lp.append(i)
        for i in range(0, len(lp), 2):
            for j in range(lp[i] + 1, lp[i + 1]):
                ss = self.script[j].split()
                datas = {}
                kw = kws[int(i / 2)]
                for jj in range(len(kw)):
                    datas[kw[jj]] = float(ss[jj])
                self.lmp_thermo.append(datas)
