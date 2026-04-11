import shutil
import subprocess

from ..utils import logger


packmol_path = shutil.which("packmol")
#WATER_XYZ = share_template.water_xyz
#CL_XYZ = share_template.cl_xyz
#NA_XYZ = share_template.na_xyz
#AION_XYZ = template_dir / "ion/AION.xyz"



class Packmol:
    def __init__(self,):
        pass

    @staticmethod
    def packmol(sections,fn,output_dir,pbc=None):
        packmol_input = f"{output_dir}/{fn}.inp"
        packmol_output = f"{fn}.xyz"
        with open(packmol_input, "w") as inf:
            inf.write("seed 0\n")
            inf.write("tolerance 2.0\n")
            inf.write("filetype xyz\n\n")
            if pbc is not None:
                inf.write(f"pbc {' '.join([str(round(ll,3)) for ll in pbc])}\n")
            inf.write(f"output {packmol_output}\n\n")
            for section in sections:
                inf.write(f"structure {section[0]}\n")
                inf.write(f"  number {section[1]}\n")
                text_box = " ".join([ss if isinstance(ss,str) else str(round(ss,3)) for ss in section[2]])
                if section[3] is not None:
                    if "center" in section[3]:
                        inf.write("  center\n")
                        text = " ".join([str(round(ss,3)) for ss in section[3]["center"]] + ["0.0","0.0","0.0"])
                        inf.write(f"    fixed {text}\n")
                    else:
                        inf.write(f"  {text_box}\n")
                else:
                    inf.write(f"  {text_box}\n")
                inf.write("end structure\n\n")

        Packmol._run_packmol(packmol_input,output_dir)
        #return Packmol._get_coordinates(packmol_output)

    @staticmethod 
    def packmol_bilayer(sections,fn,output_dir,pbc=None):
        """
        structure water.xyz
            number 100
            inside box 0. 0. 0. 40. 40. 40.
        end structure
        
        structure surfactant.xyz
            number 100
            inside box 0. 0. 42. 40. 40. 50.0
            atom 8 9 10
                below 0. 0. 1.0 42.0
            end atom
        end structure
        
        structure surfactant.xyz
            number 1000
            inside box 0. 0. 52. 40. 40. 60.0
            atom 8 9 10
                above 0. 0. 1.0 58.0
            end atom
        end structure
        
        structure water.xyz
            number 100
            inside box 0. 0. 62. 40. 40. 100.
        end structure
        """
        packmol_input = f"{output_dir}/{fn}.inp"
        packmol_output = f"{fn}.xyz"
        with open(packmol_input, "w") as inf:
            inf.write("seed 0\n")
            inf.write("tolerance 2.0\n")
            inf.write("filetype xyz\n\n")
            if pbc is not None:
                inf.write(f"pbc {' '.join([str(round(ll,3)) for ll in pbc])}\n")
            inf.write(f"output {packmol_output}\n\n")   
            for section in sections:
                inf.write(f"structure {section[0]}\n")
                inf.write(f"  number {section[1]}\n")
                text_box = " ".join([ss if isinstance(ss,str) else str(round(ss,3)) for ss in section[2]])
                inf.write(f"  {text_box}\n")
                if len(section) > 3:
                    text = " ".join([ss if isinstance(ss,str) else str(ss) for ss in section[3]])
                    inf.write(f"  {text}\n")
                    text = " ".join([ss if isinstance(ss,str) else str(ss) for ss in section[4]])
                    inf.write(f"    {text}\n")
                    inf.write("  end atom\n")
                inf.write("end structure\n\n") 

    @staticmethod
    def _run_packmol(input,output_dir):
        #logger.info(f"packmol run :{input}")
        if packmol_path is None:
            logger.error(
                "Please install packmol before running this script: "
                "See https://m3g.iqm.unicamp.br/packmol/home.shtml for more details"
            )
            return

        with open(input, "r") as f:
            # p = subprocess.run([packmol_path], stdin=f,cwd=self.output_dir)
            p = subprocess.run([packmol_path], stdin=f, stdout=subprocess.DEVNULL, cwd=output_dir)
            assert p.returncode == 0

    @staticmethod
    def _get_coordinates(fn,output_dir):
        with open(f"{output_dir}/{fn}") as inf:
            lines = inf.readlines()
        coors = []
        for line in lines[2:]:
            coors.append([float(rr) for rr in line.split()[1:]])
        return coors
