import shutil
import subprocess

from ..utils import logger


packmol_path = shutil.which("packmol")
#WATER_XYZ = share_template.water_xyz
#CL_XYZ = share_template.cl_xyz
#NA_XYZ = share_template.na_xyz
#AION_XYZ = template_dir / "ion/AION.xyz"



class Packmol:
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
    
    def __init__(self,):
        pass

    @staticmethod
    def packmol(text,input,output,output_dir,style):
        with open(input,'w') as outf:
            outf.write(text)

        Packmol._run_packmol(input,output_dir,style)
        return Packmol._get_coordinates(output)

    @staticmethod
    def _run_packmol(input,output_dir,style):
        #logger.info(f"packmol run :{input}")
        if packmol_path is None:
            logger.error(
                "Please install packmol before running this script: "
                "See https://m3g.iqm.unicamp.br/packmol/home.shtml for more details"
            )
            return

        with open(input, "r") as f:
            p = subprocess.run(
                [packmol_path],
                stdin=f,
                stdout=subprocess.PIPE,
                cwd=output_dir,
                text=True,
            )
            if p.stdout:
                print("---------- Last 200 lines of packmol stdout ----------")
                for line in p.stdout.splitlines()[-200:]:
                    print(line)
            #if style != "layer":
            #    assert p.returncode == 0

    @staticmethod
    def _get_coordinates(output):
        with open(output) as inf:
            lines = inf.readlines()
        coors = []
        for line in lines[2:]:
            coors.append([float(rr) for rr in line.split()[1:]])
        return coors
