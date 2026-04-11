import os,sys
from .pocket_analyzer import PocketAnalyzer as PA
from .docking import VinaRun as VR
from ..utils.commons import parallel_run

def get_pocket(protein,
               csv_file=None,
               cavity_file=None,
               no_dispaly=True,
               step=0.6,
               probe_in=1.4,
               probe_out=1.4,
               removal_distance=2.4,
               volume_cutoff=5,
               include_depth=True,
               include_hydropathy=False,
               verbose=False
               ):
    pock = PA(protein,
              csv_file=csv_file,
              cavity_file=cavity_file,
              no_dispaly=no_dispaly,
              step=step,
              probe_in=probe_in,
              probe_out=probe_out,
              removal_distance=removal_distance,
              volume_cutoff=volume_cutoff,
              include_depth=include_depth,
              include_hydropathy=include_hydropathy,
              verbose=verbose
              )
    results = pock.run()
    return results

class Dock:

    def __init__(self,protein,ligands,center,box_size,output_directory=None,parallel=True):
        self.protein = protein
        self.ligands = ligands
        self.center = center
        self.box_size = box_size
        self.output_directory = output_directory
        self.parallel = parallel

    def run_docking(self,ligand,idx=None):
        if self.output_directory is not None:
            if not ligand[0].startswith(self.output_directory):
                output_pdbqt = f"{self.output_directory}/{ligand[0]}_docked.pdbqt"
            else:
                output_pdbqt = f"{ligand[0]}_docked.pdbqt"
        else:
            output_pdbqt = None
        vina = VR(self.protein,ligand[1],center=self.center,box_size=self.box_size,output_pdbqt=output_pdbqt)
        score = vina.run_vina()
        if idx is None:
            return score
        else:
            return score,idx

    def docking(self):
        if self.parallel:
            results = parallel_run(self.run_docking,self.ligands)
        else:
            results = []
            for ligand in self.ligands:
                results.append(self.run_docking(ligand))
        return results
        


       
        

