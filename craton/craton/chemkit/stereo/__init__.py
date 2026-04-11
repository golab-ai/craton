from .create_coordinates import create_coordinates


class Stereo:
    def __init__(self) -> None:
        pass    
    
    @staticmethod
    def create_coor(molecule,loss_atoms):
        if loss_atoms:
            set_coor_atoms = []
            for atom in molecule.Atoms:
                if atom.element not in  ["H","F","Cl","Br","I"] and hasattr(atom,"coordinates"):
                    #loss_H = [self.protein.Atoms[an] for an in atom.connectivity if self.protein.Atoms[an].element == "H" and not hasattr(self.protein.Atoms[an],"coordinates")]
                    loss_H = [molecule.Atoms[an] for an in atom.connectivity if not hasattr(molecule.Atoms[an],"coordinates")]
                    if loss_H:
                        loss_H_ids = [hatom.ID for hatom in loss_H]
                        neight_atoms = [molecule.Atoms[an] for an in atom.connectivity if an not in loss_H_ids]
                        molecule.roots.append(atom)
                        molecule.addHs.append(loss_H_ids)
                        create_coordinates(atom,neight_atoms,loss_H)
                        set_coor_atoms.extend(loss_H_ids)
            loss_atoms = list(set(loss_atoms).difference(set(set_coor_atoms)))
            Stereo.create_coor(molecule,loss_atoms)

