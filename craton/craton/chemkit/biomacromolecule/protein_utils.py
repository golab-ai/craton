import json
from pathlib import Path
from copy import deepcopy
import itertools
from ...chem.elements import get_bonded_type_distance
import numpy as np
from simtk import openmm as mm
from simtk import unit

from ...utils import logger
from ..stereo.create_coordinates import create_coordinates


from ... import CRATON_CONFIGURE
DEFAULT_PROTEIN_FORCEFIELD = CRATON_CONFIGURE["ForceFieldSetting"]["PROTEIN_FORCE_FIELD_FILE"]
DEFAULT_COMPEN_TYPING_FILE = CRATON_CONFIGURE["ForceFieldSetting"]["COMPEN_TYPING_FILE"]
DEFAULT_COMPEN_FORCE_FIELD_FILE = CRATON_CONFIGURE["ForceFieldSetting"]["COMPEN_FORCE_FIELD_FILE"]

amino_acid_json_f = f'{CRATON_CONFIGURE["path"]["template"]}/amino_acid.json'
non_normal_amino_acid_json_f = f'{CRATON_CONFIGURE["path"]["template"]}/non_normal_amino_acid.json'

amino_acid = json.loads(open(amino_acid_json_f).read())
non_normal_amino_acid_total = json.loads(open(non_normal_amino_acid_json_f).read())
non_AA_register = non_normal_amino_acid_total["registered"]
non_normal_amino_acid = {}

for kk,vv in non_normal_amino_acid_total.items():
    if kk not in ["registered"]:
        amino_acid[kk] = vv
        non_normal_amino_acid[kk] = vv



TERMINAL_RESIDUE = ["ACE","NME","MEC","MEN","NHE","NH2"]
CUTOFF_RESIDUE = 2.2
N_TERMINAL = "ACE"
C_TERMINAL = "NME"

def template_molecule(template,residue_name):
    from ...chem.atom import Atom
    from ...chem.molecule import Molecule
    _labels = [
            "element",
            "formal_charge",
            "atom_name",
            ]
    atoms = []
    for name in template["template"][residue_name][0]:
        atoms.append(Atom("aa"))
        atoms[-1].residue = residue_name
        for attr in _labels:
            setattr(atoms[-1],attr,template[name][attr])
    rm = Molecule("residue")
    rm.Atoms = atoms
    return rm

n_terminal_molecule = template_molecule(amino_acid[N_TERMINAL],N_TERMINAL)
c_terminal_molecule = template_molecule(amino_acid[C_TERMINAL],C_TERMINAL)

def assign_force_field(protein):
    from ...force_field import MolForceField as MFF #grasp_force_field
    for atom in protein.Atoms:
        atom.nonb_atom_type = atom.atom_type_name
        atom.binc_atom_type = atom.atom_type_name
        atom.atc_atom_type = atom.atom_type_names
    for attr in ["Bonds","Angles","Dihedrals","Impropers","Pair14"]:
        for term in getattr(protein,attr,[]):
            for ii in range(1,5):
                if hasattr(term,f"a{ii}"):
                    setattr(term,f"a{ii}_atom_type",protein.Atoms[getattr(term,f"a{ii}")].atom_type_name)
                    setattr(term,f"a{ii}_atom_type_used",protein.Atoms[getattr(term,f"a{ii}")].atom_type_name)
            
    proteins = MFF.grasp_force_field(protein,None,DEFAULT_PROTEIN_FORCEFIELD,empi_ff_flag=False,use_scalevdw=False,parallel=False)
    return proteins[0]

def get_residue_protein(protein,residue_name):
    residues = {}
    atom0 = protein.Atoms[0]
    pre_name = f"{atom0.residue}_{atom0.residue_ID}_{atom0.chain_name}"
    tmp = []
    for atom in protein.Atoms:
            #if atom.atom_name[0].isdigit():
            #    atom.atom_name = atom.atom_name[1:] + atom.atom_name[0]
        res_name = f"{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
        if res_name == pre_name:
            tmp.append(atom)
        else:
            residues[pre_name] = tmp
            pre_name = res_name
            tmp = [atom]
    residues[pre_name] = tmp
    return residues[residue_name]

def protein_atom_mapping(protein1,protein2):
    atom_mapping = {}
    if hasattr(protein2,"relative_shift"):
        ii = protein2.relative_shift[0]
        jj = protein2.relative_shift[1]
    else:
        ii = 0
        jj = 0
    
    nn = min(len(protein1.Groups) - ii,len(protein2.Groups) - jj)
    
    for kk in range(nn):
        group1 = protein1.Groups[ii+kk]
        group2 = protein2.Groups[jj+kk]

        an_dict1 = {protein1.Atoms[an].atom_name:an for an in group1.atoms}
        an_dict2 = {protein2.Atoms[an].atom_name:an for an in group2.atoms}
        if group1.group_str == group2.group_str:
            for aname,an in an_dict1.items():
                if aname in an_dict2:
                    atom_mapping[an] = an_dict2[aname]
        else:
            for aname in ["N", "H", "CA", "HA","CB","C", "O","OC1", "OC2","H1", "H2", "H3","HB","HB1","HB2"]:
                if aname in an_dict1 and aname in an_dict2:
                        atom_mapping[an_dict1[aname]] = an_dict2[aname]
    return atom_mapping

def o_protein_atom_mapping(protein1,protein2):
    atom_mapping = {}
    def get_dicts(protein):
        residues = {}
        pre_name = f"{protein.Atoms[0].residue}_{protein.Atoms[0].residue_ID}_{protein.Atoms[0].chain_name}"
        tmp = {}
        for atom in protein.Atoms:
            atom.bond_type_aromatic = atom.bond_type
            res_name = f"{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if res_name == pre_name:
                tmp[atom.atom_name] = atom.ID
            else:
                residues[pre_name]=tmp
                pre_name = res_name
                tmp = {atom.atom_name:atom.ID}
        residues[pre_name] = tmp
        return residues
    
    residues_1 = get_dicts(protein1)
    residues_2 = get_dicts(protein2)
    residues_2_s = {}
    for kk,vv in residues_2.items():
        ss = kk.split("_")
        nkk = f"{ss[1]}_{ss[2]}"
        residues_2_s[nkk] = vv
    
    
    for kk,vv1 in residues_1.items():
        if kk in  residues_2:
            vv2 = residues_2[kk]
            for an,zz in vv1.items():
                if an in vv2:
                    atom_mapping[zz] = vv2[an]
        else:
            ss = kk.split("_")
            nkk = f"{ss[1]}_{ss[2]}" 
            if nkk in residues_2_s:
                vv2 = residues_2_s[nkk]
                for an in ["N", "H", "CA", "HA","CB","C", "O","OC1", "OC2","H1", "H2", "H3","HB","HB1","HB2"]:
                    if an in vv1 and an in vv2:
                        atom_mapping[vv1[an]] = vv2[an]
    return atom_mapping

def old_protein_atom_mapping(protein1,protein2):
    atom_mapping = {}
    def get_dicts(protein):
        residues = {}
        pre_name = f"{protein.Atoms[0].residue}_{protein.Atoms[0].residue_ID}_{protein.Atoms[0].chain_name}"
        tmp = {}
        for atom in protein.Atoms:
            atom.bond_type_aromatic = atom.bond_type
            res_name = f"{atom.residue}_{atom.residue_ID}_{atom.chain_name}"
            if res_name == pre_name:
                tmp[atom.atom_name] = atom.ID
            else:
                residues[pre_name]=tmp
                pre_name = res_name
                tmp = {atom.atom_name:atom.ID}
        residues[pre_name] = tmp
        return residues
    
    residues_1 = get_dicts(protein1)
    residues_2 = get_dicts(protein2)
    residues_2_s = {}
    for kk,vv in residues_2.items():
        ss = kk.split("_")
        nkk = f"{ss[1]}_{ss[2]}"
        residues_2_s[nkk] = vv
    
    
    for kk,vv1 in residues_1.items():
        if kk in  residues_2:
            vv2 = residues_2[kk]
            for an,zz in vv1.items():
                if an in vv2:
                    atom_mapping[zz] = vv2[an]
        else:
            ss = kk.split("_")
            nkk = f"{ss[1]}_{ss[2]}" 
            if nkk in residues_2_s:
                vv2 = residues_2_s[nkk]
                for an in ["N", "H", "CA", "HA","CB","C", "O","OC1", "OC2","H1", "H2", "H3","HB","HB1","HB2"]:
                    if an in vv1 and an in vv2:
                        atom_mapping[vv1[an]] = vv2[an]
    return atom_mapping

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
        create_coor(molecule,loss_atoms)
    
def optimize_position(molecule,roots,addHs):
    if roots:
        addHs_all = [a for l in addHs for a in l]
        msg = f"Adding {len(addHs_all)} hydrogen atoms based on residue templates"
        logger.warning(msg)

        system = mm.System()
        for atom in molecule.Atoms:
            system.addParticle(0)
        for lr in addHs:
            for an in lr:
                system.setParticleMass(an, 1)

        bforce = mm.HarmonicBondForce()
        for root, added in zip(roots, addHs):
            for kk,neigh in enumerate(root.connectivity):
                if neigh in added:
                    bond_dist = get_bonded_type_distance(root.elem,molecule.Atoms[neigh].elem,root.bond_type[kk])
                    bforce.addBond(root.ID, neigh, bond_dist/10.0, 100.0)
            for neigh1, neigh2 in itertools.combinations(root.connectivity, 2):
                if (molecule.Atoms[neigh1].elem == "H" and neigh2 in added) or (neigh1 in added and molecule.Atoms[neigh2].elem == "H"):
                    bforce.addBond(neigh1, neigh2, 0.16, 100.0)
        system.addForce(bforce)

        aforce = mm.HarmonicAngleForce()
        for angle in molecule.Angles:
            a1, a2, a3 = angle.a1, angle.a2, angle.a3
            theta = 120.0 if len(molecule.Atoms[a2].connectivity) == 3 else 109.5
            aforce.addAngle(a1, a2, a3, theta / 180 * np.pi, 100.0)
        system.addForce(aforce)

        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName("CPU")
        context = mm.Context(system, integrator, platform)
        #positions = [a.coordinates for a in self.protein.Atoms]
        positions = [[coor/10.0 for coor in a.coordinates]  for a in molecule.Atoms ]
            
        context.setPositions(positions)
        # logger.debug(f"Energy before optimization: {context.getState(getEnergy=True).getPotentialEnergy()}")
        mm.LocalEnergyMinimizer.minimize(context, tolerance=1e-4)
        # logger.debug(f"Energy after optimization: {context.getState(getEnergy=True).getPotentialEnergy()}")

        positions = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        for atom, pos in zip(molecule.Atoms, positions):
            atom.coordinates = [pp*10.0 for pp in pos] 
            
def get_max_sequence(origin_sequ,sequ):
    tmp = []
    nn = len(sequ)
    mm = len(origin_sequ)
    for ii in range(nn):
        for jj in range(mm):
            tmp.append([ii,jj,0,0])
            for iii in range(ii,nn):
                jjj = jj + (iii-ii)
                if jjj < mm:
                    tmp[-1][2] += 1
                    if origin_sequ[jjj] == sequ[iii]:
                        tmp[-1][3] += 1
    tmp = sorted(tmp,key=lambda x:x[3],reverse=True)
    nnn = tmp[0][3]
    _tmp_ = [rr for rr in tmp if rr[3] == nnn]
    if len(_tmp_) > 1:
        _tmp_ = sorted(_tmp_,key=lambda x:x[2],reverse=True)
        mmm = _tmp_[0][2]
        _tmp_tmp = [rr for rr in _tmp_ if rr[2] == mmm]
        if len(_tmp_tmp) > 1:
            _tmp_tmp = sorted(_tmp_tmp,key=lambda x:x[1])
            kkk = _tmp_tmp[0][1]
            _tmp_tmp_ = [rr for rr in _tmp_tmp if rr[1] == kkk]
            if len(_tmp_tmp_) > 1:
                return sorted(_tmp_tmp_,key=lambda x:x[0])[0]
            else:
                return _tmp_tmp_[0]
        else:
            return _tmp_tmp[0]
    else:
        return _tmp_[0]         
    
def create_mutation(protein,sequences):
    _label_ = {"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY","H":"HIS","I":"ILE",
              "K":"LYS","L":"LEU","M":"MET","N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER",
              "T":"THR","V":"VAL","W":"TRP","Y":"TYR","NME":"NME","ACE":"ACE","MEC":"MEA","MEN":"MEN","NHE":"NHE","NH2":"NH2"}
    _label = {vv:kk for kk,vv in _label_.items()}
    changes = []
    origin_sequence_str = [group.group_str for group in protein.Groups]
    origin_sequence = [_label[group.group_name] for group in protein.Groups]
    snn = len(origin_sequence)
    for sequence in sequences:
        tnn = len(sequence)
        rr = get_max_sequence(origin_sequence,sequence)
        mutation = []
        ii = rr[0] # target sequence start matched residue
        jj = rr[1] # reference sequence start matched residue
        kk = tnn - ii# the length of target sequence matched residue
        mm = snn - jj# the length of reference sequence matched residue
        
        if ii != jj:
            if ii > jj:
                mutation.append([[_label[sequence[iii]] for iii in range(ii-1,-1,-1)],"add",None,None,"n_terminal"])
            else:
                mutation.append([[origin_sequence_str[iii] for iii in range(0,jj)],"delete",None,None,"n_terminal"])
        tmp = []
        for iii in range(ii,min(kk,mm)):
            if sequence[iii] != origin_sequence[jj+iii]:
                tmp.append([origin_sequence_str[jj+iii],_label_[sequence[ii+iii]]])
        if tmp:
            mutation.append([tmp,"mutation"])
        if kk != mm:
            if kk > mm:
                mutation.append([[_label[sequence[iii]] for iii in range(ii+mm,tnn)],"add",None,None,"c_terminal"])
            else:
                mutation.append([[origin_sequence_str[iii] for iii in range(jj+kk,snn)],"delete",None,None,"c_terminal"])
        changes.append(mutation)
    return changes

def molecule_alignment(self):
    if self.config["AlignmentSetting"]["mutation"] is not None:
        changes = self.config["AlignmentSetting"]["mutation"]
    else:
        changes = []
    
        
    if self.config["AlignmentSetting"]["sequences"] is not None:
        changes.extend(create_mutation())     
    
    for change in changes:
        pre_pep = self.molecules["ligands"][0]
        pre_pep_name = ""
        for chag in change:
            pre_pep = MX.protein_process(pre_pep,chag)
            pre_pep_name += pre_pep.mole_name
        self.molecules["ligands"].append(pre_pep)                
            
        #self.molecules["ligands"].append(MX.protein_process(self.molecules["ligands"][0],change))
    self._before_integral_ligands = {molecule.mole_name: deepcopy(molecule) for molecule in self.molecules["ligands"]}