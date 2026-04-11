from .amino_acid import MoleculeAssembly


def create_aminoacid(n,left_cap="ACE",right_cap="NME",terminal_flag=True,templates=None):
    MoleculeAssembly.peptide_gen(n,left_cap=left_cap,right_cap=right_cap,terminal_flag=terminal_flag,templates=templates)
    #AA = AminoAcid(rc,left_cap=left_cap,right_cap=right_cap,template=template,output_dir=output_dir)
    #return AA.combine_residue()
    
def create_rnadna(n,templates=None):
    MoleculeAssembly.dnarna_gen(n,templates=templates)
    #AA = AminoAcid(rc,left_cap=left_cap,right_cap=right_cap,template=template,output_dir=output_dir)
    #return AA.combine_residue()
    
