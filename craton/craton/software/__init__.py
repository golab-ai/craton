from .amino_acid_template import read_amberff_to_ff, read_rtp_file, create_data, restructure_amino_acid_data

class AmberUtils:
    def __init__(self) -> None:
        pass
    
    @staticmethod
    def aminoacid_json(ff):
        tt = read_rtp_file(ff)
        molecules = create_data(tt)
        r_molecules = restructure_amino_acid_data(molecules)
        for rname, atoms in r_molecules.items():
            for kk,vv in atoms["template"].items():
                formal_charge = sum([atoms[aname][vv[2]] for aname in vv[0]])
                ff_charge = sum([atoms[aname][vv[3]] for aname in vv[0]])
                if abs(formal_charge - ff_charge) > 0.00000000001:
                    print(rname,kk,"formal_charge:",formal_charge," ff_charge:",ff_charge)
        import json
        with open("amino_acid_json.json",'w') as outf:
            outf.write(json.dumps(r_molecules))
    
    @staticmethod
    def amberff_to_ff(atf,nonbf,bondf):
        read_amberff_to_ff(atf,nonbf,bondf)