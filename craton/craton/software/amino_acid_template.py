from copy import deepcopy
import sys

terminals = {
    "ACE": {
        "template":{
            "ACE":[["C","O","CH3", "HH31", "HH32", "HH33"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","LT"]
        },
        "C": {
            "element": "C",
            "atom_name": "C",
            "atom_type_name": "C",
            "connectivity": [
                "O",
                "CH3",
                "L*"
            ],
            "bond_type": [
                "2",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": 0.5972
        },
        "O": {
            "element": "O",
            "atom_name": "O",
            "atom_type_name": "O",
            "connectivity": [
                "C"
            ],
            "bond_type": [
                "2"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": -0.5679
        },
        "HH31": {
            "element": "H",
            "atom_name": "HH31",
            "atom_type_name": "HC",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.1123
        },
        "CH3": {
            "element": "C",
            "atom_name": "CH3",
            "atom_type_name": "CT",
            "connectivity": [
                "C",
                "HH31",
                "HH32",
                "HH33"
            ],
            "bond_type": [
                "1",
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": -0.3662
        },
        "HH32": {
            "element": "H",
            "atom_name": "HH32",
            "atom_type_name": "HC",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.1123
        },
        "HH33": {
            "element": "H",
            "atom_name": "HH33",
            "atom_type_name": "HC",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.1123
        }
    },
    "NME": {
        "template":{
            "NME":[["N","H","CH3", "HH31", "HH32", "HH33"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","RT"]
        },
        "N": {
            "element": "N",
            "atom_name": "N",
            "atom_type_name": "N",
            "connectivity": [
                "H",
                "CH3",
                "R*"
            ],
            "bond_type": [
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": -0.4157
        },
        "H": {
            "element": "H",
            "atom_name": "H",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.2719
        },
        "HH31": {
            "element": "H",
            "atom_name": "HH31",
            "atom_type_name": "H1",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        },
        "CH3": {
            "element": "C",
            "atom_name": "CH3",
            "atom_type_name": "CT",
            "connectivity": [
                "N",
                "HH31",
                "HH32",
                "HH33"
            ],
            "bond_type": [
                "1",
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": -0.1490
        },
        "HH32": {
            "element": "H",
            "atom_name": "HH32",
            "atom_type_name": "H1",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        },
        "HH33": {
            "element": "H",
            "atom_name": "HH33",
            "atom_type_name": "H1",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        }
    },
    "NMA": {
        "template":{
            "NMA":[["N","H","CA", "HA1", "HA2", "HA3"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","RT"]
        },
        "N": {
            "element": "N",
            "atom_name": "N",
            "atom_type_name": "N",
            "connectivity": [
                "H",
                "CA",
                "R*"
            ],
            "bond_type": [
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": -0.4157
        },
        "H": {
            "element": "H",
            "atom_name": "H",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.2719
        },
        "HA1": {
            "element": "H",
            "atom_name": "HA1",
            "atom_type_name": "H1",
            "connectivity": [
                "CA"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        },
        "CA": {
            "element": "C",
            "atom_name": "CA",
            "atom_type_name": "CT",
            "connectivity": [
                "N",
                "HA1",
                "HA2",
                "HA3"
            ],
            "bond_type": [
                "1",
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": -0.1490
        },
        "HA2": {
            "element": "H",
            "atom_name": "HA2",
            "atom_type_name": "H1",
            "connectivity": [
                "CA"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        },
        "HA3": {
            "element": "H",
            "atom_name": "HA3",
            "atom_type_name": "H1",
            "connectivity": [
                "CA"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        }
    },
    "MEC": {
        "template":{
            "MEC":[["CH3", "HH31", "HH32", "HH33"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","RT"]
        },
        "HH31": {
            "element": "H",
            "atom_name": "HH31",
            "atom_type_name": "HC",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.1123
        },
        "CH3": {
            "element": "C",
            "atom_name": "CH3",
            "atom_type_name": "CT",
            "connectivity": [
                "HH31",
                "HH32",
                "HH33",
                "R*"
            ],
            "bond_type": [
                "1",
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": -0.3369
        },
        "HH32": {
            "element": "H",
            "atom_name": "HH32",
            "atom_type_name": "HC",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.1123
        },
        "HH33": {
            "element": "H",
            "atom_name": "HH33",
            "atom_type_name": "HC",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.1123
        }
    },
    "MEN": {
        "template":{
            "MEN":[["CH3", "HH31", "HH32", "HH33"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","LT"]
        },
        "CH3": {
            "element": "C",
            "atom_name": "CH3",
            "atom_type_name": "CT",
            "connectivity": [
                "HH31",
                "HH32",
                "HH33",
                "L*"
            ],
            "bond_type": [
                "1",
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": -0.2928
        },
        "HH31": {
            "element": "H",
            "atom_name": "HH31",
            "atom_type_name": "H1",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        },
        "HH32": {
            "element": "H",
            "atom_name": "HH32",
            "atom_type_name": "H1",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        },
        "HH33": {
            "element": "H",
            "atom_name": "HH33",
            "atom_type_name": "H1",
            "connectivity": [
                "CH3"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.0976
        }
    },
    "NHE": {
        "template":{
            "NHE":[["N", "H1", "H2"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","RT"]
        },
        "N": {
            "element": "N",
            "atom_name": "N",
            "atom_type_name": "N",
            "connectivity": [
                "H1",
                "H2",
                "R*"
            ],
            "bond_type": [
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": -0.463
        },
        "H1": {
            "element": "H",
            "atom_name": "H1",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.2315
        },
        "H2": {
            "element": "H",
            "atom_name": "H2",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.2315
        }
    },
    "NH2": {
        "template":{
            "NH2":[["N", "H1", "H2"],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","RT"]
        },
        "N": {
            "element": "N",
            "atom_name": "N",
            "atom_type_name": "N",
            "connectivity": [
                "H1",
                "H2",
                "R*"
            ],
            "bond_type": [
                "1",
                "1",
                "1"
            ],
            "formal_charge": 0,
            "plate": "yes",
            "ff_charge": -0.463
        },
        "H1": {
            "element": "H",
            "atom_name": "H1",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.2315
        },
        "H2": {
            "element": "H",
            "atom_name": "H2",
            "atom_type_name": "H",
            "connectivity": [
                "N"
            ],
            "bond_type": [
                "1"
            ],
            "formal_charge": 0,
            "plate": "no",
            "ff_charge": 0.2315
        }
    },
}

AA_property = {
    "ALA":["alanine","A","丙氨酸","hydrophobic",6.11,2.35,9.87,None,7.8,"-","GCU, GCC, GCA, GCG"],
    "CYS":["cysteine","C","","hydrophilic",5.05,1.92,10.70,8.37,1.9,"-","UGU, UGC"],
    "ASP":["aspartate","D","天冬氨酸","acid",2.85,1.99,9.90,3.90,5.3,"-","GAU, GAC"],
    "GLU":["glutamate","E","谷氨酸","acid",3.15,2.10,9.47,4.07,6.3,"-","GAA, GAG"],
    "PHE":["phenylalanine","F","苯丙氨酸","hydrophobic",5.49,2.20,9.31,None,3.9,"X","UUU, UUC"],
    "GLY":["glycine","G","甘氨酸","hydrophobic",6.06,2.35,9.78,None,7.2,"-","GGU, GGC, GGA, GGG"],
    "HIE":["histidine","H","组氨酸","base",7.60,1.80,9.33,6.04,2.3,"X","CAU, CAC"],
    "ILE":["isoleucine","I","异亮氨酸","hydrophobic",6.05,2.32,9.76,None,5.3,"X","AUU, AUC, AUA"],
    "LYS":["lysine","K","赖氨酸","base",9.60,2.16,9.06,10.54,5.9,"X","AAA, AAG"],
    "LEU":["luecine","L","亮氨酸","hydrophobic",6.01,2.33,9.74,None,9.1,"X","UUA, UUG, CUU, CUC, CUA, CUG"],
    "MET":["methionine","M","甲硫氨酸","hydrophobic",5.74,2.13,9.28,None,2.3,"X","AUG"],
    "ASN":["asparagine","N","天冬酰胺","hydrophilic",5.41,2.14,8.72,None,4.3,"-","AAU, AAC"],
    "PRO":["proline","P","脯氨酸","hydrophobic",6.30,1.95,10.64,None,5.2,"-","CCU, CCC, CCA, CCG"],
    "GLN":["glutamine","Q","谷氨酰胺","hydrophilic",5.65,2.17,9.13,None,4.2,"-","CAA, CAG"],
    "ARG":["arginine","R","精氨酸","base",10.76,1.82,8.99,12.48,5.1,"-","CGU, CGC, CGA, CGG, AGA, AGG"],
    "SER":["serine","S","丝氨酸","hydrophilic",5.68,2.19,9.21,None,6.8,"_","UCU, UCC, UCA, UCG, AGU, AGC"],
    "THR":["threonine","T","苏氨酸","hydrophilic",5.60,2.09,9.10,None,5.9,"X","ACU, ACC, ACA, ACG"],
    "VAL":["valine","V","缬氨酸","hydrophobic",6.00,2.39,9.74,None,6.6,"X","GUU, GUC, GUA, GUG"],
    "TRP":["tryptophan","W","色氨酸","hydrophobic",5.89,2.46,9.41,None,1.4,"X","UGG"],
    "TYR":["tyrosine","Y","酷氨酸","hydrophilic",5.64,2.20,9.21,10.46,3.2,"_","UAU, UAC"],
}


def read_rtp_file(ff):
   molecules = {}
   with open(ff) as inf:
      lines = inf.readlines()
      
   _label = []
   tmp = []
   ii = -1
   while 1:
      ii += 1
      if ii >= len(lines):
         break
      line = lines[ii]
      if line[:2] == "[ ":
         tmp.append(ii)
         while 1:
            ii += 1
            line = lines[ii]
            if line.strip() == "[ atoms ]":
               tmp.append(ii)
            elif line.strip() == "[ bonds ]":
               tmp.append(ii)
            elif line.strip() == "[ impropers ]":
               tmp.append(ii)
            elif line.strip() == "":
               tmp.append(ii)
               _label.append(tmp)
               tmp = []
               break
   
   for arr in _label:
      rname = lines[arr[0]].split(";")[0].strip().strip("[ ").strip(" ]")
      molecules[rname] = {}
      if len(arr) > 2:
         molecules[rname]["atoms"] = {lines[ii].strip().split()[0]:lines[ii].strip().split() for ii in range(arr[1]+1,arr[2])}
   
      if len(arr) > 3:
         molecules[rname]["bonds"] = [lines[ii].strip().split() for ii in range(arr[2]+1,arr[3])]
      if len(arr) > 4:
         molecules[rname]["impropers"] = [lines[ii].strip().split()[2] for ii in range(arr[3]+1,arr[4])]
   
   return molecules
               
def create_data(data):
   molecules = {}
   del data["bondedtypes"]
   for rname,dd in data.items():
      molecules[rname] = {}
      if "impropers" in dd:
         impropers = dd["impropers"]
      else:
         impropers = []
      connectivity = {kk:[] for kk in dd["atoms"].keys()}
      if "bonds" in dd:
         for rr in dd["bonds"]:
            #if rr[0] != "-C" and rr[1] != "+N":
            if rr[0] in dd["atoms"] and rr[1] in dd["atoms"]:
               connectivity[rr[0]].append(rr[1])
               connectivity[rr[1]].append(rr[0])
      
      for aname,atom in dd["atoms"].items():
         molecules[rname][aname] = {"element":aname[0],"atom_name":atom[0],"atom_type_name":atom[1],
                                    "connectivity":connectivity[aname],"bond_type":["1" for __ in connectivity[aname]],
                                    "formal_charge":0,"plate": "yes" if aname in impropers else "no", 
                                    "ff_charge":float(atom[2])}
   return molecules
         
def restructure_amino_acid_data(data):
   molecules = deepcopy(terminals)
   
   ###molecules = {"ACE":data["ACE"],"NME":data["NME"],"NHE":data["NHE"],"NH2":data["NH2"]}
   for rname in ["ALA","GLY","SER","THR","LEU","ILE","VAL","ASN","GLN","ARG",
                 "HIE","TRP","PHE","TYR","GLU","ASP","LYS","PRO","CYS","MET"]:
      molecules[rname] = deepcopy(data[rname])
      
      c_atoms = data[f"C{rname}"]
      n_atoms = data[f"N{rname}"]
      atoms = data[rname]

      n_ans = [aname for aname in n_atoms if aname not in atoms]
      c_ans = [aname for aname in c_atoms if aname not in atoms]
      try:
         _tmp_charge = round(0-atoms["H"]["ff_charge"]-atoms["N"]["ff_charge"],4)
      except:
         _tmp_charge = round(0-atoms["N"]["ff_charge"],4)
      for aname, atom in molecules[rname].items():
         atom["c_formal_charge"] = 0
         atom["n_formal_charge"] = 0
         if aname == 'N':
            atom["n_formal_charge"] = 1
         if aname == "C":
            atom["bond_type"][atom["connectivity"].index("O")] = "2"
         if aname == "H":
            atom["n_ff_charge"] = n_atoms["H1"]["ff_charge"]
            atom["c_ff_charge"] = c_atoms["H"]["ff_charge"]
            atom["h_ff_charge"] = 0.2315
         elif aname == "O":
            atom["n_ff_charge"] = n_atoms["O"]["ff_charge"]
            atom["c_ff_charge"] = c_atoms["OC1"]["ff_charge"]
            atom["h_ff_charge"] = atom["ff_charge"]
            atom["bond_type"][0] = "2"
         elif aname == "N":
            atom["n_ff_charge"] = n_atoms[aname]["ff_charge"]
            atom["c_ff_charge"] = c_atoms[aname]["ff_charge"]
            atom["h_ff_charge"] = -0.4630
         elif aname == "CA":
            atom["n_ff_charge"] = n_atoms[aname]["ff_charge"]
            atom["c_ff_charge"] = c_atoms[aname]["ff_charge"]
            atom["h_ff_charge"] = atom["ff_charge"] - _tmp_charge
         else:
            atom["n_ff_charge"] = n_atoms[aname]["ff_charge"]
            atom["c_ff_charge"] = c_atoms[aname]["ff_charge"]
            atom["h_ff_charge"] = atom["ff_charge"]
      
      
      for aname in n_ans:
         molecules[rname][aname] = deepcopy(n_atoms[aname])
         molecules[rname][aname]["n_formal_charge"] = molecules[rname][aname]["formal_charge"]
         molecules[rname][aname]["n_ff_charge"] = molecules[rname][aname]["ff_charge"]
         molecules[rname][aname]["h_ff_charge"] = 0.2315
         molecules[rname][aname]["ff_charge"] = None
         molecules[rname][aname]["c_ff_charge"] = None
         molecules[rname][aname]["formal_charge"] = 0
         molecules[rname][aname]["c_formal_charge"] = None
         for aa in n_atoms[aname]["connectivity"]:
            molecules[rname][aa]["connectivity"].append(aname)
            molecules[rname][aa]["bond_type"].append("1")
      for aname in c_ans:
         molecules[rname][aname] = deepcopy(c_atoms[aname])
         if aname == "OC1":
            molecules[rname][aname]["bond_type"][0] = "2"
         if aname == "OC2":
            molecules[rname][aname]["c_formal_charge"] = -1
         else:
            molecules[rname][aname]["c_formal_charge"] = molecules[rname][aname]["formal_charge"]
         molecules[rname][aname]["c_ff_charge"] = molecules[rname][aname]["ff_charge"]
         molecules[rname][aname]["ff_charge"] = None
         molecules[rname][aname]["n_ff_charge"] = None
         molecules[rname][aname]["formal_charge"] = None
         molecules[rname][aname]["n_formal_charge"] = None
         for aa in c_atoms[aname]["connectivity"]:
            molecules[rname][aa]["connectivity"].append(aname)
            if aname == "OC1":
               molecules[rname][aa]["bond_type"].append("2")
            else:
               molecules[rname][aa]["bond_type"].append("1")
      
   
   molecules["ASP"]["OD1"]["bond_type"][0] = "2"
   molecules["ASP"]["CG"]["bond_type"][molecules["ASP"]["CG"]["connectivity"].index("OD1")] = "2"
   molecules["ASP"]["OD2"]["formal_charge"] = -1
   molecules["ASP"]["OD2"]["c_formal_charge"] = -1
   molecules["ASP"]["OD2"]["n_formal_charge"] = -1
   
   molecules["GLU"]["OE1"]["bond_type"][0] = "2"
   molecules["GLU"]["CD"]["bond_type"][molecules["GLU"]["CD"]["connectivity"].index("OE1")] = "2"
   molecules["GLU"]["OE2"]["formal_charge"] = -1
   molecules["GLU"]["OE2"]["c_formal_charge"] = -1
   molecules["GLU"]["OE2"]["n_formal_charge"] = -1
   
   molecules["PHE"]["CG"]["bond_type"][molecules["PHE"]["CG"]["connectivity"].index("CD1")] = "2"
   molecules["PHE"]["CD1"]["bond_type"][molecules["PHE"]["CD1"]["connectivity"].index("CG")] = "2"
   molecules["PHE"]["CE1"]["bond_type"][molecules["PHE"]["CE1"]["connectivity"].index("CZ")] = "2"
   molecules["PHE"]["CZ"]["bond_type"][molecules["PHE"]["CZ"]["connectivity"].index("CE1")] = "2"
   molecules["PHE"]["CD2"]["bond_type"][molecules["PHE"]["CD2"]["connectivity"].index("CE2")] = "2"
   molecules["PHE"]["CE2"]["bond_type"][molecules["PHE"]["CE2"]["connectivity"].index("CD2")] = "2"
   
   molecules["HIE"]["CG"]["bond_type"][molecules["HIE"]["CG"]["connectivity"].index("CD2")] = "2"
   molecules["HIE"]["CD2"]["bond_type"][molecules["HIE"]["CD2"]["connectivity"].index("CG")] = "2"
   molecules["HIE"]["ND1"]["bond_type"][molecules["HIE"]["ND1"]["connectivity"].index("CE1")] = "2"
   molecules["HIE"]["CE1"]["bond_type"][molecules["HIE"]["CE1"]["connectivity"].index("ND1")] = "2"
   
   molecules["LYS"]["NZ"]["formal_charge"] = 1
   molecules["LYS"]["NZ"]["c_formal_charge"] = 1
   molecules["LYS"]["NZ"]["n_formal_charge"] = 1
   
   molecules["ASN"]["OD1"]["bond_type"][0] = "2"
   molecules["ASN"]["CG"]["bond_type"][molecules["ASN"]["CG"]["connectivity"].index("OD1")] = "2"
   
   molecules["GLN"]["OE1"]["bond_type"][0] = "2"
   molecules["GLN"]["CD"]["bond_type"][molecules["GLN"]["CD"]["connectivity"].index("OE1")] = "2"
   
   molecules["ARG"]["NH1"]["formal_charge"] = 1
   molecules["ARG"]["NH1"]["c_formal_charge"] = 1
   molecules["ARG"]["NH1"]["n_formal_charge"] = 1
   molecules["ARG"]["CZ"]["bond_type"][molecules["ARG"]["CZ"]["connectivity"].index("NH1")] = "2"
   molecules["ARG"]["NH1"]["bond_type"][molecules["ARG"]["NH1"]["connectivity"].index("CZ")] = "2"
   
   molecules["TRP"]["CG"]["bond_type"][molecules["TRP"]["CG"]["connectivity"].index("CD1")] = "2"
   molecules["TRP"]["CD1"]["bond_type"][molecules["TRP"]["CD1"]["connectivity"].index("CG")] = "2"
   molecules["TRP"]["CD2"]["bond_type"][molecules["TRP"]["CD2"]["connectivity"].index("CE2")] = "2"
   molecules["TRP"]["CE2"]["bond_type"][molecules["TRP"]["CE2"]["connectivity"].index("CD2")] = "2"
   molecules["TRP"]["CZ2"]["bond_type"][molecules["TRP"]["CZ2"]["connectivity"].index("CH2")] = "2"
   molecules["TRP"]["CH2"]["bond_type"][molecules["TRP"]["CH2"]["connectivity"].index("CZ2")] = "2"
   molecules["TRP"]["CE3"]["bond_type"][molecules["TRP"]["CE3"]["connectivity"].index("CZ3")] = "2"
   molecules["TRP"]["CZ3"]["bond_type"][molecules["TRP"]["CZ3"]["connectivity"].index("CE3")] = "2"
   
   molecules["TYR"]["CG"]["bond_type"][molecules["TYR"]["CG"]["connectivity"].index("CD1")] = "2"
   molecules["TYR"]["CD1"]["bond_type"][molecules["TYR"]["CD1"]["connectivity"].index("CG")] = "2"
   molecules["TYR"]["CE1"]["bond_type"][molecules["TYR"]["CE1"]["connectivity"].index("CZ")] = "2"
   molecules["TYR"]["CZ"]["bond_type"][molecules["TYR"]["CZ"]["connectivity"].index("CE1")] = "2"
   molecules["TYR"]["CD2"]["bond_type"][molecules["TYR"]["CD2"]["connectivity"].index("CE2")] = "2"
   molecules["TYR"]["CE2"]["bond_type"][molecules["TYR"]["CE2"]["connectivity"].index("CD2")] = "2"
   
   
   special_lable = {"HID":"HIE","HIP":"HIE","CYM":"CYS","CYX":"CYS","ASH":"ASP","GLH":"GLU","LYN":"LYS"}
   for aa,bb in special_lable.items():
      t_atoms = data[aa]
      r_atoms = molecules[bb]
      t_ans = [aname for aname in t_atoms if aname not in r_atoms]
      r_ans = [aname for aname in r_atoms if aname not in t_atoms]
      t_r = [aname for aname in t_atoms if aname in r_atoms]
      for aname in t_r:
         r_atoms[aname][f"{aa}_formal_charge"] = 0
         r_atoms[aname][f"{aa}_ff_charge"] = t_atoms[aname]["ff_charge"]
         r_atoms[aname][f"{aa}_h_ff_charge"] = t_atoms[aname]["ff_charge"]
      for aname in r_ans:
         r_atoms[aname][f"{aa}_formal_charge"] = None
         r_atoms[aname][f"{aa}_ff_charge"] = None
         r_atoms[aname][f"{aa}_h_ff_charge"] = None
      for aname in t_ans:
         r_atoms[aname] = t_atoms[aname]
         r_atoms[aname][f"{aa}_formal_charge"] = r_atoms[aname]["formal_charge"]
         r_atoms[aname][f"{aa}_ff_charge"] = r_atoms[aname]["ff_charge"]
         r_atoms[aname][f"{aa}_h_ff_charge"] = r_atoms[aname]["ff_charge"]
         r_atoms[aname]["formal_charge"] = None
         r_atoms[aname]["ff_charge"] = None
         r_atoms[aname]["n_formal_charge"] = None
         r_atoms[aname]["n_ff_charge"] = None
         r_atoms[aname]["c_formal_charge"] = None
         r_atoms[aname]["c_ff_charge"] = None
         for an in t_atoms[aname]["connectivity"]:
            r_atoms[an]["connectivity"].append(aname)
            r_atoms[an]["bond_type"].append("1")
   for rname in ["HID","HIP","CYX"]:
      n_t_atoms = data[f"N{rname}"]
      c_t_atoms = data[f"C{rname}"]
      r_atoms = molecules[special_lable[rname]]
      for aname,atom in r_atoms.items():
         if aname in n_t_atoms:
            atom[f"{rname}_n_formal_charge"] = n_t_atoms[aname]["formal_charge"]
            atom[f"{rname}_n_ff_charge"] = n_t_atoms[aname]["ff_charge"]
         else:
            atom[f"{rname}_n_formal_charge"] = None
            atom[f"{rname}_n_ff_charge"] = None   
            
         if aname in c_t_atoms:
            atom[f"{rname}_c_formal_charge"] = c_t_atoms[aname]["formal_charge"]
            atom[f"{rname}_c_ff_charge"] = c_t_atoms[aname]["ff_charge"]
         else:
            atom[f"{rname}_c_formal_charge"] = None
            atom[f"{rname}_c_ff_charge"] = None      
   
   for rname,atoms in molecules.items():
      for aname,atom in atoms.items():
         if aname == "N":
            atom["connectivity"].append("R*")
            atom["bond_type"].append("1")
         elif aname == "C":
            atom["connectivity"].append("L*")
            atom["bond_type"].append("1")
   
   molecules["CYS"]["SG"]["CYM_formal_charge"] = -1            

   molecules["HIE"]["ND1"]["HIP_formal_charge"] = 1
   molecules["HIE"]["ND1"]["HIP_c_formal_charge"] = 1
   molecules["HIE"]["ND1"]["HIP_n_formal_charge"] = 1
   
   molecules["CYS"]["N"]["CYX_n_formal_charge"] = 1
   molecules["CYS"]["OC2"]["CYX_c_formal_charge"] = -1
   
   molecules["HIE"]["N"]["HIP_n_formal_charge"] = 1
   molecules["HIE"]["N"]["HID_n_formal_charge"] = 1
   
   molecules["HIE"]["OC2"]["HIP_c_formal_charge"] = -1
   molecules["HIE"]["OC2"]["HID_c_formal_charge"] = -1
   
   for aname, atom in molecules["HIE"].items():
      atom["HID_bond_type"] = deepcopy(atom["bond_type"])
   
   molecules["HIE"]["ND1"]["HID_bond_type"][molecules["HIE"]["ND1"]["connectivity"].index("CE1")] = "1"
   molecules["HIE"]["CE1"]["HID_bond_type"][molecules["HIE"]["CE1"]["connectivity"].index("ND1")] = "1"
   molecules["HIE"]["CE1"]["HID_bond_type"][molecules["HIE"]["CE1"]["connectivity"].index("NE2")] = "2"
   molecules["HIE"]["NE2"]["HID_bond_type"][molecules["HIE"]["NE2"]["connectivity"].index("CE1")] = "2"
   
   for rname in ["ALA","GLY","SER","THR","LEU","ILE","VAL","ASN","GLN","ARG",
                 "HIE","TRP","PHE","TYR","GLU","ASP","LYS","CYS","MET"]:
      
      c_atoms = data[f"C{rname}"]
      n_atoms = data[f"N{rname}"]
      atoms = data[rname]
      molecules[rname]["template"] = {
         rname:[[aname for aname in atoms],"atom_type_name","formal_charge","ff_charge","connectivity","bond_type","IM"],
         f"C{rname}":[[aname for aname in c_atoms], "atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type","RT"],
         f"N{rname}":[[aname for aname in n_atoms], "atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type","LT"],  
         f"H{rname}":[[aname for aname in n_atoms if aname != "H3"], "atom_type_name","formal_charge","h_ff_charge","connectivity","bond_type","LT"],                        
                                      }
      molecules[rname]["property"] = AA_property[rname]   
   
   molecules["PRO"]["template"] = {}
   molecules["PRO"]["template"]["PRO"] = [["N", "CD", "HD1", "HD2", "CG", "HG1", "HG2", "CB", "HB1", "HB2", "CA", "HA", "C", "O"],
                    "atom_type_name","formal_charge","ff_charge","connectivity","bond_type","IM"]
   molecules["PRO"]["template"]["CPRO"] = [["N", "CD", "HD1", "HD2", "CG", "HG1", "HG2", "CB", "HB1", "HB2", "CA", "HA", "C", "OC1", "OC2"],
                    "atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type","RT"]
   molecules["PRO"]["template"]["NPRO"] = [["N", "H1", "H2", "CD", "HD1", "HD2", "CG", "HG1", "HG2", "CB", "HB1", "HB2", "CA", "HA", "C", "O"],
                    "atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type","LT"]
   molecules["PRO"]["property"] = AA_property["PRO"]
   
   molecules["HIE"]["template"]["HID"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "CD2", "HD2", "C", "O"],
                    "HID_atom_type_name","HID_formal_charge","HID_ff_charge","connectivity","HID_bond_type","IM"]
   molecules["HIE"]["template"]["CHID"]=[["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "CD2", "HD2", "C", "OC1", "OC2"],
                    "HID_atom_type_name","HID_c_formal_charge","HID_c_ff_charge","connectivity","HID_bond_type","RT"]
   molecules["HIE"]["template"]["NHID"]=[["N", "H1", "H2", "H3", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "CD2", "HD2", "C", "O"],
                    "HID_atom_type_name","HID_n_formal_charge","HID_n_ff_charge","connectivity","HID_bond_type","LT"]
   molecules["HIE"]["template"]["HIP"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "C", "O"],
                    "HIP_atom_type_name","HIP_formal_charge","HIP_ff_charge","connectivity","bond_type","IM"]
   molecules["HIE"]["template"]["CHIP"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "C", "OC1", "OC2"],
                    "HIP_atom_type_name","HIP_c_formal_charge","HIP_c_ff_charge","connectivity","bond_type","RT"]
   molecules["HIE"]["template"]["NHIP"]=[["N", "H1", "H2", "H3", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "C", "O"],
                    "HIP_atom_type_name","HIP_n_formal_charge","HIP_n_ff_charge","connectivity","bond_type","LT"]
   
   molecules["GLU"]["template"]["GLH"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "OE1", "OE2", "HE2", "C", "O"],
                    "GLH_atom_type_name","GLH_formal_charge","GLH_ff_charge","connectivity","bond_type","IM"]
   molecules["ASP"]["template"]["ASH"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "OD1", "OD2", "HD2", "C", "O"],
                    "ASH_atom_type_name","ASH_formal_charge","ASH_ff_charge","connectivity","bond_type","IM"]
   molecules["LYS"]["template"]["LYN"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "HG1", "HG2", "CD", "HD1", "HD2", "CE", "HE1", "HE2", "NZ", "HZ1", "HZ2", "C", "O"],
                    "atom_type_name","LYN_formal_charge","LYN_ff_charge","connectivity","bond_type","IM"]
   
   molecules["CYS"]["template"]["CYM"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "SG", "C", "O"],
                    "atom_type_name","CYM_formal_charge","CYM_ff_charge","connectivity","bond_type","IM"]
   molecules["CYS"]["template"]["CYX"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "SG", "C", "O"],
                    "CYX_atom_type_name","CYX_formal_charge","CYX_ff_charge","connectivity","bond_type","IM"]
   molecules["CYS"]["template"]["CCYX"] = [["N", "H", "CA", "HA", "CB", "HB1", "HB2", "SG", "C", "OC1", "OC2"],
                    "CYX_atom_type_name","CYX_c_formal_charge","CYX_c_ff_charge","connectivity","bond_type","RT"]
   molecules["CYS"]["template"]["NCYX"] = [["N", "H1", "H2", "H3", "CA", "HA", "CB", "HB1", "HB2", "SG", "C", "O"],
                    "CYX_atom_type_name","CYX_n_formal_charge","CYX_n_ff_charge","connectivity","bond_type","LT"]
   
   molecules["HIS"] = deepcopy(molecules["HIE"])
   del molecules["HIE"]
   molecules["HIS"]["template"]["HIS"] = deepcopy(molecules["HIS"]["template"]["HIE"])
   molecules["HIS"]["template"]["CHIS"] = deepcopy(molecules["HIS"]["template"]["CHIE"])
   molecules["HIS"]["template"]["NHIS"] = deepcopy(molecules["HIS"]["template"]["NHIE"])
   molecules["HIS"]["template"]["HHIS"] = deepcopy(molecules["HIS"]["template"]["HHIE"])

   for atom_name in ["N", "H", "CA", "HA", "CB", "HB1", "HB2", "CG", "CE1", "HE1","HE2", "HD2", "C", "O", "H1","H2","H3", "OC1", "OC2","HD1"]:
      molecules["HIS"][atom_name]["HIP_atom_type_name"] = molecules["HIS"][atom_name]["atom_type_name"]
      molecules["HIS"][atom_name]["HID_atom_type_name"] = molecules["HIS"][atom_name]["atom_type_name"]
   molecules["HIS"]["ND1"]["HIP_atom_type_name"] = "NA"
   molecules["HIS"]["ND1"]["HID_atom_type_name"] = "NA"
   molecules["HIS"]["NE2"]["HIP_atom_type_name"] = "NA"
   molecules["HIS"]["NE2"]["HID_atom_type_name"] = "NB"
   molecules["HIS"]["CD2"]["HIP_atom_type_name"] = "CW"
   molecules["HIS"]["CD2"]["HID_atom_type_name"] = "CV"
    
   
   
   
   #molecules["HIE"]["template"]["HHID"]=[["N", "H1", "H2", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "CD2", "HD2", "C", "O"],
   #                 "HID_atom_type_name","HID_formal_charge","HID_h_ff_charge","connectivity","HID_bond_type","LT"]
   #molecules["HIE"]["template"]["HHIP"]=[["N", "H1", "H2", "H3", "CA", "HA", "CB", "HB1", "HB2", "CG", "ND1", "HD1", "CE1", "HE1", "NE2", "HE2", "CD2", "HD2", "C", "O"],
   #                 "HIP_atom_type_name","HIP_formal_charge","HIP_h_ff_charge","connectivity","bond_type","LT"]
   #molecules["CYS"]["template"]["HCYX"] = [["N", "H1", "H2", "CA", "HA", "CB", "HB1", "HB2", "SG", "C", "O"],
   #                 "CYX_atom_type_name","CYX_formal_charge","CYX_h_ff_charge","connectivity","bond_type","LT"]
   
   
   return molecules

def restructure_rna_data(data):
   molecules = {}
   for rname in ["RA","RU","RG","RC","DA","DT","DG","DC"]:
      if rname in data:
      
         c_atoms = data[f"{rname}3"]
         n_atoms = data[f"{rname}5"]
         end_atoms = data[f"{rname}N"]


         atoms = data[rname]
         ans = list(atoms.keys())

         n_ans = [aname for aname in n_atoms if aname not in ans]
         c_ans = [aname for aname in c_atoms if aname not in ans + n_ans]
         end_ans = [aname for aname in end_atoms if aname not in ans+n_ans+c_ans]
         add_atoms = [[aname,n_atoms[aname]] for aname in n_ans] + [[aname,c_atoms[aname]] for aname in c_ans] + [[aname,end_atoms[aname]] for aname in end_ans]

         molecules[rname] = deepcopy(data[rname])
         for rr in add_atoms:
            molecules[rname][rr[0]] = deepcopy(rr[1])

         
         for aname, atom in molecules[rname].items():
            atom["c_formal_charge"] = 0
            atom["n_formal_charge"] = 0
            atom["end_formal_charge"] = 0
            if aname == 'O1P':
               atom["formal_charge"] = -1
               atom["n_formal_charge"] = -1
               atom["c_formal_charge"] = -1
               atom["end_formal_charge"] = -1
            if aname == "O2P":
               atom["bond_type"][atom["connectivity"].index("P")] = "2"
            if aname == "P":
               atom["bond_type"][atom["connectivity"].index("O2P")] = "2"
               atom["bond_type"].append("1")
               atom["connectivity"].append("R*")
            if aname == "O3'":
               atom["bond_type"].append("1")
               atom["connectivity"].append("L*")
               atom["bond_type"].append("1")
               atom["connectivity"].append("H3T")
            if aname == "O5'":
               atom["bond_type"].append("1")
               atom["connectivity"].append("H5T")

            atom["ff_charge"] = atoms[aname]["ff_charge"] if aname in atoms else 0.000
            atom["c_ff_charge"] = c_atoms[aname]["ff_charge"] if aname in c_atoms else 0.000
            atom["n_ff_charge"] = n_atoms[aname]["ff_charge"] if aname in n_atoms else 0.000
            atom["end_ff_charge"] = end_atoms[aname]["ff_charge"] if aname in end_atoms else 0.000
            atom["atom_type_name"] = atoms[aname]["atom_type_name"] if aname in atoms else None
            atom["c_atom_type_name"] = c_atoms[aname]["atom_type_name"] if aname in c_atoms else None
            atom["n_atom_type_name"] = n_atoms[aname]["atom_type_name"] if aname in n_atoms else None
            atom["end_atom_type_name"] = end_atoms[aname]["atom_type_name"] if aname in end_atoms else None
         molecules[rname]["template"] = {
                                  rname:[ans,"atom_type_name","formal_charge","ff_charge","connectivity","bond_type"],
                                  f"{rname}3":[list(c_atoms.keys()),"c_atom_type_name","c_formal_charge","c_ff_charge","connectivity","bond_type"],
                                  f"{rname}5":[list(n_atoms.keys()),"n_atom_type_name","n_formal_charge","n_ff_charge","connectivity","bond_type"],
                                  f"{rname}N":[list(end_atoms.keys()),"end_atom_type_name","end_formal_charge","end_ff_charge","connectivity","bond_type"]
                                  }
   
   assign_special = {
                     "RA":[["C8","N7"],["C5","C4"],["C6","N1"],["C2","N3"]],
                     "DA":[["C8","N7"],["C5","C4"],["C6","N1"],["C2","N3"]],
                     "RG":[["C8","N7"],["C5","C4"],["C6","O6"],["C2","N3"]],
                     "DG":[["C8","N7"],["C5","C4"],["C6","O6"],["C2","N3"]],
                     "RU":[["C5","C6"],["C4","O4"],["C2","O2"]],
                     "DT":[["C5","C6"],["C4","O4"],["C2","O2"]],
                     "RC":[["C5","C6"],["C4","N3"],["C2","O2"]],
                     }
   for attr in assign_special:
      if attr in molecules:
         tmp_ = assign_special[attr]
         for pair in tmp_:
            molecules[attr][pair[0]]["bond_type"][molecules[attr][pair[0]]["connectivity"].index(pair[1])] = "2"
            molecules[attr][pair[1]]["bond_type"][molecules[attr][pair[1]]["connectivity"].index(pair[0])] = "2"
            molecules[attr][pair[0]]["plate"] = "yes"
            molecules[attr][pair[1]]["plate"] = "yes"
         
   return molecules


def read_amberff_to_ff(atf,nonbf,bondf):
   with open(atf) as inf:
      ats = inf.readlines()
   with open(nonbf) as inf:
      nonbs = inf.readlines()
   with open(bondf) as inf:
      bonds = inf.readlines()
   
   
   nns = []
   for ii,line in enumerate(bonds):
      if line.strip() in ["[ bondtypes ]","[ constrainttypes ]","[ angletypes ]","[ dihedraltypes ]"]:
         nns.append(ii)
   bondterms = bonds[nns[0]+1:nns[1]]
   angleterms = bonds[nns[2]+1:nns[3]]
   improperterms = bonds[nns[3]+1:nns[4]]
   dihedralterms = bonds[nns[4]+1:]
   
   
   at_dict = {}
   for line in ats:
      ss = line.strip().split()
      if len(ss) > 0:
         at_dict[ss[0]] = ss
   
   force_field_dict = {"atomtype":[],"bondterm":[],"angleterm":[],"dihedralterm":[],"improperterm":[]}
   for line in nonbs:
      line = line.strip()
      if line != "" and line[0] not in ["[",";"]:
         ss = line.split()
         force_field_dict["atomtype"].append(["atomtype",ss[0],"LJ12_6",round(float(ss[2]),4),round(float(ss[5])*10,5),round(float(ss[6])/4.184,5),"amber99sb",80,80])
   
   for line in bondterms:
      line = line.strip()
      if line != "" and line[0] not in [";"]:
         ss = line.split()
         force_field_dict["bondterm"].append(["bondterm",ss[0],ss[1],"harmonic",round(float(ss[3])*10,4),round(float(ss[4])/836.8,4),"amber99sb",80,80])
   
   for line in angleterms:
      line = line.strip()
      if line != "" and line[0] not in [";"]:
         ss = line.split()
         force_field_dict["angleterm"].append(["angleterm",ss[0],ss[1],ss[2],"harmonic",round(float(ss[4]),4),round(float(ss[5])/8.368,4),"amber99sb",80,80])
   
   for line in improperterms:
      line = line.strip()
      if line != "" and line[0] not in [";"]:
         ss = line.split()
         for iii in range(4):
            if ss[iii] == "X":
               ss[iii] = "Y"
         force_field_dict["improperterm"].append(["improperterm",ss[2],ss[0],ss[1],ss[3],"amber",round(float(ss[6])/4.184,4),"amber99sb",80,80])
   
   dihe_name = "Y$Y$Y$Y"
   ats = ["Y","Y","Y","Y"]
   paras =[0.0000, 0.0000, 0.0000, 180.0000, 0.0000, 0.0000, 0.0000, 180.0000] 
   for line in dihedralterms:
      line = line.strip()
      if line != "" and line[0] not in [";"]:
         ss = line.split()
         for iii in range(4):
            if ss[iii] == "X":
               ss[iii] = "Y"
         _dihe_name = "$".join(ss[:4])
         if _dihe_name != dihe_name:
            force_field_dict["dihedralterm"].append(["dihedralterm",ats[0],ats[1],ats[2],ats[3],"amber"])
            force_field_dict["dihedralterm"][-1] += paras
            force_field_dict["dihedralterm"][-1] += ["amber99sb",80,80]
            ats = [ss[0],ss[1],ss[2],ss[3]]
            dihe_name = _dihe_name
            paras = [0.0000, 0.0000, 0.0000, 180.0000, 0.0000, 0.0000, 0.0000, 180.0000]
            pn = int(ss[7])
            pp = float(ss[6])/4.184
            pa = float(ss[5])
            if pn > 0:
               paras[(pn-1)*2] = round(pp,4)
               paras[(pn-1)*2+1] = pa
         else:
            pn = int(ss[7])
            pp = float(ss[6])/4.184
            pa = float(ss[5])
            if pn > 0:
               paras[(pn-1)*2] = round(pp,4)
               paras[(pn-1)*2+1] = pa
               

   force_field_dict["dihedralterm"].append(["dihedralterm",ats[0],ats[1],ats[2],ats[3],"amber"])
   force_field_dict["dihedralterm"][-1] += paras
   force_field_dict["dihedralterm"][-1] += ["amber99sb",80,80]   
   
   with open("tmp.ff",'w') as outf:
      outf.write("combination_rule LB\nspecial_bond None None 0.8333 None None 0.5\nqmodel none\nequ_table no\n\n\n")
      for attr,terms in force_field_dict.items():
         for term in terms:
            outf.write("     ".join([str(tt) for tt in term]))
            outf.write("\n")
      
         outf.write("\n\n")

def create_amino_acid(ff):
   tt = read_rtp_file(ff)
   molecules = create_data(tt)
   r_molecules = restructure_amino_acid_data(molecules)
   import json
   with open("amino_acid_json.json",'w') as outf:
      outf.write(json.dumps(r_molecules))



if __name__ == "__main__":
   ff = sys.argv[1]
   tt = read_rtp_file(ff)
   molecules = create_data(tt)
   molecules_r = restructure_rna_data(molecules)
   import json
   with open("rna.json",'w') as outf:
      outf.write(json.dumps(molecules_r))