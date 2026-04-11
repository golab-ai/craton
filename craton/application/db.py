import os
from craton import molxpert as MX
from bson.binary import Binary

class WithDB:
    def __init__(self,inputs,config=None):
        self.inputs = inputs
        if config is None:
            self.config = {}
        else:
            self.config = config
        
        self.config = MX.update_configure(self.config)

    def molecule_to_vcompound(self):
        molecules = MX.molecule_create(self.inputs)
        for molecule in molecules:
            imgs = MX.molecule_show(molecule,save_file=True, opath="./IMG",return_flag=True)

            molecule.structure_img = f"{molecule.mole_name}.png"
        MX.insert_to_db("prj_vcompound",molecules)

    def molecule_to_db(self):
        if not isinstance(self.inputs,list):
            self.inputs = [self.inputs]
        for input in self.inputs:
            if os.path.isdir(input):
                ffs = [ff for ff in os.listdir(input)]
                for ff in ffs:
                    ff_path = f"{input}/{ff}"
                    #try:
                    molecules = MX.molecule_create(ff_path,show_figure=False)
                    molecules = MX.molecule_structure(molecules)
                    molecules = MX.molecule_function_group(molecules)
                    molecules = MX.molecule_torsion(molecules)

                        #MX.molecule_to_db(molecules_1,style=self.style)
                    MX.insert_to_db("compound",molecules,config=self.config)
                    #except:
                    #    print("error:",ff_path)
            else:
                #try:
                molecules = MX.molecule_create(input,show_figure=False)
                molecules = MX.molecule_structure(molecules)
                    #molecules = MX.molecule_function_group(molecules)
                molecules = MX.molecule_torsion(molecules,parallel=False)
                MX.insert_to_db("compound",molecules,config=self.config)
        
                    #MX.molecule_to_db(molecules_1,style=self.style)
                #except:
                #    print("error:",input)

    def qmdata_to_db(self):
        extra_config = {"read_optimizing":True}
        if not isinstance(self.inputs,list):
            self.inputs = [self.inputs]
        for input in self.inputs:
            if os.path.isdir(input):
                ffs = [ff for ff in os.listdir(input) if ff.find(".log") != -1]
                for ff in ffs:
                    ff_path = f"{input}/{ff}"
                    try:
                        molecules = MX.molecule_create(ff_path,show_figure=False,extra_var=extra_config)
                        MX.insert_to_db("qmdata",molecules,config=self.config)
                    except:
                        print("error:",ff_path)
            else:
                try:
                    molecules = MX.molecule_create(input,show_figure=False,extra_var=extra_config)
                    MX.insert_to_db("qmdata",molecules,config=self.config)
                except:
                    print("error:",input)
                    
    def get_from_db(self):
        return MX.get_from_db(config=self.config)
