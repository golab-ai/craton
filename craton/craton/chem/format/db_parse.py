import datetime
from copy import deepcopy
from pymongo import ReplaceOne
import time
import sys
from pathlib import Path

#from ..database.mongodb import MongoDB
CRATON_DIR = Path(__file__).parent
ROOT_DIR = CRATON_DIR.parent.parent.parent.parent
DATABASE_EXISTS = Path(f"{ROOT_DIR}/database").exists()


def covert_frag_to_dict(molecule):
    elem_frags = []
    for index, frag in getattr(molecule, "elem_frag", {}).items():
        elem_frags.append(
            {
                "frag_type": "elem_frag",
                "index": index,
                "components": frag["components"],
                "offsets": frag["offsets"],
                "connects": frag["connects"],
                "smiles": frag["smiles"],
                "inchi": frag["inchi"],
                "inchi_key": frag["inchi_key"],
                "label": frag["label"],
                "ef_type": frag["ef_type"],
            }
        )
    # Prepare other fragments
    other_frags = {}
    for frag_type in ["csf_frag", "rsf_frag", "sssr_frag", "scaffold_frag", "tf_frag"]:
        
        other_frags[frag_type] = []
        for index, frag in getattr(molecule, frag_type, {}).items():
            # frag.mol
            other_frags[frag_type].append(
                {
                    "frag_type": frag_type,
                    "index": index,
                    "components": frag["components"],
                    "offsets": frag["offsets"],
                    "smiles": frag["smiles"],
                    "inchi": frag["inchi"],
                    "inchi_key": frag["inchi_key"],
                    "label": frag["label"],
                }
            )
    fragments = {"elem_frag": elem_frags}
    fragments.update(other_frags)
    return fragments

class DBData:
    def __init__(self):
        if DATABASE_EXISTS:
            from ....database import DataDB
            self.dbobj = DataDB()
        else:
            sys.exit("数据库模块不存在")

    def _parse(self,inchi_keys,extra_var=None):
        
        datasearch={}
        if extra_var is not None and "datasearch" in extra_var:
            if extra_var["datasearch"] is not None:
                if not isinstance(inchi_keys,list):
                    inchi_keys = [inchi_keys]
                for kk,vv in extra_var["datasearch"].items():
                    datasearch[kk] = vv
                if inchi_keys is not None:
                    datasearch["db_molecules"] = inchi_keys
        else:
            if inchi_keys is not None:
                if not isinstance(inchi_keys,list):
                    inchi_keys = [inchi_keys]
                datasearch = {"db_molecules":inchi_keys,"data_type":"compound","compound_style":"molecule"}
            else:
                datasearch = {"data_type":"compound","compound_style":"molecule"}
        
        docs = self.dbobj.db_get(config={"datasearch":datasearch})
        return docs    
 
    def _convert(self,molecules,extra_var=None):
        pass
 
    def get_mole_from_db(self, inchi_keys, conformer_setting=None):
        """
        get molecule data from compound or fragment database
        also can select the conformation from conformation database
        """
        selector = {"inchi_key": {"$in": inchi_keys}}
        docs = list(self.mongodb.compound_coll.find(selector,{"_id":0,"updated":0}))
        _found = [doc["inchi_key"] for doc in docs]
        _not_found = set(inchi_keys) - set(_found)
        if len(_not_found) > 0:
            docs.extend(
                self.mongodb.compound_pep_coll.find({"inchi_key": {"$in": list(_not_found)},},{"_id":0,"updated":0}))
        
        
        if conformer_setting:
            docs_dict = {doc["inchi_key"]:deepcopy(doc) for doc in docs}
            
            d_inchi_confs = {doc["inchi_key"]: [] for doc in docs}
            selector.update(conformer_setting["conformer_selector"])
            docs = []

            for doc_conf in self.mongodb.conformation_coll.find(selector):
                if doc_conf["inchi_key"] in d_inchi_confs.keys():
                    d_inchi_confs[doc_conf["inchi_key"]].append(doc_conf)

            for doc_confs in d_inchi_confs.values():
                doc_confs.sort(
                    key=lambda x: (x["constrain_term"], x["constrain_value"])
                    if x["search_method"] == "r1"
                    else ([-1, -1, -1, -1], -1)
                )
                for doc in doc_confs:
                    docs.append(deepcopy(docs_dict[doc["inchi_key"]]))
                    for ii in range(len(doc["coordinates"])):
                        docs[-1]["coordinates"][ii] = doc["coordinates"][ii]
                    for attr,value in doc.items():
                        if attr not in ["coordinates","_id","remark"]:
                            docs[-1][attr] = value
        return docs

    def get_mole_from_db(self, inchi_keys, conformer_setting=None):
        """
        get molecule data from compound or fragment database
        also can select the conformation from conformation database
        """
        selector = {"inchi_key": {"$in": inchi_keys}}
        
        docs = list(self.mongodb.fragment_coll.find(selector,{"_id":0,"updated":0}))
        _found = [doc["inchi_key"] for doc in docs]
        _not_found = set(inchi_keys) - set(_found)
        if len(_not_found) > 0:
            docs.extend(
                self.mongodb.compound_coll.find({"inchi_key": {"$in": list(_not_found)},},{"_id":0,"updated":0}))
        if conformer_setting:
            docs_dict = {doc["inchi_key"]:deepcopy(doc) for doc in docs}
            
            d_inchi_confs = {doc["inchi_key"]: [] for doc in docs}
            selector.update(conformer_setting["conformer_selector"])
            docs = []

            for doc_conf in self.mongodb.conformation_coll.find(selector):
                if doc_conf["inchi_key"] in d_inchi_confs.keys():
                    d_inchi_confs[doc_conf["inchi_key"]].append(doc_conf)

            for doc_confs in d_inchi_confs.values():
                doc_confs.sort(
                    key=lambda x: (x["constrain_term"], x["constrain_value"])
                    if x["search_method"] == "r1"
                    else ([-1, -1, -1, -1], -1)
                )
                for doc in doc_confs:
                    docs.append(deepcopy(docs_dict[doc["inchi_key"]]))
                    for ii in range(len(doc["coordinates"])):
                        docs[-1]["coordinates"][ii] = doc["coordinates"][ii]
                    for attr,value in doc.items():
                        if attr not in ["coordinates","_id","remark"]:
                            docs[-1][attr] = value
        return docs

    def _convert(self,molecules,extra_var=None):
        if not isinstance(molecules,list):
            molecules = [molecules]
        fragment_flag = extra_var["fragment_flag"] if extra_var is not None and "fragment_flag" in extra_var \
                        else False
        if fragment_flag:
            coll = self.mongodb.fragment_coll
        else:
            coll = self.mongodb.compound_coll
        batch_num, cnt = 50, 0
        upsert_compound_reqs = []
        for molecule in molecules:
            if fragment_flag:
                resultor = {"source":1,"count":1}
            else:
                resultor = {"source":1}
            op = self._create_insert_data(molecule,fragment_flag=fragment_flag)
            doc = coll.find_one({"inchi_key": molecule.inchi_key, "deleted": False}, resultor)
            if doc is not None and "source" in doc and getattr(molecule,"source",None) is not None:
                op = self._create_update_data(op,molecule,doc,fragment_flag=fragment_flag)


            upsert_compound_reqs.append(ReplaceOne(
                {"inchi_key": molecule.inchi_key, "deleted": False},
                op,
                upsert=True,
            ))
            cnt += 1
            if cnt >= batch_num:
                coll.bulk_write(upsert_compound_reqs)
                upsert_compound_reqs = []
                cnt = 0
        if len(upsert_compound_reqs) > 0:
            coll.bulk_write(upsert_compound_reqs)

    def _create_update_data(self,op,molecule,doc,fragment_flag=False):
        op["source"] = list(set(molecule.source + doc["source"]))
        if fragment_flag:
            if getattr(molecule,"count",None) is not None:
                #op["count"] = doc["count"].update(molecule.count)
                doc["count"].update(molecule.count)
                op["count"] = doc["count"]
            else:
                #op["count"] = doc["count"].update({getattr(molecule,"source")[ii]:1  
                               #for ii in range(len(getattr(molecule,"source")))})
                doc["count"].update({getattr(molecule,"source")[ii]:1  
                               for ii in range(len(getattr(molecule,"source")))})
                op["count"] = doc["count"]
        return op

    def _create_insert_data(self,molecule,fragment_flag=False):
        _attrs = {"ID":"","inchi_key":"","inchi":"","smiles":"",
                "inpac_name":"","nick_name":"","drug_name":"","internal_name":"",
                "function_group":[],"function_group_label":"","formula":"","mass":0.0,"net_charge":-1000,
                "elements":[],"connectivity":[],"bond_type":[],"formal_charge":[],"coordinates":[],
                "heavy_atoms":-1,"torsion_number":-1,
                "ring_number":-1,"ring_size":[],"ring_property":[],
                "topol_label":"","heterocycle":[],"chirality":[],"cis_trans":[],
                "frag_type":"small_molecule",
                "source":[],"deleted":False,"update":datetime.datetime.now()
        }
        speical = ["frag_type","topol_label"]
        op = {attr:getattr(molecule,attr,value) for attr,value in _attrs.items()}
        op["create_user"] = "default"
        ss=["fragments","count","create_user"]
        if not fragment_flag:
            op["fragments"] = covert_frag_to_dict(molecule)
        else:
            if getattr(molecule,"count",None) is not None:
                op["count"] = molecule.count
            else:
                op["count"] = {getattr(molecule,"source")[ii]:1  
                               for ii in range(len(getattr(molecule,"count"))) }
            for attr in speical:
                op[attr] = getattr(molecule,attr) if hasattr(molecule,attr) else None
        return op



class DBData_OLD:
    def __init__(self):
        self.mongodb = MongoDB()

    def _parse(self,inchi_keys,extra_var=None):
        if not isinstance(inchi_keys,list):
            inchi_keys = [inchi_keys]
        conformer_setting = extra_var["conformer_setting"] \
                if extra_var is not None and "conformer_setting" in extra_var \
                else None
        docs = self.get_mole_from_db(inchi_keys,conformer_setting=conformer_setting)
        return docs        

    def get_mole_from_db(self, inchi_keys, conformer_setting=None):
        """
        get molecule data from compound or fragment database
        also can select the conformation from conformation database
        """
        selector = {"inchi_key": {"$in": inchi_keys}}
        docs = list(self.mongodb.compound_coll.find(selector,{"_id":0,"updated":0}))
        _found = [doc["inchi_key"] for doc in docs]
        _not_found = set(inchi_keys) - set(_found)
        if len(_not_found) > 0:
            docs.extend(
                self.mongodb.compound_pep_coll.find({"inchi_key": {"$in": list(_not_found)},},{"_id":0,"updated":0}))
        
        
        if conformer_setting:
            docs_dict = {doc["inchi_key"]:deepcopy(doc) for doc in docs}
            
            d_inchi_confs = {doc["inchi_key"]: [] for doc in docs}
            selector.update(conformer_setting["conformer_selector"])
            docs = []

            for doc_conf in self.mongodb.conformation_coll.find(selector):
                if doc_conf["inchi_key"] in d_inchi_confs.keys():
                    d_inchi_confs[doc_conf["inchi_key"]].append(doc_conf)

            for doc_confs in d_inchi_confs.values():
                doc_confs.sort(
                    key=lambda x: (x["constrain_term"], x["constrain_value"])
                    if x["search_method"] == "r1"
                    else ([-1, -1, -1, -1], -1)
                )
                for doc in doc_confs:
                    docs.append(deepcopy(docs_dict[doc["inchi_key"]]))
                    for ii in range(len(doc["coordinates"])):
                        docs[-1]["coordinates"][ii] = doc["coordinates"][ii]
                    for attr,value in doc.items():
                        if attr not in ["coordinates","_id","remark"]:
                            docs[-1][attr] = value
        return docs

    def get_mole_from_db(self, inchi_keys, conformer_setting=None):
        """
        get molecule data from compound or fragment database
        also can select the conformation from conformation database
        """
        selector = {"inchi_key": {"$in": inchi_keys}}
        
        docs = list(self.mongodb.fragment_coll.find(selector,{"_id":0,"updated":0}))
        _found = [doc["inchi_key"] for doc in docs]
        _not_found = set(inchi_keys) - set(_found)
        if len(_not_found) > 0:
            docs.extend(
                self.mongodb.compound_coll.find({"inchi_key": {"$in": list(_not_found)},},{"_id":0,"updated":0}))
        if conformer_setting:
            docs_dict = {doc["inchi_key"]:deepcopy(doc) for doc in docs}
            
            d_inchi_confs = {doc["inchi_key"]: [] for doc in docs}
            selector.update(conformer_setting["conformer_selector"])
            docs = []

            for doc_conf in self.mongodb.conformation_coll.find(selector):
                if doc_conf["inchi_key"] in d_inchi_confs.keys():
                    d_inchi_confs[doc_conf["inchi_key"]].append(doc_conf)

            for doc_confs in d_inchi_confs.values():
                doc_confs.sort(
                    key=lambda x: (x["constrain_term"], x["constrain_value"])
                    if x["search_method"] == "r1"
                    else ([-1, -1, -1, -1], -1)
                )
                for doc in doc_confs:
                    docs.append(deepcopy(docs_dict[doc["inchi_key"]]))
                    for ii in range(len(doc["coordinates"])):
                        docs[-1]["coordinates"][ii] = doc["coordinates"][ii]
                    for attr,value in doc.items():
                        if attr not in ["coordinates","_id","remark"]:
                            docs[-1][attr] = value
        return docs

    def _convert(self,molecules,extra_var=None):
        if not isinstance(molecules,list):
            molecules = [molecules]
        fragment_flag = extra_var["fragment_flag"] if extra_var is not None and "fragment_flag" in extra_var \
                        else False
        if fragment_flag:
            coll = self.mongodb.fragment_coll
        else:
            coll = self.mongodb.compound_coll
        batch_num, cnt = 50, 0
        upsert_compound_reqs = []
        for molecule in molecules:
            if fragment_flag:
                resultor = {"source":1,"count":1}
            else:
                resultor = {"source":1}
            op = self._create_insert_data(molecule,fragment_flag=fragment_flag)
            doc = coll.find_one({"inchi_key": molecule.inchi_key, "deleted": False}, resultor)
            if doc is not None and "source" in doc and getattr(molecule,"source",None) is not None:
                op = self._create_update_data(op,molecule,doc,fragment_flag=fragment_flag)


            upsert_compound_reqs.append(ReplaceOne(
                {"inchi_key": molecule.inchi_key, "deleted": False},
                op,
                upsert=True,
            ))
            cnt += 1
            if cnt >= batch_num:
                coll.bulk_write(upsert_compound_reqs)
                upsert_compound_reqs = []
                cnt = 0
        if len(upsert_compound_reqs) > 0:
            coll.bulk_write(upsert_compound_reqs)

    def _create_update_data(self,op,molecule,doc,fragment_flag=False):
        op["source"] = list(set(molecule.source + doc["source"]))
        if fragment_flag:
            if getattr(molecule,"count",None) is not None:
                #op["count"] = doc["count"].update(molecule.count)
                doc["count"].update(molecule.count)
                op["count"] = doc["count"]
            else:
                #op["count"] = doc["count"].update({getattr(molecule,"source")[ii]:1  
                               #for ii in range(len(getattr(molecule,"source")))})
                doc["count"].update({getattr(molecule,"source")[ii]:1  
                               for ii in range(len(getattr(molecule,"source")))})
                op["count"] = doc["count"]
        return op

    def _create_insert_data(self,molecule,fragment_flag=False):
        _attrs = {"ID":"","inchi_key":"","inchi":"","smiles":"",
                "inpac_name":"","nick_name":"","drug_name":"","internal_name":"",
                "function_group":[],"function_group_label":"","formula":"","mass":0.0,"net_charge":-1000,
                "elements":[],"connectivity":[],"bond_type":[],"formal_charge":[],"coordinates":[],
                "heavy_atoms":-1,"torsion_number":-1,
                "ring_number":-1,"ring_size":[],"ring_property":[],
                "topol_label":"","heterocycle":[],"chirality":[],"cis_trans":[],
                "frag_type":"small_molecule",
                "source":[],"deleted":False,"update":datetime.datetime.now()
        }
        speical = ["frag_type","topol_label"]
        op = {attr:getattr(molecule,attr,value) for attr,value in _attrs.items()}
        op["create_user"] = "default"
        ss=["fragments","count","create_user"]
        if not fragment_flag:
            op["fragments"] = covert_frag_to_dict(molecule)
        else:
            if getattr(molecule,"count",None) is not None:
                op["count"] = molecule.count
            else:
                op["count"] = {getattr(molecule,"source")[ii]:1  
                               for ii in range(len(getattr(molecule,"count"))) }
            for attr in speical:
                op[attr] = getattr(molecule,attr) if hasattr(molecule,attr) else None
        return op


