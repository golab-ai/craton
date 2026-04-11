import textwrap
from pathlib import Path

import click
#import pandas as pd
import yaml
import json
#import simplejson

from craton.craton.utils import logger
from craton.application.db import WithDB


@click.group(name="data")
def data():
    """insert or get data from database:

    \b
    * insert, insert data to database
    * get, get data from database
    """

@data.command("insert")
@click.argument(
    "data_type",
    type=click.Choice(["compound","qmdata","prj_vcompound"]),
    )
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-t","--molecule_type",help="molecule type, such as small molecule, peptide, protein, ......",default="molecule",show_default=True)
@click.option("-table","--table_set",help="set mongodb table",default=None,show_default=True)
@click.option("-elem","--element_flag",help="write element symbol for each record, ......",default="F",show_default=True)
def insert_data_to_db(data_type,inputs,molecule_type,table_set,element_flag):
    """
    insert data to db

    \b
    compound: insert molecule to COMPOUND
    qmdata: insert qmdata to QMDATA
    """
    if element_flag == "T":
        flag = True
    else:
        flag = False
    configure = {"compound_style":molecule_type,"element_flag":flag}
    if table_set is not None:
        ss = table_set.split("-")
        for pp in ss:
            configure[pp.split(":")[0]] = pp.split(":")[1]
        
    DBP = WithDB(inputs,config=configure)
    if data_type == "compound":
        DBP.molecule_to_db()
    elif data_type == "qmdata":
        DBP.qmdata_to_db()
    elif data_type == 'prj_vcompound':
        DBP.molecule_to_vcompound()

@data.command("get")
@click.argument(
    "data_type",
    type=click.Choice(["compound","qmdata"]),
    )
@click.option("-i","--inputs",help="file or directory",default=".",show_default=True,)
@click.option("-f","--yaml_files",type=click.File("r"),help="search selector for qmdata",default=None,show_default=True,)
def get_data_from_db(data_type,inputs,yaml_files):
    """
    insert data to db

    \b
    compound: insert molecule to COMPOUND
    qmdata: insert qmdata to QMDATA
    """
    
    if yaml_files != "no":
        selector = yaml.safe_load(yaml_files.read())
    else:
        selector = {"search_method":"r0"}
    DBP = WithDB(inputs,config=selector)
    js = DBP.get_from_db()
    with open("result.json",'w') as outf:
        outf.write(json.dumps(js))
    
