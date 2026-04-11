#!/usr/bin/env python3
import requests
import sys
import json
import os
from urllib.parse import quote

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

def get_molecule_info(input_str, input_type='smiles'):
    input_encoded = quote(input_str)
    
    properties_url = f"{PUBCHEM_BASE_URL}/compound/{input_type}/{input_encoded}/property/IUPACName,MolecularFormula,MolecularWeight,XLogP,CanonicalSMILES,IsomericSMILES,InChI,InChIKey/JSON"
    synonyms_url = f"{PUBCHEM_BASE_URL}/compound/{input_type}/{input_encoded}/synonyms/JSON"
    
    try:
        props_response = requests.get(properties_url, timeout=30)
        synonyms_response = requests.get(synonyms_url, timeout=30)
        
        props_response.raise_for_status()
        synonyms_response.raise_for_status()
        
        props_data = props_response.json()
        synonyms_data = synonyms_response.json()
        
        if not props_data.get('PropertyTable', {}).get('Properties'):
            print(f"未找到{input_type.upper()}: {input_str}")
            return None
        
        prop = props_data['PropertyTable']['Properties'][0]
        synonyms = synonyms_data.get('InformationList', {}).get('Information', [{}])[0].get('Synonym', [])
        
        cas_numbers = [syn for syn in synonyms if syn.startswith('CAS')]
        
        result = {
            'smiles': prop.get('CanonicalSMILES', 'N/A'),
            'iupac_name': prop.get('IUPACName', 'N/A'),
            'molecular_formula': prop.get('MolecularFormula', 'N/A'),
            'molecular_weight': prop.get('MolecularWeight', 'N/A'),
            'xlogp': prop.get('XLogP', 'N/A'),
            'inchi': prop.get('InChI', 'N/A'),
            'inchi_key': prop.get('InChIKey', 'N/A'),
            'cas_no': cas_numbers[0] if cas_numbers else 'N/A',
            'synonyms': synonyms[:10] if synonyms else [],
            'cid': prop.get('CID', 'N/A')
        }
        
        return result
        
    except requests.exceptions.RequestException as e:
        print(f"请求失败: {e}")
        return None

def print_molecule_info(info):
    if not info:
        return
    
    print("\n" + "="*50)
    print("分子信息")
    print("="*50)
    print(f"SMILES:           {info['smiles']}")
    print(f"IUPAC Name:       {info['iupac_name']}")
    print(f"分子式:          {info['molecular_formula']}")
    print(f"分子量:          {info['molecular_weight']}")
    print(f"XLogP:           {info['xlogp']}")
    print(f"CAS号:           {info['cas_no']}")
    print(f"InChI:           {info['inchi']}")
    print(f"InChI Key:       {info['inchi_key']}")
    print(f"CID:             {info['cid']}")
    
    if info['synonyms']:
        print(f"\n别名 (前10个):")
        for i, syn in enumerate(info['synonyms'], 1):
            print(f"  {i}. {syn}")
    print("="*50 + "\n")

def run(str,typ,print_flag=True,opath="."):
    type_mapping = {
        'smiles': 'smiles',
        'smi': 'smiles',
        'cas': 'name',
        'iupac': 'name',
        'name': 'name'
    }
    
    if typ not in type_mapping:
        print(f"不支持的查询类型: {typ}")
        print("支持的类型: smiles, cas, iupac, name")
        sys.exit(1)
    
    pubchem_type = type_mapping[typ]
    info = get_molecule_info(str, pubchem_type)

    if info:
        if print_flag:
            print_molecule_info(info)
            if not os.path.exists(opath):
                os.makedirs(opath)
            output_file = f"{opath}/molecule_info_{info.get('cid', 'unknown')}.json"
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(info, f, indent=2, ensure_ascii=False)
            print(f"信息已保存到: {output_file}")
    
        return info
    return None


def main():
    if len(sys.argv) < 3:
        print("用法: python pubchem_query.py <查询类型> <查询值>")
        print("\n查询类型:")
        print("  smiles    - 通过SMILES查询")
        print("  cas       - 通过CAS号查询")
        print("  iupac     - 通过IUPAC名称查询")
        print("  name      - 通过别名/通用名查询")
        print("\n示例:")
        print("  python pubchem_query.py smiles CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
        print("  python pubchem_query.py cas 50-78-2")
        print("  python pubchem_query.py iupac 'ibuprofen'")
        print("  python pubchem_query.py name 'aspirin'")
        sys.exit(1)
    
    query_type = sys.argv[1].lower()
    query_value = sys.argv[2]
    
    run(query_value,query_type,print_flag=True)

if __name__ == "__main__":
    main()
