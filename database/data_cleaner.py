# -- NUMBER 1 --
from database.sorting import * 
from utility.smiles import * 
from utility.services import *
from openbabel import openbabel as ob
from pathlib import Path
from rdkit.Chem import MolFromSmiles 
from rdkit.Chem.rdMolDescriptors import CalcMolFormula 
import pandas as pd
from typing import Optional, List 
from pathlib import Path 
import json 
import sqlite3 
import re 
import os
#|%%--%%| <cig1EKuGeP|iCzaKnuydH>
# -- NUMBER 2 --

base = Path.cwd() 
print(base)
rpath = str("relative_to(os.getcwd())") 

smi_data_path = base / "./qchem_data/smi/log_to_smi_out"
smi_data_file = smi_data_path / "smi_290.txt"

log_data_path = base / "./qchem_data/log"
log_files = list(log_data_path.glob("*.log")) 

csv_data_path = base / "./qchem_data/csv"
csv_data_file = csv_data_path / "nasa7_parms_final.csv" # must be this file path
PFAS_290_csv  = csv_data_path / "nasa7_290.csv"

txt_data_path = base / "./qchem_data/txt" 
nasa7_290 = txt_data_path / "nasa7_290.txt" 

with open(csv_data_file, 'r') as f:
    nasa7_csv_arr = f.read() 

def df_generator(sort: Optional[bool]=False):
    pfas_data_df = pd.merge(csvdf, logcat1, on=['SMILES', 'Log Files'], how='left', sort=sort, suffixes=['_csv', '_txt']) 
    return pfas_data_df

def csv_generator(df, fname: str, index: Optional[bool]=False):
    filename = f"{fname}.csv" 
    csvdf = df.to_csv(filename, index=index) # include index positional argument for to_csv() 
    return csvdf

def txt_generator(arr: list, fname: Optional[str]=None):
    with open("missing_log.txt", "w") as txtf: 
        for elm in arr: 
            txtf.write(f"{elm}\n") 

#|%%--%%| <iCzaKnuydH|D8DWRZPK7T>
# -- NUMBER 3 --

smi_log_290_dict = {
        "SMILES": [],
        "Molecule": [], 
        "Veri. SMILES": [], 
        "Log Name": [],
        "Log Path": [] }
log_smis = []
log_names = []
log_files = []
log_smi_duds = [] 

with open(smi_data_file) as f:
    for line in f: 
        if not line.strip():
            continue 
        parts = line.split() 
        log_smis.append(parts[0])
        log_files.append(parts[1]) 

# Array containing smiles, file, and fname created. 
for idx, (smi, file) in enumerate(zip(log_smis, log_files)):
    fname = Path(file).stem
    file = Path(file) 
    try: 
        vsmi = InternalValid.validator(smi) 
    except ValueError: 
        log_abbrv = None
        continue
    mol_abbrv = CalcMolFormula(MolFromSmiles(vsmi)) 
    smi_log_290_dict["SMILES"].append(smi) 
    smi_log_290_dict["Veri. SMILES"].append(vsmi) 
    smi_log_290_dict["Log Name"].append(fname) 
    smi_log_290_dict["Log Path"].append((file))
    smi_log_290_dict["Molecule"].append(mol_abbrv)  

#|%%--%%| <D8DWRZPK7T|QsOzaZYYL6>
# -- NUMBER 4 --
## Parsing .log files form ./qchem_data/log into a single array and df 
#def logf_dict_pop(): 
#    logf_dict = {
#            "frpath_col": [],
#            "frpath_stem_col": []}
#    for fpath in log_files:
#        frpath = f"{fpath.relative_to(os.getcwd())}" 
#        frpath_stem = fpath.name
#        logf_dict['frpath_stem_col'].append(frpath_stem) 
#        logf_dict['frpath_col'].append(frpath) 
#    return logf_dict 
#


#logf_dict = logf_dict_pop()
#logf_dict 
#logfdf = pd.DataFrame(logf_dict)
smi_290_df = pd.DataFrame(smi_log_290_dict) 
smi_221_df = smi_290_df.drop_duplicates(subset=["Veri. SMILES"]) 
#total_pfas_df_clean = total_pfas_df.drop(columns=['Log Files (Rel. Path)']) 
#|%%--%%| <QsOzaZYYL6|sQYTnalZ3a>
shomate_221 = smi_221_df.drop(columns=['Veri. SMILES', 'Log Path', 'Molecule'])  
csv_generator(shomate_221, "shomate_221") 
print(shomate_221.head(20).to_string(index=False)) 
#|%%--%%| <sQYTnalZ3a|bsu9wHOMlV>
len(smi_221_df)
smi_221_df.keys()

csv_generator(smi_221_df, "smi_221_df") 
#|%%--%%| <bsu9wHOMlV|xg6ouOpjlK>
|%%--%%| <xg6ouOpjlK|E1I76lVlbT>
# Parsing nasa7 parameter csv file for smiles
csv = pd.read_csv(PFAS_290_csv, dtype={"big_id": "Int64"}) 
csvdf = csv.rename(columns={"S_300K ": "S(300K)", "Log_file": "Log Files"}) 
#csvdf.keys()
#csvdf["SMILES"] = csvdf["SMILES"].str.strip()
csvdf = csvdf.map(lambda x: x.strip() if isinstance(x, str) else x) 
csvdf.keys()
#csvdf[['a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)']] = csvdf[['a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)']].round()
print(csvdf.head(290).to_string(index=False))
print(csvdf.keys())
#|%%--%%| <E1I76lVlbT|6C4y7p8lWU>
# -- NUMBER 7 --
# Concatenating a0, a1 ... S(300K), SMILES, Abbreviations, .log names 
log_csv_df = pd.merge(csvdf, logfdf, on="Log Files", how="outer", sort=False
log_csv_df = logcat1.map(lambda x: x.strip() if isinstance(x, str) else x) 
log_csv_df.keys() 
print(log_csv_df.head(100).to_string(index=False))
len(log_csv_df)
#|%%--%%| <6C4y7p8lWU|7qvyDuOXKN>
## -- NUMBER 8 --
merge_part1 = logfdf.merge(csvdf[["Abbreviation", "SMILES", "Log Files", "a0", "a1", "a2", "a3", "a4", "H_f_0K", "S(300K)"]], on=["Log Files"], how="right")
merge_part1 = merge_part1.rename(columns={'Abbreviation': 'Molecule'})
print(merge_part1.head(10).to_string(index=False))
csv_generator(merge_part1, fname='PFAS_data_130_part1')
#|%%--%%| <7qvyDuOXKN|5OwoIpBpEc>
# -- NUMBER 9 --
# part2 contains extra info
merge_part2 = pd.read_csv("PFAS_data_130_part2.csv") 
merge_part2 
## Parsing nasa7 parameter txt file for abbreviations 
#txtdf = pd.DataFrame(nasa7_txt_arr) 
#txtdf
##txtdf.columns = ['Abbreviation', 'a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)']
##txtdf = txtdf.map(lambda x: x.strip() if isinstance(x, str) else x) 
#txtdf
#print(txtdf.head(20).to_string(index=False))
#logcat3 
#print(logcat3.head(20).to_string(index=False))
#pfas_130_df = pd.read_csv('PFAS_data_130.csv') 
#pfas_130_df
##|%%--%%| <DjojCOhCLc|yzpfas_df = df_generator(sort=False)
total_pfas_df = pd.merge(merge_part1, merge_part2, on=['Molecule'], how='outer', sort=False)
total_pfas_df_clean = total_pfas_df.drop(columns=['Log Files (Rel. Path)']) 
total_pfas_df_clean.keys()
##|%%--%%| <DjojCOhCLc|yzpfas_df = df_generator(sort=False)
csv_generator(total_pfas_df_clean, fname='PFAS_data_130')
csv_generator(total_pfas_df, fname='PFAS_data_130_personal')
#tot_df = pd.merge(pfas_df, pfas_130_df, on=['Molecule'], how='outer', sort=False)l_pfas_df = pd.merge(pfas_df, pfas_130_df, on=['Molecule'], how='outer', sort=False)

#|%%--%%| <|9ahirKdJXe>
import csv as csv_mod
def path_matcher(
        csv_path: Optional[str]='./qchem_data/csv/nasa7_parms_final.csv', 
        dir1: Optional[str]='/mnt/d/J', 
        dir2: Optional[str]='/mnt/d/K', 
        col_idx: Optional[int]=0
):
    dirs = [Path(dir1),Path(dir2)]
    csv_path = Path(csv_path)
#    all_files = set() 
#    for d in dirs: 
#        for p in d.rglob("*.log"):
#            if p.is_file():
#                all_files.add(p.name)
    all_files = { 
                p.name
                for d in dirs 
                for p in d.rglob("*.log") 
                if p.is_file() and p.name != ".DS_Store" and not p.name.startswith("._") 
    } 
    with csv_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv_mod.reader(f) 
        next(reader, None)
        results = {} 
        missing_results = [] 
        found_results = [] 
        for row in reader: 
#            print(row)
            if len(row) <= col_idx:
                continue
            name = row[col_idx].strip()
            if not name:
                continue
            status = "FOUND" if name in all_files else "MISSING" 
            if name in all_files: 
                found_results.append(name) 
            else: 
                missing_results.append(name) 
            results = {"Found": found_results,
                       "Missing": missing_results
            } 
    return (results, all_files)

    
#|%%--%%| <9ahirKdJXe|JeDJzfEcwi>
#csv_data_file = str(csv_data_file)
csv_data_file
result, files = path_matcher(col_idx=1) 

#|%%--%%| <JeDJzfEcwi|l0JiY4r8BC>
print(len(result["Missing"]))
print(len(result["Found"]))
print(len(set(all_files)))

miss_arr = result["Missing"]

txt_generator(miss_arr)

#|%%--%%| <l0JiY4r8BC|p32xSXam5j>
dir1= '/mnt/d/J' 
dir2= '/mnt/d/K'
dirs = [Path(dir1), Path(dir2)]
all_files = set() 
all_files = { 
            p.name
            for d in dirs 
            for p in d.rglob("*.log") 
            if p.is_file() and p.name != ".DS_Store" and not p.name.startswith("._") 
} 
all_files

with csv_path.open(newline="") as f:
    reader = csv_mod.reader(f) 
    for row in reader: 
        if len(row) > 1:
            print(row[1].strip()) 




#|%%--%%| <p32xSXam5j|LinBL5yqln>
# Cleaning up txt files 
nasa_df = pd.read_csv(nasa7_290, sep=r"\s+", header=None, names=["Molecule", "a0", "a1", "a2", "a3", "a4", "H_f_0K", "S(300K)"]) 
nasa_df = nasa_df.drop_duplicates() 
nasa_df_mols = nasa_df["Molecule"]

print(nasa_df.head(290).to_string(index=False)) 

#j#|%%--%%| <LinBL5yqln|fYhdO5RQ1e>
## -- NUMBER 5 --
#
## Dataframe generated from .log file directory containing SMILES and .log names
#smidf = pd.DataFrame(smi_nlog_flog, columns=["SMILES", "Log Names", "Log Files"]) 
#smidf = smidf.drop_duplicates(keep=False) 
#smidf = smidf.map(lambda x: x.strip() if isinstance(x, str) else x) 
##smidf.keys() 
##pd.set_option("display.max_columns", None) 
#print(smidf.head(100).to_string(index=False))  


