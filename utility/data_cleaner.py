from utility.sorting import * 
from utility.smiles import * 
from openbabel import openbabel as ob
from pathlib import Path
import pandas as pd
from typing import Optional, List 
from pathlib import Path 
import json 
import sqlite3 
import re 
import os
import csv 
#|%%--%%| <zNgtESxaFm|QsOzaZYYL6>
base = Path.cwd() 
print(base)
rpath = str("relative_to(os.getcwd())") 
print(rpath)

smi_data_path = base / "./qchem_data/smi"
smi_data_file = smi_data_path / "pfas_smi_log_153.txt"
print(smi_data_file.stem)

log_data_path = base / "./qchem_data/log"
log_files = list(log_data_path.glob("*.log")) 
print(log_data_path)

csv_data_path = base / "./qchem_data/csv"
csv_data_file = csv_data_path / "nasa7_parms_final.csv"
with open(csv_data_file, 'r') as f:
    nasa7_csv_arr = f.read() 

#txt_data_path = base / "./qchem_data/txt"
#txt_data_file = txt_data_path / "nasa7_parms.txt"
#with open(txt_data_file, 'r') as f:
#    content = f.read()
#    nasa7_txt_arr = content.split('\n')
#arr = [line.strip() for line in nasa7_txt_arr]
#|%%--%%| <QsOzaZYYL6|E1I76lVlbT>
smi_nlog_flog = []
log_smis = []
log_names = []
log_files = []

with open(f"{smi_data_file.parent / smi_data_file.stem}.txt") as f:
    for line in f: 
        parts = line.split() 
        log_names.append(parts[0])
        log_smis.append(parts[1])
        log_files.append(parts[2]) 

for idx, (smi, name, file) in enumerate(zip(log_smis, log_names, log_files)):
    smi_nlog_flog.append([smi, name, file])

#|%%--%%| <E1I76lVlbT|fYhdO5RQ1e>
smidf = pd.DataFrame(smi_nlog_flog, columns=["SMILES", "Log Names", "Log Files"]) 
smidf = smidf.drop_duplicates(keep=False) 
smidf = smidf.map(lambda x: x.strip() if isinstance(x, str) else x) 
#smidf.keys() 
#pd.set_option("display.max_columns", None) 
print(smidf.head(20).to_string(index=False))  
|%%--%%| <fYhdO5RQ1e|5OwoIpBpEc>
# Logdf and Smidf are the same thing
logdf = pd.DataFrame({
    "SMILES": log_smis,
    "Log Names": log_names,
    "Log Paths": log_files,
}) 

logdf = logdf.map(lambda x: x.strip() if isinstance(x, str) else x) 
print(logdf.head(20).to_string(index=False))
##|%%--%%| <5OwoIpBpEc|E43WgvZIQC>
## Parsing nasa7 parameter txt file for abbreviations 
#txtdf = pd.DataFrame(nasa7_txt_arr) 
#txtdf
##txtdf.columns = ['Abbreviation', 'a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)']
##txtdf = txtdf.map(lambda x: x.strip() if isinstance(x, str) else x) 
#txtdf
#print(txtdf.head(20).to_string(index=False))
j#|%%--%%| <E43WgvZIQC|HDUs1aDqWw>
# Parsing nasa7 parameter csv file for smiles
csv = pd.read_csv(csv_data_file) 
csvdf = csv.rename(columns={"S_300K ": "S(300K)", "Log_file": "Log Files"}) 
csvdf.keys()
csvdf["SMILES"] = csvdf["SMILES"].str.strip()
csvdf = csvdf.map(lambda x: x.strip() if isinstance(x, str) else x) 
csvdf.keys()
csvdf[['a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)']] = csvdf[['a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)']].round(6)
print(csvdf.head(20).to_string(index=False))
#|%%--%%| <HDUs1aDqWw|D8DWRZPK7T>
# Parsing .log files form ./qchem_data/log into a single array and df 
def logf_dict_pop(): 
    logf_dict = {
            "frpath_col": [],
            "frpath_stem_col": []}
    for fpath in log_files:
        frpath = f"{fpath.relative_to(os.getcwd())}" 
        frpath_stem = fpath.name
        logf_dict['frpath_stem_col'].append(frpath_stem) 
        logf_dict['frpath_col'].append(frpath) 
    return logf_dict 


logf_dict = logf_dict_pop()
logf_dict 
logfdf = pd.DataFrame(logf_dict)
logfdf = logfdf.rename(columns={"frpath_col": "Log Files (Rel. Path)", "frpath_stem_col": "Log Files"})
logfdfA

#logfdf = logfdf.map(lambda x: x.strip() if isinstance(x, str) else x)
print(logfdf.head(20).to_string(index=False))
#|%%--%%| <D8DWRZPK7T|6C4y7p8lWU>
# Concatenating a0, a1 ... S(300K), SMILES, Abbreviations, .log names 
logcat1 = logfdf.merge(smidf, on="Log Files", how="inner", sort=False) 
logcat1 = logcat1.map(lambda x: x.strip() if isinstance(x, str) else x) 
logcat1.keys() 
print(logcat1.head(20).to_string(index=False))
#|%%--%%| <6C4y7p8lWU|xIFeuvCqnz>
#final_merge = logcat1.merge(csvdf[["Abbreviation", "SMILES", "Log Files", "a0", "a1", "a2", "a3", "a4", "H_f_0K", "S(300K)"]], on=["SMILES", "Log Files"], how="outer")
#logcat3 
#print(logcat3.head(20).to_string(index=False))
#|%%--%%| <xIFeuvCqnz|DjojCOhCLc>
def df_generator(sort: Optional[bool]=False):
    pfas_data_df = pd.merge(csvdf, logcat1, on='SMILES', how='left', sort=sort, suffixes=['_csv', '_txt']) 
    print(pfas_data_df.head(20).to_string(index=False)) 
    return pfas_data_df

def csv_generator(df, fname: str, index: Optional[bool]=False):
    filename = f"{fname}.csv" 
    csvdf = df.to_csv(filename, index=index) # include index positional argument for to_csv() 
    return csvdf

def txt_generator(arr: list, fname: Optional[str]=None):
    with open("missing_log.txt", "w") as txtf: 
        for elm in arr: 
            txtf.write(f"{elm}\n") 

#|%%--%%| <DjojCOhCLc|yzWnLVRcEJ>
pfas_df = df_generator(sort=False)
csv_generator(pfas_df, fname='nasa7_total') 
csv_generator(logfdf, fname='imported_logs')
csv_generator(smidf, fname='imported_smi_nlog_flog')
len(pfas_df)

#|%%--%%| <yzWnLVRcEJ|9ahirKdJXe>
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



