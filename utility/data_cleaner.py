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
import csv 
#|%%--%%| <zNgtESxaFm|QsOzaZYYL6>
base = Path.cwd() 
smi_data_path = base / "./qchem_data/smi"
smi_data_file = smi_data_path / "pfas_smi_log_153.txt"
#smi_data_file.stem

log_data_path = base / "./qchem_data/log"
log_files = list(log_data_path.glob("*.log")) 

csv_data_path = base / "./qchem_data/csv"
csv_data_file = csv_data_path / "nasa7_parms_final.csv"
with open(csv_data_file, 'r') as f:
    nasa7_csv_arr = f.read() 

txt_data_path = base / "./qchem_data/txt"
txt_data_file = txt_data_path / "nasa7_parms.txt"
with open(txt_data_file, 'r') as f:
    nasa7_txt_arr = f.read() 
#|%%--%%| <QsOzaZYYL6|E1I76lVlbT>
smi_nlog_flog = []
log_smis = []
log_names = []
log_files = []

with open(f"{smi_data_file.parent / smi_data_file.stem}.txt") as f:
    for line in f: 
        parts = line.split() 
        log_smis.append(parts[0])
        log_names.append(parts[1])
        log_files.append(parts[2]) 

for idx, (smi, name, file) in enumerate(zip(log_smis, log_names, log_files)):
    smi_nlog_flog.append([idx, smi, name, file])

smidf = pd.DataFrame(smi_nlog_flog, columns=["ID", "SMILES", "Log Names", "Log Files"]) 
smidf = df.drop_duplicates(keep=False) 
pd.set_option("display.max_columns", None) 
print(smidf.to_string(index=False))  
|%%--%%| <E1I76lVlbT|5OwoIpBpEc>
logdf = pd.DataFrame({
    "SMILES": log_smis,
    "Log Names": log_names,
    "Log Paths": log_files,
}) 
print(logdf.to_string(index=False))
#|%%--%%| <5OwoIpBpEc|E43WgvZIQC>
# Parsing nasa7 parameter txt file for abbreviations 
txtdf = pd.read_csv(txt_data_file) 
txtdf.columns
j#|%%--%%| <E43WgvZIQC|HDUs1aDqWw>
# Parsing nasa7 parameter csv file for smiles
csvdf = pd.read_csv(csv_data_file) 
#|%%--%%| <HDUs1aDqWw|D8DWRZPK7T>
# Concatenating a0, a1 ... S(300K), SMILES, Abbreviations, .log names 

key1 = set(csvdf[

