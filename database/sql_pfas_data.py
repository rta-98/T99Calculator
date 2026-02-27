from utility.sorting import * 
from utility.smiles import * 
from openbabel import openbabel as ob
from pathlib import Path
import pandas as pd
from typing import Optional, List 
from pathlib import Path 
import csv
import json 
import sqlite3 
import re 
import os
#|%%--%%| <UCbKgaEgCp|VAdfGzcLGV>
base = Path.cwd() 
print(base)
rpath = str("relative_to(os.getcwd())") 
print(rpath)

#|%%--%%| <VAdfGzcLGV|0BdP3srebF>
base = Path.cwd() 
smi_data_path = base / "./qchem_data/smi"
smi_data_file = smi_data_path / "pfas_smi_log_153.txt"
print(smi_data_file.stem)

#|%%--%%| <0BdP3srebF|haqcvukBen>
base = Path.cwd() 
csv_data_path = base / "./qchem_data/csv"
csv_data_file = csv_data_path / "nasa7_total.csv"
with open(csv_data_file, 'r') as f:
    nasa7_csv_arr = f.read() 

#|%%--%%| <haqcvukBen|qsozazyyl6>
base = Path.cwd() 
txt_data_path = base / "./qchem_data/txt"
txt_data_file = txt_data_path / "nasa7_parms.txt"
with open(txt_data_file, 'r') as f:
   content = f.read()
   nasa7_txt_arr = content.split('\n')
arr = [line.strip() for line in nasa7_txt_arr]


#|%%--%%| <qsozazyyl6|Z5izUkmsat>
base = Path.cwd() 
log_data_path = base / "./qchem_data/log"
log_files = list(log_data_path.glob("*.log")) 
print(log_data_path)


#|%%--%%| <Z5izUkmsat|bGOK5d7RiK>



