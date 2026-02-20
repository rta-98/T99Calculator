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
smi_data_file.stem
#|%%--%%| <QsOzaZYYL6|BOvmqhbdzQ>iu\asssk
cols = []
with open(f"{smi_data_file.parent / smi_data_file.stem}.txt") as f:
    for line in f: 
        parts = line.split() 
        cols.append((parts[0], parts[1], parts[2])) 

