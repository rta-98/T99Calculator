from utility.sorting import * 
import pandas as pd
import json 
import sqlite3 
import re 
import csv 
#|%%--%%| <2hkp2sVB7w|FHoHj5Arg7>
#conn = sqlite3.connect("PFAS_thermal_data.db") 
#df.to_sql("pfas", conn, if_exists="replace", index=False)
#result = pd.read_sql("SELECT * FROM people", conn) 
#print(result) 
#conn.close() 
#
#
#df = pd.read_json("PFAS_data.json") 
#for col in df.columns:
#    df[col] = df[col].apply(lambda x: json.dumps(x) if isinstance(x, (dict, list, tuple)) else x)
#conn = sqlite3.connect("PFAS_thermal_data.db") 
#df.to_sql("pfas", conn, if_exists="replace", index=False) 
#conn.close()
#|%%--%%| <FHoHj5Arg7|TOY63qQAWN>
data = get_pfas_data() 
def flatten(obj, parent_key="", sep="__"):
    items = {} 
    if isinstance(obj, dict):
        for k, v in obj.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else str(k) 
            items.update(flatten(v, new_key, sep=sep)) 
    elif isinstance(obj, (list, tuple)): 
        for i, v in enumerate(obj): 
            new_key = f"{parent_key}{sep}{i}" if parent_key else str(i) 
            items.update(flatten(v, new_key, sep=sep)) 
    else: 
        items[parent_key] = obj
    return items
#|%%--%%| <TOY63qQAWN|8tJXQAAfnQ>
rows = [] 
for mol_id, payload in data.items():
    row = {"Molecule": mol_id, **flatten(payload, sep=" ")}
    rows.append(row) 

#|%%--%%| <8tJXQAAfnQ|iDvk6pEcQj>
def canon_col(c):
    m = re.match(r"^(Motif\s+\w+\s+Pair Count)\s+([A-Za-z]+)-([A-Za-z]+)$", c)
    if not m:
        return c
    prefix, a, b = m.groups() 
    a, b = sorted([a, b])
    return f"{prefix} {a}-{b}"
def rename_torsion_cols(c): 
    if "Torsion Counts" in c:
        tail = c.split()[-1]
        return f"{tail} Torsions" 
    return c 
def rename_motif_cols(c):
    if c.startswith("Motif "):
        return c.replace("Motif ", "", 1) 
    return c 
#|%%--%%| <iDvk6pEcQj|5sgAiiqsae>
df = pd.DataFrame(rows).fillna(0)
cols = ["Molecule"] + sorted(c for c in df.columns if c != "Molecule")
df = df[cols]
df = df.rename(columns=canon_col) 
df = df.rename(columns=rename_torsion_cols)
df = df.rename(columns=rename_motif_cols) 
df = df.groupby(df.columns, axis=1, sort=False).sum() 
df.to_csv("PFAS_data.csv", index=False) 
#|%%--%%| <5sgAiiqsae|mD754zr3PI>
# Reading in new csv file from Tony

df = pd.read_csv("/home/tau/projects/t99_calc/v1/static/storage/spreadsheet/nasa7_parms_final.csv") 
abbrv_col = df['Abbreviation'].tolist() 
smiles_col = df['SMILES'].tolist() 
#|%%--%%| <mD754zr3PI|EmirmaL8Et>j
data_path 
inst = BytesPDB(abbrv=abbrv_col, smiles=smiles_col)
tmp = MoleculeSorter(inst)
tmp.analyze_all() 
##|%%--%%| <emirmal8et|5ncx47evrh>
from openbabel import openbabel as ob
from pathlib import Path
def log_to_pdb(path: Path, fmt: str) -> Path:
   conv = ob.OBConversion()
   conv.SetInAndOutFormats(fmt, "pdb") 
   mol = ob.OBMol() 
   if not conv.ReadFile



