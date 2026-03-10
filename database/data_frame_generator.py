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
#|%%--%%| <2hkp2sVB7w|FHoHj5Arg7>
class Flattener: 

    def __init__(self, dict_inst: Optional[dict] = None): 
        self.dict_inst = dict_inst 
        self.rows = [] 
        self.final_df: pd.DataFrame = pd.DataFrame() 

    def amalgam(self): 
        for mol_id, payload in self.dict_inst.items():
            row = {"Molecule": mol_id, **flatten(payload, sep=" ")}
            self.rows.append(row) 

        final_df = pd.DataFrame(rows).fillna(0)
        cols = ["Molecule"] + sorted(c for c in final_df.columns if c != "Molecule")
        final_df = final_df[cols]
        final_df = final_df.rename(columns=self.canon_col) 
        final_df = final_df.rename(columns=self.rename_torsion_cols)
        final_df = final_df.rename(columns=self.rename_motif_cols) 
        final_df = final_df.T.groupby(final_df.columns, sort=False).sum().T
        final_df = final_df.T.groupby(final_df.columns, sort=False).first().T
        return final_df 

    def flatten(self, obj, parent_key="", sep="__"):
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

    def canon_col(self, c):
        m = re.match(r"^(Motif\s+\w+\s+Pair Count)\s+([A-Za-z]+)-([A-Za-z]+)$", c)
        if not m:
            return c
        prefix, a, b = m.groups() 
        a, b = sorted([a, b])
        return f"{prefix} {a}-{b}"

    def rename_torsion_cols(self, c): 
        if "Torsion Counts" in c:
            tail = c.split()[-1]
            return f"{tail} Torsions" 
        return c 

    def rename_motif_cols(self, c):
        if c.startswith("Motif "):
            return c.replace("Motif ", "", 1) 
        return c 

#|%%--%%| <FHoHj5Arg7|mD754zr3PI>
df = pd.read_csv("/home/tau/projects/t99_calc/v1/qchem_data/csv/merged_smi_221_mol_221.csv")
df.keys()
df.drop(columns=["Unnamed: 0"]) 
abbrv_col = df['Molecule'].tolist() 
log_col = df['Log Name'].tolist()
log_mols = df['Log Mol. Objects'].tolist() 
smiles_col = df['SMILES'].tolist() 
#|%%--%%| <mD754zr3PI|LhPKZvIDLK>
zed = BytesPDB(abbrv=abbrv_col, log_abbrvs=log_col, smiles=smiles_col, log_mols=log_mols) 
inst = MoleculeSorter(zed)
inst.analyze_all()j
dict_inst = inst.mol_sorted_dict
#|%%--%%| <LhPKZvIDLK|8tJXQAAfnQ>

###|%%--%%| <8tJXQAAfnQ|EmirmaL8Et>
# --FINAL STEP---
df.to_csv("PFAS_data_130_part2.csv", index=False) 
#|%%--%%| <EmirmaL8Et|pfbkhAiokF>
#data_path 
#inst = BytesPDB(abbrv=abbrv_col, smiles=smiles_col)
#tmp = MoleculeSorter(inst)
#tmp.analyze_all() 
##|%%--%%| <pfbkhAiokF|5ncx47evrh>
#def log_to_pdb(path: Optional[Path] = None, fmt: Optional[str] = None) -> Path:
#   conv = ob.OBConversion()
#   conv.SetInAndOutFormats(fmt, "pdb") 
#   mol = ob.OBMol() 
#   if not conv.ReadFile(mol, str(path)):
#       raise RuntimeError(f"Failed to read {path}") 
#   out = path.with_suffix(".pdb") 
#   conv.WriteFile(mol, str(out)) 
#   return out
#
## log_to_pdb
##|%%--%%| <5ncx47evrh|QVKUWv3fEd>
#results = {"orca": [], "g16": [], "unk": []} 
#for idx, i in enumerate(files): 
#    head = i.read_text(errors="ignore", encoding="utf-8")[:20000]
#    if "O R C A" in head or "ORCA" in head: 
#        results["orca"].append(tuple([idx,i])) 
#    elif "Entering Gaussian" in head or "G16" in head or "Gaussian Inc." in head:
#        results["g16"].append(tuple([idx,i])) 
#    else: 
#        results["unk"].append(tuple([idx,i])) 
#print(len(results["g16"]))
#print(len(files))
##|%%--%%| <QVKUWv3fEd|KR0AMgJf4y>
#logs_dir = Path("/home/tau/projects/t99_calc/v1/static/storage/logs") 
#files = list(logs_dir.glob("*.log")) 
#files 
##|%%--%%| <KR0AMgJf4y|1GSm9aLx4H>
#inst1 = SmileFileParser(logs_dir)
#inst1.log_parse() 
##|%%--%%| <1GSm9aLx4H|ZhIBoaDdkH>
##class LogConverter:
##    def __init__(self, logs_dir: Optional[list] = None):
##        self.logs_dir = logs_dir
##        self.log_files  = log_dir.stem() 
##        self.log_sdfs = []
##        self.sdf_mols = [] 
##        self.
##
##|%%--%%| <ZhIBoaDdkH|agsbb5VZiu>
#logs_dir = Path("./static/storage/logs") 
#files = list(logs_dir.glob("*.log")) 
#log_parser = SmileFileParser(logs_dir)
#temp = log_parser.log_parse()
##|%%--%%| <agsbb5VZiu|3059aPYzyj>
#temp
##|%%--%%| <3059aPYzyj|rmLbzznwNZ>
#mol_col = temp["mol"]
#log_col = temp[".log path"] 
#log_col
##|%%--%%| <rmLbzznwNZ|tMQLbCQO6K>
#zed = BytesPDB() 
#inst = MoleculeSorter(zed)
#inst.pdb_mols
##|%%--%%| <tMQLbCQO6K|Tk4Rq9bmeg>
#base = Path.cwd() 
#data_path = base / "./smi_pdb_data"
#data_path
##|%%--%%| <Tk4Rq9bmeg|LUNbaD5uR7>
##             data_path: Optional[str] = None, 
##             abbrv: Optional[list[str]] = None, # .csv mol abbreviations 
##             log_abbrvs: Optional[list[str]] = None, 
##             smiles: Optional[list[str]] = None, 
#df = pd.read_csv('/home/tau/projects/t99_calc/v1/qchem_data/csv/nasa7_total.csv') 
#log_col = df["Log Files"]
#abbrv_col = df["Abbreviation"]
#smi_col = df["SMILES"]
##|%%--%%| <LUNbaD5uR7|3zYSIqsYe5>
