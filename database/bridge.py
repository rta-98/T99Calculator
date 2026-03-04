from database.sorting import * 
from utility.smiles import * 
from utility.services import * 
from utility.display import *
from openbabel import openbabel as ob
from pathlib import Path
import pandas as pd
from typing import Optional, List 
from pathlib import Path 
import json 
import sqlite3 
import re 
import os
from typing import Optional, List 
from dataclasses import dataclass 
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors, rdchem, PandasTools 
from rdkit.Chem.rdchem import Mol
from pathlib import Path 
from itertools import zip_longest, product 
from collections import Counter 
import subprocess
from rdkit.Chem.MolStandardize import rdMolStandardize 
#|%%--%%| <VS1Fzw2SX5|xPKJXsl1PN>
class LogToMol: 
    def __init__(self, log_col_fpath: Optional[list]=None): 
        self.log_col_fpath = log_col_fpath
#        self.log_col = log_col 
        self.log_mols = {
                'Log Files (Rel. Path)': [],
                'Log Mol. Objects': []
        }
        self.log_mols_df: pd.DataFrame = pd.DataFrame() 
    
    def log_mol_list_gen(self) -> dict: 
        fmt = "g09"
        for log_file in self.log_col_fpath:  
            log_file = Path(log_file)
#            print(log_file.with_suffix('.sdf'))
            sdf = log_file.with_suffix('.sdf')
            subprocess.run(
                ["obabel", f"-i{fmt}", str(log_file), "-osdf", "-O", str(sdf)],
                check=True
            )
            mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
            self.log_mols['Log Files (Rel. Path)'].append(log_file) 
            self.log_mols['Log Mol. Objects'].append(mol) 
        return self.log_mols 

    def log_mol_df_gen(self):
        mols = self.log_mols
        log_mols_df = pd.DataFrame(mols) 
        return log_mols_df

#|%%--%%| <xPKJXsl1PN|kJDVW0mpA2>
class AppendToCSV: 
    def __init__(self, 
                 csv_path: Optional[Path] = None,
                 imp_df: Optional[pd.DataFrame] = None,
                 merge_key: Optional[str] = None):
        self.csv_path = csv_path
        self.imp_df = imp_df 
        self.merge_key = merge_key
        self.csv_df: pd.DataFrame = pd.DataFrame() 
        self.merged_df: pd.DataFrame = pd.DataFrame() 

    def converter(self):
        self.csv_df = pd.read_csv(self.csv_path) 
        df_left = df1.set_index('Log Files (Rel. Path)')
        df_right = self.csv_df.set_index('Log Files (Rel. Path)') 
#        self.merged_df = pd.merge(self.csv_df, self.imp_df, on=[f"{self.merge_key}"], how="left") 
        self.merged_df = self.csv_df.join(df_right, how="outer", sort=False).reset_index() 
        return self.merged_df

def csv_generator(df, fname: str, index: Optional[bool]=False):
    filename = f"{fname}.csv" 
    csvdf = df.to_csv(filename, index=index) # include index positional argument for to_csv() 
    return csvdf
#|%%--%%| <kJDVW0mpA2|gKORvpbmVp>
df = pd.read_csv('/mnt/d/academic/tmp/PFAS_data_130_personal.csv') 
MolAbbrv = df['Molecule']
SmilesAbbrv = df['SMILES']
log_col_rfpath = df['Log Files (Rel. Path)']
rfpath_list = log_col_rfpath.tolist()
#|%%--%%| <gKORvpbmVp|orETn1G454>
log_mol_inst = LogToMol(rfpath_list) 
list1 = log_mol_inst.log_mol_list_gen()
list1
#|%%--%%| <orETn1G454|MJaYpbuhf9>
df1 = log_mol_inst.log_mol_df_gen()
df1 # this is working
#csv_generator(mol_merge_df, fname='PFAS_data_130_plus_mol') 
#|%%--%%| <MJaYpbuhf9|ecrzYxG8p2>
df_left = df1.set_index('Log Files (Rel. Path)')
df_right = csv
csv_path = Path('/mnt/d/academic/tmp/PFAS_data_130_personal.csv') 
csv = pd.read_csv(csv_path)
#self.merged_df = pd.merge(self.csv_df, self.imp_df, on="Log Files (Rel. Path), how="inner") 
append_inst = AppendToCSV(csv_path=csv_path, imp_df=df1, merge_key='Log Files (Rel. Path)')
mol_merge_df = append_inst.converter()
mol_merge_df['Log Mol. Objects']
#|%%--%%| <ecrzYxG8p2|9smJxrBgHU>
keys = df.keys().tolist()
with open("pfas_csv_keys_requested.txt", "w", encoding="utf-8") as f: 
    f.writelines(f"{item}\n" for item in keys) 

#|%%--%%| <9smJxrBgHU|etdZtaVQG3>
with open("pfas_SMILES_for InChIKey.txt", "w", encoding="utf-8") as f: 
    f.writelines(f"{item}\n" for item in keys) 
keys
