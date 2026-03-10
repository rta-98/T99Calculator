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
#|%%--%%| <zWzZr5LIkF|fw77UvbU0w>
base = Path.cwd() 
log_data_path = base / "./qchem_data/log"
log_data_path 
#|%%--%%| <fw77UvbU0w|Z3cJvKMS9u>
class Flattener: 

    def __init__(self, dict_inst: Optional[dict] = None): 
        self.dict_inst = dict_inst 
        self.rows = [] 
        self.items = {}  
        self.final_df: pd.DataFrame = pd.DataFrame() 

    def amalgam(self): 
        for mol_id, payload in self.dict_inst.items():
            row = {"Molecule": mol_id, **self.flatten(payload, sep=" ")}
            self.rows.append(row) 

        final_df = pd.DataFrame(self.rows).fillna(0)
        cols = ["Molecule"] + sorted(c for c in final_df.columns if c != "Molecule")
        final_df = final_df[cols]
        final_df = final_df.rename(columns=self.canon_col) 
        final_df = final_df.rename(columns=self.rename_torsion_cols)
        final_df = final_df.rename(columns=self.rename_motif_cols) 
        final_df = final_df.T.groupby(final_df.columns, sort=False).sum().T
        final_df = final_df.T.groupby(final_df.columns, sort=False).first().T

        return final_df 

    def flatten(self, obj, parent_key="", sep="__"):
        self.items = {} 
        if isinstance(obj, dict):
            for k, v in obj.items():
                new_key = f"{parent_key}{sep}{k}" if parent_key else str(k) 
                self.items.update(self.flatten(v, new_key, sep=sep)) 
        elif isinstance(obj, (list, tuple)): 
            for i, v in enumerate(obj): 
                new_key = f"{parent_key}{sep}{i}" if parent_key else str(i) 
                self.items.update(self.flatten(v, new_key, sep=sep)) 
        else: 
            self.items[parent_key] = obj

        return self.items

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

#|%%--%%| <Z3cJvKMS9u|MjpSa4aA3S>
class LogToMol: 
    def __init__(self, log_list: Optional[list] = None): 
        self.log_list = log_list
        self.log_mols = {
                'Log Files (Rel. Path)': [],
                'Log Mol. Objects': []
        }
        self.log_mols_df: pd.DataFrame = pd.DataFrame() 
    
    def log_mol_list_gen(self) -> dict: 
        fmt = "g09"
        out_dir = Path("./sdf_out")
        out_dir.mkdir(exist_ok=True)
        for log in self.log_list:
            log_file = log_data_path / log
            sdf = out_dir / f"{log_file.stem}.sdf"
            subprocess.run(
                ["obabel", f"-i{fmt}", str(log_file), "-osdf", "-O", str(sdf)],
                check=True
            )
            mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
            self.log_mols['Log Files (Rel. Path)'].append(log_file.stem) 
            self.log_mols['Log Mol. Objects'].append(mol) 
        return self.log_mols 

    def log_mol_df_gen(self):
        mols = self.log_mols
        log_mols_df = pd.DataFrame(mols) 
        return log_mols_df

#|%%--%%| <MjpSa4aA3S|VPRHdGIqUU>
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
#        df_right = self.csv_df.set_index('Log Files (Rel. Path)') 
        self.merged_df = pd.merge(self.csv_df, self.imp_df, on=[f"{self.merge_key}"], how="left") 
        self.merged_df = self.csv_df.join(df_right, how="left", sort=False).reset_index() 
        return self.merged_df

#|%%--%%| <VPRHdGIqUU|IMACPFeKJo>
# Mol objects generated in-situ saved to a single SDF file for reuse purposes 
def save_mols_sdf(mols, path): 
    writer = Chem.SDWriter(str(path)) 
    try: 
        for m in mols: 
            if m is not None:
                writer.write(m) 
    finally: 
        writer.close() 

#|%%--%%| <IMACPFeKJo|2b5Iwe9IVf>
# Mol objects loaded from the single sdf. 
def load_mols_sdf(path): 
    suppl = Chem.SDMolSupplier(str(path), sanitize=True, removeHs=False) 
    return [m for m in suppl if m is not None] 
#|%%--%%| <2b5Iwe9IVf|kJDVW0mpA2>
def csv_generator(df, fname: str, index: Optional[bool]=False):
    filename = f"{fname}.csv" 
    csvdf = df.to_csv(filename, index=index) # include index positional argument for to_csv() 
    return csvdf
#|%%--%%| <kJDVW0mpA2|Q2IXbtVEa9>
# smi_22_df.csv contains smiles, molecule abbreviations, log files and names generated from openbabel for all 221 unique molecules 
smi_221_df = pd.read_csv('~/projects/t99_calc/v1/smi_221_df.csv')
smi_221_df.keys() 
# Convert each df column to a list 
MolAbbrv = smi_221_df['Molecule']
abbrv = MolAbbrv.tolist() 
abbrv

SmilesAbbrv = smi_221_df['Veri. SMILES']
smiles = SmilesAbbrv.tolist() 
smiles 

log_col_rfpath = smi_221_df['Log Path']
log_abbrvs = []
for i in log_col_rfpath:
    log_abbrvs.append(Path(i).name) 
#log_abbrvs = log_col_rfpath.tolist()  
log_abbrvs 

# This is used to generate a list of log fpath names which is used to correlate a list of mol objects. 
#df = pd.read_csv('/mnt/d/academic/tmp/PFAS_data_130_personal.csv') 
rfpath_list = log_abbrvs 
path_names = [] 
for i in rfpath_list: 
    path_name = Path(i).name 
    path_names.append(path_name) 

log_mol_inst = LogToMol(path_names) 
smi_221_mol_dict = log_mol_inst.log_mol_list_gen()
smi_221_mol_dict.keys()
log_mols = smi_221_mol_dict.get('Log Mol. Objects') 
log_mols

#|%%--%%| <Q2IXbtVEa9|MJaYpbuhf9>
# Saving mols to a single sdf file
save_mols_sdf(mols=log_mols, path=(str("./qchem_data/cache/221_log_mol.sdf"))) 

# Loading mols from the single sdf file 
log_mols_sdf = load_mols_sdf("./qchem_data/cache/221_log_mol.sdf") 
#|%%--%%| <MJaYpbuhf9|lD3N60D8bz>
zed = BytesPDB(log_abbrvs=log_abbrvs, smiles=smiles, log_mols=log_mols) 
#|%%--%%| <lD3N60D8bz|Ueh57pSU6k>
inst = MoleculeSorter(zed)
#|%%--%%| <Ueh57pSU6k|5iYrY7LQ8Y>
dict_inst = inst.analyze_all()
print(
    len(dict_inst[1]),
    len(dict_inst[0]) 
    ) 

no_sn = dict_inst[0]
yes_sn = dict_inst[1]

yes_sn
#|%%--%%| <5iYrY7LQ8Y|VtJOwyfZHF>
out_tmp = []
out_smis_list = [] 
for mol, smile in zip(log_mols, smiles):
    out = inst.count_dihedral(mol, smiles)
    out_tmp.append(out['Number of Rot. Bonds']) 
    print(out, '\n') 

dict_inst = inst.mol_sorted_dict
len(dict_inst)
list_inst = inst.dud_list 
#|%%--%%| <VtJOwyfZHF|hZCm21zB9V>
flat_dict_inst = Flattener(dict_inst) 
df3 = flat_dict_inst.amalgam() 
csv_generator(df3, "tmp_df3") 
#|%%--%%| <hZCm21zB9V|UW43u8j1ir>
df2 = pd.DataFrame(dict_inst)

# WORKING  
len(dict_inst)

#|%%--%%| <UW43u8j1ir|orETn1G454>
smi_221_mol_df = log_mol_inst.log_mol_df_gen()
smi_221_mol_df = smi_221_mol_df.rename(columns={"Log Files (Rel. Path)": "Log Name"}) 
smi_221_mol_df.keys()  # this is working

df1_mol_csv = smi_221_mol_df.to_csv('mol_221_df.csv') 
#csv_generator(mol_merge_df, fname='PFAS_data_130_plus_mol') 

#|%%--%%| <orETn1G454|Qdt8pYRChn>
#df_left = df1.set_index('Log Files (Rel. Path)')
#df_right = csv
csv_path = Path('/mnt/d/academic/tmp/PFAS_data_130_personal.csv') 
df2_mol = pd.read_csv(csv_path)
df2_mol

print(df2_mol.head(100).to_string(index=False))

#|%%--%%| <Qdt8pYRChn|NJZceCMubo>
merged_df = pd.merge(smi_221_df, smi_221_mol_df, on="Log Name", how="right") 
merged_df.keys() 
print(merged_df.head(20).to_string(index=False)) 
#|%%--%%| <NJZceCMubo|ecrzYxG8p2>
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

#|%%--%%| <etdZtaVQG3|5uqXfId71G>
# Concatenating Shomate polynomial coefs.  ---------------------------------
fmt = str('g09') 
log_file = Path(str('/qchem_data/log_to_pdb_in/o.log')) 
sdf = log_file.with_suffix('.sdf')
result = subprocess.run(
    ["obabel", f"-i{fmt}", str('/qchem_data/log_to_pdb_in/o.log'), "-osdf", "-O", str(sdf)],
    capture_output=True,
    text=True
)
print("returncode:", result.returncode)
print("stdout:", result.stdout)
print("stderr:", result.stderr)
print("sdf exists:", sdf.exists())
print("sdf size:", sdf.stat().st_size if sdf.exists() else "missing")
