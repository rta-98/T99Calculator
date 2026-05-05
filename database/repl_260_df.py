import pandas as pd
import os
from pathlib import Path 
from rdkit import Chem 
from rdkit.Chem import MolFromSmiles 

#|%%--%%| <X5YfMo6nSA|YtzL7VqY3s>j
# File paths ---------------------------------
base = Path.cwd() 

csv = base / "./qchem_data/csv" 
dup_290_csv = base / "./nasa7_290.csv" 

no_rot = csv / "./no_rot.csv" 
no_sn = csv / "./no_sn.csv" 
yes_sn = csv / "./yes_sn.csv" 

original_nasa7_csv = csv / "./nasa7_parms_final.csv"
rdkit_mol_obab_smi_csv = csv / "./merged_smi_221_mol_221.csv" 

# Junk ---------------------------------


#|%%--%%| <YtzL7VqY3s|RftDWjuiQ4>
# Df generation ---------------------------------

dup_290_df = pd.read_csv(dup_290_csv, dtype={"big_id": "Int64"}) 
dedup_290_df = dup_290_df.drop_duplicates(subset=['S(300K)'], keep='first') 

no_rot = pd.read_csv(no_rot, dtype={"big_id": "Int64"}) 
no_sn = pd.read_csv(no_sn, dtype={"big_id": "Int64"}) 
yes_sn = pd.read_csv(yes_sn, dtype={"big_id": "Int64"}) 

org_nasa_df = pd.read_csv(original_nasa7_csv, dtype={"big_id": "Int64"})
rd_mol_df = pd.read_csv(rdkit_mol_obab_smi_csv, dtype={"big_id": "Int64"})

# Junk ---------------------------------


#|%%--%%| <RftDWjuiQ4|Onesi0p6GT>
# Df modification ---------------------------------

# Renaming Log_file from org_nasa_df to Log File; removing .log from file 
mod_df0 = org_nasa_df["Log_file"].tolist()
new_list0 = [] 

for i in mod_df: 
    fname = Path(i).stem
    new_list0.append(fname) 

mod_df1 = org_nasa_df.drop(columns=['Log_file']) 
mod_df1.insert(int(1), "Log Name", new_list0) 


# Canonicalizing the SMILES strings in nasa7_parms_final.csv for merging 
mod_df2 = org_nasa_df["Veri. SMILES"].tolist()
new_list1 = [] 

for i in mod_df2: 
    mol = Chem.MolFromSmiles(i)
    smile = Chem.MolToSmiles(mol) 
    new_list1.append(smile) 

mod_df1 = org_nasa_df.drop(columns=['Veri. SMILES']) 
mod_df1.insert(int(2), 'Veri. SMILES', new_list1) 

org_nasa_df = mod_df1 

# Junk ---------------------------------


#|%%--%%| <Onesi0p6GT|YjYSETzf3P>
# Merge ---------------------------------

orgin_removed_df = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="inner", sort=False) 

m = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
only_right = m.loc[m["_merge"] == "right_only"].drop(columns="_merge") 

# Junk ---------------------------------
# log_csv_df = pd.merge(csvdf, logfdf, on="Log Files", how="outer", sort=False)

#|%%--%%| <YjYSETzf3P|6MLrIVQF6W>
# Print ---------------------------------

print(dedup_290_df.head(100).to_string(index=False))
len(dedup_290_df)
dup_290_df.keys()

print(no_rot.head(100).to_string(index=False)) 
no_rot.keys()

org_nasa_df.keys()                 
len(org_nasa_df)
rd_mol_df.keys()
len(rd_mol_df) 

orgin_removed_df.keys()
len(orgin_removed_df)
print(orgin_removed_df.head(100).to_string(index=False))

only_right.keys() 
len(only_right) 
print(only_right.head(100).to_string(index=False)) 

print(mod_df1.head(100).to_string(index=False))

# Junk ---------------------------------
