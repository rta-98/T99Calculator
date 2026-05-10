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

unmatch_molecules_csv = csv / "./unmatch_molecules.csv" 

pfas_data_130_csv = csv / "./PFAS_data_130_strip.csv"

Molecule_unmatch_log_csv = csv / "./Molecule_unmatch_log.csv" 

#nasagen_fit_csv = csv / "./nasagen_fit_results.csv" 
nasagen_fit_csv = csv / "./nasagen_fit_130.csv" 

nasa7_242_parms_csv = csv / "./nasa7_242_parms_df.csv"

nasagen_112_clean_csv = csv / "./nasagen_plus_112_clean.csv"

nasa7_202_clean_csv = csv / "./nasa7_202_clean.csv"

# Junk ---------------------------------
#nasagen_fit_csv = csv / "./nasagen_fit_results.csv" 
#|%%--%%| <YtzL7VqY3s|RftDWjuiQ4>
# df generation ---------------------------------
#nasa7_242_parms_df = pd.read_csv(nasa7_242_parms_csv, dtype={"big_id": "Int64"})
#
#no_rot = pd.read_csv(no_rot, dtype={"big_id": "Int64"}) 
#no_sn = pd.read_csv(no_sn, dtype={"big_id": "Int64"}) 
#yes_sn = pd.read_csv(yes_sn, dtype={"big_id": "Int64"}) 
#
#org_nasa_df = pd.read_csv(original_nasa7_csv, dtype={"big_id": "Int64"})
#org_nasa_df_clean = org_nasa_df.drop(columns=['Log_file', 'Veri. SMILES']) 
#rd_mol_df = pd.read_csv(rdkit_mol_obab_smi_csv, dtype={"big_id": "Int64"})
#
#unmatch_molecules_df = pd.read_csv(unmatch_molecules_csv, dtype={"big_id": "Int64"}) 
#
#pfas_data_130_df = pd.read_csv(pfas_data_130_csv, dtype={"big_id": "Int64"})
#
#Molecule_unmatch_log_df = pd.read_csv(Molecule_unmatch_log_csv, dtype={"big_id": "Int64"}) 
#
#nasagen_fit_df = pd.read_csv(nasagen_fit_csv, dtype={"big_id": "Int64"})
#
#nasagen_112_clean_df = pd.read_csv(nasagen_112_clean_csv, dtype={"big_id": "Int64"})
#
nasa7_202_clean_df = pd.read_csv(nasa7_202_clean_csv, dtype={"big_id": "Int64"})
# Junk ---------------------------------
#dup_290_df = pd.read_csv(dup_290_csv, dtype={"big_id": "Int64"}) 
#nasa7_242_parms_df = dedup_290_df
#|%%--%%| <RftDWjuiQ4|Onesi0p6GT>
# Df modification ---------------------------------


# Junk ---------------------------------
#''' Duplicate Removal ''' 
#df_diff_account = df_unmatch_smi_log_parm.drop_duplicates(subset=['Abbreviation'], keep='first') 
#
#''' 
#1. Renaming Log_file from org_nasa_df to Log Name; removing .log from file
#2. Canonicalizing the SMILES column from Tony's 130 PFAS data 
#''' 
#df_log_col = org_nasa_df["Log_file"].tolist()
#df_smile_col = org_nasa_df["Veri. SMILES"].tolist()
#
#new_list0 = [] 
#new_list1 = [] 
#
#for log, smile in zip(df_log_col, df_smile_col): 
#    fname = Path(log).stem
#    new_list0.append(fname) 
#
#    mol = Chem.MolFromSmiles(smile)
#    can_smile = Chem.MolToSmiles(mol
#    new_list1.append(can_smile) 
#
#mod_df = org_nasa_df.drop(columns=['Log_file', 'Veri. SMILES']) 
#mod_df.insert(int(1), "Log Name", new_list0)
#mod_df.insert(int(2), 'Veri. SMILES', new_list1) 
#
#org_nasa_df = mod_df

#nasa_Molecule = nasagen_fit_df['Molecule'].to_list()
#
#new_names = [] 
#
#for old_name in nasa_Molecule:
#    artifact = "_2000_" 
#    new_name = old_name.replace(artifact, "")  
#    new_names.append(new_name)
#   
#old_names = nasagen_fit_df.pop('Molecule') 
#nasagen_fit_df.insert(0, 'Molecule', new_names) 
#
#Log_Name_unmatch_strip = Log_Name_unmatch_df.drop(columns=['a0', 'a1', 'a2', 'a3', 'a4', 'H_f_0K', 'S(300K)'])
#
#nasagen_fit_df.rename(columns={"S": "S(300K)"}) 


#suffix = "_x" 
#unmatch_nasa7_112 = unmatch_nasa7_242.drop(columns=[c for c in unmatch_nasa7_242.columns if c.endswith(suffix)]) 
#unmatch_nasa7_112.rename(columns={c: c[:-2] for c in unmatch_nasa7_112.columns if c.endswith("_y")}, inplace=True)

#unmatch_nasa7_112.drop(columns=['Log Files', 'SMILES']) 
#
#col_molecule = unmatch_molecules_df.pop("Molecule")
#unmatch_molecules_df.insert(0, "Molecule", col_molecule)

#col_log_name = unmatch_molecules_df.pop("Log Name")
#unmatch_molecules_df.insert(0, "Log Name", col_log_name)
#
#i = 3
#key = left.columns[i]
#assert key == right.columns[i]
#pd.merge(left, right, on=key, how="inner")
#
#nasagen_fit_df.head(100)
#nasa7_242_parms_df.head(5)
#
#entropy_col_nasagen = nasagen_fit_df["S(300K)"].to_list() 
#entropy_col_242 = nasa7_242_parms_df["S(300K)"].to_list() 
#entropy_col_112 = unmatch_nasa7_112["S(300K)"].to_list() 
#
#ngen_rounded = [] 
#n7_rounded = [] 
#n7_112_rounded = []
#
#for nasagen in entropy_col_nasagen:
#    ngrnd = round(nasagen, 2)
#    ngen_rounded.append(ngrnd) 
#
#for nasa7 in entropy_col_242: 
#    n7rnd = round(nasa7, 2) 
#    n7_rounded.append(n7rnd) 

#for nasa7 in entropy_col_112: 
#    n7rnd = round(nasa7, 2) 
#    n7_112_rounded.append(n7rnd) 
#
##entropy_col_nasagen = nasagen_fit_df.pop("S(300K)") 
#nasagen_fit_df.insert(7,"S(300K)_rnd", ngen_rounded)
#
##entropy_col_242 = nasa7_242_parms_df.pop("S(300K)")
#nasa7_242_parms_df.insert(8,"S(300K)_rnd", n7_rounded) 
#
#col_112 = unmatch_nasa7_112.pop("S(300K)_rnd") 
#unmatch_nasa7_112.insert(8, "S(300K)_rnd", n7_112_rounded) 
#
#col = org_nasa_df.pop('Abbreviation') 
#org_nasa_df.insert(0, 'Molecule', col) 
#
#col_abbrv = org_nasa_df_clean.pop('Abbreviation')
#org_nasa_df_clean.insert(0, 'Molecule', col_abbrv) 
#
#entropy_unmatch_df["Molecule"] 
#
#nasagen_112_clean_df.drop(columns=['Unnamed: 0'])
#nasagen_112_clean_df.keys()
#|%%--%%| <Onesi0p6GT|sMVNkOfdM0>
# Merge 0 ---------------------------------
org_nasa_df.head(5)
nasa7_242_parms.head(5)

m2_Molecule = pd.merge(org_nasa_df, nasa7_242_parms, on="Molecule", how="right", indicator=True, sort=False) 
unmatch1_nasa7_242 = m2_Molecule.loc[m2_Molecule["_merge"] == "right_only"].drop(columns="_merge") 

unmatch1_nasa7_242.head(100)

# Junk ---------------------------------
# log_csv_df = pd.merge(csvdf, logfdf, on="Log Files", how="outer", sort=False)
#m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
#m_log = pd.merge(org_nasa_df, rd_mol_df, on="Log Name", how="right", indicator=True, sort=False) 
#df_unmatch_smi_parm = m_smi.loc[m_smi["_merge"] == "right_only"].drop(columns="_merge") # 97/221 
#df_unmatch_log_parm = m_log.loc[m_log["_merge"] == "right_only"].drop(columns="_merge") # 97/221 
#df_unmatch_smi_log_parm = pd.concat([df_unmatch_log_parm, df_unmatch_smi_parm], axis=0, ignore_index=True) 

#orgin_removed_df_m_log = pd.merge(org_nasa_df, rd_mol_df, on="Log Name", how="inner", sort=False) # 104/130 
#orgin_removed_df_m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="inner", sort=False) # 127/130 
#
#m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
#df_unmatch_smi_log_parm = m_smi.loc[m_smi["_merge"] == "right_only"].drop(columns="_merge") # 97/221 

# m_entropy = pd.merge(pfas_data_130_df, nasa7_242_parms_df, on="Molecule", how="right", indicator=True, sort=False)
# unmatch_nasa7_242 = m_entropy.loc[m_entropy["_merge"] == "right_only"].drop(columns="_merge") 

#i = 0
#key = unmatch_nasa7_112.columns[i] 
#assert key == unmatch_molecules_df.columns[i]
#m_Molecule = pd.merge(unmatch_nasa7_112, unmatch_molecules_df, on=key, how='right', indicator=True, sort=False) 
#Molecule_unmatch = m_Molecule.loc[m_Molecule['_merge'] == "right_only"].drop(columns='_merge')

#i = 0
#key = unmatch_nasa7_112.columns[i] 
#assert key == Molecule_unmatch_log_df.columns[i]
#m_Log_Name = pd.merge(unmatch_nasa7_112, Molecule_unmatch_log_df, on=key, how='right', indicator=True, sort=False) 
#Log_Name_unmatch_df = m_Log_Name.loc[m_Log_Name['_merge'] == "right_only"].drop(columns='_merge')  
#|%%--%%| <sMVNkOfdM0|YjYSETzf3P>
# Correlating matches; merging nasagen fits, and the 60 remaining unmatches ---------------------------------
len(org_nasa_df)
len(nasagen_112_clean_df)
org_nasa_df.head(100)
nasagen_112_clean_df.head(100)

org_nasa_df_clean.keys()
nasagen_112_clean_df.keys()

nasa7_202_df = pd.concat([org_nasa_df_clean, nasagen_112_clean_df], axis=0, ignore_index=True) 

# Junk ---------------------------------
## 1. nasagen_fit_df + Log_Name_unmatch_strip merged on right (unmatched should be the df with log paths to see what remains) 
##    and key='Molecule'
#
#m1_Molecule = pd.merge(nasagen_fit_df, Log_Name_unmatch_strip, on="Molecule", how="right", indicator=True, sort=False) 
#Log_Name_unmatch_df1 = m1_Molecule.loc[m1_Molecule["_merge"] == "right_only"].drop(columns="_merge") 
#len(Log_Name_unmatch_df1)
#
#Log_Name_match_df1 = pd.merge(nasagen_fit_df, Log_Name_unmatch_strip, on="Molecule", how="inner", sort=False)
#
#m2_entropy = pd.merge(nasagen_fit_df, unmatch_nasa7_112, on="S(300K)_rnd", how="right", indicator=True, sort=False)
#entropy_unmatch_df = m2_entropy.loc[m2_entropy["_merge"] == "right_only"].drop(columns="_merge") 
#len(entropy_unmatch_df) # 0
#
#nasagen_unmatch_nasa7_112_merge_df = pd.merge(nasagen_fit_df, unmatch_nasa7_112, on="S(300K)_rnd", sort=False) 

|%%--%%| <YjYSETzf3P|6MLrIVQF6W>
# Conversion ---------------------------------
#df_unmatch_smi_log_parm.to_csv("df_unmatch_smi_log_parm.csv", index=False) 
#unmatch_molecules_df.to_csv("unmatch_molecules_df.csv", index=False) 
#unmatch_nasa7_112.to_csv("unmatch_nasa7_112.csv", index=False)
#Molecule_unmatch.to_csv("Molecule_unmatch.csv", index=False)
#Log_Name_unmatch_df.to_csv("Log_Name_unmatch.csv", index=False) 
Log_Name_match_df1.to_csv("./qchem_data/csv/Log_Name_match.csv", index=False)
nasa7_242_parms_df.to_csv("./qchem_data/csv/nasa7_242_parms_df.csv", index=False)

# Junk ---------------------------------
#|%%--%%| <6MLrIVQF6W|DIpRDtjC7d>
# Correlating unmatches; merge 
''' 
Merging unmatch_molecules_df + nasa7_242_parms_df on "Log Name" and "Molecule" 
Keep "Log Path", "SMILES" and "Veri. SMILES" in the final frame.

0. Cut molecules with known data from nasa7_242_parms_df; step 1 is complete, 
   "Molecule" column was used to make a perfect cut
1. Merge on "Log Name" first, using the same logic in # Merge cell. 
2. The unmatched molecules, then will be matched once more on "Molecule"  
'''

orgin_removed_df_m_log = pd.merge(org_nasa_df, rd_mol_df, on="Log Name", how="inner", sort=False) # 104/130 
orgin_removed_df_m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="inner", sort=False) # 127/130 

m_smi = pd.merge(org_nasa_df, rd_mol_df, on="Veri. SMILES", how="right", indicator=True, sort=False) 
df_unmatch_smi_log_parm = m_smi.loc[m_smi["_merge"] == "right_only"].drop(columns="_merge") # 97/221 


#|%%--%%| <DIpRDtjC7d|0lRSAsXXvY>
# Correlating SMILES with log files. 


#|%%--%%| <0lRSAsXXvY|6IguWlK26Q>
# Print ---------------------------------
print(unmatch_molecules_df.head(10).to_string(index=False))  
print(nasa7_242_parms_df.head(10).to_string(index=False))

print(pfas_data_130_df.head(10).to_string(index=False))

len(unmatch_nasa7_242) 
print(unmatch_nasa7_242.head(10).to_string(index=False))

nasagen_fit_df
# Junk ---------------------------------
#
#print(dedup_290_df.head(10).to_string(index=False))
#len(dedup_290_df)
#dup_290_df.keys()
#
#print(no_rot.head(10).to_string(index=False)) 
#no_rot.keys()
#
#org_nasa_df.keys()                 
#len(org_nasa_df)
#print(org_nasa_df.head(10).to_string(index=False))
#
#rd_mol_df.keys()
#len(rd_mol_df) # 221
#print(rd_mol_df.head(10).to_string(index=False))
#
#orgin_removed_df_m_log.keys()
#len(orgin_removed_df_m_log)
#print(orgin_removed_df_m_log.head(10).to_string(index=False))
#
#orgin_removed_df_m_smi.keys()
#len(orgin_removed_df_m_smi)
#print(orgin_removed_df_m_smi.head(10).to_string(index=False))
#
#df_unmatch_smi_log_parm.keys() 
#len(df_unmatch_smi_log_parm) 
#print(df_unmatch_smi_log_parm.head(10).to_string(index=False)) 
#
#m.keys() 
#len(m)
#print(m.head(10).to_string(index=False)) 
#
#df_diff_account # frd_903_cof_OH 
