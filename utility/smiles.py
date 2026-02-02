from pathlib import Path 
from utility.services import InternalValid
from rdkit import Chem 
from rdkit.Chem import rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize 

# valid_smi_inst = InternalValid()

class SmileFileParser(InternalValid): 

    def __init__(self, dir_path):  
        self.dir = dir_path
        self.smiles_list = []
        self.mols_list = []
        self.smi_fname = []
        self.smiles_dud_list = []
        self.pdb_files = []

    def smi_populate(self):
        for smi_file in Path(self.dir).rglob('*.smi'): 
            pdb_file = smi_file.with_suffix('.pdb') 
            with open(smi_file) as f: 
                for line in f: 
                    tmp_smile = line.split()
                    if not tmp_smile: 
                        continue 
                    smiles = tmp_smile[0]
                    try:
                        canon_smi = InternalValid.validator(smiles)  
                        if canon_smi is None: 
                            self.smiles_dud_list.append(smiles)
                            continue 
                        mol = Chem.MolFromSmiles(canon_smi)
                        mol = Chem.AddHs(mol)
                        if mol is not None:
                            self.smiles_list.append(canon_smi)
                            self.mols_list.append(mol)
                            self.smi_fname.append(smi_file) 
                            self.pdb_files.append(pdb_file) 
                        else: 
                            self.smiles_dud_list.append(smiles)
#                         for i in self.smiles_dud_list:
#                             print(self.smiles_dud_list[i])
                    except Exception as e: 
                        self.smiles_dud_list.append(smiles)
                        print(f"Error in smi_populate() {smiles}: {e}") 
                        continue 

# #                     if not parts: 
# #                         continue 
# #                     smiles = parts[0]
# #                     if smiles is not None: 
#                     try:
#                         canon = InternalValid.validator(smiles) 
#                         mol = Chem.MolFromSmiles(canon)
#                         if mol is None: 
#                             self.smiles_dud_list.append(mol)
#                             continue 
#                         self.smiles_list.append(canon) 
#                         self.mols_list.append(mol) 
#                         self.smi_fname.append(smi_file)
#                     except Exception as e: 
#                         self.smiles_dud_list.append(smiles) 
#                         continue 
 #                             print(InternalValid.validator(smiles))  
