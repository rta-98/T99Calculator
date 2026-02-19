import subprocess
from pathlib import Path 
from utility.services import InternalValid
from rdkit import Chem 
from rdkit.Chem import rdchem
from rdkit.Chem.MolStandardize import rdMolStandardize 
#|%%--%%| <MH2NKWQ4ll|tvdFZzanew>

# valid_smi_inst = InternalValid()

class SmileFileParser(InternalValid): 
    def __init__(self, dir_path):  
        self.dir = dir_path
        self.smiles_list = []
        self.mols_list = []
        self.smi_fname = []
        self.smiles_dud_list = []
        self.pdb_files = []
        self.log_files = []
        self.log_mols = [] 
        self.log_files_mols_dict = {
                ".log path": [],
                "mol": []
        }

    def log_parse(self) -> dict: 
        fmt = "g09" 
        for log_file in Path(self.dir).glob('*.log'): 
            sdf = log_file.with_suffix('.sdf')
            subprocess.run(
                    ["obabel", f"-i{fmt}", str(log_file), "-osdf", "-O", str(sdf)],
                    check=True
            )
            mol = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
            self.log_files_mols_dict[".log path"].append(log_file.stem) 
            self.log_files_mols_dict["mol"].append(mol) 

        return self.log_files_mols_dict 

#|%%--%%| <tvdFZzanew|fQAtR2jmMw>
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



