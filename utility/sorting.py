from utility.services import * 
from utility.display import *
from utility.smiles import * 
from rdkit import Chem
from pathlib import Path 
from itertools import zip_longest 
#|%%--%%| <sbUXtDskGX|02x1QgN63r>
# Only works in script file 
# root_dir = Path(__file__).resolve().parent 

# txt_path = root_dir / "../smi_pdb_data

# --- Imported CLass Ref. ---
#         self.dir = dir_path
#         self.smiles_list = []
#         self.mols_list = []
#         self.smi_fname = []
#         self.smiles_dud_list = []
#         self.pdb_files = []

base = Path.cwd() 
data_path = base / "./smi_pdb_data"


print(data_path)

# for i in data_path.rglob('*.smi'):
#     print(i)


base = Path.cwd() 
data_path = base / "./smi_pdb_data"

class BytesPDB:

    def __init__(self, data_path, fparse: SmileFileParser | None = None): 
        self.data_path = data_path 
        self.fparse = fparse or SmileFileParser(data_path) 
        self.pdb_posix = []
        self.pdb_bytes = []
        self.pdb_IUPAC = []
        self.pdb_smarts = []
        self.pdb_smiles = [] 
        self.smi_pdb_dict = {
                "SMILES": [],
                "SMARTS": [], 
                "PDB Bytes": [],
                "PDB Files": [],
                "PDB Posix": [],
                "IUPAC": [],
                "PDB Bytes": [],
                "IDX": []
        }

    def bookeeper(self):
        self.fparse = SmileFileParser(str(data_path)) 
        self.fparse.smi_populate()
        for idx, (smiles, pdbs) in enumerate(zip_longest(self.fparse.smiles_list, self.fparse.pdb_files)):
            self.smi_pdb_dict["IDX"].append(idx) 
            self.smi_pdb_dict["SMILES"].append(smiles) 
            self.smi_pdb_dict["PDB Files"].append(str(pdbs)) 
            self.smi_pdb_dict["PDB Posix"].append(pdbs) 
            self.smi_pdb_dict["IUPAC"].append(pdbs.stem) 
            
#             self.smi_pdb_dict["PDB Bytes"].append(pdbs.read_bytes()) 
#             self.smi_pdb_dict[i] = {
#                 "SMILES": smiles,
#                 "PDB Files": str(pdbs),
#                 "PDB Posix": [pdbs],
#                 "None Idxs": [
#                     i for i, (smiles, pdbs) in enumerate(zip_longest(self.fparse.smiles_list, self.fparse.pdb_files)) if smiles is None or pdbs is None
#                 ]
#             }
#            self.pdb_posix.append(pdbs) 

    def correlator(self): 
        pass 
        
    
  
#|%%--%%| <02x1QgN63r|xeQWkYYp1c>


#|%%--%%| <xeQWkYYp1c|9OGAofAeBD>

#|%%--%%| <9OGAofAeBD|7wgqOXTKKq>
zed = BytesPDB(data_path) 
zed.bookeeper()
d = zed.smi_pdb_dict
print(d["IDX"], d["SMILES"])
