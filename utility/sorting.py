from utility.services import * 
from utility.display import *
from utility.smiles import * 
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from pathlib import Path 
from itertools import zip_longest 
from collections import Counter 
#|%%--%%| <xeQWkYYp1c|zBTodN4fRQ>
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
                "MOLS": [],
                "PDB Bytes": [],
                "PDB Files": [],
                "PDB MOLS": [],
                "PDB Posix": [],
                "IUPAC": [],
                "PDB Bytes": [],
                "IDX": []
        }

    def mol_pdb(self):
        for pdb in self.fparse.pdb_files:
            try: 
                mol_pdb = Chem.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
            except OSError: 
                mol_pdb = None 
            self.smi_pdb_dict["PDB MOLS"].append(mol_pdb)  

    def bookeeper(self) -> dict:
        self.fparse = SmileFileParser(str(data_path)) 
        self.fparse.smi_populate()
        for idx, (smiles, pdbs) in enumerate(zip_longest(self.fparse.smiles_list, self.fparse.pdb_files)):
            self.smi_pdb_dict["IDX"].append(idx) 
            self.smi_pdb_dict["SMILES"].append(smiles) 
            self.smi_pdb_dict["PDB Files"].append(str(pdbs)) 
#            self.smi_pdb_dict["PDB MOLS"].append(pdbs) 
            self.smi_pdb_dict["PDB Posix"].append(pdbs) 
            self.smi_pdb_dict["IUPAC"].append(pdbs.stem) 
            self.smi_pdb_dict["MOLS"].append(Chem.MolFromSmiles(smiles)) 

        return self.smi_pdb_dict 

#|%%--%%| <zBTodN4fRQ|G7fwsrep6J>
class MoleculeSorter: 
    def __init__(self, molecule_sorter: BytesPDB): 
        self.molecule_sorter: BytesPDB = molecule_sorter 
        self.imported_mol_data: dict = molecule_sorter.bookeeper() 
        self.iupacs = self.imported_mol_data["IUPAC"]
        self.mols = self.imported_mol_data["MOLS"]
        self.pdb_mols = self.imported_mol_data["PDB MOLS"]
        self.mol_sorted_dict: dict = {}

    def analyze_all(self): 
        for iupac, mol in zip(self.iupacs, self.mols):
            self.mol_sorted_dict[iupac] = self.analyzer(mol) 
        return self.mol_sorted_dict 

    def analyzer(self, mol):
        return {
            "Num. Atoms": self.count_atoms(mol),
            "Connectivity": self.count_motif(mol),
            "Torsions": self.count_dihedrals(mol),
        }

    def count_atoms(self, mol): 
        cnt = Counter() 
        for atom in mol.GetAtoms():
            Z = atom.GetAtomicNum()
            if Z > 0: 
                E = Chem.GetPeriodicTable().GetElementSymbol(Z)
            else: 
                E = f"query({atom.GetSmarts()})" 
            cnt[E] += 1
        return dict(cnt) 

    def count_motif(self, mol): 
        # creates a map
        bond_map = {
            Chem.BondType.SINGLE: "sp3", 
            Chem.BondType.DOUBLE: "sp2", 
            Chem.BondType.TRIPLE: "sp",
            Chem.BondType.AROMATIC: "aromatic",
        }
        bonds_in_mol = {
            "sp3": [],
            "sp2": [],
            "sp" : [],
            "Aromatic": [],
            "Unk": []
        } 
        for bond_obj in mol.GetBonds():
            # -- Two Steps: --
            # 1. bond.GetBondType() returns RDKit enum
            # 2. .get(...) searchs that enum in the map -- if the bond type isnt in the map then return "other" 
            bond_obj_type = bond_map.get(bond_obj.GetBondType(), "Unk") 
            # defining two points
            atom1 = bond_obj.GetBeginAtom()
            atom2 = bond_obj.GetEndAtom() 
            # convert two points to atomic numbers (Z)
            z1 = atom1.GetAtomicNum() 
            z2 = atom2.GetAtomicNum() 
            # convert atomic numbers (Z) to atomic symbols (E) 
            e1 = Chem.GetPeriodicTable().GetElementSymbol(z1) if z1 > 0 else atom1.GetSmarts() 
            e2 = Chem.GetPeriodicTable().GetElementSymbol(z2) if z2 > 0 else atom2.GetSmarts() 
            # append atomic symbols (E) to the bonds_in_mol dict.
            bonds_in_mol[bond_obj_type].append({
                "Pair": f"{e1}-{e2}",
                "Mol. Obj Idxs": (bond_obj.GetBeginAtomIdx(), bond_obj.GetEndAtomIdx())
            }) 

        result = {} 

        for hybridization, bond_idx in bonds_in_mol.items():
            if bond_idx: 
                result[hybridization] = {
                    "Total Count": len(bond_idx),
                    "Pair Count": dict(Counter(bond["Pair"] for bond in bond_idx)), 
                    "More Details": bond_idx,
                }
        return result 

    def count_dihedrals(self, mol): 
        pass 


#|%%--%%| <G7fwsrep6J|Tk4Rq9bmeg>
zed = BytesPDB(data_path) 
inst = MoleculeSorter(zed)
inst.analyze_all() 
#|%%--%%| <Tk4Rq9bmeg|Ct5CBZ9ecE>
zed = BytesPDB(data_path) 
zed.bookeeper()
zed.mol_pdb() 
xyz = zed.smi_pdb_dict
xyz.keys()
xyz["PDB MOLS"]

