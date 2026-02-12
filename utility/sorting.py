from utility.services import * 
from utility.display import *
from utility.smiles import * 
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from pathlib import Path 
from itertools import zip_longest, product 
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
##
#    def mol_pdb(self):
#        for pdb in self.fparse.pdb_files:
#            try: 
#                mol_pdb = Chem.MolFromPDBFile(pdb, sanitize=False, removeHs=False)
#            except OSError: 
#                mol_pdb = None 

    
    def bookeeper(self) -> dict:
        self.fparse = SmileFileParser(str(data_path)) 
        self.fparse.smi_populate()
        for idx, (smiles, pdbs) in enumerate(zip_longest(self.fparse.smiles_list, self.fparse.pdb_files)):
            try: 
                mol_pdb = Chem.MolFromPDBFile(str(pdbs))
            except OSError: 
                mol_pdb = None 
            self.smi_pdb_dict["PDB MOLS"].append(mol_pdb)  
            self.smi_pdb_dict["IDX"].append(idx) 
            self.smi_pdb_dict["SMILES"].append(smiles) 
            self.smi_pdb_dict["PDB Files"].append(str(pdbs)) 
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
#
    def analyze_all(self): 
       for iupac, mol, pdb_mol in zip(self.iupacs, self.mols, self.pdb_mols):
           self.mol_sorted_dict[iupac] = self.analyzer(mol, pdb_mol) 
       return self.mol_sorted_dict 

    def analyzer(self, mol, pdb_mol):
        return {
            "Num. Atoms": self.count_atoms(mol),
            "Connectivity": self.count_motif(mol),
            "Torsions": self.count_dihedral(pdb_mol),
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

    def count_dihedral(self, pdb_mol=None):
        ROT_BONDS = str('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
        DUMMY_MOL = self.pdb_mols[0]
        rot_bonds_smarts = Chem.MolFromSmarts(ROT_BONDS)
        try: 
            confs = pdb_mol.GetConformer() 
            matches = pdb_mol.GetSubstructMatches(rot_bonds_smarts)
            imp_mol = pdb_mol
        except Exception: 
            confs = DUMMY_MOL.GetConformer() 
            matches = DUMMY_MOL.GetSubstructMatches(rot_bonds_smarts) 
            imp_mol = DUMMY_MOL
        traversed = set() 
        unique_matches = []
        for j, k in matches: 
            bond = (min(j, k), max(j, k))
            if bond not in traversed:
                traversed.add(bond)
                unique_matches.append(bond) 
        torsions = []
        for j, k in unique_matches:
            atom_j = imp_mol.GetAtomWithIdx(j)
            atom_k = imp_mol.GetAtomWithIdx(k)
            i = [n.GetIdx() for n in atom_j.GetNeighbors() if n.GetIdx() != k]
            l = [n.GetIdx() for n in atom_k.GetNeighbors() if n.GetIdx() != j]
            num_tbond = 0 
            for m, n in product(i, l):
                if m != n:
                    num_tbond += 1 
            # All possible combinations via cartesian product 
            for m, n in product(i, l):
                if m == n: 
                    continue
                phi = rdMolTransforms.GetDihedralDeg(confs, m, j, k, n)
                atoms = [imp_mol.GetAtomWithIdx(idx) for idx in (m, j, k, n)]
                torsions.append({
                    "Atoms": tuple(a.GetSymbol() for a in atoms),
                    "Phi": round(phi, 4),
                    "Torsion Indices": (m, j, k, n),
                    "Torsions Per Bond": num_tbond
                }) 
        return {"Number of Rot. Bonds": len(unique_matches), "Torsions Info": torsions}

#|%%--%%| <G7fwsrep6J|tMQLbCQO6K>
zed = BytesPDB(data_path) 
inst = MoleculeSorter(zed)
inst.pdb_mols
#|%%--%%| <tMQLbCQO6K|Tk4Rq9bmeg>
zed = BytesPDB(data_path) 
inst = MoleculeSorter(zed)
inst.analyze_all() 
inst.mol_sorted_dict
#|%%--%%| <Tk4Rq9bmeg|Ct5CBZ9ecE>
zed = BytesPDB(data_path) 
zed.bookeeper()
zed.mol_pdb()
inst.analyze_all()
#xyz = zed.smi_pdb_dict
#xyz.keys()
#xyz["PDB MOLS"]


#                "Torsion Idxs": (i,
#            print(f"Rot. Bond: ({j}-{k}), Torsion: ({i},{j},{k},{l}): {angle:.2f}°") 
