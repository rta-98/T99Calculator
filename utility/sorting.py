from utility.services import * 
from utility.display import *
from utility.smiles import * 
from typing import Optional, List 
from dataclasses import dataclass 
from rdkit import Chem
from rdkit.Chem import rdMolTransforms, rdMolDescriptors
from rdkit.Chem.rdchem import Mol
from pathlib import Path 
from itertools import zip_longest, product 
from collections import Counter 

base = Path.cwd() 
data_path = base / "./qchem_data/csv"

class BytesPDB:
    def __init__(self, 
                 data_path: Optional[str] = None, 
                 abbrv: Optional[list[str]] = None, # .csv mol abbreviations 
                 log_abbrvs: Optional[list[str]] = None, 
                 smiles: Optional[list[str]] = None, 
                 fparse: SmileFileParser | None = None): 
        self.data_path = data_path
        self.fparse = fparse or SmileFileParser(data_path) 
        self.abbrv = [] if abbrv is None else list(abbrv) 
        self.log_abbrvs = [] if log_abbrvs is None else list(log_abbrvs) 
        self.smiles = [] if smiles is None else list(smiles) 
        self.pdb_posix = []
        self.pdb_bytes = []
        self.pdb_IUPAC = []
        self.pdb_smarts = []
        self.pdb_smiles = [] 
        self.smi_pdb_dict = {
                "SMILES": [],
                "MOLS": [],
                "Log Files": [],
                "Abbrv. | IUPAC": [],
        }

    def bookeeper(self) -> dict:
        self.fparse = SmileFileParser(str(self.data_path)) 
        self.fparse.smi_populate()
        for idx, (smiles, abbrv, log_abbrvs) in enumerate(zip_longest(self.smiles, self.abbrv, self.log_abbrvs)):
            self.smi_pdb_dict["SMILES"].append(smiles) 
            self.smi_pdb_dict["Abbrv. | IUPAC"].append(abbrv) 
            self.smi_pdb_dict["MOLS"].append(Chem.MolFromSmiles(InternalValid.validator(smiles))) 
            self.smi_pdb_dict["Log Files"].append(log_abbrvs) 
        return self.smi_pdb_dict 

class MoleculeSorter: 
    def __init__(self, molecule_sorter: BytesPDB): 
        self.molecule_sorter: BytesPDB = molecule_sorter 
        self.imported_mol_data: dict = molecule_sorter.bookeeper() 
        self.iupacs = self.imported_mol_data["Abbrv. | IUPAC"]
        self.mols = self.imported_mol_data["MOLS"]
        self.mol_sorted_dict: dict = {}
#
    def analyze_all(self): 
       for i, (iupac, mol) in enumerate(zip_longest(self.iupacs, self.mols)):
           self.mol_sorted_dict[iupac] = self.analyzer(mol) 
       return self.mol_sorted_dict 

    def analyzer(self, mol):
        return {
            "Num. Atoms": self.count_atoms(mol),
            "Motif": self.count_motif(mol),
            "Torsions": self.count_dihedral(mol),
#            "Torsions": self.count_dihedral(mol)
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
#                    "More Details": bond_idx,
                }
        return result 

    def count_dihedral(self, mol=None):
        ROT_BONDS = str('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
        rot_bonds_smarts = Chem.MolFromSmarts(ROT_BONDS)
        dud_list = []
        try: 
#            confs = mol.GetConformer() 
            matches = mol.GetSubstructMatches(rot_bonds_smarts)
            imp_mol = mol
        except ValueError: 
            dud_list.append(mol)
            return 
#            confs = DUMMY_MOL.GetConformer() 
#            matches = DUMMY_MOL.GetSubstructMatches(rot_bonds_smarts) 
#            imp_mol = DUMMY_MOL
        traversed = set() 
        unique_matches = []
        for j, k in matches: 
            bond = (min(j, k), max(j, k))
            if bond not in traversed:
                traversed.add(bond)
                unique_matches.append(bond) 
        torsions = []
        bonds = []
        torsion_counts = {} 
        for j, k in unique_matches:
            atom_j = imp_mol.GetAtomWithIdx(j)
            atom_k = imp_mol.GetAtomWithIdx(k)
            i = [n.GetIdx() for n in atom_j.GetNeighbors() if n.GetIdx() != k] # Neighbor atom that is not k = left side
            l = [n.GetIdx() for n in atom_k.GetNeighbors() if n.GetIdx() != j] # Neighbor atom that is not j = right side 
            num_tbond = 0 # counter for num torsions around central bond
            for m, n in product(i, l): # iterates all neighbor pairs across the bond.
                if m != n: # skips identical indices 
                    num_tbond += 1 
            # All possible combinations via cartesian product 
            for m, n in product(i, l):
                if m == n: # skips degenerate pairs
                    continue
#                phi = rdMolTransforms.GetDihedralDeg(confs, m, j, k, n) 
                atoms = [imp_mol.GetAtomWithIdx(idx) for idx in (m, j, k, n)]
                atom_symbols = tuple(a.GetSymbol() for a in atoms) 
                torsions.append({
                    "Atoms": tuple(a.GetSymbol() for a in atoms),
#                    "Phi": round(phi, 4),
                    "Torsion Indices": (m, j, k, n),
                }) 
                label = "-".join(atom_symbols) 
                torsion_counts[label] = torsion_counts.get(label, 0) + 1 
            bonds.append({
                "Central Bond": f"{imp_mol.GetAtomWithIdx(j).GetSymbol()}-{imp_mol.GetAtomWithIdx(k).GetSymbol()}",
#                "Central Bond Idx": (j, k),
                "Torsions Per Bond": num_tbond, 
#                "Torsions": torsions
            }) 

        return {
                "Number of Rot. Bonds": len(unique_matches), 
                "Torsion Counts": torsion_counts, 
#                "Torsions Info": bonds,
        }





