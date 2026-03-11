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
from collections.abc import Collection 

FORBIDDEN = frozenset({ "S", "N" }) 

class BytesPDB:
    def __init__(self, 
                 data_path: Optional[str] = None, 
                 abbrv: Optional[list[str]] = None, # .csv mol abbreviations 
                 log_abbrvs: Optional[list[str]] = None, 
                 log_mols: Optional[list] = None, 
                 smiles: Optional[list[str]] = None, 
                 fparse: SmileFileParser | None = None): 
        self.data_path = data_path
        self.fparse = fparse or SmileFileParser(data_path) 
        # --------------------------------- 
        self.abbrv = [] if abbrv is None else list(abbrv) 
        self.log_abbrvs = [] if log_abbrvs is None else list(log_abbrvs) 
        self.smiles = [] if smiles is None else list(smiles) 
        self.log_mols = [] if log_mols is None else list(log_mols) 
        # ---------------------------------
        self.pdb_posix = []
        self.pdb_bytes = []
        self.pdb_IUPAC = []
        self.pdb_smarts = []
        self.pdb_smiles = [] 
        self.smi_pdb_dict = {
                "SMILES": [],
                "Mol. Object": [],
                ".log": [],
                "Mol. Formula": [],
        }

    def bookeeper(self) -> dict:
        self.fparse = SmileFileParser(str(self.data_path)) 
        self.fparse.smi_populate()
        for idx, (smiles, abbrv, log_abbrvs, log_mols) in enumerate(zip_longest(self.smiles, self.abbrv, self.log_abbrvs, self.log_mols)):
            self.smi_pdb_dict["SMILES"].append(smiles) 
            self.smi_pdb_dict["Mol. Formula"].append(abbrv) 
            self.smi_pdb_dict["Mol. Object"].append(log_mols) 
            self.smi_pdb_dict[".log"].append(log_abbrvs) 
        return self.smi_pdb_dict 

class MoleculeSorter: 

    def __init__(self, molecule_sorter: BytesPDB): 
        self.molecule_sorter: BytesPDB = molecule_sorter 
        self.imported_mol_data: dict = molecule_sorter.bookeeper() 
        self.iupacs = self.imported_mol_data["Mol. Formula"]
        self.smiles = self.imported_mol_data["SMILES"] 
        self.mols = self.imported_mol_data["Mol. Object"]
        self.logs = self.imported_mol_data[".log"] 
        self.mol_sorted_dict: dict = {} 
        self.forbidden_dict: dict = {} 
        self.non_rot_dict: dict = {} 
        self.dud_list = [] 

    def debug(self):
        for i in self.logs:
            if i is not None: 
                zed = i 
        return self.logs 

    def analyze_all(self): 
        for i, (log, mol, smiles, iupac) in enumerate(
                zip(self.logs, self.mols, self.smiles, self.iupacs, strict=True)
        ):
            if not self.has_rot(mol):
                self.non_rot_dict[log] = self.analyzer(mol, smiles, iupac) 
                continue 
            if not self.has_atom(mol):
                self.mol_sorted_dict[log] = self.analyzer(mol, smiles, iupac) 
            elif self.has_atom(mol):  
                self.forbidden_dict[log] = self.analyzer(mol, smiles, iupac)
        return (self.mol_sorted_dict, self.forbidden_dict, self.non_rot_dict) 
    
    def has_atom(
            self,
            mol: Optional[Mol] = None, 
            smiles: Optional[str] = None, 
            forbidden: Collection[str] = FORBIDDEN) -> bool: 
        return any(atom.GetSymbol() in forbidden for atom in mol.GetAtoms()) 

    def has_rot(
            self, 
            mol: Optional[Mol] = None) -> bool: 
        try: 
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(
                mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
            ) 
        except AttributeError:
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol, strict=True) 
        return n_rot > 0 

#        if mol is None:
#            return False
#        if rdMolDescriptors.CalcNumRotatableBonds(mol) == 0: 
#            return False
#        if require_conformer and mol.GetNumConformers() == 0: 
#            return False 
#        return True 

    def analyzer(self, mol, smiles, iupac):
        return {
            "SMILES": smiles,  
            "IUPAC": iupac,  
            "Num. Atoms": self.count_atoms(mol),
            "Motif": self.count_motif(mol),
            "Torsions": self.count_dihedral(mol)
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
            Chem.BondType.AROMATIC: "Aromatic",
        }
        bonds_in_mol = {
            "sp3": [],
            "sp2": [],
            "sp" : [],
            "Aromatic": [],
            "Unk": []
        } 

        mol_h = Chem.AddHs(mol) 
        for bond_obj in mol_h.GetBonds():
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
                "Pair": "-".join(sorted((e1, e2))),
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

    def count_dihedral(self, mol, smiles: Optional[str] = None):
        ROT_BONDS = str('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
        rot_bonds_smarts = Chem.MolFromSmarts(ROT_BONDS)
        confs = mol.GetConformer() 
        matches = mol.GetSubstructMatches(rot_bonds_smarts)
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
            atom_j = mol.GetAtomWithIdx(j)
            atom_k = mol.GetAtomWithIdx(k)
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
                phi = rdMolTransforms.GetDihedralDeg(confs, m, j, k, n) 
                atoms = [mol.GetAtomWithIdx(idx) for idx in (m, j, k, n)]
                atom_symbols = tuple(a.GetSymbol() for a in atoms) 
                torsions.append({
                    "Atoms": tuple(a.GetSymbol() for a in atoms),
                    "Phi": round(phi, 4),
                    "Torsion Indices": (m, j, k, n),
                }) 
                label = "-".join(atom_symbols) 
                torsion_counts[label] = torsion_counts.get(label, 0) + 1 
            bonds.append({
                "Central Bond": f"{mol.GetAtomWithIdx(j).GetSymbol()}-{mol.GetAtomWithIdx(k).GetSymbol()}",
#                "Central Bond Idx": (j, k),
                "Torsions Per Bond": num_tbond, 
#                "Torsions": torsions
            }) 
        is_rot = self.has_rot(mol) 
        return {
                "Number of Rot. Bonds": len(unique_matches) if is_rot else 0, 
                "Torsion Counts": torsion_counts if is_rot else {}
#                "Torsions Info": bonds,
        }








