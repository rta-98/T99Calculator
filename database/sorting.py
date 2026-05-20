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

FORBIDDEN = frozenset({ "S", "N" }) # molecules to exclude 

"""
    INPUT: file name, mol objects, and SMILES as separate, equally sized lists
    OUTPUT: flat dictionary of molecular info.
"""
class BytesPDB:

    def __init__(self, 
                 name: list, # Usually derived from Path(<file>).stem and contained in "Molecules"  
                 smiles: list, # Accepts non-canonical SMILES
                 mol: list): # Mol objects stored in RAM  
        # --------------------------------- 
        self.name = name
        self.smiles = smiles
        self.mol = mol
        # ---------------------------------
        self.smiles_names_mols_dict = {
                "Name": [],
                "SMILES": [],
                "Mol Object": []
        }

    def bookeeper(self) -> dict:
        for name, smiles, mol in zip(self.name, self.smiles, self.mol):

            self.smiles_names_mols_dict["Name"].append(name) 
            self.smiles_names_mols_dict["SMILES"].append(InternalValid.validator(smiles)) 
            self.smiles_names_mols_dict["Mol Object"].append(mol) 

        return self.smiles_names_mols_dict 

class MoleculeSorter: 
    def __init__(self, molecule_sorter: BytesPDB): 
        self.molecule_sorter: BytesPDB = molecule_sorter 
        self.imported_mol_data: dict = molecule_sorter.bookeeper() 
        self.name = self.imported_mol_data["Name"]
        self.smiles = self.imported_mol_data["SMILES"] 
        self.mol = self.imported_mol_data["Mol Object"]
        self.mol_sorted_dict: dict = {} # rotatable bonds; S and N devoid
        self.forbidden_dict: dict = {} # rotatable bonds; S and N containing 
        self.non_rot_dict: dict = {} # non-rotatable bonds; S and N devoid
        self.custom_dict: dict = {} # 5/19: Viet wishes to combine non-rotatable bonds, with Sulfur & Nitrogen devoid dataset
        self.dud_list = [] # err 

    def analyze_all(self): 
        for i, (name, smiles, mol) in enumerate(
                zip(self.name, self.smiles, self.mol, strict=True)
        ):
            if not self.has_rot(mol):
                self.non_rot_dict[name] = self.analyzer(name, smiles, mol) 
                continue 
            if not self.has_atom(mol):
                self.mol_sorted_dict[name] = self.analyzer(name, smiles, mol) 
            elif self.has_atom(mol):  
                self.forbidden_dict[name] = self.analyzer(name, smiles, mol)
            if not self.has_rot(mol) and not self.has_atom(mol):
                self.custom_dict[name] = self.anaylzer(name, smiles, mol)

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

    def analyzer(self, name, smiles, mol):
        analyzer_dict = {}
        analyzer_dict["Molecule"] = name
        analyzer_dict["SMILES"] = smiles

        # parsing count_atoms() dict output; appending to analyzer dict
        num_atoms_dict = self.count_atoms(mol) 
        for atom, atom_count in num_atoms_dict.items():
            analyzer_dict[f"{atom} Count"] = atom_count

        # parsing count_motif() dict output; appending to analyzer dict
        motif_dict = self.count_motif(mol) 
        for hybrid, hybrid_vals in motif_dict.items():
            analyzer_dict[hybrid] = hybrid_vals
       
        # parsing count_dihedral() dict output; appending to analyzer dict
        torsions_dict = self.count_dihedral(mol)
        for tor, tor_vals in torsions_dict.items():
            analyzer_dict[tor] = tor_vals

        return analyzer_dict

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
        # creates first dict. (the inner-most, ultimately set to nest inside result)

        bonds_in_mol = {
            "sp3": [],
            "sp2": [],
            "sp" : [],
            "Aromatic": [],
            "Unk": []
        } 

        mol_h = Chem.AddHs(mol) # mol but with added H's

        for mol_bond in mol_h.GetBonds():

            # GetBondType() also produces: GetBeginAtom() and GetEndAtom() as options
            bond_type_mapped = bond_map.get(mol_bond.GetBondType(), "Unk") 

            # defining two points
            atom1 = mol_bond.GetBeginAtom()
            atom2 = mol_bond.GetEndAtom() 

            # convert two points to atomic numbers (Z)
            z1 = atom1.GetAtomicNum() 
            z2 = atom2.GetAtomicNum() 

            # convert atomic numbers (Z) to atomic symbols (E) 
            e1 = Chem.GetPeriodicTable().GetElementSymbol(z1) if z1 > 0 else atom1.GetSmarts() 
            e2 = Chem.GetPeriodicTable().GetElementSymbol(z2) if z2 > 0 else atom2.GetSmarts() 
            
            # append atomic symbols (E) to the bonds_in_mol dict.
            bonds_in_mol[bond_type_mapped].append({
                "Bond Pair ID": "-".join(sorted((e1, e2))),
            }) 
        
        # The nested dictionary result is initialized
        motif_result = {} 

        # Loop over the key (hybrid) and the value (bond_idx) in bonds_in_mol dict.  
        # first for loop simply takes the list associated with one of the 4 hybridization 
        # values, and iterates over the total number of entries, i.e., 
        # for a given sp3 you get 
        # [Bond Pair ID: C-C, Bond Pair ID: C-C, Bond Pair ID: C-N], which has length 3.

        for hybrid_key, bond_list in bonds_in_mol.items():
            if bond_list:
                motif_result[f"Total {hybrid_key} Count"] = len(bond_list) 
                pair_list = []
                pair_count_dict = {}
                for bonds_in_mol_dict in bond_list: 
                    name = bonds_in_mol_dict["Bond Pair ID"] 
                    pair_list.append(name)
                    if name in pair_list:
                        pair_count_dict[name] = pair_count_dict.get(name, 0) + 1 # alternative to dict[key] += 1 
                    else:
                        pair_count_dict[name] = 1
                        
                for pair_name, pair_count in pair_count_dict.items():
                    motif_result[f"{hybrid_key} {pair_name} Count"] = pair_count

        return motif_result 

                #pair_count_dict = dict(Counter(bond["Bond Pair ID"] for bond in bond_list))

    def count_dihedral(self, mol, smiles: Optional[str] = None):
        # Template for rotatable bonds
        ROT_BONDS = str('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
        ROT_BONDS_SMARTS = Chem.MolFromSmarts(ROT_BONDS)

        rot_matches = mol.GetSubstructMatches(ROT_BONDS_SMARTS)
        confs = mol.GetConformer() 
        traversed = set() 
        unique_rot_matches = []

        for j, k in rot_matches: 
            bond = (min(j, k), max(j, k)) # e.g., min(7, 3), max(7, 3) -> (3, 7) 
            if bond not in traversed:
                traversed.add(bond) # object of type set naturally removes duplicates
                unique_rot_matches.append(bond) # now append to list  

        torsions = []
        bonds = []
        torsion_counts = {} 

        for j, k in unique_rot_matches:
            atom_j = mol.GetAtomWithIdx(j)
            atom_k = mol.GetAtomWithIdx(k)

            # Neighbor atoms not including k (left side)
            i = [n.GetIdx() for n in atom_j.GetNeighbors() if n.GetIdx() != k]

            # Neighbor atoms not including j (right side)
            l = [n.GetIdx() for n in atom_k.GetNeighbors() if n.GetIdx() != j] 
            num_tbond = 0 # counter for num. torsions around central bond

            for m, n in product(i, l): # All possible combinations via cartesian product 
                if m != n: # skips identical indices 
                    num_tbond += 1 

                phi = rdMolTransforms.GetDihedralDeg(confs, m, j, k, n) 
                atoms = [mol.GetAtomWithIdx(idx) for idx in (m, j, k, n)]
                atom_symbols = tuple(a.GetSymbol() for a in atoms) 

                torsions.append({
                    "Atoms": atom_symbols,
                    "Phi": round(phi, 4),
                    "Torsion Indices": (m, j, k, n),
                }) 
                
                # label being X-X-X-X, a key for torsion_counts dict. 
                label = "-".join(atom_symbols) 
                
                # if label/key exists, append increment 
                torsion_counts[label] = torsion_counts.get(label, 0) + 1 
        
        is_rot = self.has_rot(mol) 
        rotatable_bonds = len(unique_rot_matches) if is_rot else 0

        # final return dictionary 
        torsions_result = {}

        for label, val in torsion_counts.items():
            torsions_result[f"Torsions {label}"] = val 
        torsions_result["Torsional Axes Count"] = rotatable_bonds 

        return torsions_result
