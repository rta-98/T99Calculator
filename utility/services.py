from fastapi import FastAPI, Form, Request, HTTPException
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles 
from pydantic import BaseModel, ConfigDict, Field, field_validator 
from typing import Optional, List  
import io 
import csv
import os
import rdkit
import requests
import glob
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdchem, Draw, PandasTools
import subprocess
import re
import pubchempy as pubpy
from pathlib import Path
from openpyxl import Workbook
from openpyxl.drawing.image import Image as pyxl_img
from openpyxl.utils import get_column_letter as col_char
from PIL import Image as pil_img
#|%%--%%| <KnOGX2r6W3|DyftkKJNz8>
class InputValid(BaseModel):
    smi_in: str = Field(..., min_length=1)

    @field_validator('smi_in')
    @classmethod 
    def veri_smi(cls, v):
        v = v.strip()
        if not v: 
            raise TypeError('SMILES are of type str, not None') 
        smi_2_mol = Chem.MolFromSmiles(v)
        if smi_2_mol is not None:
            cannon_smi = Chem.MolToSmiles(smi_2_mol) 
            return cannon_smi
        else: 
            raise ValueError(f'Invalid SMILES string: {v}')
            
#|%%--%%| <DyftkKJNz8|MrLEK7p93w>
class MolData:
    mol_data = []
    
    def img_2_bytes(self, smi_in):
        mol_obj = Chem.MolFromSmiles(smi_in) 
        if mol_obj is None: 
            raise ValueError(f'Failed to parse SMILES: {smi_in}')
        try: 
            mol_img = Draw.MolToImage(mol_obj, size=(100,100)) 
            img_bytes = io.BytesIO()
            mol_img.save(img_bytes, format="PNG")
            img_bytes.seek(0)
            return img_bytes.getvalue()
        except Exception as e:
            raise RuntimeError(f'General Error: {e}')

    def get_iupac(self, smi_in) -> Optional[str]:
        try:
            smi_iupac = pubpy.get_compounds(smi_in, 'smiles')
            return smi_iupac[0].iupac_name
        except Exception as e:
            return RuntimeError(f'Failed to determine IUPAC from SMILES: {e}') 

    def append_mol_data(self, smi_in: str) -> dict:
        col_iupac = self.get_iupac(smi_in)
        col_img_bytes = self.img_2_bytes(smi_in)
        
        dict_mol_data = {
            "smi_in": smi_in,
            "iupac": col_iupac or "Unk",
            "img_bytes": col_img_bytes or "Unk"
        }

        self.mol_data.append(dict_mol_data)
        return dict_mol_data 
    def remove_mol_data(self):
        self.mol_data = []
#|%%--%%| <MrLEK7p93w|8EIFGRuXUj>
valid_smi = InputValid(smi_in="OC(C)C") 
print(valid_smi)
mol_obj = MolData()
mol_obj.img_2_bytes(valid_smi.smi_in)
mol_obj.get_iupac(valid_smi.smi_in)
mol_obj.append_mol_data(valid_smi.smi_in)
