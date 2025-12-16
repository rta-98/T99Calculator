from fastapi import FastAPI, Form, Request
from fastapi.responses import HTMLResponse, RedirectResponse
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles 
# from starlette.middleware.sessions import SessionMiddleware
from pydantic import BaseModel, ConfigDict 
from typing import Optional 
import numpy as np
import os
import rdkit
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit import DataStructs
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
import subprocess
import re
import pubchempy as pubpy

app = FastAPI()
# app.add_middleware() 
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates") 

# class UserProfie(BaseModel):
#     def __init__(self, user_id=None, is_auth=False):
#         self.user_id = user_id
#         self.is_auth = is_auth
#     def save_result     

#|%%--%%| <KnOGX2r6W3|DyftkKJNz8>
class User(BaseModel):
    user_id: Optional[str] = None 
    smi_in: Optional[str] = None
    smi_canon: Optional[str] = None
    mol_in: Optional[rdchem.Mol] = None
    mol_out: Optional[rdchem.Mol] = None

    model_config = ConfigDict(arbitrary_types_allowed=True) 

    class UserInput:
        def __init__(self, smi_in: str):
            self.smi_in = smi_in
            self.user_id = user_id

        def udir(self):
            directory = f"user_{user_id}"
            os.makedirs(directory, 
        def veri_smi(self):
            self.smi_in.strip() == "":
            if not self.smi_in: 
                print("Blank field")
                return
            try:
                smi_mol  = Chem.CannonSmiles(self.smi_in)
                return smi_mol
            except Exception as e: 
                print({f"Invalid SMILES: {e}"}) 
                return None 
        
        def smi_2_mol(self):
            verid_smi = self.veri_smi()
            if verid_smi is None:
                return None
            else:
                mol_obj = Chem.MolFromSmiles(verid_smi) 
                if mol_obj is not None:
                    return verid_smi, mol_obj
                else: 
                    return None

        def mol_2_img(self):
            verid_smi = self.veri_smi()

            try:
                smi_2_iupac = pubpy.get_compounds(self.smi_in, 'smiles')
            except Exception as e:
                return None:
          
           
            
            
#|%%--%%| <DyftkKJNz8|f3uF9mS9gd>


#|%%--%%| <f3uF9mS9gd|ngfGKgl8C5>
obj = User.UserInput("test") 
p1 = obj.veri_smi("test")
print(p1)
#|%%--%%| <ngfGKgl8C5|tKnZ0ck33L>
# Serve "Home" as page 0 
@app.get("/", response_class=HTMLResponse)
async def root(request: Request): 
    return templates.TemplateResponse(
        "index.html", {
        "request": request
        "page": {0: "Home"}
    })

@app.post("/input-smiles")
async def grab_smi(
        request: Request,
        in_smi: string = Form(...),
        '''
        in_smi_csv: string = Form(...),
        ''' 
        note_smile: bool = Form(False)
): 
    '''/input-smiles'''
        
        molStruct = Chem.MolFromSmiles(
        smiles = Chem.MolToSmiles(       
#|%%--%%| <tKnZ0ck33L|UJDrNroodF>



#|%%--%%| <UJDrNroodF|NU92Cne6FE>


#|%%--%%| <NU92Cne6FE|h6D92NwjoB>



