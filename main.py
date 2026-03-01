from fastapi import FastAPI, Form, Request, Depends, HTTPException
from fastapi.responses import Response, HTMLResponse, RedirectResponse
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
from utility.services import InputValid, MolData, SubstMatch 
from utility.display import FileParser  
from rdkit.Chem.Draw import MolDrawOptions
from utility.smiles import SmileFileParser

OUTPUT_DIR = Path("static/storage/imgs") 
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SIZE = (200,200)
BG_COLOR = (.29, .31, .33)

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates") 
drawOptions = MolDrawOptions() 
drawOptions.setBackgroundColour(BG_COLOR) 

async def user_valid(smi_in = Form(...)) -> InputValid:
    return InputValid.model_validate({"smi_in": smi_in})

# Serve "Home" as page 0 by accessing {{ page[0] }}  
@app.get("/", response_class=HTMLResponse)
async def root(request: Request): 
    context = {"result": None, "error": None} 
    return templates.TemplateResponse(
        "base.html", 
        {
            "request": request,
            "context": context, 
            "page": 1 
        })

# # Canonicalize the smiles string and return to the user in real time as they type
@app.post("/canonical-smiles")
async def note_smile(usrStr: InputValid):
    try:
        result = usrStr.smi_in 
        return {"result": result, "error": None}
    except Exception as e:
        return {"result": None, "error": str(e)} 

# Return a mol image from the SMILES string in real time as the user types
@app.post("/smiles-structure")
async def display_image(usrStr: InputValid):
    mol_data = MolData()  
    try:
        image_bytes = mol_data.img_2_bytes(usrStr.smi_in) 
        return Response(content=image_bytes, media_type='image/png') # gets picked up by input-form.js
    except HTTPException:
        raise # re-raise error so FastAPI can handle it
    except Exception as e:
        raise HTTPException(status_code=500, detail=f'ERROR: {e}') 

# Return a grid of mol pngs 
@app.get("/grid", response_class=HTMLResponse) 
async def grid_page(request: Request): 
    try: 
#         data_dir = Path(__file__).parent / "data" # equivalent to os.path.join(parent, "data") 
        data_dir = Path('./data') # name the path 
        parser = FileParser(str(data_dir)) # assign   
        parser.parse() # instance of FileParser created  
        grid_png = parser.smi_to_png() # a single png of the different mols created 
        grid_png.save("./static/storage/imgs/tmp_grid.png") # save the dir  
    
        return templates.TemplateResponse(
            "display.html",
            {
                "request": request,
                "grid_png": "/static/storage/imgs/tmp_grid.png",
                "count": len(parser.mols), 
                "page": 2
            })

    except Exception as e: 
        raise RuntimeError(f'Error in grid_page(): {e}') 

@app.get("/filter-grid", response_class=HTMLResponse) 
async def filter_grid_page(request: Request): 
    data_dir = Path('./qchem_data/log') # name the path 
    parser = SmileFileParser(str(data_dir)) # assign   
    parser.smi_populate() # instance of FileParser created  
    matcher = SubstMatch() 
    items = []
    categories = set() 
    for zed, (mol, smi, pdb) in enumerate(zip(parser_1.mols_list, parser_1.smiles_list, parser_1.pdb_files)): 
        if mol is None: 
            continue
        cats = matcher.classify(mol) 
        if not cats:
            cats = ["Fluorocarbon Chain"]
        categories.update(cats) 
        smi_fnames = parser_1.smi_fname[zed] 
        mol_fname = f"{Path(smi_fnames).stem}.png" 
        pdb_path = parser_1.pdb_files[zed]
        pdb_url = f"/static/storage/pdbs/{Path(pdb_path).name}"
        path = OUTPUT_DIR / mol_fname
        Draw.MolToImage(mol, size=SIZE, options=drawOptions).save(path) 
        items.append({
            "filename": f"/static/storage/imgs/{mol_fname}",
            "label": smi,
            "pdb_url": pdb_url,
            "category": " ".join(cats), 
        }) 
    return templates.TemplateResponse(
        "display.html",
        {
            "request": request,
            "items": items, 
            "categories": sorted(categories),
            "count": len(items),
            "show_viewer": False,
        } 
    ) 

@app.post("/append-molecule")
async def add_molecule(
    request: Request,
    data: InputValid
): 
    mol_data = MolData()  
    mol_data.molecules = request.session.get('mols_data', []) # creates 'mols_data' key and grabs values 
    output = mol_data.add_molecule(data.smi_in) 
    output.session['mols_data'] = mol_data.molecules 
    return output 


