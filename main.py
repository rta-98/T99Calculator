from fastapi import FastAPI, Form, Request, HTTPException
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
from utility.services import InputValid, MolData 

#|%%--%%| <tKnZ0ck33L|NLVW3NX2Sg>
app = FastAPI()
# app.add_middleware() 
app.mount("/static", StaticFiles(directory="static"), name="static")
templates = Jinja2Templates(directory="templates") 

# Serve "Home" as page 0 
@app.get("/", response_class=HTMLResponse)
async def root(request: Request): 
    return templates.TemplateResponse(
        "index.html", {
        "request": request,
        "page": {0: "Home"}
    })

@app.post("/canonical-smiles"):
async def 

@app.post("/smiles-structure")
async def display_image(data: InputValid):
    try:
        image_bytes = MolData.img_2_bytes(data.smi_in) 
        return Response(content=image_bytes, media_type='image/png')
    except HTTPException:
        raise # re-raise error so FastAPI can handle it
    except Exception as e:
        raise HTTPException(status_code=500, detail=f'ERROR: {e}'}

#|%%--%%| <NLVW3NX2Sg|SVBOFjN41o>
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
#|%%--%%| <SVBOFjN41o|UJDrNroodF>



#|%%--%%| <UJDrNroodF|NU92Cne6FE>


#|%%--%%| <NU92Cne6FE|h6D92NwjoB>



