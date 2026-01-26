from fastapi import FastAPI, Form, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, field_validator
from typing import Optional, List
import io
import json
import hashlib
import re
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw
import pubchempy as pubpy

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "images"
OUTPUT_DIR.mkdir(exist_ok=True)
ITEMS_PATH = BASE_DIR / "items.json"
SIZE = (200, 200)


class InputValid(BaseModel):
    smi_in: str = Field(..., min_length=1)

    @field_validator("smi_in")
    @classmethod
    def veri_smi(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise TypeError("SMILES are of type str, not None")
        smi_2_mol = Chem.MolFromSmiles(v)
        if smi_2_mol is not None:
            return Chem.MolToSmiles(smi_2_mol)
        raise ValueError(f"Invalid SMILES string: {v}")


class MolData:
    def __init__(self) -> None:
        self.mol_data: List[dict] = []

    def img_2_bytes(self, smi_in: str) -> bytes:
        mol_obj = Chem.MolFromSmiles(smi_in)
        if mol_obj is None:
            raise ValueError(f"Failed to parse SMILES: {smi_in}")
        try:
            mol_img = Draw.MolToImage(mol_obj, size=(100, 100))
            img_bytes = io.BytesIO()
            mol_img.save(img_bytes, format="PNG")
            img_bytes.seek(0)
            return img_bytes.getvalue()
        except Exception as e:
            raise RuntimeError(f"General Error: {e}")

    def get_iupac(self, smi_in: str) -> Optional[str]:
        try:
            smi_iupac = pubpy.get_compounds(smi_in, "smiles")
            return smi_iupac[0].iupac_name
        except Exception:
            return None

    def append_mol_data(self, smi_in: str) -> dict:
        try:
            col_iupac = self.get_iupac(smi_in)
            col_img_bytes = self.img_2_bytes(smi_in)
            dict_mol_data = {
                "smi_in": smi_in,
                "iupac": col_iupac or "Unk",
                "img_bytes": col_img_bytes or b"",
            }
            self.mol_data.append(dict_mol_data)
            return dict_mol_data
        except Exception as e:
            raise RuntimeError(f"Error in append_mol_data: {e}")

    def clear_mol_data(self) -> None:
        self.mol_data = []


class SubstMatch:
    sub_mol_pfeca = Chem.MolFromSmarts("[O]C(C(O)=O)F")
    sub_mol_pfsa = Chem.MolFromSmarts("O=[S](O)=O")
    sub_mol_ftoh = Chem.MolFromSmarts("OC[CH2]")
    sub_mol_mefasaa = Chem.MolFromSmarts("CN(CC(O)=O)[S](=O)=O")
    sub_mol_ftca = Chem.MolFromSmarts("[CH2]CC(O)=O")
    sub_mol_fts = Chem.MolFromSmarts("[CH2]CS(=O)(O)=O")
    sub_mol_pfca = Chem.MolFromSmarts("O=[C]O")
    sub_mol_fasa = Chem.MolFromSmarts("N[S](=O)=O")

    def __init__(self) -> None:
        self.pfeca: List[tuple] = []
        self.fasa: List[tuple] = []
        self.pfca: List[tuple] = []
        self.fts: List[tuple] = []
        self.ftca: List[tuple] = []
        self.mefasaa: List[tuple] = []
        self.ftoh: List[tuple] = []
        self.pfsa: List[tuple] = []

    def match_pfeca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfeca):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfeca.append((matched_smiles, mol_in))
            return True
        return False

    def match_fasa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fasa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.fasa.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfca):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfca.append((matched_smiles, mol_in))
            return True
        return False

    def match_fts(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_fts):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.fts.append((matched_smiles, mol_in))
            return True
        return False

    def match_ftca(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftca):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.ftca.append((matched_smiles, mol_in))
            return True
        return False

    def match_mefasaa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_mefasaa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.mefasaa.append((matched_smiles, mol_in))
            return True
        return False

    def match_ftoh(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_ftoh):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.ftoh.append((matched_smiles, mol_in))
            return True
        return False

    def match_pfsa(self, mol_in) -> bool:
        if mol_in and mol_in.HasSubstructMatch(self.sub_mol_pfsa):
            matched_smiles = Chem.MolToSmiles(mol_in)
            self.pfsa.append((matched_smiles, mol_in))
            return True
        return False

    def classify(self, mol_in) -> List[str]:
        mol_cats: List[str] = []
        if self.match_pfeca(mol_in):
            mol_cats.append("pfeca")
        if self.match_fasa(mol_in):
            mol_cats.append("fasa")
        if self.match_pfca(mol_in):
            mol_cats.append("pfca")
        if self.match_fts(mol_in):
            mol_cats.append("fts")
        if self.match_ftca(mol_in):
            mol_cats.append("ftca")
        if self.match_mefasaa(mol_in):
            mol_cats.append("mefasaa")
        if self.match_ftoh(mol_in):
            mol_cats.append("ftoh")
        if self.match_pfsa(mol_in):
            mol_cats.append("pfsa")
        return mol_cats

    def create_images(self, mol_array: List[tuple]) -> None:
        try:
            for smi, mol in mol_array:
                fname_match = safe_filename(smi)
                path = OUTPUT_DIR / fname_match
                img_match = Draw.MolToImage(mol, size=SIZE)
                img_match.save(path)
        except Exception as e:
            raise RuntimeError(f"Error in create_images(): {e}")


def safe_filename(smiles: str) -> str:
    base = re.sub(r"[^A-Za-z0-9_-]+", "_", smiles).strip("_") or "mol"
    digest = hashlib.sha1(smiles.encode("utf-8")).hexdigest()[:8]
    return f"{base}_{digest}.png"


def load_items() -> List[dict]:
    if ITEMS_PATH.exists():
        return json.loads(ITEMS_PATH.read_text())
    return []


def save_items(items: List[dict]) -> None:
    ITEMS_PATH.write_text(json.dumps(items, indent=2))


def build_item(smiles: str) -> dict:
    canonical = InputValid(smi_in=smiles).smi_in
    mol = Chem.MolFromSmiles(canonical)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {canonical}")

    matcher = SubstMatch()
    categories = matcher.classify(mol)
    if not categories:
        categories = ["unknown"]

    filename = safe_filename(canonical)
    Draw.MolToImage(mol, size=SIZE).save(OUTPUT_DIR / filename)

    iupac = MolData().get_iupac(canonical) or "Unk"

    return {
        "smi_in": canonical,
        "iupac": iupac,
        "filename": filename,
        "category": categories[0],
        "categories": categories,
    }


def build_items_from_smiles(smiles_list: List[str]) -> List[dict]:
    items: List[dict] = []
    for smi in smiles_list:
        items.append(build_item(smi))
    return items


app = FastAPI()
app.mount("/images", StaticFiles(directory=str(OUTPUT_DIR)), name="images")


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    items = load_items()
    rows = []
    for item in items:
        rows.append(
            f'<div style="display:inline-block;margin:6px;text-align:center;">'
            f'<img src="/images/{item["filename"]}" width="120" height="120">'
            f'<div>{item["category"]}</div>'
            f"</div>"
        )
    return "<html><body>" + "".join(rows) + "</body></html>"


@app.get("/items")
def get_items() -> List[dict]:
    return load_items()


@app.post("/smiles")
def add_smiles(smi_in: str = Form(...)) -> dict:
    try:
        item = build_item(smi_in)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))
    items = load_items()
    items.append(item)
    save_items(items)
    return item


if __name__ == "__main__":
    sample_smiles = [
        "O=C(O)c1ccccc1",
        "O=C(O)C(=O)O",
        "CC(=O)O",
        "CN(CCO)S(=O)(=O)O",
        "O=S(=O)(O)O",
    ]
    items = build_items_from_smiles(sample_smiles)
    save_items(items)
    print("Wrote items.json with", len(items), "items")
