from pathlib import Path
from datetime import datetime
import re
import time

import pandas as pd
import requests
import streamlit as st
from rdkit import Chem
from streamlit_ketcher import st_ketcher

st.set_page_config(page_title="Chemical Inventory Locator", layout="wide")

DATA_DIR = Path(".inventory")
DATA_DIR.mkdir(exist_ok=True)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


def clean_user_id(user_text):
    user_text = str(user_text).strip().lower()
    return re.sub(r"[^a-z0-9_-]+", "_", user_text) or "default"


def user_file(user_key):
    return DATA_DIR / f"{user_key}.csv"


def meta_file(user_key):
    return DATA_DIR / f"{user_key}.txt"


def remove_unwanted_columns(df):
    return df.drop(columns=["Location (space)"], errors="ignore")


def infer_inventory_date(filename, upload_time):
    match = re.search(r"(20\d{2})[-_](\d{1,2})[-_](\d{1,2})", filename)
    if match:
        y, m, d = match.groups()
        return f"{int(y):04d}-{int(m):02d}-{int(d):02d}"
    return upload_time


@st.cache_data(show_spinner=False)
def lookup_smiles(identifier):
    if identifier is None:
        return ""

    identifier = str(identifier).strip()
    if not identifier or identifier.lower() == "nan":
        return ""

    url = (
        f"{PUBCHEM}/name/"
        f"{requests.utils.quote(identifier)}"
        "/property/CanonicalSMILES/JSON"
    )

    try:
        r = requests.get(url, timeout=12)

        if r.status_code != 200:
            return ""

        props = r.json().get("PropertyTable", {}).get("Properties", [])

        if not props:
            return ""

        return props[0].get("CanonicalSMILES", "") or ""

    except Exception:
        return ""


def mol_from_smiles(smiles):
    try:
        if smiles is None:
            return None

        smiles = str(smiles).strip()

        if not smiles or smiles.lower() == "nan":
            return None

        return Chem.MolFromSmiles(smiles)

    except Exception:
        return None


def has_substructure(target_smiles, query_mol):
    mol = mol_from_smiles(target_smiles)

    if mol is None or query_mol is None:
        return False

    return mol.HasSubstructMatch(query_mol)


def query_from_drawn_or_smarts(drawn_smiles, smarts):
    smarts = str(smarts or "").strip()

    if smarts:
        return Chem.MolFromSmarts(smarts), smarts, "SMARTS"

    if drawn_smiles:
        return Chem.MolFromSmarts(Chem.MolToSmarts(Chem.MolFromSmiles(drawn_smiles))), drawn_smiles, "Drawn structure"

    return None, "", ""


st.title("Chemical Inventory Locator")

st.caption(
    "Upload a chemical inventory, look up SMILES from PubChem, then draw a structure "
    "or enter SMARTS to find matching substructures and locations."
)

user_input = st.text_input("User name or initials", value="default")
user_key = clean_user_id(user_input)

uploaded = st.file_uploader("Upload inventory CSV", type=["csv"])

if uploaded is not None:
    df_uploaded = pd.read_csv(uploaded)
    df_uploaded = remove_unwanted_columns(df_uploaded)

    upload_time = datetime.now().strftime("%Y-%m-%d %H:%M")
    inventory_date = infer_inventory_date(uploaded.name, upload_time)

    df_uploaded.to_csv(user_file(user_key), index=False)

    meta_file(user_key).write_text(
        f"original_filename={uploaded.name}\n"
        f"uploaded_at={upload_time}\n"
        f"inventory_date={inventory_date}\n"
    )

    st.success("Inventory saved for this user.")

path = user_file(user_key)

if not path.exists():
    st.warning("No saved inventory found for this user. Upload a CSV to begin.")
    st.stop()

df = pd.read_csv(path)
df = remove_unwanted_columns(df)
df.to_csv(path, index=False)

metadata = {}

if meta_file(user_key).exists():
    for line in meta_file(user_key).read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            metadata[k] = v

st.info(
    f"Current inventory date: **{metadata.get('inventory_date', 'Unknown')}** | "
    f"Uploaded at: **{metadata.get('uploaded_at', 'Unknown')}** | "
    f"Original file: **{metadata.get('original_filename', 'Unknown')}**"
)

st.subheader("Inventory preview")
st.dataframe(df.head(25), use_container_width=True)

columns = list(df.columns)

st.subheader("Column setup")

name_col = st.selectbox("Chemical name column", columns)

cas_col = st.selectbox("CAS column, if available", ["None"] + columns)

if "SMILES" not in df.columns:
    df["SMILES"] = ""

smiles_col = st.selectbox(
    "SMILES column",
    list(df.columns),
    index=list(df.columns).index("SMILES") if "SMILES" in df.columns else 0,
)

if smiles_col != "SMILES":
    df["SMILES"] = df[smiles_col].fillna("").astype(str)

st.subheader("SMILES lookup")

missing_count = (
    df["SMILES"]
    .fillna("")
    .astype(str)
    .str.strip()
    .eq("")
    .sum()
)

st.write(f"Missing SMILES: {missing_count}")

if st.button("Look up only missing SMILES from PubChem"):
    progress = st.progress(0)
    updated = 0

    for i, row in df.iterrows():
        current = str(row.get("SMILES", "")).strip()

        if current:
            progress.progress((i + 1) / len(df))
            continue

        found = ""

        if cas_col != "None":
            found = lookup_smiles(row.get(cas_col, ""))

        if not found:
            found = lookup_smiles(row.get(name_col, ""))

        if found:
            df.at[i, "SMILES"] = found
            updated += 1

        progress.progress((i + 1) / len(df))
        time.sleep(0.12)

    df = remove_unwanted_columns(df)
    df.to_csv(path, index=False)

    st.success(f"SMILES lookup complete. Added {updated} new SMILES.")

st.download_button(
    "Download inventory with saved SMILES",
    df.to_csv(index=False),
    file_name="inventory_with_smiles.csv",
    mime="text/csv",
)

st.subheader("Substructure search")

drawn_smiles = st_ketcher()

smarts_input = st.text_input(
    "Optional SMARTS query instead of drawn structure",
    placeholder="Examples: C=O, c1ccccc1, [OH], [N+](=O)[O-]",
)

query_mol, query_text, query_type = query_from_drawn_or_smarts(
    drawn_smiles,
    smarts_input,
)

st.write("Drawn SMILES:", drawn_smiles)
st.write("Query text:", query_text)
st.write("Valid query:", query_mol is not None)
st.write("Rows with SMILES:", df["SMILES"].astype(str).str.strip().ne("").sum())

location_cols = [
    c for c in df.columns
    if (
        "location" in c.lower()
        or c.lower()
        in [
            "bench",
            "shelf",
            "room",
            "cabinet",
            "box",
            "specific location note",
        ]
    )
    and c != "Location (space)"
]

if query_text:
    st.write(f"{query_type}: `{query_text}`")

    if query_mol is None:
        st.error("Could not parse the drawn structure or SMARTS query.")
    else:
        matches = df[
            df["SMILES"].apply(lambda s: has_substructure(s, query_mol))
        ].copy()

        st.write(f"Found **{len(matches)}** matching chemicals.")

        display_cols = [name_col, "SMILES"] + location_cols
        display_cols = list(
            dict.fromkeys(
                [c for c in display_cols if c in matches.columns]
            )
        )

        st.dataframe(matches[display_cols], use_container_width=True)

        st.download_button(
            "Download matching results",
            matches.to_csv(index=False),
            file_name="substructure_matches.csv",
            mime="text/csv",
        )

st.subheader("Text search")

text_query = st.text_input("Search all inventory fields")

if text_query:
    mask = pd.Series(False, index=df.index)

    for col in df.columns:
        mask |= df[col].astype(str).str.contains(
            text_query,
            case=False,
            na=False,
        )

    text_matches = df[mask]

    st.write(f"Found **{len(text_matches)}** text matches.")
    st.dataframe(text_matches, use_container_width=True)
