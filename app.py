import re
import time
from pathlib import Path
from datetime import datetime

import pandas as pd
import requests
import streamlit as st
from rdkit import Chem
from streamlit_ketcher import st_ketcher


st.set_page_config(page_title="Chemical Inventory Locator", layout="wide")

DATA_DIR = Path(".user_uploads")
DATA_DIR.mkdir(exist_ok=True)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


st.title("Chemical Inventory Locator")


def safe_user_key(user_text):
    user_text = user_text.strip().lower()
    return re.sub(r"[^a-z0-9_-]+", "_", user_text) or "default_user"


def user_file_path(user_key):
    return DATA_DIR / f"{user_key}_latest_inventory.csv"


def user_meta_path(user_key):
    return DATA_DIR / f"{user_key}_metadata.txt"


def infer_inventory_date(filename, upload_time=None):
    match = re.search(r"(20\d{2})[-_](\d{1,2})[-_](\d{1,2})", filename)
    if match:
        y, m, d = match.groups()
        return f"{int(y):04d}-{int(m):02d}-{int(d):02d}"
    return upload_time or "Unknown"


def save_user_inventory(user_key, uploaded_file):
    path = user_file_path(user_key)
    meta_path = user_meta_path(user_key)

    bytes_data = uploaded_file.getvalue()
    path.write_bytes(bytes_data)

    upload_time = datetime.now().strftime("%Y-%m-%d %H:%M")
    inventory_date = infer_inventory_date(uploaded_file.name, upload_time)

    meta_path.write_text(
        f"original_filename={uploaded_file.name}\n"
        f"uploaded_at={upload_time}\n"
        f"inventory_date={inventory_date}\n"
    )

    return path


def load_metadata(user_key):
    meta_path = user_meta_path(user_key)
    if not meta_path.exists():
        return {}

    meta = {}
    for line in meta_path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            meta[k] = v
    return meta


@st.cache_data(show_spinner=False)
def pubchem_smiles_lookup(identifier):
    if not identifier or str(identifier).lower() == "nan":
        return ""

    identifier = str(identifier).strip()
    url = (
        f"{PUBCHEM_BASE}/name/"
        f"{requests.utils.quote(identifier)}"
        f"/property/CanonicalSMILES/JSON"
    )

    try:
        r = requests.get(url, timeout=12)
        if r.status_code != 200:
            return ""

        data = r.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return ""

        return props[0].get("CanonicalSMILES", "")
    except Exception:
        return ""


def mol_from_smiles(smiles):
    try:
        if not smiles or str(smiles).lower() == "nan":
            return None
        return Chem.MolFromSmiles(str(smiles))
    except Exception:
        return None


def has_substructure(target_smiles, query_mol):
    mol = mol_from_smiles(target_smiles)
    return mol is not None and query_mol is not None and mol.HasSubstructMatch(query_mol)


user_id = st.text_input(
    "User name or initials",
    value="default_user",
    help="Used to remember the last uploaded inventory for this user.",
)

user_key = safe_user_key(user_id)

uploaded = st.file_uploader("Upload chemical inventory CSV", type=["csv"])

if uploaded is not None:
    save_user_inventory(user_key, uploaded)
    st.success("Uploaded file saved as this user's latest inventory.")

latest_path = user_file_path(user_key)
meta = load_metadata(user_key)

if latest_path.exists():
    st.info(
        f"Loaded latest inventory for `{user_id}`. "
        f"Inventory date: **{meta.get('inventory_date', 'Unknown')}**. "
        f"Uploaded at: **{meta.get('uploaded_at', 'Unknown')}**. "
        f"Original file: **{meta.get('original_filename', 'Unknown')}**."
    )

    df = pd.read_csv(latest_path)

    st.subheader("Inventory preview")
    st.dataframe(df.head(), use_container_width=True)

    columns = list(df.columns)

    st.subheader("Column setup")

    name_col = st.selectbox("Chemical name column", columns)

    cas_col = st.selectbox("CAS column, if available", ["None"] + columns)

    smiles_existing = st.selectbox(
        "Existing SMILES column, if available",
        ["None"] + columns,
    )

    if smiles_existing != "None":
        df["SMILES"] = df[smiles_existing].fillna("").astype(str)
    elif "SMILES" not in df.columns:
        df["SMILES"] = ""

    st.subheader("SMILES lookup")

    if st.button("Look up missing SMILES from PubChem"):
        progress = st.progress(0)
        smiles_values = []

        for i, row in df.iterrows():
            current_smiles = str(row.get("SMILES", "")).strip()

            if current_smiles:
                smiles_values.append(current_smiles)
            else:
                found = ""

                if cas_col != "None":
                    found = pubchem_smiles_lookup(row.get(cas_col, ""))

                if not found:
                    found = pubchem_smiles_lookup(row.get(name_col, ""))

                smiles_values.append(found)
                time.sleep(0.15)

            progress.progress((i + 1) / len(df))

        df["SMILES"] = smiles_values
        df.to_csv(latest_path, index=False)
        st.success("SMILES lookup complete and saved to this user's latest inventory.")

    st.download_button(
        "Download inventory with SMILES",
        df.to_csv(index=False),
        file_name="inventory_with_smiles.csv",
        mime="text/csv",
    )

    st.subheader("Draw or enter substructure")

    drawn_smiles = st_ketcher()

    smarts_input = st.text_input(
        "Optional SMARTS query",
        placeholder="Examples: C=O, c1ccccc1, [OH]",
    )

    if smarts_input.strip():
        query_mol = Chem.MolFromSmarts(smarts_input.strip())
        query_text = smarts_input.strip()
    elif drawn_smiles:
        query_mol = Chem.MolFromSmiles(drawn_smiles)
        query_text = drawn_smiles
    else:
        query_mol = None
        query_text = ""

    location_cols = [
        c for c in df.columns
        if "location" in c.lower()
        or c.lower() in ["bench", "shelf", "room", "cabinet", "box"]
    ]

    if query_text:
        st.write("Query:", query_text)

        if query_mol is None:
            st.error("Could not parse the drawn structure or SMARTS query.")
        else:
            matches = df[df["SMILES"].apply(lambda s: has_substructure(s, query_mol))]

            st.subheader("Matching chemicals")
            st.write(f"Found {len(matches)} matching chemicals.")

            display_cols = [name_col, "SMILES"] + location_cols
            display_cols = list(dict.fromkeys([c for c in display_cols if c in df.columns]))

            st.dataframe(matches[display_cols], use_container_width=True)

            st.download_button(
                "Download matching results",
                matches.to_csv(index=False),
                file_name="substructure_matches.csv",
                mime="text/csv",
            )

else:
    st.warning("No saved inventory found for this user. Upload a CSV to begin.")
