import re
from pathlib import Path
from datetime import datetime

import pandas as pd
import requests
import streamlit as st
from streamlit_ketcher import st_ketcher

st.set_page_config(page_title="Chemical Inventory Locator", layout="wide")

DATA_DIR = Path(".inventory")
DATA_DIR.mkdir(exist_ok=True)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


def user_file(user):
    return DATA_DIR / f"{user}.csv"


def meta_file(user):
    return DATA_DIR / f"{user}.txt"


@st.cache_data
def lookup_smiles(identifier):
    if not identifier:
        return ""

    try:
        url = (
            f"{PUBCHEM}/name/"
            f"{requests.utils.quote(str(identifier))}"
            "/property/CanonicalSMILES/JSON"
        )

        r = requests.get(url, timeout=10)

        if r.status_code != 200:
            return ""

        props = (
            r.json()
            .get("PropertyTable", {})
            .get("Properties", [])
        )

        if not props:
            return ""

        return props[0]["CanonicalSMILES"]

    except Exception:
        return ""


st.title("Chemical Inventory Locator")

user = st.text_input("User", value="default")

upload = st.file_uploader("Upload inventory CSV", type=["csv"])

if upload:
    df = pd.read_csv(upload)
    user_file(user).write_bytes(upload.getvalue())
    meta_file(user).write_text(datetime.now().strftime("%Y-%m-%d %H:%M"))
    st.success("Inventory saved")

if user_file(user).exists():
    df = pd.read_csv(user_file(user))
    uploaded_date = meta_file(user).read_text()

    st.info(f"Current inventory date: {uploaded_date}")
    st.dataframe(df, width="stretch")

    name_col = st.selectbox("Chemical Name Column", df.columns)
    cas_col = st.selectbox("CAS Column", ["None"] + list(df.columns))

    if st.button("Lookup Missing SMILES"):
        smiles = []
        for _, row in df.iterrows():
            value = ""
            if cas_col != "None":
                value = lookup_smiles(row[cas_col])
            if not value:
                value = lookup_smiles(row[name_col])
            smiles.append(value)

        df["SMILES"] = smiles
        df.to_csv(user_file(user), index=False)
        st.success("SMILES Added")

    st.subheader("Draw Molecule")
    drawn = st_ketcher()

    search = st.text_input(
        "Search Chemical Name / CAS / SMILES / Functional Group"
    )

    filtered = df.copy()

    if search:
        mask = pd.Series(False, index=df.index)
        for col in df.columns:
            mask |= (
                df[col].astype(str).str.contains(
                    search, case=False, na=False
                )
            )
        filtered = df[mask]

    st.write(f"{len(filtered)} matching records")
    st.dataframe(filtered, width="stretch")

    st.download_button(
        "Download Current Inventory",
        df.to_csv(index=False),
        file_name="inventory.csv",
    )
else:
    st.warning("Upload an inventory file first.")
