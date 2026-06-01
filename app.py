import time
import requests
import pandas as pd
import streamlit as st
from rdkit import Chem
from streamlit_ketcher import st_ketcher

st.set_page_config(page_title="Chemical Substructure Locator", layout="wide")

st.title("Chemical Substructure Locator")

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


@st.cache_data(show_spinner=False)
def pubchem_smiles_lookup(identifier: str):
    """Look up canonical SMILES from PubChem by chemical name or CAS."""
    if not identifier or str(identifier).lower() == "nan":
        return None

    identifier = str(identifier).strip()
    url = f"{PUBCHEM_BASE}/name/{requests.utils.quote(identifier)}/property/CanonicalSMILES/JSON"

    try:
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            return None

        data = r.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None

        return props[0].get("CanonicalSMILES")
    except Exception:
        return None


def mol_from_smiles(smiles):
    try:
        return Chem.MolFromSmiles(str(smiles))
    except Exception:
        return None


def has_substructure(target_smiles, query_mol):
    mol = mol_from_smiles(target_smiles)
    if mol is None or query_mol is None:
        return False
    return mol.HasSubstructMatch(query_mol)


uploaded = st.file_uploader("Upload chemical inventory CSV", type=["csv"])

if uploaded:
    df = pd.read_csv(uploaded)
    st.write("Preview:")
    st.dataframe(df.head())

    st.subheader("Column setup")

    columns = list(df.columns)

    name_col = st.selectbox(
        "Chemical name column",
        columns,
        index=0,
    )

    cas_col = st.selectbox(
        "CAS column, if available",
        ["None"] + columns,
    )

    existing_smiles_col = st.selectbox(
        "Existing SMILES column, if available",
        ["None"] + columns,
    )

    if existing_smiles_col != "None":
        df["SMILES"] = df[existing_smiles_col]
    else:
        df["SMILES"] = ""

    st.subheader("Look up missing SMILES")

    if st.button("Look up SMILES from PubChem"):
        progress = st.progress(0)
        results = []

        for i, row in df.iterrows():
            current = row.get("SMILES", "")

            if pd.notna(current) and str(current).strip():
                results.append(current)
            else:
                smiles = None

                if cas_col != "None":
                    smiles = pubchem_smiles_lookup(row.get(cas_col))

                if not smiles:
                    smiles = pubchem_smiles_lookup(row.get(name_col))

                results.append(smiles or "")

                # Be polite to PubChem
                time.sleep(0.15)

            progress.progress((i + 1) / len(df))

        df["SMILES"] = results
        st.success("SMILES lookup complete.")

    st.download_button(
        "Download CSV with SMILES",
        df.to_csv(index=False),
        file_name="chemical_inventory_with_smiles.csv",
        mime="text/csv",
    )

    st.subheader("Draw substructure")

    drawn_smiles = st_ketcher()
    smarts_input = st.text_input(
        "Optional SMARTS query instead of drawn structure",
        placeholder="Example: C=O, c1ccccc1, [OH]",
    )

    if smarts_input.strip():
        query_mol = Chem.MolFromSmarts(smarts_input.strip())
        query_label = smarts_input.strip()
    elif drawn_smiles:
        query_mol = Chem.MolFromSmiles(drawn_smiles)
        query_label = drawn_smiles
    else:
        query_mol = None
        query_label = None

    if query_label:
        st.write("Query:", query_label)

    location_cols = [
        c for c in df.columns
        if "location" in c.lower()
        or c.lower() in ["bench", "shelf", "room", "cabinet"]
    ]

    if query_mol is not None and "SMILES" in df.columns:
        matches = df[df["SMILES"].apply(lambda s: has_substructure(s, query_mol))]

        st.subheader("Matches")
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

    elif query_label:
        st.error("Could not parse the drawn structure or SMARTS query.")
else:
    st.info("Upload your chemical inventory CSV to begin.")
