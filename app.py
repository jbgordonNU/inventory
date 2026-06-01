import pandas as pd
import streamlit as st

st.set_page_config(page_title="Chemical Locator", layout="wide")
st.title("Chemical Locator")
st.write("Upload a chemical inventory CSV, then search by chemical name, CAS, functional group keyword, or location fields.")

uploaded = st.file_uploader("Upload inventory CSV", type=["csv"])

FUNCTIONAL_GROUPS = {
    "Aromatic / benzene": ["benzene", "phenyl", "toluene", "xylene", "anisole", "pyridine"],
    "Alcohol": ["ol", "alcohol", "methanol", "ethanol", "propanol", "butanol"],
    "Amine": ["amine", "aniline", "amino"],
    "Carboxylic acid": ["acid", "carboxylic"],
    "Ester": ["ester", "acetate", "benzoate"],
    "Ether": ["ether", "methoxy", "ethoxy"],
    "Halide": ["fluoro", "chloro", "bromo", "iodo", "chloride", "bromide", "iodide"],
    "Ketone / aldehyde": ["one", "aldehyde", "benzaldehyde"],
    "Nitrile": ["nitrile", "cyano"],
    "Nitro": ["nitro"],
}

if uploaded is None:
    st.info("Upload your Chemical Container export CSV to begin.")
    st.stop()

df = pd.read_csv(uploaded)
st.subheader("Inventory preview")
st.dataframe(df.head(50), use_container_width=True)

cols = df.columns.tolist()
name_col = st.selectbox("Chemical name column", cols, index=cols.index("Chemical Name") if "Chemical Name" in cols else 0)
cas_col = st.selectbox("CAS column", ["(none)"] + cols, index=(["(none)"] + cols).index("CAS Number") if "CAS Number" in cols else 0)

mode = st.radio("Search mode", ["Functional group / keyword", "Chemical name or CAS", "Location contains"], horizontal=True)

result = df.copy()

if mode == "Functional group / keyword":
    group = st.selectbox("Pick a group", list(FUNCTIONAL_GROUPS.keys()))
    extra = st.text_input("Optional extra keyword")
    terms = FUNCTIONAL_GROUPS[group] + ([extra] if extra else [])
    pattern = "|".join([t for t in terms if t])
    mask = result[name_col].fillna("").str.contains(pattern, case=False, regex=True)
    result = result[mask]
    st.caption("This MVP uses name-based matching. For true substructure search, add SMILES/InChI and use RDKit.")

elif mode == "Chemical name or CAS":
    q = st.text_input("Search text")
    if q:
        mask = result[name_col].fillna("").str.contains(q, case=False, regex=False)
        if cas_col != "(none)":
            mask = mask | result[cas_col].fillna("").astype(str).str.contains(q, case=False, regex=False)
        result = result[mask]

else:
    location_cols = [c for c in ["Location (space)", "Bench", "Shelf", "Specific Location Note"] if c in cols]
    chosen = st.multiselect("Location columns", location_cols or cols, default=location_cols)
    q = st.text_input("Location text")
    if q and chosen:
        mask = pd.Series(False, index=result.index)
        for c in chosen:
            mask = mask | result[c].fillna("").astype(str).str.contains(q, case=False, regex=False)
        result = result[mask]

st.subheader(f"Matches: {len(result)}")
show_cols = [c for c in ["Chemical Name", "CAS Number", "Location (space)", "Bench", "Shelf", "Specific Location Note", "Amount", "Units", "Container ID", "Chemical Owner", "Storage Group Category"] if c in result.columns]
st.dataframe(result[show_cols] if show_cols else result, use_container_width=True)

csv = result.to_csv(index=False).encode("utf-8")
st.download_button("Download matches as CSV", csv, "chemical_location_matches.csv", "text/csv")
