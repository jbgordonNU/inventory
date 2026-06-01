# Chemical Locator Streamlit App

## Run locally
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy
Push these files to a GitHub repository, then deploy with Streamlit Community Cloud.

## Notes
This MVP uses name/keyword matching because the sample CSV does not include SMILES or InChI structures. For true molecule/substructure search, add a SMILES or InChI column and use RDKit.
