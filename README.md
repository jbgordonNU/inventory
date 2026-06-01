# Chemical Inventory Locator

Docker-ready Streamlit app for chemical inventory lookup.

Features:
- Upload CSV inventory
- Saves latest inventory per user name/initials
- Removes the `Location (space)` column from saved/displayed inventory
- Looks up missing SMILES from PubChem
- Lets users draw a structure with Ketcher
- Uses RDKit for true substructure matching
- Shows matching chemicals and location columns

## Files

- `app.py`
- `requirements.txt`
- `Dockerfile`
- `.dockerignore`

## Run locally with Docker

```bash
docker build -t chemical-inventory .
docker run -p 8501:8501 chemical-inventory
```

Open:

```text
http://localhost:8501
```

## Deploy on Render

1. Push these files to GitHub.
2. Go to Render.
3. New > Web Service.
4. Connect the GitHub repo.
5. Environment: Docker.
6. Deploy.

Render should use the included Dockerfile.
