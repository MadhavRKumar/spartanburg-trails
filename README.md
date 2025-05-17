# About
A python script to calculate the distribution of trails within each Spartanburg
City Council District. 

It then creates a geojson file of the City Council Districts
with the trail data included. That is used to create [this map](https://umap.openstreetmap.fr/en/map/spartanburg-city-district-trails_1225351#13/34.9422/-81.9294)
which shows the distribution of trails within each district.

# Usage
Use `venv` to create a virtual environment and install the required packages:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows use .venv\Scripts\activate
pip install -r requirements.txt
```

Copy `.env.example` to `.env` and fill in the required variables. The values can
be remote URLs or local paths. 

Then run the script with:
```bash
python main.py
```

