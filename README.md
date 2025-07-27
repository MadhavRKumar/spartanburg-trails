# About
A repository with various scripts to analyze and visualize data about Spartanburg.

# Usage
Use `venv` to create a virtual environment and install the required packages:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows use .venv\Scripts\activate
pip install -r requirements.txt
```

Copy `.env.example` to `.env` and fill in the required variables. The values can
be remote URLs or local paths. 

Then run whatever script you want. For example, to run the main script:
```bash
python main.py
```

