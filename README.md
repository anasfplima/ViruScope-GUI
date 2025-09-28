# ViruScope (Shiny GUI version) - Short Documentation

## Overview
ViruScope is a Python Shiny web app for viral primer retrieval from literature, in-silico primer generation, scoring, annotation and combinations. The UI is organized in tabs:
- Home
- AROLit (fetch articles, scrape primers, validate primers)
- iSOP (in-silico primers generation and BLAST filtering)
- Primer Combinations (conservation scores, loci mapping, primer combinations, structure stability insights)
- Utilities (file tools, FASTA operations, Venn diagrams)

The main file is `app.py` which constructs the UI and server. Supporting modules (in same folder) are imported by the app automatically:
- `Conservation_Scores_V2.py` (iSOP and Primer Combinations logic)
- `arolit.py` (AROLit logic)
- `safedir_generator.py` (custom clean up module due to specific tool-related requirements)
- `utilities_fasta.py` (extra data processing tools)
- `www/` (static assets: CSS, JS, logos, favicon)

## Prerequisites
- Python 3.8+
- BLAST+ command-line tools (if you want BLAST functionality)
- Recommended Python packages (see `requirements.txt`):
  - shiny
  - biopython
  - pandas
  - numpy
  - matplotlib
  - matplotlib-venn
  - others used by helper modules
- Windows Subsystem for Linux is required for automatic article retrieval.

## Installing

```bash
git clone https://github.com/anasfplima/ViruScope-GUI.git
cd ViruScope-GUI
python -m venv .venv
.venv\Scripts\activate     # Windows
# or source .venv/bin/activate   # Mac/Linux
pip install -r requirements.txt
```
 
## Running the app

```bash
shiny run --reload app.py
```

Open http://localhost:8000 in your browser.
