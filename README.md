# Protein Classification with Deep Reinforcement Learning – Phase 1

## Project Overview

This is the first phase of the project: a Django web app that allows users to submit a protein sequence in FASTA format and get predicted GO terms (Gene Ontology functional annotations). Currently, the GO terms are predicted using the DeepGO API. The app is structured to integrate additional models or APIs as needed.

## How to Run

1. Clone the repository:

```bash
git clone <REPO_URL>
cd Protein-Classification-with-Deep-Reinforcement-Learning
```

2. Create and activate a virtual environment:

Mac / Linux:

```bash
python3 -m venv venv
source venv/bin/activate
```

Windows (PowerShell):

```powershell
venv\\Scripts\\Activate.ps1
```

3. Install dependencies (Django):

```bash
pip install django
```

4. Run the development server:

```bash
python manage.py runserver
```

5. Open your browser and go to:

http://127.0.0.1:8000/

Submit a protein sequence in FASTA format and press "Predict GO Terms".

## Example FASTA Input

```
>TestProtein1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ
```

## Example Output

```
GO:0008150
GO:0003674
```

## Next Steps

- Validate predictions with real protein sequences and benchmark performance.
- Store input sequences and predicted GO terms in the database for analysis (model already implemented).

## Recent Changes

- Translated all Italian comments, docstrings, error messages, and user interface text to English across the entire project.
- Updated `views.py` with English comments and strings.
- Modified `forms.py` labels and help text to English.
- Translated HTML templates (`predict.html` and `structure_2d.html`) to English for consistency.

## Notes

This is Phase 1 of the full pipeline for protein classification using Deep Reinforcement Learning. Phase 2 will involve HP folding in 2D and reinforcement learning to refine functional predictions.