# Protein Classification with Deep Reinforcement Learning

## Project Overview

This project is a Django-based web platform for protein sequence analysis. It accepts one or more FASTA sequences and provides three main workflows:

- Functional annotation with Gene Ontology (GO) terms through the DeepGO API.
- 2D HP lattice structure prediction with Hill Climbing, Simulated Annealing, Monte Carlo, Replica Exchange Monte Carlo, Tabular Q-Learning experiments, and trained DQN inference.
- 3D structure retrieval through AlphaFold DB, with ESMFold used as a fallback for sequences without an AlphaFold entry.

The 2D folding module compares classical heuristic optimization methods with reinforcement learning approaches. The final DQN model uses a constructive formulation: it builds the protein one residue at a time using relative actions (Forward, Left, Right) and a local rotation-invariant 2 x 21 x 21 grid representation.

## Repository Structure

- `protein_project/` - Django web application.
- `protein_project/go_predictor/` - Main app with forms, views, templates, static files, and API integration.
- `tests/` - Algorithm implementations, benchmark scripts, reports, and trained DQN inference module.
- `model/` - Trained DQN weights used by the web interface.
- `colab_export/` - Colab training scripts for the final and previous DQN models.
- `reports/` - Project reports and presentation files.

## Requirements

Use Python 3.12 or newer. The project was generated with Django 6.0.2.

Install dependencies from the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

## How to Run the Web Platform

From the repository root:

```bash
cd protein_project
python manage.py migrate
python manage.py runserver
```

Open:

```text
http://127.0.0.1:8000/
```

## Example FASTA Input

```text
>sp|P04637|P53_HUMAN Cellular tumor antigen p53
MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
```

## Main Features

- Parses multiple FASTA sequences and validates sequence lines.
- Converts amino acid sequences into HP strings using the hydrophobic set `ACFILMVWY`.
- Runs 2D HP folding algorithms with configurable hyperparameters.
- Shows the final lattice conformation, HP string, energy, length, H-count, and step-by-step folding frames.
- Calls DeepGO for GO prediction and groups terms by Biological Process, Molecular Function, and Cellular Component.
- Retrieves 3D PDB structures from AlphaFold DB or predicts them through ESMFold when needed.

## Algorithm Notes

The web interface exposes interactive defaults that may differ from controlled benchmark settings. For example, the Monte Carlo benchmark uses the script-level baseline temperature `T = 2.0`, while the web interface allows user configuration and defaults to `T = 5.0`. Similarly, REMC benchmark tables use 10 replicas for controlled runtime comparison, while the web default is 20 replicas.

## Outputs and Reports

Benchmark outputs and DQN training summaries are stored in `tests/`, including:

- `convergence_summary.txt`
- `remc_convergence_summary.txt`
- `ql_convergence_summary.txt`
- `dqn_report.txt`
- `dqn_summary.csv`
- `dqn_episode_log.csv`

The final thesis document is maintained separately in the top-level `thesis/` directory.
