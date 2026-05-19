"""
views.py — Django backend for the protein analysis web platform.

Handles three types of user requests:
  1. 'predict'   → call DeepGO API to get GO term predictions
  2. 'structure' → run one of 6 algorithms (HC/SA/MC/REMC/QL/DQN) to fold the protein
  3. 'alphafold' → fetch 3D structure from AlphaFold DB or ESMFold API
"""

from django.shortcuts import render
from .forms import ProteinSearchForm
import requests
import random
import math
import sys
import os

# ── Path setup: make the tests/ directory importable ─────────────────────────
# Django's working directory is protein_project/, so we must add tests/ manually
current_dir = os.path.dirname(os.path.abspath(__file__))        # .../go_predictor/
project_dir = os.path.dirname(os.path.dirname(current_dir))     # .../ProjetoFinal/
tests_dir   = os.path.join(project_dir, 'tests')                # .../ProjetoFinal/tests/
if tests_dir not in sys.path:
    sys.path.append(tests_dir)  # allow 'import test_hc' etc. from anywhere

# ── Import algorithm functions from tests/ ────────────────────────────────────
try:
    from test_hc   import generate_2d_structure          as run_hc    # Hill Climbing
    from test_sa   import generate_2d_structure_sa       as run_sa    # Simulated Annealing
    from test_mc   import generate_2d_structure_mc       as run_mc    # Monte Carlo
    from test_remc import generate_2d_structure_remc     as run_remc  # Replica Exchange MC
    from test_ql   import generate_2d_structure_ql       as _run_ql_raw  # Tabular Q-Learning (raw)
    from agent_dqn import run_dqn_inference                           # DQN inference

    # Absolute path to the pre-trained DQN weights file
    dqn_weights_path = os.path.join(tests_dir, 'dqn_weights.pth')

    def run_ql(hp_string, **kwargs):
        """
        Wrapper that replaces slow tabular Q-Learning with fast DQN inference
        for real-time web use. Runs only 200 greedy inference steps so the
        web response is near-instant.
        """
        positions, energy = run_dqn_inference(hp_string, dqn_weights_path, max_steps=200)
        return positions, energy

except ImportError as e:
    # If imports fail (e.g., missing PyTorch), print a warning but don't crash
    print(f"Warning: Could not import algorithms from tests folder. {e}")


# ── Utility: extract UniProt ID from a FASTA header ──────────────────────────

def get_uniprot_id_from_fasta(header):
    """
    Parses a standard UniProt FASTA header to extract the accession ID.
    Standard format: >sp|P04637|P53_HUMAN Description...
    The accession (e.g. 'P04637') is the second pipe-delimited field.

    Returns None if the header doesn't match the expected format.
    """
    if not header:
        return None
    parts = header.split('|')
    if len(parts) >= 2:
        return parts[1]   # return the accession code
    return None


# ── 3D structure fetchers ─────────────────────────────────────────────────────

def fetch_alphafold_pdb(uniprot_id):
    """
    Queries the EMBL-EBI AlphaFold API for a given UniProt accession ID
    and returns the URL of the PDB file if found.

    API endpoint: https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}
    Returns: str (URL) or None if not found / request failed.
    """
    if not uniprot_id:
        return None
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()   # raise exception for 4xx/5xx status codes
        data = response.json()
        if data and isinstance(data, list) and len(data) > 0:
            return data[0].get("pdbUrl")  # extract the PDB download URL
    except requests.RequestException as e:
        print(f"AlphaFold API Error for {uniprot_id}: {e}")
    return None


def fetch_esmfold_pdb(sequence):
    """
    Sends a protein sequence to the ESMFold REST API and returns the
    predicted 3D structure in PDB format as a raw string.

    Used as a FALLBACK when AlphaFold DB doesn't have the requested structure.
    ESMFold can predict structures for any sequence without needing a UniProt ID.

    Returns: str (PDB content) or None if request failed.
    """
    if not sequence:
        return None
    # Remove all whitespace/newlines — ESMFold API rejects them (HTTP 422 error)
    clean_seq = "".join(sequence.split())
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=clean_seq)
        response.raise_for_status()
        return response.text   # raw PDB-format string
    except requests.RequestException as e:
        print(f"ESMFold API Error: {e}")
    return None


# ── Main view: handles all user actions ───────────────────────────────────────

def predict_go(request):
    """
    Main Django view function. Handles POST requests from the prediction form.

    Three actions are dispatched based on request.POST['action']:
      'predict'   → DeepGO GO term prediction
      'structure' → 2D HP folding with the selected algorithm
      'alphafold' → 3D structure retrieval (AlphaFold DB or ESMFold fallback)

    Returns an HTTP response by rendering the 'predict.html' template with
    all result data passed as context variables.
    """
    # Initialise all template context variables to safe defaults
    result = None
    results = []
    structure = None
    structures = []
    alphafold_results = []
    hp_string = None
    energy = 0
    error_message = None
    selected_algorithm = 'sa'   # default algorithm shown in the UI

    if request.method == 'POST':
        form = ProteinSearchForm(request.POST)
        if form.is_valid():
            # Extract validated form fields
            fasta_str          = form.cleaned_data['fasta_sequence'].strip()
            action             = request.POST.get('action')              # predict / structure / alphafold
            selected_algorithm = request.POST.get('algorithm', 'sa')     # hc / sa / mc / remc / ql
            sequences          = parse_fasta(fasta_str)                  # list of {header, sequence}

            # ── Branch 1: GO term prediction ─────────────────────────────
            if action == 'predict':
                if sequences:
                    results = []
                    for seq_data in sequences:
                        sequence = seq_data['sequence']
                        result = call_deepgo_api(sequence)    # call DeepGO REST API
                        if result:
                            # Use the FASTA header as the protein name in the UI
                            result['protein']['name'] = seq_data['header'] or f'Sequence {len(results)+1}'
                            results.append(result)
                    if not results:
                        error_message = "Error in GO prediction. Check the FASTA sequences."
                else:
                    error_message = "Enter valid FASTA sequences."

            # ── Branch 2: 2D structure prediction ─────────────────────────
            elif action == 'structure':
                structures = []
                algorithm  = request.POST.get('algorithm', 'sa')

                # Parse all algorithm hyperparameters from the POST form.
                # Use try/except to handle invalid (non-numeric) input gracefully.
                try:
                    hc_iterations   = int(request.POST.get('hc_iterations', 50000))
                    sa_iterations   = int(request.POST.get('sa_iterations', 50000))
                    sa_initial_t    = float(request.POST.get('sa_initial_t', 30.0))
                    sa_final_t      = float(request.POST.get('sa_final_t', 0.001))
                    mc_iterations   = int(request.POST.get('mc_iterations', 50000))
                    mc_temperature  = float(request.POST.get('mc_temperature', 2.0))
                    remc_iterations = int(request.POST.get('remc_iterations', 50000))
                    remc_replicas   = int(request.POST.get('remc_replicas', 5))
                    remc_swap       = int(request.POST.get('remc_swap', 200))
                    remc_tmin       = float(request.POST.get('remc_tmin', 0.1))
                    remc_tmax       = float(request.POST.get('remc_tmax', 30.0))
                    ql_episodes     = int(request.POST.get('ql_episodes', 100))
                    ql_steps        = int(request.POST.get('ql_steps', 200))
                except (ValueError, TypeError):
                    # If any value is invalid, fall back to safe defaults
                    error_message = "Invalid hyperparameter values. Using defaults."
                    hc_iterations = sa_iterations = mc_iterations = remc_iterations = 50000
                    sa_initial_t = 30.0;  sa_final_t = 0.001
                    mc_temperature = 2.0; remc_replicas = 5; remc_swap = 200
                    remc_tmin = 0.1;      remc_tmax = 30.0
                    ql_episodes = 100;    ql_steps = 200

                if sequences:
                    # Process at most 5 sequences per request (performance limit)
                    for seq_data in sequences[:5]:
                        sequence  = seq_data['sequence']
                        hp_string = sequence_to_hp(sequence)   # convert AA → HP string

                        # ── Dispatch to the selected algorithm ────────────
                        if algorithm == 'hc':
                            positions, energy = run_hc(hp_string, iterations=hc_iterations)
                        elif algorithm == 'mc':
                            positions, energy = run_mc(hp_string, iterations=mc_iterations, temperature=mc_temperature)
                        elif algorithm == 'remc':
                            positions, energy = run_remc(hp_string, iterations=remc_iterations,
                                                         num_replicas=remc_replicas,
                                                         t_min=remc_tmin, t_max=remc_tmax,
                                                         swap_interval=remc_swap)
                        elif algorithm == 'ql':
                            positions, energy = run_ql(hp_string, episodes=ql_episodes,
                                                       max_steps_per_episode=ql_steps)
                        else:  # 'sa' is the default
                            positions, energy = run_sa(hp_string, iterations=sa_iterations,
                                                       initial_t=sa_initial_t, final_t=sa_final_t)

                        # ── Compute SVG grid bounds ────────────────────────
                        # Find the bounding box of all residue positions so
                        # the SVG renderer can auto-scale the lattice to fit
                        xs = [x for x, y in positions]
                        ys = [y for x, y in positions]
                        min_x, max_x = min(xs), max(xs)
                        min_y, max_y = min(ys), max(ys)
                        width  = max_x - min_x + 1   # lattice width in cells
                        height = max_y - min_y + 1   # lattice height in cells

                        # Build a list of residue dicts for the template.
                        # grid_x / grid_y are 1-indexed positions inside the bounding box.
                        structure = [
                            {
                                'x': x, 'y': y,
                                'grid_x': x - min_x + 1,   # normalise to 1..width
                                'grid_y': y - min_y + 1,   # normalise to 1..height
                                'type': hp_string[i]        # 'H' or 'P'
                            }
                            for i, (x, y) in enumerate(positions)
                        ]

                        # Build a dict of the hyperparameters actually used
                        used_params = {
                            'Iterations':        hc_iterations if algorithm == 'hc' else sa_iterations if algorithm == 'sa' else mc_iterations if algorithm == 'mc' else remc_iterations if algorithm == 'remc' else f'{ql_episodes} episodes',
                            'Temperature':       f'{sa_initial_t}' if algorithm == 'sa' else f'{mc_temperature}' if algorithm == 'mc' else None,
                            'Final Temperature': f'{sa_final_t}'   if algorithm == 'sa' else None,
                            'Replicas':          remc_replicas      if algorithm == 'remc' else None,
                            'Swap Interval':     remc_swap          if algorithm == 'remc' else None,
                            'T Min':             f'{remc_tmin}'     if algorithm == 'remc' else None,
                            'T Max':             f'{remc_tmax}'     if algorithm == 'remc' else None,
                            'Steps per Episode': ql_steps           if algorithm == 'ql'   else None,
                        }

                        # Build a dict of the recommended "best" parameters per algorithm
                        best_params = {}
                        if   algorithm == 'hc':   best_params = {'Iterations': 50000}
                        elif algorithm == 'sa':   best_params = {'Iterations': 50000, 'Temperature': '30.0', 'Final Temperature': '0.001'}
                        elif algorithm == 'mc':   best_params = {'Iterations': 50000, 'Temperature': '2.0'}
                        elif algorithm == 'remc': best_params = {'Iterations': 50000, 'Replicas': 5, 'Swap Interval': 200, 'T Min': '0.1', 'T Max': '30.0'}
                        elif algorithm == 'ql':   best_params = {'Episodes': 100, 'Steps per Episode': 200}

                        structures.append({
                            'name':        seq_data['header'] or f'Sequence {len(structures)+1}',
                            'structure':   structure,    # list of residue dicts for SVG
                            'hp_string':   hp_string,
                            'energy':      energy,       # best HP energy found
                            'grid_width':  width,
                            'grid_height': height,
                            'algorithm':   algorithm,
                            'used_params': used_params,
                            'best_params': best_params
                        })
                else:
                    error_message = "Enter valid FASTA sequences."

            # ── Branch 3: 3D structure retrieval ──────────────────────────
            elif action == 'alphafold':
                alphafold_results = []
                if sequences:
                    for seq_data in sequences[:5]:   # max 5 sequences at once
                        header    = seq_data['header']
                        sequence  = seq_data['sequence']
                        uniprot_id = get_uniprot_id_from_fasta(header)

                        pdb_url  = None
                        pdb_data = None
                        source   = None

                        # Phase 1: try AlphaFold DB (requires a known UniProt ID)
                        if uniprot_id:
                            pdb_url = fetch_alphafold_pdb(uniprot_id)
                            if pdb_url:
                                source = "AlphaFold DB"

                        # Phase 2 fallback: use ESMFold to predict structure de-novo
                        # if AlphaFold doesn't have this protein
                        if not pdb_url:
                            pdb_data = fetch_esmfold_pdb(sequence)
                            if pdb_data:
                                source = "ESMFold API"

                        alphafold_results.append({
                            'name':       header or f'Sequence {len(alphafold_results)+1}',
                            'uniprot_id': uniprot_id,
                            'pdb_url':    pdb_url,    # URL to PDB file (AlphaFold)
                            'pdb_data':   pdb_data,   # raw PDB string (ESMFold)
                            'source':     source      # which API provided the result
                        })

                    # If BOTH AlphaFold and ESMFold failed for all sequences
                    if not any(res['pdb_url'] or res['pdb_data'] for res in alphafold_results):
                        error_message = "Could not generate 3D structures. The sequences might be too long or malformatted."
                else:
                    error_message = "Please enter valid FASTA sequences."

    else:
        # GET request: display the empty form
        form = ProteinSearchForm()

    # Inject the current year for the footer copyright notice
    from datetime import datetime
    year = datetime.utcnow().year

    # Render the template with all context variables
    return render(request, 'go_predictor/predict.html', {
        'form':               form,
        'results':            results,            # GO prediction results
        'structures':         structures,          # 2D folding results
        'alphafold_results':  alphafold_results,   # 3D structure results
        'error_message':      error_message,
        'selected_algorithm': selected_algorithm,
        'year':               year
    })


# ── Legacy demo function (not used in production) ─────────────────────────────

def fake_search(query):
    """
    Returns hardcoded fake GO data for UI testing purposes.
    Not used in the live application — kept for development reference.
    """
    return {
        'protein': {'name': 'Tumor protein p53', 'description': 'Tumor suppressor'},
        'categories': [
            {'title': 'Biological Processes', 'subtitle': 'Biological Process',
             'terms': [{'go': 'GO:0006915', 'label': 'apoptotic process'},
                       {'go': 'GO:0006974', 'label': 'cellular response to DNA damage stimulus'},
                       {'go': 'GO:0007050', 'label': 'cell cycle arrest'}]},
            {'title': 'Molecular Functions', 'subtitle': 'Molecular Function',
             'terms': [{'go': 'GO:0003677', 'label': 'DNA binding'},
                       {'go': 'GO:0003700', 'label': 'DNA-binding transcription factor activity'}]},
            {'title': 'Cellular Components', 'subtitle': 'Cellular Component',
             'terms': [{'go': 'GO:0005634', 'label': 'nucleus'},
                       {'go': 'GO:0005737', 'label': 'cytoplasm'}]}
        ]
    }


# ── FASTA parser ──────────────────────────────────────────────────────────────

def parse_fasta(fasta_str):
    """
    Parses a multi-sequence FASTA string into a list of dicts.

    Handles:
      - Multiple sequences separated by '>' headers
      - Mixed-case sequences (converted to uppercase)
      - Empty lines (skipped)
      - Sequences without headers (header stored as empty string)

    Parameters
    ----------
    fasta_str : str  Raw FASTA text from the textarea

    Returns
    -------
    list of {'header': str, 'sequence': str}
    """
    sequences = []
    current_header = ''
    current_seq = ''

    for line in fasta_str.strip().splitlines():
        line = line.strip()
        if not line:
            continue   # skip blank lines
        if line.startswith('>'):
            if current_seq:
                # Save the previous sequence before starting a new one
                sequences.append({'header': current_header, 'sequence': current_seq})
            current_header = line[1:]   # strip the '>' character from the header
            current_seq = ''
        else:
            current_seq += line.upper()  # accumulate residues in uppercase

    if current_seq:
        # Save the last sequence (no trailing '>' to trigger it)
        sequences.append({'header': current_header, 'sequence': current_seq})

    return sequences


# ── HP translation ────────────────────────────────────────────────────────────

def sequence_to_hp(sequence):
    """
    Converts an amino acid sequence to a binary HP string.

    Classification:
      'H' (Hydrophobic) : A, C, F, I, L, M, V, W, Y
      'P' (Polar)       : all other amino acids

    Parameters
    ----------
    sequence : str  One-letter amino acid sequence (uppercase)
    Returns  : str  HP binary string of the same length
    """
    hydrophobic = set('ACFILMVWY')  # standard hydrophobic residue set
    hp_string = ''
    for aa in sequence:
        if aa in hydrophobic:
            hp_string += 'H'
        else:
            hp_string += 'P'
    return hp_string


# ── Secondary view: structure_2d (standalone page) ────────────────────────────

def structure_2d(request):
    """
    Alternative view for a dedicated 2D structure prediction page.
    Functionally identical to the 'structure' branch of predict_go(),
    but renders a different template ('structure_2d.html').
    Useful if the UI separates GO prediction and structure prediction
    into different pages.
    """
    structures = []
    alphafold_results = []
    error_message = None
    selected_algorithm = 'sa'

    if request.method == 'POST':
        form = ProteinSearchForm(request.POST)
        if form.is_valid():
            fasta_str          = form.cleaned_data['fasta_sequence'].strip()
            action             = request.POST.get('action')
            selected_algorithm = request.POST.get('algorithm', 'sa')
            sequences          = parse_fasta(fasta_str)

            if action == 'structure' or not action:
                algorithm = request.POST.get('algorithm', 'sa')
                if sequences:
                    for seq_data in sequences[:5]:
                        sequence  = seq_data['sequence']
                        hp_string = sequence_to_hp(sequence)

                        # Run the selected algorithm with default parameters
                        if   algorithm == 'hc':   positions, energy = run_hc(hp_string, iterations=50000)
                        elif algorithm == 'mc':   positions, energy = run_mc(hp_string, iterations=50000, temperature=2.0)
                        elif algorithm == 'remc': positions, energy = run_remc(hp_string, iterations=50000, num_replicas=5)
                        elif algorithm == 'ql':   positions, energy = run_ql(hp_string, episodes=500, max_steps_per_episode=200)
                        else:                     positions, energy = run_sa(hp_string, iterations=50000)

                        # Compute bounding box for SVG scaling
                        xs = [x for x, y in positions]; ys = [y for x, y in positions]
                        min_x, max_x = min(xs), max(xs); min_y, max_y = min(ys), max(ys)
                        width  = max_x - min_x + 1
                        height = max_y - min_y + 1

                        # Build residue list for the template
                        structure = [
                            {'x': x, 'y': y,
                             'grid_x': x - min_x + 1, 'grid_y': y - min_y + 1,
                             'type': hp_string[i]}
                            for i, (x, y) in enumerate(positions)
                        ]
                        structures.append({
                            'name': seq_data['header'] or f'Sequence {len(structures)+1}',
                            'structure': structure, 'hp_string': hp_string,
                            'energy': energy, 'grid_width': width, 'grid_height': height
                        })
                else:
                    error_message = "Enter valid FASTA sequences."

            elif action == 'alphafold':
                if sequences:
                    for seq_data in sequences[:5]:
                        header    = seq_data['header']; sequence = seq_data['sequence']
                        uniprot_id = get_uniprot_id_from_fasta(header)
                        pdb_url = pdb_data = source = None

                        # Try AlphaFold DB first, fall back to ESMFold
                        if uniprot_id:
                            pdb_url = fetch_alphafold_pdb(uniprot_id)
                            if pdb_url: source = "AlphaFold DB"
                        if not pdb_url:
                            pdb_data = fetch_esmfold_pdb(sequence)
                            if pdb_data: source = "ESMFold API"

                        alphafold_results.append({
                            'name': header or f'Sequence {len(alphafold_results)+1}',
                            'uniprot_id': uniprot_id, 'pdb_url': pdb_url,
                            'pdb_data': pdb_data, 'source': source
                        })
                    if not any(r['pdb_url'] or r['pdb_data'] for r in alphafold_results):
                        error_message = "Could not generate 3D structures."
                else:
                    error_message = "Please enter valid FASTA sequences."
    else:
        form = ProteinSearchForm()

    from datetime import datetime
    year = datetime.utcnow().year
    return render(request, 'go_predictor/structure_2d.html', {
        'form': form, 'structures': structures,
        'alphafold_results': alphafold_results,
        'error_message': error_message,
        'selected_algorithm': selected_algorithm, 'year': year
    })


# ── DeepGO API caller ─────────────────────────────────────────────────────────

def call_deepgo_api(sequence):
    """
    Sends an amino acid sequence to the DeepGO REST API and returns
    predicted Gene Ontology (GO) terms with confidence scores.

    DeepGO predicts terms across three GO domains:
      - Biological Process (BP)
      - Molecular Function (MF)
      - Cellular Component (CC)

    Parameters
    ----------
    sequence : str  Raw amino acid sequence (one-letter codes)

    Returns
    -------
    dict (parsed result) or None if the API call failed
    """
    url = 'https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create'
    # Payload format required by the DeepGO v1.0.3 API
    json_data = {
        'data_format': 'enter',    # 'enter' means raw sequence (not file upload)
        'data': sequence,
        'threshold': 0.3,          # minimum confidence score to include a GO term
        'version': '1.0.3'
    }
    try:
        response = requests.post(url, json=json_data)
        response.raise_for_status()                 # raise on HTTP errors
        api_result = response.json()
        return parse_deepgo_result(api_result)      # parse and restructure the response
    except requests.RequestException as e:
        print(f"API Error: {e}")
        return None


def parse_deepgo_result(api_result):
    """
    Transforms the raw DeepGO API JSON response into a structured dict
    suitable for rendering in the Django template.

    Input structure (simplified):
      api_result['predictions'][0]['functions'] = [
          {'name': 'Biological Process', 'functions': [(go_id, label, score), ...]},
          {'name': 'Molecular Function', ...},
          {'name': 'Cellular Component', ...},
      ]

    Output structure:
      {
        'protein': {'name': str, 'description': str},
        'categories': [{'title': str, 'subtitle': str, 'terms': [...]}]
      }

    Returns None if the response is empty or malformed.
    """
    if not api_result or 'predictions' not in api_result or not api_result['predictions']:
        return None   # malformed or empty response

    pred = api_result['predictions'][0]   # first prediction in the list

    # Extract protein name from protein_info field
    protein_info = pred.get('protein_info')
    if protein_info:
        name = protein_info.split()[0] if '|' in protein_info else protein_info
    else:
        name = 'Unknown Protein'

    # Map DeepGO category names to human-readable display titles
    categories = []
    for func_cat in pred.get('functions', []):
        cat_name = func_cat['name']
        if   cat_name == 'Biological Process': title, subtitle = 'Biological Processes', 'Biological Process'
        elif cat_name == 'Molecular Function':  title, subtitle = 'Molecular Functions',  'Molecular Function'
        elif cat_name == 'Cellular Component':  title, subtitle = 'Cellular Components',  'Cellular Component'
        else:
            continue   # skip unknown category names

        # Build the list of GO terms: each is (go_id, label, score)
        terms = [{'go': go, 'label': label, 'score': score}
                 for go, label, score in func_cat['functions']]
        categories.append({'title': title, 'subtitle': subtitle, 'terms': terms})

    return {
        'protein':    {'name': name, 'description': 'Predicted by DeepGO'},
        'categories': categories
    }