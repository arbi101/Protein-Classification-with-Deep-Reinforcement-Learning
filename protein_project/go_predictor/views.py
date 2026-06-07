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
import math
import sys
import os

ALGORITHM_DEFAULTS = {
    'hc_iterations': 100000,
    'sa_iterations': 100000,
    'sa_initial_t': 30.0,
    'sa_final_t': 0.001,
    'mc_iterations': 100000,
    'mc_temperature': 5.0,
    'remc_iterations': 100000,
    'remc_replicas': 20,
    'remc_swap': 200,
    'remc_tmin': 0.1,
    'remc_tmax': 30.0,
    'dqn_rollouts': 200,
    'trace_interval': 500,
}

BEST_ALGORITHM_PARAMS = {
    'hc': {'Iterations': ALGORITHM_DEFAULTS['hc_iterations']},
    'sa': {
        'Iterations': ALGORITHM_DEFAULTS['sa_iterations'],
        'Initial Temperature': f"{ALGORITHM_DEFAULTS['sa_initial_t']}",
        'Final Temperature': f"{ALGORITHM_DEFAULTS['sa_final_t']}",
    },
    'mc': {
        'Iterations': ALGORITHM_DEFAULTS['mc_iterations'],
        'Temperature': f"{ALGORITHM_DEFAULTS['mc_temperature']}",
    },
    'remc': {
        'Iterations': ALGORITHM_DEFAULTS['remc_iterations'],
        'Replicas': ALGORITHM_DEFAULTS['remc_replicas'],
        'Swap Interval': ALGORITHM_DEFAULTS['remc_swap'],
        'T Min': f"{ALGORITHM_DEFAULTS['remc_tmin']}",
        'T Max': f"{ALGORITHM_DEFAULTS['remc_tmax']}",
    },
    'ql': {
        'Model': 'DQN',
        'Rollouts': ALGORITHM_DEFAULTS['dqn_rollouts'],
        'Construction Step Interval': 1,
    },
}

ALGORITHM_LABELS = {
    'hc': 'Hill Climbing',
    'sa': 'Simulated Annealing',
    'mc': 'Monte Carlo',
    'remc': 'Replica Exchange MC',
    'ql': 'Deep Q-Learning',
}


def build_used_params(algorithm, params):
    return {
        'Model': 'DQN' if algorithm == 'ql' else None,
        'Iterations': (
            params['hc_iterations'] if algorithm == 'hc'
            else params['sa_iterations'] if algorithm == 'sa'
            else params['mc_iterations'] if algorithm == 'mc'
            else params['remc_iterations'] if algorithm == 'remc'
            else None
        ),
        'Initial Temperature': f"{params['sa_initial_t']}" if algorithm == 'sa' else None,
        'Temperature': f"{params['mc_temperature']}" if algorithm == 'mc' else None,
        'Final Temperature': f"{params['sa_final_t']}" if algorithm == 'sa' else None,
        'Replicas': params['remc_replicas'] if algorithm == 'remc' else None,
        'Swap Interval': params['remc_swap'] if algorithm == 'remc' else None,
        'T Min': f"{params['remc_tmin']}" if algorithm == 'remc' else None,
        'T Max': f"{params['remc_tmax']}" if algorithm == 'remc' else None,
        'Rollouts': params['dqn_rollouts'] if algorithm == 'ql' else None,
        'Step Interval': params.get('trace_interval'),
    }


def get_algorithm_budget(algorithm, params, hp_string):
    if algorithm == 'hc':
        return params['hc_iterations']
    if algorithm == 'sa':
        return params['sa_iterations']
    if algorithm == 'mc':
        return params['mc_iterations']
    if algorithm == 'remc':
        return params['remc_iterations']
    return len(hp_string)


def clamp_trace_interval(algorithm, requested_interval, total_steps, max_frames=200):
    min_interval = 1 if algorithm == 'ql' else 500
    requested_interval = max(min_interval, int(requested_interval))
    frame_limited_interval = max(1, math.ceil(max(1, total_steps) / max_frames))
    return max(requested_interval, frame_limited_interval), max_frames


def parse_int_param(post_data, name, default, minimum=None, maximum=None):
    value = int(post_data.get(name, default))
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be at most {maximum}.")
    return value


def parse_float_param(post_data, name, default, minimum=None, maximum=None):
    value = float(post_data.get(name, default))
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise ValueError(f"{name} must be at most {maximum}.")
    return value


def parse_structure_params(post_data):
    return {
        'hc_iterations': parse_int_param(
            post_data, 'hc_iterations', ALGORITHM_DEFAULTS['hc_iterations'],
            minimum=1000,
        ),
        'sa_iterations': parse_int_param(
            post_data, 'sa_iterations', ALGORITHM_DEFAULTS['sa_iterations'],
            minimum=1000,
        ),
        'sa_initial_t': parse_float_param(
            post_data, 'sa_initial_t', ALGORITHM_DEFAULTS['sa_initial_t'],
            minimum=0.1,
        ),
        'sa_final_t': parse_float_param(
            post_data, 'sa_final_t', ALGORITHM_DEFAULTS['sa_final_t'],
            minimum=0.0001,
        ),
        'mc_iterations': parse_int_param(
            post_data, 'mc_iterations', ALGORITHM_DEFAULTS['mc_iterations'],
            minimum=1000,
        ),
        'mc_temperature': parse_float_param(
            post_data, 'mc_temperature', ALGORITHM_DEFAULTS['mc_temperature'],
            minimum=0.1,
        ),
        'remc_iterations': parse_int_param(
            post_data, 'remc_iterations', ALGORITHM_DEFAULTS['remc_iterations'],
            minimum=1000,
        ),
        'remc_replicas': parse_int_param(
            post_data, 'remc_replicas', ALGORITHM_DEFAULTS['remc_replicas'],
            minimum=2, maximum=30,
        ),
        'remc_swap': parse_int_param(
            post_data, 'remc_swap', ALGORITHM_DEFAULTS['remc_swap'],
            minimum=50,
        ),
        'remc_tmin': parse_float_param(
            post_data, 'remc_tmin', ALGORITHM_DEFAULTS['remc_tmin'],
            minimum=0.01,
        ),
        'remc_tmax': parse_float_param(
            post_data, 'remc_tmax', ALGORITHM_DEFAULTS['remc_tmax'],
            minimum=5,
        ),
        'dqn_rollouts': parse_int_param(
            post_data, 'ql_steps', ALGORITHM_DEFAULTS['dqn_rollouts'],
            minimum=1, maximum=300,
        ),
        'trace_interval': parse_int_param(
            post_data, 'trace_interval', ALGORITHM_DEFAULTS['trace_interval'],
            minimum=1,
        ),
    }


def positions_to_structure(positions, hp_string):
    if not positions:
        return [], 0, 0

    xs = [x for x, y in positions]
    ys = [y for x, y in positions]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    width = max_x - min_x + 1
    height = max_y - min_y + 1

    structure = [
        {
            'x': x,
            'y': y,
            'grid_x': x - min_x + 1,
            'grid_y': y - min_y + 1,
            'type': hp_string[i],
        }
        for i, (x, y) in enumerate(positions)
    ]
    return structure, width, height


def build_step_frames(trace_frames, hp_string, algorithm):
    frame_label = 'Residue' if algorithm == 'ql' else 'Iteration'
    step_frames = []
    for frame in trace_frames:
        frame_structure, _, _ = positions_to_structure(frame['positions'], hp_string)
        step_frames.append({
            'step': frame['step'],
            'energy': frame['energy'],
            'structure': frame_structure,
            'label': frame_label,
        })
    return step_frames


# ── Path setup: make the tests/ directory importable ─────────────────────────
# Django's working directory is protein_project/, so we must add tests/ manually
current_dir = os.path.dirname(os.path.abspath(__file__))        # .../go_predictor/
project_dir = os.path.dirname(os.path.dirname(current_dir))     # .../Protein-Classification.../
tests_dir   = os.path.join(project_dir, 'tests')                # .../Protein-Classification.../tests/
model_dir   = os.path.join(project_dir, 'model')                # .../Protein-Classification.../model/
if tests_dir not in sys.path:
    sys.path.append(tests_dir)  # allow 'import test_hc' etc. from anywhere

# ── Import algorithm functions from tests/ ────────────────────────────────────
try:
    from test_hc   import generate_2d_structure          as run_hc    # Hill Climbing
    from test_sa   import generate_2d_structure_sa       as run_sa    # Simulated Annealing
    from test_mc   import generate_2d_structure_mc       as run_mc    # Monte Carlo
    from test_remc import generate_2d_structure_remc     as run_remc  # Replica Exchange MC
    from agent_dqn_model import run_dqn_inference                    # DQN inference

    # Absolute path to the Colab-trained DQN weights file
    dqn_weights_path = os.path.join(model_dir, 'dqn_weights.pth')

    def run_ql(hp_string, **kwargs):
        """
        Wrapper that replaces slow tabular Q-Learning with DQN
        inference for real-time web use.
        """
        max_rollouts = kwargs.get(
            'max_steps_per_episode',
            kwargs.get('steps', ALGORITHM_DEFAULTS['dqn_rollouts']),
        )
        result = run_dqn_inference(
            hp_string,
            dqn_weights_path,
            max_rollouts=max_rollouts,
            return_trace=kwargs.get('return_trace', False),
            trace_interval=kwargs.get('trace_interval', 1),
            max_trace_frames=kwargs.get('max_trace_frames', 200),
        )
        return result

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
    algorithm_params = ALGORITHM_DEFAULTS.copy()

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
                    algorithm_params = parse_structure_params(request.POST)
                except (ValueError, TypeError):
                    # If any value is invalid, fall back to safe defaults
                    error_message = "Invalid parameter values. Using defaults."
                    algorithm_params = ALGORITHM_DEFAULTS.copy()

                hc_iterations   = algorithm_params['hc_iterations']
                sa_iterations   = algorithm_params['sa_iterations']
                sa_initial_t    = algorithm_params['sa_initial_t']
                sa_final_t      = algorithm_params['sa_final_t']
                mc_iterations   = algorithm_params['mc_iterations']
                mc_temperature  = algorithm_params['mc_temperature']
                remc_iterations = algorithm_params['remc_iterations']
                remc_replicas   = algorithm_params['remc_replicas']
                remc_swap       = algorithm_params['remc_swap']
                remc_tmin       = algorithm_params['remc_tmin']
                remc_tmax       = algorithm_params['remc_tmax']
                ql_steps        = algorithm_params['dqn_rollouts']
                trace_interval  = algorithm_params['trace_interval']

                if sequences:
                    # Process at most 5 sequences per request (performance limit)
                    for seq_data in sequences[:5]:
                        sequence  = seq_data['sequence']
                        hp_string = sequence_to_hp(sequence)   # convert AA → HP string
                        total_steps = get_algorithm_budget(algorithm, algorithm_params, hp_string)
                        effective_trace_interval, max_trace_frames = clamp_trace_interval(
                            algorithm,
                            trace_interval,
                            total_steps,
                        )
                        algorithm_params['trace_interval'] = effective_trace_interval

                        # ── Dispatch to the selected algorithm ────────────
                        if algorithm == 'hc':
                            positions, energy, trace_frames = run_hc(
                                hp_string,
                                iterations=hc_iterations,
                                return_trace=True,
                                trace_interval=effective_trace_interval,
                                max_trace_frames=max_trace_frames,
                            )
                        elif algorithm == 'mc':
                            positions, energy, trace_frames = run_mc(
                                hp_string,
                                iterations=mc_iterations,
                                temperature=mc_temperature,
                                return_trace=True,
                                trace_interval=effective_trace_interval,
                                max_trace_frames=max_trace_frames,
                            )
                        elif algorithm == 'remc':
                            positions, energy, trace_frames = run_remc(
                                hp_string,
                                iterations=remc_iterations,
                                num_replicas=remc_replicas,
                                t_min=remc_tmin,
                                t_max=remc_tmax,
                                swap_interval=remc_swap,
                                return_trace=True,
                                trace_interval=effective_trace_interval,
                                max_trace_frames=max_trace_frames,
                            )
                        elif algorithm == 'ql':
                            positions, energy, trace_frames = run_ql(
                                hp_string,
                                max_steps_per_episode=ql_steps,
                                return_trace=True,
                                trace_interval=effective_trace_interval,
                                max_trace_frames=max_trace_frames,
                            )
                        else:  # 'sa' is the default
                            positions, energy, trace_frames = run_sa(
                                hp_string,
                                iterations=sa_iterations,
                                initial_t=sa_initial_t,
                                final_t=sa_final_t,
                                return_trace=True,
                                trace_interval=effective_trace_interval,
                                max_trace_frames=max_trace_frames,
                            )

                        # ── Compute SVG grid bounds ────────────────────────
                        structure, width, height = positions_to_structure(positions, hp_string)
                        step_frames = build_step_frames(trace_frames, hp_string, algorithm)

                        used_params = build_used_params(algorithm, algorithm_params)
                        best_params = BEST_ALGORITHM_PARAMS.get(algorithm, {})

                        structures.append({
                            'name':        seq_data['header'] or f'Sequence {len(structures)+1}',
                            'structure':   structure,    # list of residue dicts for SVG
                            'hp_string':   hp_string,
                            'energy':      energy,       # best HP energy found
                            'grid_width':  width,
                            'grid_height': height,
                            'algorithm':   algorithm,
                            'algorithm_label': ALGORITHM_LABELS.get(algorithm, algorithm.upper()),
                            'used_params': used_params,
                            'best_params': best_params,
                            'step_frames': step_frames,
                            'trace_interval': effective_trace_interval,
                            'total_steps': total_steps,
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
                        error_message = "Could not generate 3D structures. The sequences might be too long or incorrectly formatted."
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
        'algorithm_params':   algorithm_params,
        'year':               year
    })


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
