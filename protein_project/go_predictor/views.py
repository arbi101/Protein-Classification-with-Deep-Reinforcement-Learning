from django.shortcuts import render
from .forms import ProteinSearchForm
import requests
import random
import math

def get_uniprot_id_from_fasta(header):
    """
    Attempts to extract a UniProt ID from a standard FASTA header.
    Format is typically >sp|P04637|P53_HUMAN
    """
    if not header:
        return None
    parts = header.split('|')
    if len(parts) >= 2:
        return parts[1]
    return None

def fetch_alphafold_pdb(uniprot_id):
    """
    Queries the EMBL-EBI AlphaFold API for the given UniProt ID
    and returns the PDB file URL if available.
    """
    if not uniprot_id:
        return None
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        if data and isinstance(data, list) and len(data) > 0:
            return data[0].get("pdbUrl")
    except requests.RequestException as e:
        print(f"AlphaFold API Error for {uniprot_id}: {e}")
    return None

def fetch_esmfold_pdb(sequence):
    """
    Predicts the 3D structure solely from the amino acid sequence
    using the ESMFold API. Returns raw PDB string format.
    """
    if not sequence:
        return None
    # Clean the sequence of any whitespace/newlines that causes 422 errors
    clean_seq = "".join(sequence.split())
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=clean_seq)
        response.raise_for_status() # Raises an exception on bad status
        return response.text # Raw PDB string output
    except requests.RequestException as e:
        print(f"ESMFold API Error: {e}")
    return None

def predict_go(request):
    result = None
    results = []
    structure = None
    structures = []
    alphafold_results = []
    hp_string = None
    energy = 0
    error_message = None
    if request.method == 'POST':
        form = ProteinSearchForm(request.POST)
        if form.is_valid():
            fasta_str = form.cleaned_data['fasta_sequence'].strip()
            action = request.POST.get('action')
            sequences = parse_fasta(fasta_str)

            if action == 'predict':
                if sequences:
                    results = []
                    for seq_data in sequences:
                        sequence = seq_data['sequence']
                        result = call_deepgo_api(sequence)
                        if result:
                            result['protein']['name'] = seq_data['header'] or f'Sequence {len(results)+1}'
                            results.append(result)
                    if not results:
                        error_message = "Error in GO prediction. Check the FASTA sequences."
                else:
                    error_message = "Enter valid FASTA sequences."

            elif action == 'structure':
                structures = []
                if sequences:
                    for seq_data in sequences[:5]:  # Limit to 5 for performance
                        sequence = seq_data['sequence']
                        hp_string = sequence_to_hp(sequence)
                        positions, energy = generate_2d_structure_sa(hp_string)
                        structure = [{'x': x, 'y': y, 'type': hp_string[i]} for i, (x, y) in enumerate(positions)]
                        structures.append({
                            'name': seq_data['header'] or f'Sequence {len(structures)+1}',
                            'structure': structure,
                            'hp_string': hp_string,
                            'energy': energy
                        })
                else:
                    error_message = "Enter valid FASTA sequences."
                    
            elif action == 'alphafold':
                alphafold_results = []
                if sequences:
                    for seq_data in sequences[:5]: # Limit to 5 max
                        header = seq_data['header']
                        sequence = seq_data['sequence']
                        uniprot_id = get_uniprot_id_from_fasta(header)
                        
                        pdb_url = None
                        pdb_data = None
                        source = None
                        
                        # Phase 1: Try AlphaFold BD via UniProt
                        if uniprot_id:
                            pdb_url = fetch_alphafold_pdb(uniprot_id)
                            if pdb_url:
                                source = "AlphaFold DB"

                        # Phase 2 (Fallback): Try ESMFold via raw sequence if AlphaFold failed
                        if not pdb_url:
                            pdb_data = fetch_esmfold_pdb(sequence)
                            if pdb_data:
                                source = "ESMFold API"
                        
                        alphafold_results.append({
                            'name': header or f'Sequence {len(alphafold_results)+1}',
                            'uniprot_id': uniprot_id,
                            'pdb_url': pdb_url,
                            'pdb_data': pdb_data,
                            'source': source
                        })
                    
                    # If all methods failed to fetch
                    if not any(res['pdb_url'] or res['pdb_data'] for res in alphafold_results):
                        error_message = "Could not generate 3D structures. The sequences might be too long or malformatted."
                else:
                    error_message = "Please enter valid FASTA sequences."

    else:
        form = ProteinSearchForm()

    # add current year for footer
    from datetime import datetime
    year = datetime.utcnow().year
    return render(request, 'go_predictor/predict.html', {
        'form': form, 
        'results': results, 
        'structures': structures, 
        'alphafold_results': alphafold_results,
        'error_message': error_message, 
        'year': year
    })


def fake_search(query):
    """Returns fake data structure based on the entered name."""
    # for demonstration return hardcoded data regardless of `query`
    return {
        'protein': {
            'name': 'Tumor protein p53',
            'description': 'Tumor suppressor protein that regulates the cell cycle'
        },
        'categories': [
            {
                'title': 'Biological Processes',
                'subtitle': 'Biological Process',
                'terms': [
                    {'go': 'GO:0006915', 'label': 'apoptotic process'},
                    {'go': 'GO:0006974', 'label': 'cellular response to DNA damage stimulus'},
                    {'go': 'GO:0007050', 'label': 'cell cycle arrest'},
                ]
            },
            {
                'title': 'Molecular Functions',
                'subtitle': 'Molecular Function',
                'terms': [
                    {'go': 'GO:0003677', 'label': 'DNA binding'},
                    {'go': 'GO:0003700', 'label': 'DNA-binding transcription factor activity'},
                ]
            },
            {
                'title': 'Cellular Components',
                'subtitle': 'Cellular Component',
                'terms': [
                    {'go': 'GO:0005634', 'label': 'nucleus'},
                    {'go': 'GO:0005737', 'label': 'cytoplasm'},
                ]
            }
        ]
    }


def parse_fasta(fasta_str):
    """
    Parses FASTA format and returns a list of sequences.
    Each sequence is a dict with 'header' and 'sequence'.
    """
    sequences = []
    current_header = ''
    current_seq = ''
    lines = fasta_str.strip().splitlines()
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_seq:
                sequences.append({'header': current_header, 'sequence': current_seq})
            current_header = line[1:]  # remove >
            current_seq = ''
        else:
            current_seq += line.upper()
    if current_seq:
        sequences.append({'header': current_header, 'sequence': current_seq})
    return sequences


def sequence_to_hp(sequence):
    """
    Converts an amino acid sequence to HP string
    H: Hydrophobic, P: Polar
    """
    hydrophobic = set('ACFILMVWY')  # Common hydrophobic amino acids
    hp_string = ''
    for aa in sequence:
        if aa in hydrophobic:
            hp_string += 'H'
        else:
            hp_string += 'P'
    return hp_string


def generate_2d_structure_sa(hp_string, iterations=50000, initial_t=10.0, final_t=0.01):
    """
    Generates a 2D structure using Simulated Annealing with Pivot Moves
    Returns list of positions (x,y) and energy
    """
    if not hp_string:
        return [], 0

    def calculate_energy(pos_list):
        e = 0
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':
                x, y = pos_list[i]
                for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        return -(e // 2)

    current_positions = [(i, 0) for i in range(len(hp_string))]
    current_energy = calculate_energy(current_positions)
    
    best_positions = list(current_positions)
    best_energy = current_energy

    for i in range(iterations):
        if len(hp_string) <= 2:
            break
            
        if iterations > 1:
            fraction = i / float(iterations - 1)
            T = initial_t * ((final_t / initial_t) ** fraction)
        else:
            T = final_t
            
        pivot_idx = random.randint(1, len(hp_string) - 2)
        angle = random.choice([90, -90, 180])
        
        cx, cy = current_positions[pivot_idx]
        new_positions = list(current_positions)
        
        if angle == 90:
            cos_a, sin_a = 0, 1
        elif angle == -90:
            cos_a, sin_a = 0, -1
        else:
            cos_a, sin_a = -1, 0
            
        for k in range(pivot_idx + 1, len(hp_string)):
            x, y = current_positions[k]
            tx, ty = x - cx, y - cy
            rx = tx * cos_a - ty * sin_a
            ry = tx * sin_a + ty * cos_a
            new_positions[k] = (rx + cx, ry + cy)
            
        if len(set(new_positions)) == len(new_positions):
            new_energy = calculate_energy(new_positions)
            
            if new_energy <= current_energy:
                current_positions = new_positions
                current_energy = new_energy
                
                if current_energy < best_energy:
                    best_positions = list(current_positions)
                    best_energy = current_energy
            else:
                delta_e = new_energy - current_energy
                probability = math.exp(-delta_e / T) 
                
                if random.random() < probability:
                    current_positions = new_positions
                    current_energy = new_energy

    return best_positions, best_energy


def structure_2d(request):
    structures = []
    alphafold_results = []
    error_message = None
    if request.method == 'POST':
        form = ProteinSearchForm(request.POST)
        if form.is_valid():
            fasta_str = form.cleaned_data['fasta_sequence'].strip()
            action = request.POST.get('action')
            sequences = parse_fasta(fasta_str)
            
            if action == 'structure' or not action:
                if sequences:
                    for seq_data in sequences[:5]:  # Limit to 5
                        sequence = seq_data['sequence']
                        hp_string = sequence_to_hp(sequence)
                        positions, energy = generate_2d_structure_sa(hp_string)
                        # Convert positions to list for template
                        structure = [{'x': x, 'y': y, 'type': hp_string[i]} for i, (x, y) in enumerate(positions)]
                        structures.append({
                            'name': seq_data['header'] or f'Sequence {len(structures)+1}',
                            'structure': structure,
                            'hp_string': hp_string,
                            'energy': energy
                        })
                else:
                    error_message = "Enter valid FASTA sequences."
                    
            elif action == 'alphafold':
                if sequences:
                    for seq_data in sequences[:5]: # Limit to 5 max
                        header = seq_data['header']
                        sequence = seq_data['sequence']
                        uniprot_id = get_uniprot_id_from_fasta(header)
                        
                        pdb_url = None
                        pdb_data = None
                        source = None
                        
                        # Phase 1: Try AlphaFold BD via UniProt
                        if uniprot_id:
                            pdb_url = fetch_alphafold_pdb(uniprot_id)
                            if pdb_url:
                                source = "AlphaFold DB"

                        # Phase 2 (Fallback): Try ESMFold via raw sequence if AlphaFold failed
                        if not pdb_url:
                            pdb_data = fetch_esmfold_pdb(sequence)
                            if pdb_data:
                                source = "ESMFold API"
                        
                        alphafold_results.append({
                            'name': header or f'Sequence {len(alphafold_results)+1}',
                            'uniprot_id': uniprot_id,
                            'pdb_url': pdb_url,
                            'pdb_data': pdb_data,
                            'source': source
                        })
                    
                    # If all methods failed to fetch
                    if not any(res['pdb_url'] or res['pdb_data'] for res in alphafold_results):
                        error_message = "Could not generate 3D structures. The sequences might be too long or malformatted."
                else:
                    error_message = "Please enter valid FASTA sequences."
    else:
        form = ProteinSearchForm()

    from datetime import datetime
    year = datetime.utcnow().year
    return render(request, 'go_predictor/structure_2d.html', {
        'form': form, 
        'structures': structures, 
        'alphafold_results': alphafold_results,
        'error_message': error_message,
        'year': year
    })


def call_deepgo_api(sequence):
    url = 'https://deepgo.cbrc.kaust.edu.sa/deepgo/api/create'
    json_data = {
        'data_format': 'enter',
        'data': sequence,
        'threshold': 0.3,
        'version': '1.0.3'
    }
    try:
        response = requests.post(url, json=json_data)
        response.raise_for_status()
        api_result = response.json()
        return parse_deepgo_result(api_result)
    except requests.RequestException as e:
        print(f"API Error: {e}")  # for debugging
        return None


def parse_deepgo_result(api_result):
    if not api_result or 'predictions' not in api_result or not api_result['predictions']:
        return None
    pred = api_result['predictions'][0]
    protein_info = pred.get('protein_info')
    if protein_info:
        name = protein_info.split()[0] if '|' in protein_info else protein_info
    else:
        name = 'Unknown Protein'
    
    categories = []
    for func_cat in pred.get('functions', []):
        cat_name = func_cat['name']
        if cat_name == 'Biological Process':
            title = 'Biological Processes'
            subtitle = 'Biological Process'
        elif cat_name == 'Molecular Function':
            title = 'Molecular Functions'
            subtitle = 'Molecular Function'
        elif cat_name == 'Cellular Component':
            title = 'Cellular Components'
            subtitle = 'Cellular Component'
            
        else:
            continue
        terms = [{'go': go, 'label': label, 'score': score} for go, label, score in func_cat['functions']]
        categories.append({
            'title': title,
            'subtitle': subtitle,
            'terms': terms
        })
    
    return {
        'protein': {
            'name': name,
            'description': 'Predicted by DeepGO'
        },
        'categories': categories
    }