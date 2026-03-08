from django.shortcuts import render
from .forms import ProteinSearchForm
import requests
import random

def predict_go(request):
    result = None
    structure = None
    hp_string = None
    energy = 0
    error_message = None
    if request.method == 'POST':
        form = ProteinSearchForm(request.POST)
        if form.is_valid():
            fasta_str = form.cleaned_data['fasta_sequence'].strip()
            action = request.POST.get('action')
            if action == 'predict':
                if fasta_str:
                    sequence = parse_fasta(fasta_str)
                    if sequence:
                        result = call_deepgo_api(sequence)
                        if result is None:
                            error_message = "Error in GO prediction. Check the FASTA sequence."
                    else:
                        error_message = "Invalid FASTA sequence. Make sure it contains an amino acid sequence."
                else:
                    error_message = "Enter a valid FASTA sequence."
            elif action == 'structure':
                if fasta_str:
                    sequence = parse_fasta(fasta_str)
                    if sequence:
                        hp_string = sequence_to_hp(sequence)
                        positions, energy = generate_2d_structure(hp_string)
                        structure = [{'x': x, 'y': y, 'type': hp_string[i]} for i, (x, y) in enumerate(positions)]
                    else:
                        error_message = "Invalid FASTA sequence."
                else:
                    error_message = "Enter a valid FASTA sequence."
    else:
        form = ProteinSearchForm()

    # add current year for footer
    from datetime import datetime
    year = datetime.utcnow().year
    return render(request, 'go_predictor/predict.html', {
        'form': form, 
        'result': result, 
        'structure': structure, 
        'hp_string': hp_string, 
        'energy': energy, 
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
    Removes the header and concatenates any subsequent lines
    Returns only the amino acid sequence in uppercase
    """
    lines = fasta_str.strip().splitlines()
    seq = ''
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            continue  # skip header
        seq += line.upper()
    return seq


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


def generate_2d_structure(hp_string):
    """
    Generates a simple 2D structure using self-avoiding walk
    Returns list of positions (x,y) and energy
    """
    if not hp_string:
        return [], 0
    
    positions = [(0, 0)]  # Start from (0,0)
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]  # Up, Right, Down, Left
    
    for i in range(1, len(hp_string)):
        # Try random directions until finding a free one
        random.shuffle(directions)
        placed = False
        for dx, dy in directions:
            nx, ny = positions[-1][0] + dx, positions[-1][1] + dy
            if (nx, ny) not in positions:
                positions.append((nx, ny))
                placed = True
                break
        if not placed:
            # If not found, stop (for simplicity)
            break
    
    # Calculate energy: number of non-consecutive HH contacts
    energy = 0
    for i in range(len(positions)):
        if hp_string[i] == 'H':
            x, y = positions[i]
            # Check non-consecutive neighbors
            for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
                nx, ny = x + dx, y + dy
                if (nx, ny) in positions:
                    j = positions.index((nx, ny))
                    if abs(i - j) > 1 and hp_string[j] == 'H':
                        energy += 1
    energy //= 2  # Each contact counted twice
    
    return positions, energy


def structure_2d(request):
    structure = None
    hp_string = None
    energy = 0
    if request.method == 'POST':
        form = ProteinSearchForm(request.POST)
        if form.is_valid():
            fasta_str = form.cleaned_data['fasta_sequence'].strip()
            if fasta_str:
                sequence = parse_fasta(fasta_str)
                if sequence:
                    hp_string = sequence_to_hp(sequence)
                    positions, energy = generate_2d_structure(hp_string)
                    # Convert positions to list for template
                    structure = [{'x': x, 'y': y, 'type': hp_string[i]} for i, (x, y) in enumerate(positions)]
    else:
        form = ProteinSearchForm()

    from datetime import datetime
    year = datetime.utcnow().year
    return render(request, 'go_predictor/structure_2d.html', {
        'form': form, 
        'structure': structure, 
        'hp_string': hp_string, 
        'energy': energy, 
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