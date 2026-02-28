from django.shortcuts import render
from .forms import FastaForm

def predict_go(request):
    go_terms = None
    if request.method == 'POST':
        form = FastaForm(request.POST)
        if form.is_valid():
            fasta_input = form.cleaned_data['fasta_sequence']
            sequence = parse_fasta(fasta_input)
            # qui chiamerai il modello reale DeepGO
            go_terms = fake_predict_go(sequence)  # al momento placeholder
    else:
        form = FastaForm()
    
    return render(request, 'go_predictor/predict.html', {'form': form, 'go_terms': go_terms})


def parse_fasta(fasta_str):
    """
    Rimuove l'header e concatena eventuali righe successive
    Restituisce solo la sequenza di amminoacidi in maiuscolo
    """
    lines = fasta_str.strip().splitlines()
    seq = ''
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            continue  # salta header
        seq += line.upper()
    return seq


def fake_predict_go(sequence):
    # placeholder temporaneo
    return ['GO:0008150', 'GO:0003674']