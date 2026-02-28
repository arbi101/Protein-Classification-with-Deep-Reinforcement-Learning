# go_predictor/forms.py
from django import forms

class FastaForm(forms.Form):
    name = forms.CharField(label='Protein Name', max_length=200)
    fasta_sequence = forms.CharField(label='FASTA Sequence', widget=forms.Textarea)