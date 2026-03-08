# go_predictor/forms.py
from django import forms

class ProteinSearchForm(forms.Form):
    fasta_sequence = forms.CharField(
        label='Protein Sequence in FASTA format',
        widget=forms.Textarea(attrs={
            'placeholder': '>ProteinName\nSEQUENCE...',
            'class': 'form-control mt-2',
            'rows': 10,
            'autocomplete': 'off'
        }),
        help_text='Enter one or more sequences in FASTA format. Maximum 10 sequences.'
    )
