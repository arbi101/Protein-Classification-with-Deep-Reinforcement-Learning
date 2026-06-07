# go_predictor/forms.py
from django import forms

class ProteinSearchForm(forms.Form):
    fasta_sequence = forms.CharField(
        label='Protein Sequence in FASTA format',
        error_messages={
            'required': 'Enter at least one FASTA sequence.',
        },
        widget=forms.Textarea(attrs={
            'placeholder': '>ProteinName\nSEQUENCE...',
            'class': 'form-control mt-2',
            'rows': 10,
            'autocomplete': 'off',
            'required': True,
        }),
        help_text='Enter one or more sequences in FASTA format. Maximum 10 sequences.'
    )

    def clean_fasta_sequence(self):
        fasta = self.cleaned_data['fasta_sequence'].strip()
        if not fasta:
            raise forms.ValidationError('Enter at least one FASTA sequence.')

        sequence_count = 0
        current_seq = []

        for line in fasta.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequence_count += 1
                    current_seq = []
                continue
            if not line.isalpha():
                raise forms.ValidationError(
                    'Use only amino-acid letters in sequence lines.'
                )
            current_seq.append(line)

        if current_seq:
            sequence_count += 1

        if sequence_count == 0:
            raise forms.ValidationError(
                'Enter a FASTA header followed by an amino-acid sequence.'
            )
        if sequence_count > 10:
            raise forms.ValidationError(
                'Enter no more than 10 FASTA sequences.'
            )

        return fasta
