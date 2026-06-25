"""Forms and server-side validation for protein sequence input."""

from django import forms


class ProteinSearchForm(forms.Form):
    """Collect and validate one or more protein sequences in FASTA form."""

    # This is a regular Form rather than a ModelForm because the current
    # workflows process the input immediately and do not save it to the DB.
    fasta_sequence = forms.CharField(
        label='Protein Sequence in FASTA format',
        # Override Django's default message for an empty required field.
        error_messages={
            'required': 'Enter at least one FASTA sequence.',
        },
        # Render the string field as a multiline HTML textarea.
        widget=forms.Textarea(attrs={
            # These attributes become attributes of the generated HTML element.
            'placeholder': '>ProteinName\nSEQUENCE...',
            'class': 'form-control mt-2',
            'rows': 10,
            'autocomplete': 'off',
            'required': True,
        }),
        help_text='Enter one or more sequences in FASTA format. Maximum 10 sequences.'
    )

    def clean_fasta_sequence(self):
        """Validate FASTA-like contents and return the normalized input."""

        # cleaned_data contains Django's initial conversion of the submitted
        # field. Strip outer whitespace before applying custom validation.
        fasta = self.cleaned_data['fasta_sequence'].strip()
        if not fasta:
            raise forms.ValidationError('Enter at least one FASTA sequence.')

        # Count completed sequences and collect the lines of the current one.
        sequence_count = 0
        current_seq = []

        for line in fasta.splitlines():
            # Whitespace-only lines are irrelevant in FASTA input.
            line = line.strip()
            if not line:
                continue

            # A line beginning with '>' starts a new FASTA record.
            if line.startswith('>'):
                # If residues were already collected, the previous record has
                # ended and can now be counted.
                if current_seq:
                    sequence_count += 1
                    current_seq = []
                continue

            # Sequence lines must contain letters only. This is a syntactic
            # check; it does not restrict input to the 20 standard amino acids.
            if not line.isalpha():
                raise forms.ValidationError(
                    'Use only amino-acid letters in sequence lines.'
                )
            # Keep track of the sequence even when it spans multiple lines.
            current_seq.append(line)

        # The final FASTA record has no following header to trigger its count.
        if current_seq:
            sequence_count += 1

        # Reject input containing headers only (or otherwise no residues).
        if sequence_count == 0:
            raise forms.ValidationError(
                'Enter a FASTA header followed by an amino-acid sequence.'
            )
        # Bound work performed by API calls and prediction algorithms.
        if sequence_count > 10:
            raise forms.ValidationError(
                'Enter no more than 10 FASTA sequences.'
            )

        # The returned value becomes form.cleaned_data['fasta_sequence'].
        return fasta
