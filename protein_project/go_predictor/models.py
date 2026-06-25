from django.db import models


class ProteinSequence(models.Model):
    """Optional database record for a protein and its GO predictions.

    The current view does not save these records yet, but the model prepares
    the schema for a future prediction history feature.
    """

    # Human-readable protein or FASTA record name.
    name = models.CharField(max_length=200)
    # TextField supports sequences longer than a short CharField limit.
    fasta_sequence = models.TextField()
    # Predictions may be empty before an analysis has completed.
    predicted_go_terms = models.TextField(blank=True)
    # Set once, automatically, when the database row is first created.
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        """Use the protein name in the admin site and Python displays."""
        return self.name
