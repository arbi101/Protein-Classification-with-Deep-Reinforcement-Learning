from django.db import models

class ProteinSequence(models.Model):
    name = models.CharField(max_length=200)
    fasta_sequence = models.TextField()
    predicted_go_terms = models.TextField(blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.name
