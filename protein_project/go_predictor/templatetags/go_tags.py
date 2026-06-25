"""Custom presentation filters used by the protein prediction templates."""

from django import template

# Django discovers filters through this Library instance when a template uses
# ``{% load go_tags %}``.
register = template.Library()


@register.filter
def get_score_color(score):
    """Return Tailwind classes representing a DeepGO confidence score.

    Higher-confidence GO terms receive stronger colors. A neutral style is
    returned when the input is missing or cannot be converted to a number.
    """

    try:
        # API/template values may arrive as strings, integers or floats.
        if isinstance(score, str):
            score = float(score)
        elif not isinstance(score, (int, float)):
            score = float(score)

        # Keep score-to-style mapping out of the HTML template so the same
        # visual rule can be reused wherever GO terms are displayed.
        if score >= 0.7:
            return 'bg-emerald-100 text-emerald-800 border-emerald-200'
        elif score >= 0.6:
            return 'bg-blue-100 text-blue-800 border-blue-200'
        elif score >= 0.5:
            return 'bg-amber-100 text-amber-800 border-amber-200'
        else:
            return 'bg-slate-100 text-slate-800 border-slate-200'
    except (ValueError, TypeError):
        # Invalid scores must not break the rendering of the entire page.
        return 'bg-slate-100 text-slate-800 border-slate-200'
