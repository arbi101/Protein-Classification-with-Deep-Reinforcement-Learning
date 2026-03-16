from django import template

register = template.Library()

@register.filter
def get_score_color(score):
    try:
        # Handle both string and numeric inputs
        if isinstance(score, str):
            score = float(score)
        elif not isinstance(score, (int, float)):
            score = float(score)
        
        if score >= 0.7:
            return 'bg-emerald-100 text-emerald-800 border-emerald-200'
        elif score >= 0.6:
            return 'bg-blue-100 text-blue-800 border-blue-200'
        elif score >= 0.5:
            return 'bg-amber-100 text-amber-800 border-amber-200'
        else:
            return 'bg-slate-100 text-slate-800 border-slate-200'
    except (ValueError, TypeError):
        return 'bg-slate-100 text-slate-800 border-slate-200'