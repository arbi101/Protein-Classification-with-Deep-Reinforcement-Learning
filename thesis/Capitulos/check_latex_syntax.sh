#!/bin/bash

FILES=("capitulo2.tex" "capitulo4.tex" "conclusao.tex")
ERRORS_FOUND=0

echo "=== LaTeX Syntax Check ==="
echo ""

for file in "${FILES[@]}"; do
    if [ ! -f "$file" ]; then
        echo "❌ File not found: $file"
        ERRORS_FOUND=$((ERRORS_FOUND + 1))
        continue
    fi
    
    echo "Checking: $file"
    echo "---"
    
    # Check for unmatched braces
    open_braces=$(grep -o '{' "$file" | wc -l)
    close_braces=$(grep -o '}' "$file" | wc -l)
    if [ "$open_braces" -ne "$close_braces" ]; then
        echo "⚠️  Unmatched braces: $open_braces '{' vs $close_braces '}'"
        ERRORS_FOUND=$((ERRORS_FOUND + 1))
    else
        echo "✅ Braces balanced: $open_braces pairs"
    fi
    
    # Check for unmatched \begin and \end
    begins=$(grep -o '\\begin{' "$file" | wc -l)
    ends=$(grep -o '\\end{' "$file" | wc -l)
    if [ "$begins" -ne "$ends" ]; then
        echo "⚠️  Unmatched \begin/\end: $begins \\begin vs $ends \\end"
        ERRORS_FOUND=$((ERRORS_FOUND + 1))
    else
        echo "✅ Begin/End balanced: $begins pairs"
    fi
    
    # Check for table syntax issues
    if grep -q '\\begin{table' "$file"; then
        table_count=$(grep -c '\\begin{table' "$file")
        table_end_count=$(grep -c '\\end{table' "$file")
        if [ "$table_count" -ne "$table_end_count" ]; then
            echo "⚠️  Unmatched table environments: $table_count \\begin{table} vs $table_end_count \\end{table}"
            ERRORS_FOUND=$((ERRORS_FOUND + 1))
        else
            echo "✅ Tables balanced: $table_count pairs"
        fi
    fi
    
    # Check for missing labels in important sections
    if grep -q '\\section\|\\subsection\|\\chapter' "$file"; then
        sections=$(grep -c '\\section\|\\subsection\|\\chapter' "$file")
        labels=$(grep -c '\\label{' "$file")
        if [ "$labels" -lt "$sections" ]; then
            echo "⚠️  Possible missing labels: $sections sections/subsections vs $labels labels"
        else
            echo "✅ Labels present for sections"
        fi
    fi
    
    # Check for broken command patterns
    broken_patterns=$(grep -cE '\\\\[a-zA-Z]+\s*$' "$file" || true)
    if [ "$broken_patterns" -gt 0 ]; then
        echo "⚠️  Found incomplete commands at line endings:"
        grep -nE '\\\\[a-zA-Z]+\s*$' "$file" | head -3
    fi
    
    echo ""
done

echo "=== Summary ==="
if [ $ERRORS_FOUND -eq 0 ]; then
    echo "✅ All syntax checks passed!"
else
    echo "⚠️  Found $ERRORS_FOUND issues - review recommended"
fi
