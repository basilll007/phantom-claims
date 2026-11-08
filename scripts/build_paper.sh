#!/usr/bin/env bash
# One command from logged results to the submission PDF:
#   numbers.tex (from results.jsonl) -> figures -> pre-flight check -> single-file source -> pdflatex/bibtex.
# Needs PC_DATA for the overview figure (it reads the processed tissue). pdflatex/bibtex are optional: without them the
# script stops after writing paper/build/, which can be uploaded to Overleaf as is.
set -euo pipefail
PY=${PYTHON:-.venv/Scripts/python.exe}
"$PY" scripts/make_paper_numbers.py
"$PY" -m scripts.paper_figures
"$PY" scripts/check_paper.py
"$PY" scripts/flatten_paper.py
if command -v pdflatex >/dev/null 2>&1 && command -v bibtex >/dev/null 2>&1; then
  cd paper/build
  pdflatex -interaction=nonstopmode -halt-on-error phantom_claims.tex >/dev/null
  bibtex phantom_claims >/dev/null
  pdflatex -interaction=nonstopmode -halt-on-error phantom_claims.tex >/dev/null
  pdflatex -interaction=nonstopmode -halt-on-error phantom_claims.tex >/dev/null
  grep -E "Overfull|Undefined|Citation .* undefined|Reference .* undefined" phantom_claims.log || true
  echo "paper/build/phantom_claims.pdf"
else
  echo "pdflatex not found: upload paper/build/ to Overleaf (main file phantom_claims.tex) to compile."
fi
