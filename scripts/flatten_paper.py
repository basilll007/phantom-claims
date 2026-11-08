"""Write a self-contained submission folder: paper/build/ with ONE .tex file (the AAAI-27 kit requires a single source
file without \\input for camera-ready; the .bib may stay separate), the figures it uses, refs.bib and the kit's
aaai2027.sty / aaai2027.bst. Upload the folder to Overleaf, or run scripts/build_paper.sh where pdflatex exists.

    python scripts/flatten_paper.py
"""
import os, re, shutil

SRC, OUT = "paper", "paper/build"
tex = open(f"{SRC}/main.tex", encoding="utf-8").read()


def inline(m):
    name = m.group(1)
    path = f"{SRC}/{name if name.endswith('.tex') else name + '.tex'}"
    body = open(path, encoding="utf-8").read().rstrip("\n")
    return f"% ---- begin inlined {name} ----\n{body}\n% ---- end inlined {name} ----"


flat = re.sub(r"\\input\{([^}]+)\}", inline, tex)
assert "\\input{" not in flat
if os.path.isdir(OUT):
    shutil.rmtree(OUT)
os.makedirs(f"{OUT}/figures")
open(f"{OUT}/phantom_claims.tex", "w", encoding="utf-8", newline="\n").write(flat)
figs = re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", flat)
for f in figs:
    shutil.copy(f"{SRC}/{f}", f"{OUT}/{f}")
for f in ("refs.bib", "aaai2027.sty", "aaai2027.bst"):
    shutil.copy(f"{SRC}/{f}", f"{OUT}/{f}")
print(f"{OUT}/phantom_claims.tex  (+ {len(figs)} figures, refs.bib, aaai2027.sty, aaai2027.bst)")
