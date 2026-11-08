"""Pre-flight check of paper/main.tex against the AAAI-27 author kit rules and the project's number rules.

    python scripts/check_paper.py          # exit 1 on any problem

Checks: every \\macro{} is defined in numbers.tex and none resolves to TODO-VERIFY; every citation key is in refs.bib
(and lists bib entries never cited); every \\ref has a \\label; commands and packages that the kit forbids
(AnonymousSubmission2027.tex, tables "Commands that must not be used" / "LaTeX style packages that must not be used");
every \\includegraphics file exists in paper/ and is a .png/.jpg/.pdf; no \\input of anything but numbers.tex and tables.
"""
import os, re, sys

P = "paper"
tex = open(f"{P}/main.tex", encoding="utf-8").read()
body = re.sub(r"(?m)(?<!\\)%.*$", "", tex)                       # drop comments
inputs = re.findall(r"\\input\{([^}]+)\}", body)
full = body + "".join(open(f"{P}/{i if i.endswith('.tex') else i + '.tex'}", encoding="utf-8").read() for i in inputs)
nums = open(f"{P}/numbers.tex", encoding="utf-8").read()
bib = open(f"{P}/refs.bib", encoding="utf-8").read()
bad = []

defined = dict(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}", nums))
local = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", body))
used = set(re.findall(r"\\([A-Za-z]+)\{\}", body))
builtin = {"ldots", "maketitle", "midrule", "toprule", "bottomrule", "centering", "small", "appendix"}
undef = sorted(u for u in used if u not in defined and u not in builtin and u not in local)
todo = sorted(u for u in used if "TODOVERIFY" in defined.get(u, ""))
if undef: bad.append(f"undefined macros: {undef}")
if todo: bad.append(f"macros resolving to TODO-VERIFY: {todo}")

cites = {k.strip() for g in re.findall(r"\\(?:cite|citeauthor|shortcite|citeyear|citep|citet)\{([^}]*)\}", body)
         for k in g.split(",")}
keys = set(re.findall(r"@\w+\{([^,\s]+),", bib))
if cites - keys: bad.append(f"cited but not in refs.bib: {sorted(cites - keys)}")
uncited = sorted(keys - cites)

labels = set(re.findall(r"\\label\{([^}]+)\}", full))
refs = set(re.findall(r"\\ref\{([^}]+)\}", full))
if refs - labels: bad.append(f"\\ref without \\label: {sorted(refs - labels)}")

FORBIDDEN_CMDS = [r"\\abovecaption", r"\\abovedisplay", r"\\addevensidemargin", r"\\addsidemargin", r"\\addtolength",
                  r"\\baselinestretch", r"\\belowcaption", r"\\belowdisplay", r"\\break\b", r"\\clearpage", r"\\clip\b",
                  r"\\columnsep", r"\\float\b", r"\\linespread", r"\\newpage", r"\\pagebreak", r"\\textheight",
                  r"\\tiny", r"\\topmargin", r"\\trim\b", r"\\vskip\s*\{?\s*-", r"\\vspace\*?\{\s*-", r"\\pagestyle",
                  r"\\nocopyright", r"\\resizebox", r"\\textcolor", r"\\color\{"]
for pat in FORBIDDEN_CMDS:
    if re.search(pat, full):
        bad.append(f"forbidden command: {pat}")
for m in re.finditer(r"\\setlength\{\\([A-Za-z]+)\}", full):
    if m.group(1) != "tabcolsep":
        bad.append(f"forbidden \\setlength{{\\{m.group(1)}}}")
FORBIDDEN_PKGS = {"authblk", "babel", "balance", "cjk", "epsf", "epsfig", "euler", "float", "fullpage", "geometry",
                  "graphics", "hyperref", "layout", "lmodern", "navigator", "pdfcomment", "pgfplots", "psfig",
                  "pstricks", "t1enc", "times", "titlesec", "tocbibind", "ulem", "flushend", "indentfirst", "multicol",
                  "nameref", "savetrees", "setspace", "stfloats", "tabu", "wrapfig", "helvet", "courier"}
for g in re.findall(r"\\usepackage(?:\[[^]]*\])?\{([^}]+)\}", body):
    for p in g.split(","):
        if p.strip().lower() in FORBIDDEN_PKGS:
            bad.append(f"forbidden package: {p.strip()}")

for f in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", full):
    if not os.path.exists(f"{P}/{f}"):
        bad.append(f"missing figure: {f}")
    elif not f.lower().endswith((".png", ".jpg", ".pdf")):
        bad.append(f"figure format not allowed by the kit: {f}")
for i in inputs:
    if not (i.startswith("numbers") or i.startswith("tables/")):
        bad.append(f"unexpected \\input: {i}")

print(f"macros used {len(used)}, citations {len(cites)}, bib entries {len(keys)}, labels {len(labels)}, figures "
      f"{len(re.findall(r'includegraphics', full))}")
if uncited:
    print("note: bib entries not cited (BibTeX drops them):", uncited)
tv = re.findall(r"\\TODOVERIFY\{([^}]*)\}", body)
if tv:
    print(f"note: {len(tv)} TODO-VERIFY placeholder(s) in the text (print in bold in the PDF): {tv}")
print("\n".join("PROBLEM: " + b for b in bad) if bad else "OK: no problems found")
sys.exit(1 if bad else 0)
