"""Check that every \\macro{} used in paper/main.tex is defined in paper/numbers.tex (or is a LaTeX built-in),
and list macros that still resolve to TODO-VERIFY."""
import re
tex = open("paper/main.tex", encoding="utf-8").read()
nums = open("paper/numbers.tex", encoding="utf-8").read()
defined = dict(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}", nums))
used = set(re.findall(r"\\([A-Za-z]+)\{\}", tex))
builtin = {"maketitle", "todo", "linewidth", "midrule", "toprule", "bottomrule", "centering", "small"}
undef = sorted(u for u in used if u not in defined and u not in builtin)
todo = sorted(u for u in used if u in defined and "TODOVERIFY" in defined[u])
print("used:", len(used), "| undefined:", undef, "| resolving to TODO-VERIFY:", todo)
