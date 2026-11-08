"""Automated registry check of paper/refs.bib (tooling for lane F's citation audit; it does not replace lane F's sign-off).

For every entry, the DOI is read from the '% doi:<doi>' comment line above it, or the arXiv id from its eprint field.
The registry record is fetched (Crossref API for DOIs, arXiv API for eprints) and compared with the bib entry:
normalised title, first author's family name, year, and for journal articles volume and first page / article number.
Year: Crossref print-issue year, falling back to online / issued (a bib year equal to either counts as a match).

  python scripts/citation_audit.py      # -> results/tables/citation_audit.csv ; exit 1 if any field mismatches
"""
import csv, json, re, sys, time, unicodedata, urllib.parse, urllib.request

UA = {"User-Agent": "phantom-claims-citation-audit/0.1"}


def norm(s):
    s = unicodedata.normalize("NFKD", s or "")
    s = re.sub(r"\\[a-zA-Z]+|[{}\\'\"`^~]", "", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def parse_bib(path):
    txt = open(path, encoding="utf-8").read()
    starts = [m.start() for m in re.finditer(r"^@\w+\{", txt, re.M)] + [len(txt)]
    out = []
    for a, b in zip(starts[:-1], starts[1:]):
        chunk = txt[a:b]
        typ, key = re.match(r"@(\w+)\{([^,]+),", chunk).groups()
        prev = txt[max(0, a - 300):a].strip().splitlines()
        cm = re.match(r"%\s*(doi|arXiv):(\S+)", prev[-1]) if prev else None
        f = dict((k.lower(), v.strip()) for k, v in
                 re.findall(r"(\w+)\s*=\s*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}", chunk))   # 2 nesting levels
        f.setdefault("eprint", "")
        out.append(dict(f, key=key, type=typ, doi=cm.group(2) if cm and cm.group(1) == "doi" else ""))
    return out


def crossref(doi):
    m = json.load(urllib.request.urlopen(urllib.request.Request("https://api.crossref.org/works/" + doi, headers=UA),
                                         timeout=40))["message"]
    years = {(m.get(k) or {}).get("date-parts", [[None]])[0][0] for k in ("published-print", "published-online", "issued")}
    au = m.get("author", [])
    return dict(title=(m.get("title") or [""])[0], first=au[0].get("family", au[0].get("name", "")) if au else "",
                years={y for y in years if y}, volume=m.get("volume", ""), page=m.get("page") or m.get("article-number", ""),
                venue=(m.get("container-title") or [""])[0])


def arxiv(eid):
    x = urllib.request.urlopen("http://export.arxiv.org/api/query?id_list=" + eid, timeout=60).read().decode()
    e = x.split("<entry>")[1]
    t = re.sub(r"\s+", " ", re.search(r"<title>(.*?)</title>", e, re.S).group(1)).strip()
    au = re.findall(r"<name>(.*?)</name>", e)
    return dict(title=t, first=au[0].split()[-1] if au else "", years={int(re.search(r"<published>(\d{4})", e).group(1))},
                volume="", page="", venue="arXiv")


def main():
    rows, bad = [], 0
    for b in parse_bib("paper/refs.bib"):
        try:
            r = crossref(b["doi"]) if b["doi"] else arxiv(b["eprint"]) if b["eprint"] else None
        except Exception as e:
            r = None; err = repr(e)
        else:
            err = ""
        if r is None:
            rows.append(dict(key=b["key"], ok=False, error=err or "no DOI/eprint to check")); bad += 1; continue
        first_bib = re.split(r"\s+and\s+", b.get("author", ""))[0].split(",")[0]
        chk = dict(title=norm(b.get("title")) == norm(r["title"]),
                   first_author=norm(first_bib) == norm(r["first"]) or norm(r["first"]) in norm(first_bib),
                   year=int(b.get("year", 0)) in r["years"])
        if b["type"] == "article":
            chk["volume"] = str(b.get("volume", "")) == str(r["volume"])
            p = str(b.get("pages", "")).split("-")[0].replace("e", "").strip()
            chk["page"] = p != "" and p in str(r["page"]).replace("e", "").replace("MSB2020", "")
        ok = all(chk.values()); bad += not ok
        rows.append(dict(key=b["key"], id=b["doi"] or "arXiv:" + b["eprint"], ok=ok, **{f"match_{k}": v for k, v in chk.items()},
                         registry_title=r["title"], registry_first_author=r["first"],
                         registry_years=";".join(map(str, sorted(r["years"]))), registry_volume=r["volume"],
                         registry_page=r["page"], registry_venue=r["venue"]))
        time.sleep(0.4)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "key", k != "id", k != "ok", k))
    with open("results/tables/citation_audit.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    for r in rows:
        if not r["ok"]:
            print("MISMATCH", {k: v for k, v in r.items() if k.startswith(("key", "match", "registry", "error"))})
    print(f"{len(rows)} entries checked, {len(rows) - bad} match the registry, {bad} need attention")
    return bad


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
