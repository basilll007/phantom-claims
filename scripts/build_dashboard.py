#!/usr/bin/env python3
"""Build the project dashboard: one self-contained HTML file (data embedded, no external libraries).
Standard library only, so it runs locally and in GitHub Actions.

Usage:
  python scripts/build_dashboard.py                 # -> dashboard/index.html (links to ../figures, ../results)
  python scripts/build_dashboard.py --site _site    # -> _site/index.html + copied assets (GitHub Pages)

Every number shown comes from files under results/ (results.jsonl, tables/, agent_logs/, claims/)."""
import argparse, csv, glob, html, json, os, re, shutil, subprocess, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def p(*a):
    return os.path.join(ROOT, *a)


def read(rel, default=""):
    try:
        with open(p(rel), encoding="utf-8", errors="replace") as f:
            return f.read()
    except OSError:
        return default


# ---------------------------------------------------------------- minimal markdown -> HTML
def _inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<![*\w])\*([^*\n]+)\*(?![*\w])", r"<em>\1</em>", s)
    s = re.sub(r"\[([^\]]+)\]\((https?://[^)\s]+)\)", r'<a href="\2" target="_blank" rel="noopener">\1</a>', s)
    return s


def md_to_html(md):
    out, lines, i = [], md.splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if not ln.strip():
            i += 1; continue
        m = re.match(r"^(#{1,6})\s+(.*)", ln)
        if m:
            lvl = min(len(m.group(1)) + 1, 6)
            out.append(f"<h{lvl}>{_inline(m.group(2))}</h{lvl}>"); i += 1; continue
        if ln.lstrip().startswith("|") and i + 1 < len(lines) and re.match(r"^\s*\|[\s:\-|]+\|\s*$", lines[i + 1]):
            head = [c.strip() for c in ln.strip().strip("|").split("|")]
            i += 2; body = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                body.append([c.strip() for c in lines[i].strip().strip("|").split("|")]); i += 1
            t = "<div class='tw'><table><thead><tr>" + "".join(f"<th>{_inline(c)}</th>" for c in head) + "</tr></thead><tbody>"
            t += "".join("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in r) + "</tr>" for r in body)
            out.append(t + "</tbody></table></div>"); continue
        if re.match(r"^\s*([-*]|\d+\.)\s+", ln):
            ordered = bool(re.match(r"^\s*\d+\.", ln)); items = []
            while i < len(lines) and re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]):
                items.append(re.sub(r"^\s*([-*]|\d+\.)\s+", "", lines[i])); i += 1
                while i < len(lines) and lines[i].startswith("   ") and lines[i].strip():
                    items[-1] += " " + lines[i].strip(); i += 1
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>" + "".join(f"<li>{_inline(x)}</li>" for x in items) + f"</{tag}>"); continue
        para = [ln]; i += 1
        while i < len(lines) and lines[i].strip() and not re.match(r"^(#{1,6}\s|\s*\||\s*([-*]|\d+\.)\s)", lines[i]):
            para.append(lines[i]); i += 1
        out.append("<p>" + _inline(" ".join(x.strip() for x in para)) + "</p>")
    return "\n".join(out)


# ---------------------------------------------------------------- data collection
def load_results():
    recs = []
    for ln in read("results/results.jsonl").splitlines():
        ln = ln.strip()
        if ln:
            try:
                recs.append(json.loads(ln))
            except json.JSONDecodeError:
                pass
    latest, versions = {}, {}
    for r in recs:
        latest[r["id"]] = r; versions[r["id"]] = versions.get(r["id"], 0) + 1
    out = []
    for rid, r in latest.items():
        m = r.get("metrics") or {}
        flat = {k: v for k, v in m.items() if isinstance(v, (int, float)) and not isinstance(v, bool)}
        out.append(dict(id=rid, phase=r.get("phase"), split=r.get("split"), tissue=r.get("tissue"),
                        condition=r.get("condition"), K=r.get("K"), hypothesis=r.get("hypothesis"),
                        ablation=r.get("ablation"), family=r.get("family"), exploratory=bool(r.get("exploratory")),
                        timestamp=r.get("timestamp"), git=r.get("git_hash"), versions=versions[rid],
                        metrics=flat, ci=r.get("ci") or {}, raw=json.dumps(r, default=str)[:4000]))
    out.sort(key=lambda x: x.get("timestamp") or "")
    return out, len(recs), {r["id"]: r for r in latest.values()}


def load_claims():
    rows = {}
    for f in glob.glob(p("results", "claims", "agent_*_claims.csv")):
        with open(f, encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("rank") in ("0", "", None):
                    continue
                rows.setdefault(r["run_id"], []).append(r)
    return rows


def load_agents():
    claims = load_claims()
    runs = []
    for f in sorted(glob.glob(p("results", "agent_logs", "*.json"))):
        name = os.path.basename(f)[:-5]
        if name == "budget":
            continue
        try:
            l = json.load(open(f, encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        m = re.match(r"^(dev|test)_(.+)_s(\d+)_(gemini|deepseek)_(unguided|gate_tool)$", name)
        if not m:
            continue
        seq = [tc["name"] for s in l.get("steps", []) for tc in s.get("tool_calls", [])]
        cl = []
        for r in claims.get(name, []):
            cl.append(dict(rank=int(r["rank"]), sender=r["sender"], receiver=r["receiver"], ligand=r["ligand"],
                           receptor=r["receptor"], label=int(float(r["label"])) if r.get("label") not in ("", None) else None,
                           neg_type=r.get("neg_type") or "", gate=r.get("gate_pass") or "",
                           evidence=(r.get("evidence") or "")[:300]))
        cl.sort(key=lambda c: c["rank"])
        n = len(cl); nt = sum(1 for c in cl if c["label"] == 1)
        u = l.get("usage", {})
        runs.append(dict(id=name, split=m.group(1), tissue=m.group(2), seed=int(m.group(3)), family=m.group(4),
                         model=l.get("model"), cond=m.group(5), status=l.get("status"),
                         prompt=u.get("prompt", 0), completion=u.get("completion", 0),
                         tool_calls=l.get("n_tool_calls"), runtime=l.get("runtime_s"), started=l.get("started"),
                         seq=seq, claims=cl, n_claims=n, n_true=nt,
                         phantom=(None if n == 0 else round(1 - nt / n, 4)),
                         n_gate_pass=sum(1 for c in cl if c["gate"] == "True")))
    return runs


def load_progress():
    out = []
    for ln in read("PROGRESS.md").splitlines():
        m = re.match(r"^- (\d{4}-\d{2}-\d{2} \d{2}:\d{2}(?: [A-Z]+)?) — (.*)$", ln)
        if m:
            out.append(dict(ts=m.group(1), text=m.group(2)))
    return out


def load_commits():
    try:
        r = subprocess.run(["git", "log", "--date=format:%Y-%m-%d %H:%M", "--pretty=format:%h\t%ad\t%s"],
                           cwd=ROOT, capture_output=True, text=True, timeout=30)
        return [dict(zip(("hash", "date", "msg"), x.split("\t", 2))) for x in r.stdout.splitlines() if x.count("\t") >= 2]
    except Exception:
        return []


def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    except Exception:
        return ""


def read_csv(rel):
    try:
        with open(p(rel), encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except OSError:
        return []


def collect():
    results, n_lines, by_id = load_results()
    verd = {}
    try:
        verd = json.load(open(p("results", "tables", "verdicts.json"), encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        pass
    try:
        budget = json.load(open(p("results", "agent_logs", "budget.json"), encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        budget = {}
    hard = (by_id.get("dev_hard_list_true_vs_niche") or {}).get("metrics", {})
    prereg = read("PREREG.md")
    amend = prereg[prereg.find("## Amendments"):] if "## Amendments" in prereg else ""
    figs = sorted(os.path.basename(x) for x in glob.glob(p("figures", "*.png")))
    return dict(
        generated=time.strftime("%Y-%m-%d %H:%M %Z"), head=git_head(),
        state=md_to_html(read("STATE.md")), morning=md_to_html(read("MORNING_REPORT.md")),
        amendments=md_to_html(amend), verdicts=verd, results=results, n_result_lines=n_lines,
        agents=load_agents(), budget=budget, progress=load_progress(), commits=load_commits(),
        planted=read_csv("results/tables/planted_pass_rates.csv"), hard=hard, figures=figs,
        budget_cap=dict(tokens=4_000_000, usd=5.0))


# ---------------------------------------------------------------- page
PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Phantom Claims — project dashboard</title>
<style>
.viz-root{color-scheme:light;--page:#f9f9f7;--surface-1:#fcfcfb;--text-primary:#0b0b0b;--text-secondary:#52514e;
 --muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;--border:rgba(11,11,11,.10);--series-1:#2a78d6;--series-2:#eb6834;
 --good:#0ca30c;--warning:#fab219;--serious:#ec835a;--critical:#d03b3b;--wash:rgba(42,120,214,.07);
 --seq-100:#cde2fb;--seq-200:#9ec5f4;--seq-300:#6da7ec;--seq-400:#3987e5;--seq-500:#256abf;--seq-600:#184f95}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])) .viz-root{color-scheme:dark;--page:#0d0d0d;
 --surface-1:#1a1a19;--text-primary:#f0efec;--text-secondary:#c3c2b7;--grid:#2c2c2a;--axis:#383835;
 --border:rgba(255,255,255,.10);--series-1:#3987e5;--series-2:#d95926;--wash:rgba(57,135,229,.12)}}
:root[data-theme="dark"] .viz-root{color-scheme:dark;--page:#0d0d0d;--surface-1:#1a1a19;--text-primary:#f0efec;
 --text-secondary:#c3c2b7;--grid:#2c2c2a;--axis:#383835;--border:rgba(255,255,255,.10);--series-1:#3987e5;
 --series-2:#d95926;--wash:rgba(57,135,229,.12)}
*{box-sizing:border-box}
body{margin:0}
.viz-root{background:var(--page);color:var(--text-primary);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif;min-height:100vh}
header{display:flex;flex-wrap:wrap;gap:12px;align-items:baseline;justify-content:space-between;padding:20px 28px 8px}
header h1{font-size:20px;margin:0;font-weight:650}
.sub{color:var(--text-secondary);font-size:13px}
nav{display:flex;flex-wrap:wrap;gap:6px;padding:4px 28px 12px;position:sticky;top:0;background:var(--page);z-index:5;border-bottom:1px solid var(--border)}
nav a{color:var(--text-secondary);text-decoration:none;font-size:13px;padding:4px 10px;border-radius:999px}
nav a:hover{background:var(--wash);color:var(--text-primary)}
main{padding:16px 28px 60px;max-width:1400px}
section{margin:22px 0 30px}
h2{font-size:16px;margin:0 0 4px;font-weight:650}
h3{font-size:14px;margin:16px 0 6px}
.note{color:var(--text-secondary);font-size:12.5px;margin:0 0 10px;max-width:980px}
.card{background:var(--surface-1);border:1px solid var(--border);border-radius:10px;padding:14px 16px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px}
.tile .k{color:var(--text-secondary);font-size:12px}
.tile .v{font-size:26px;font-weight:650;margin-top:2px}
.tile .s{color:var(--muted);font-size:12px}
.meter{height:6px;background:var(--grid);border-radius:3px;margin-top:8px;overflow:hidden}
.meter i{display:block;height:100%;background:var(--series-1);border-radius:3px}
.badge{display:inline-flex;align-items:center;gap:6px;font-weight:600;font-size:12.5px;white-space:nowrap}
.badge .ic{width:16px;height:16px;border-radius:50%;display:inline-grid;place-items:center;color:#fff;font-size:11px;line-height:1}
.badge.good .ic{background:var(--good)}.badge.warning .ic{background:var(--warning);color:#0b0b0b}
.badge.critical .ic{background:var(--critical)}.badge.none .ic{background:var(--muted)}
.tw{overflow-x:auto}
table{border-collapse:collapse;width:100%;font-size:12.5px}
th{text-align:left;color:var(--text-secondary);font-weight:600;border-bottom:1px solid var(--axis);padding:6px 8px;white-space:nowrap;position:sticky;top:0;background:var(--surface-1)}
td{border-bottom:1px solid var(--grid);padding:5px 8px;vertical-align:top}
td.n{font-variant-numeric:tabular-nums;text-align:right;white-space:nowrap}
tr.click{cursor:pointer}tr.click:hover td{background:var(--wash)}
.scroll{max-height:520px;overflow:auto}
.filters{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:0 0 10px}
.filters label{font-size:12px;color:var(--text-secondary);display:flex;gap:6px;align-items:center}
select,input[type=search]{font:inherit;font-size:12.5px;padding:4px 8px;border-radius:6px;border:1px solid var(--axis);background:var(--surface-1);color:var(--text-primary)}
button{font:inherit;font-size:12.5px;padding:4px 10px;border-radius:6px;border:1px solid var(--axis);background:var(--surface-1);color:var(--text-primary);cursor:pointer}
.legend{display:flex;gap:16px;font-size:12px;color:var(--text-secondary);margin:2px 0 6px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:10px;height:10px;border-radius:50%;display:inline-block}
svg text{fill:var(--text-secondary);font-size:11px;font-family:inherit}
svg .lbl{fill:var(--text-primary)}
#tip{position:fixed;pointer-events:none;background:var(--surface-1);color:var(--text-primary);border:1px solid var(--border);
 box-shadow:0 4px 16px rgba(0,0,0,.15);border-radius:8px;padding:8px 10px;font-size:12px;max-width:360px;display:none;z-index:50}
#tip b{font-weight:650}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:16px}
@media (max-width:900px){.grid2{grid-template-columns:1fr}}
.md h2,.md h3,.md h4{font-size:14px;margin:14px 0 6px}.md p,.md li{max-width:1000px}
.md code,code{font-family:ui-monospace,Consolas,monospace;font-size:12px;background:var(--wash);padding:0 4px;border-radius:4px}
details>summary{cursor:pointer;color:var(--text-secondary);font-size:13px;margin:4px 0}
.timeline{list-style:none;margin:0;padding:0}
.timeline li{border-left:2px solid var(--grid);padding:4px 0 8px 12px;margin-left:4px;font-size:12.5px}
.timeline .ts{color:var(--muted);font-variant-numeric:tabular-nums;display:block;font-size:11.5px}
.gallery{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:12px}
.gallery figure{margin:0}.gallery img{width:100%;border-radius:8px;border:1px solid var(--border);background:#fff}
.gallery figcaption{font-size:12px;color:var(--text-secondary);margin-top:4px}
.pill{display:inline-block;font-size:11px;padding:0 6px;border-radius:999px;border:1px solid var(--axis);color:var(--text-secondary)}
.seqtxt-dark{color:#fff}
a{color:var(--series-1)}
</style></head>
<body><div class="viz-root">
<header><div><h1>Phantom Claims — claim-level governance for agentic spatial biology</h1>
<div class="sub" id="meta"></div></div>
<div><button id="theme" aria-label="Toggle light/dark theme">Theme: auto</button></div></header>
<nav><a href="#overview">Overview</a><a href="#verdicts">Verdicts</a><a href="#ablations">Gate ablations</a>
<a href="#planted">Planted pass rates</a><a href="#agents">Agent runs</a><a href="#experiments">All experiments</a>
<a href="#timeline">Timeline</a><a href="#report">Morning report</a><a href="#figures">Figures</a></nav>
<main>
<section id="overview"><h2>Overview</h2><p class="note">Pre-registered hypotheses with verdicts computed by the
pre-registered rules (PREREG.md). All numbers on this page are read from <code>results/</code>.</p>
<div class="tiles" id="tiles"></div></section>

<section id="verdicts"><h2>Verdicts by split and condition (K = 10)</h2>
<p class="note">Phantom rate = fraction of reported top-10 claims that are not planted-true. Blue: as reported (95% CI);
orange: after the claim gate. H1 needs ≥ 0.30 with CI-low ≥ 0.20; H2 needs ≥ 50% relative reduction with ≥ 80% retention.</p>
<div class="card"><div class="legend"><span><i class="sw" style="background:var(--series-1)"></i>reported (no gate), 95% CI</span>
<span><i class="sw" style="background:var(--series-2)"></i>after claim gate</span></div><div id="dumbbell"></div></div>
<h3>Table view</h3><div class="card tw" id="verdictTable"></div></section>

<section id="ablations"><h2>Hard negatives: gate ablations (exploratory, dev)</h2>
<p class="note">Per dev seed, a list of all 10 TRUE + all 10 NICHE-confounded planted claims (phantom rate 0.50 before
gating). Bars: phantom rate after gating with each component removed (95% CI). Lower is better.</p>
<div class="grid2"><div class="card"><div id="hardBars"></div></div><div class="card tw" id="hardTable"></div></div></section>

<section id="planted"><h2>Gate pass rates on all planted claims</h2>
<p class="note">Share of planted claims that pass the gate, by category. TRUE should be high; PHANTOM, NICHE and RANDOM
should be low. Shading = value (one hue, light → dark).</p><div class="card tw" id="plantedTable"></div></section>

<section id="agents"><h2>LLM agent runs</h2>
<p class="note">Every agent run with its model, condition, token use, tool calls and scored claims. Click a row for its
tool-call sequence and claims; full transcripts are linked.</p>
<div class="filters" id="agentFilters"></div>
<div class="card"><div class="legend" id="tokLegend"></div><div id="tokChart"></div></div>
<h3>Runs</h3><div class="card scroll" id="agentTable"></div></section>

<section id="experiments"><h2>All experiments (results.jsonl)</h2>
<p class="note" id="expNote"></p>
<div class="filters"><label>Search <input type="search" id="expSearch" placeholder="id, tissue, condition…"></label>
<label>Split <select id="expSplit"><option value="">all</option></select></label></div>
<div class="card scroll" id="expTable"></div></section>

<section id="timeline"><h2>Timeline</h2><div class="grid2">
<div class="card"><h3 style="margin-top:0">Progress log (PROGRESS.md)</h3><div class="scroll"><ul class="timeline" id="progress"></ul></div></div>
<div class="card"><h3 style="margin-top:0">Commits</h3><div class="scroll"><ul class="timeline" id="commits"></ul></div></div></div></section>

<section id="report"><h2>Morning report</h2><div class="card md" id="morning"></div>
<details><summary>STATE.md</summary><div class="card md" id="state"></div></details>
<details><summary>PREREG amendments</summary><div class="card md" id="amend"></div></details></section>

<section id="figures"><h2>Figures</h2><div class="gallery" id="gallery"></div></section>
</main><div id="tip" role="tooltip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById('data').textContent);
const PRE = "__PREFIX__";
const $ = s => document.querySelector(s);
const esc = s => String(s ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const f2 = v => (v === null || v === undefined || Number.isNaN(v)) ? '—' : Number(v).toFixed(2);
const fint = v => (v ?? 0).toLocaleString('en-US');
const NS = 'http://www.w3.org/2000/svg';
function el(tag, attrs, parent){ const e = document.createElementNS(NS, tag); for (const k in attrs) e.setAttribute(k, attrs[k]); if (parent) parent.appendChild(e); return e; }
const css = v => getComputedStyle($('.viz-root')).getPropertyValue(v).trim();

// ---------- theme toggle (auto -> light -> dark)
(function(){ const modes=['auto','light','dark']; const b=$('#theme');
  const q = new URLSearchParams(location.search).get('theme'); let i = Math.max(0, modes.indexOf(q));
  if (i) { document.documentElement.setAttribute('data-theme', modes[i]); b.textContent = 'Theme: ' + modes[i]; }
  b.onclick=()=>{ i=(i+1)%3; const m=modes[i]; if(m==='auto') document.documentElement.removeAttribute('data-theme');
    else document.documentElement.setAttribute('data-theme',m); b.textContent='Theme: '+m; renderCharts(); }; })();

// ---------- tooltip
const tip = $('#tip');
function showTip(ev, html){ tip.innerHTML = html; tip.style.display='block';
  const x = Math.min(ev.clientX + 14, window.innerWidth - tip.offsetWidth - 8); const y = Math.min(ev.clientY + 14, window.innerHeight - tip.offsetHeight - 8);
  tip.style.left = x+'px'; tip.style.top = y+'px'; }
function hideTip(){ tip.style.display='none'; }
function bindTip(node, html){ node.addEventListener('mousemove', e => showTip(e, html)); node.addEventListener('mouseleave', hideTip); }

// ---------- status badge (icon + label; never colour alone)
function badge(v){ v = v || 'NOT YET RUN';
  const cls = v==='SUPPORTED' ? 'good' : v==='INCONCLUSIVE' ? 'warning' : (v==='FALSIFIED'||v==='NOT SUPPORTED') ? 'critical' : 'none';
  const ic = cls==='good' ? '✓' : cls==='warning' ? '~' : cls==='critical' ? '✕' : '○';
  return `<span class="badge ${cls}"><span class="ic" aria-hidden="true">${ic}</span>${esc(v)}</span>`; }

// ---------- meta + tiles
$('#meta').textContent = `Generated ${D.generated} · commit ${D.head || 'n/a'} · ${D.n_result_lines} result lines · ${D.agents.length} agent runs`;
function tiles(){
  const V = D.verdicts || {}; const T = V.table || [];
  const worst = (h, filt) => { const xs = T.filter(filt).map(r => r[h]); if (!xs.length) return 'NOT YET RUN';
    for (const k of ['FALSIFIED','NOT SUPPORTED','INCONCLUSIVE','SUPPORTED']) if (xs.includes(k)) return k; return xs[0]; };
  const H1t = worst('H1', r => r.split==='test' && r.condition.startsWith('expr'));
  const H1a = (() => { const r = D.results.find(x => x.id === 'agent_test_all_both_unguided_K10'); if (!r) return 'NOT YET RUN';
    const pr = r.metrics.phantom_rate, lo = (r.ci.phantom_rate||[0])[0]; return pr >= 0.30 && lo >= 0.20 ? 'SUPPORTED' : pr < 0.10 ? 'FALSIFIED' : 'INCONCLUSIVE'; })();
  const H2 = worst('H2', r => r.split==='test');
  const tile = (k, v, s, extra='') => `<div class="card tile"><div class="k">${k}</div><div class="v">${v}</div><div class="s">${s}</div>${extra}</div>`;
  const b = D.budget || {}; const cap = D.budget_cap.tokens;
  const meter = fam => { const x = b[fam] || {prompt:0,completion:0,usd_est:0}; const t = x.prompt + x.completion;
    return tile(`Tokens · ${fam}`, `${(t/1e6).toFixed(2)}M`, `of ${(cap/1e6).toFixed(0)}M cap · est. $${(x.usd_est||0).toFixed(2)} · ${x.runs||0} runs`,
      `<div class="meter" role="img" aria-label="${Math.round(100*t/cap)}% of token cap"><i style="width:${Math.min(100,100*t/cap)}%"></i></div>`); };
  $('#tiles').innerHTML =
    tile('H1 · expression tools (test)', badge(H1t), 'phantoms ≥ 30% of top-10?') +
    tile('H1 · LLM agents (test, pooled)', badge(H1a), 'phantoms ≥ 30% of top-10?') +
    tile('H2 · gate works (test)', badge(H2), '≥ 50% reduction, ≥ 80% retention') +
    tile('H3 · transfer', badge(V.H3), 'held-out breast + pancreas + lung') +
    tile('Experiments logged', fint(D.results.length), `${fint(D.n_result_lines)} lines in results.jsonl`) +
    tile('Agent runs', fint(D.agents.length), `${D.agents.filter(a=>a.status==='submitted').length} submitted`) +
    meter('gemini') + meter('deepseek');
}

// ---------- verdict dumbbell + table
function dumbbell(){
  const host = $('#dumbbell'); host.innerHTML = '';
  const rows = ((D.verdicts||{}).table || []).filter(r => r.pr !== null && r.pr !== undefined);
  const W = Math.max(640, host.clientWidth || 900), L = 230, R = 24, rh = 22, top = 26, H = top + rows.length*rh + 30;
  const svg = el('svg', {width: W, height: H, role:'img', 'aria-label':'Phantom rate before and after the gate by split and condition'}, host);
  const x = v => L + v*(W-L-R);
  for (const t of [0, .25, .5, .75, 1]) { el('line',{x1:x(t),x2:x(t),y1:top-6,y2:H-24,stroke:css('--grid'),'stroke-width':1},svg);
    const tx = el('text',{x:x(t),y:H-8,'text-anchor':'middle'},svg); tx.textContent = t.toFixed(2); }
  el('line',{x1:x(.30),x2:x(.30),y1:top-14,y2:H-24,stroke:css('--muted'),'stroke-width':1.5},svg);
  const th = el('text',{x:x(.30)+4,y:top-8},svg); th.textContent = 'H1 threshold 0.30';
  rows.forEach((r,i) => { const y = top + i*rh + rh/2;
    const lab = el('text',{x:L-10,y:y+4,'text-anchor':'end',class:'lbl'},svg);
    lab.textContent = `${r.split} · ${r.tissue.replace('breast_','breast ')} · ${r.condition.replace('expr_','tool: ').replace('agent_','LLM: ')}`;
    const ci = Array.isArray(r.pr_ci) ? r.pr_ci : null;
    if (ci) el('line',{x1:x(ci[0]),x2:x(ci[1]),y1:y,y2:y,stroke:css('--series-1'),'stroke-width':2,'stroke-linecap':'round',opacity:.45},svg);
    if (r.pr_after !== null && r.pr_after !== undefined) el('line',{x1:x(r.pr_after),x2:x(r.pr),y1:y,y2:y,stroke:css('--axis'),'stroke-width':1},svg);
    el('circle',{cx:x(r.pr),cy:y,r:5,fill:css('--series-1'),stroke:css('--surface-1'),'stroke-width':2},svg);
    if (r.pr_after !== null && r.pr_after !== undefined) el('circle',{cx:x(r.pr_after),cy:y,r:5,fill:css('--series-2'),stroke:css('--surface-1'),'stroke-width':2},svg);
    const hit = el('rect',{x:0,y:y-rh/2,width:W,height:rh,fill:'transparent'},svg);
    bindTip(hit, `<b>${esc(r.split)} · ${esc(r.tissue)} · ${esc(r.condition)}</b><br>reported ${f2(r.pr)} ${ci?`[${f2(ci[0])}, ${f2(ci[1])}]`:''}
      <br>after gate ${f2(r.pr_after)} · retention ${f2(r.ret)}<br>H1 ${esc(r.H1)} · H2 ${esc(r.H2)}`); });
}
function verdictTable(){
  const rows = (D.verdicts||{}).table || [];
  let h = '<table><thead><tr><th>split</th><th>tissue</th><th>condition</th><th>phantom@10</th><th>95% CI</th><th>H1</th><th>after gate</th><th>retention</th><th>H2</th><th>runs</th></tr></thead><tbody>';
  for (const r of rows) { if (r.pr === null || r.pr === undefined) continue;
    h += `<tr><td>${esc(r.split)}</td><td>${esc(r.tissue)}</td><td>${esc(r.condition)}</td><td class="n">${f2(r.pr)}</td>
      <td class="n">${Array.isArray(r.pr_ci)?`[${f2(r.pr_ci[0])}, ${f2(r.pr_ci[1])}]`:''}</td><td>${badge(r.H1)}</td><td class="n">${f2(r.pr_after)}</td>
      <td class="n">${f2(r.ret)}</td><td>${badge(r.H2)}</td><td class="n">${r.n_runs ?? ''}</td></tr>`; }
  $('#verdictTable').innerHTML = h + '</tbody></table>' + `<p class="note" style="margin-top:8px">H3: ${badge((D.verdicts||{}).H3)}</p>`;
}

// ---------- hard-list ablation bars
function hardBars(){
  const host = $('#hardBars'); host.innerHTML='';
  const order = [['full','full gate'],['no_split','− split consistency'],['no_perm','− permutation test'],['no_tiles','− tile restriction'],['no_min_evidence','− minimum evidence']];
  const rows = order.filter(([k]) => D.hard[k]).map(([k,l]) => ({k, l, v: D.hard[k].phantom_rate_after, ci: (D.hard[k].ci||{}).phantom_rate_after, ret: D.hard[k].retention}));
  if (!rows.length) { host.innerHTML = '<p class="note">Not yet run.</p>'; return; }
  const W = Math.max(420, host.clientWidth || 520), L = 150, R = 50, rh = 30, top = 10, H = top + rows.length*rh + 28;
  const svg = el('svg',{width:W,height:H,role:'img','aria-label':'Phantom rate after gating by ablation'},host);
  const x = v => L + v/0.5*(W-L-R);
  for (const t of [0,.1,.2,.3,.4,.5]) { el('line',{x1:x(t),x2:x(t),y1:top,y2:H-22,stroke:css('--grid'),'stroke-width':1},svg);
    const tx = el('text',{x:x(t),y:H-6,'text-anchor':'middle'},svg); tx.textContent = t.toFixed(1); }
  rows.forEach((r,i) => { const y = top + i*rh; const bh = 14, yc = y + rh/2;
    const lab = el('text',{x:L-10,y:yc+4,'text-anchor':'end',class:'lbl'},svg); lab.textContent = r.l;
    const w = Math.max(0, x(r.v)-x(0));
    el('path',{d:`M${x(0)},${yc-bh/2} h${Math.max(0,w-4)} a4,4 0 0 1 4,4 v${bh-8} a4,4 0 0 1 -4,4 h-${Math.max(0,w-4)} z`,fill:css('--series-1')},svg);
    if (r.ci) el('line',{x1:x(r.ci[0]),x2:x(r.ci[1]),y1:yc,y2:yc,stroke:css('--text-secondary'),'stroke-width':1.5},svg);
    const vt = el('text',{x:x(Math.max(r.v, r.ci?r.ci[1]:r.v))+6,y:yc+4,class:'lbl'},svg); vt.textContent = f2(r.v);
    const hit = el('rect',{x:0,y:y,width:W,height:rh,fill:'transparent'},svg);
    bindTip(hit, `<b>${esc(r.l)}</b><br>phantom rate after ${f2(r.v)}${r.ci?` [${f2(r.ci[0])}, ${f2(r.ci[1])}]`:''}<br>retention ${f2(r.ret)}`); });
  let h = '<table><thead><tr><th>variant</th><th>phantom after</th><th>95% CI</th><th>retention</th><th>rel. reduction</th></tr></thead><tbody>';
  for (const r of rows) { const m = D.hard[r.k]; h += `<tr><td>${esc(r.l)}</td><td class="n">${f2(m.phantom_rate_after)}</td><td class="n">${r.ci?`[${f2(r.ci[0])}, ${f2(r.ci[1])}]`:''}</td><td class="n">${f2(m.retention)}</td><td class="n">${f2(m.rel_reduction)}</td></tr>`; }
  $('#hardTable').innerHTML = h + '</tbody></table><p class="note" style="margin-top:8px">Source: results.jsonl id <code>dev_hard_list_true_vs_niche</code>.</p>';
}

// ---------- planted pass-rate table with sequential shading
function plantedTable(){
  const rows = D.planted || []; if (!rows.length) { $('#plantedTable').innerHTML = '<p class="note">Not yet run.</p>'; return; }
  const steps = ['--seq-100','--seq-200','--seq-300','--seq-400','--seq-500','--seq-600'];
  const cell = v => { const x = Number(v); if (Number.isNaN(x)) return '<td class="n">—</td>';
    const i = Math.min(5, Math.floor(x*6)); return `<td class="n ${i>=3?'seqtxt-dark':''}" style="background:var(${steps[i]})">${x.toFixed(2)}</td>`; };
  let h = '<table><thead><tr><th>setting</th><th>true</th><th>phantom</th><th>niche</th><th>random</th></tr></thead><tbody>';
  for (const r of rows) h += `<tr><td>${esc(r.setting)}</td>${cell(r.true)}${cell(r.phantom)}${cell(r.niche)}${cell(r.random)}</tr>`;
  $('#plantedTable').innerHTML = h + '</tbody></table>';
}

// ---------- agent runs: filters, token chart, table
const AF = {split:'', family:'', cond:'', tissue:''};
function agentFilters(){
  const opts = k => [...new Set(D.agents.map(a => a[k]))].sort();
  const mk = (k, label) => `<label>${label} <select data-k="${k}"><option value="">all</option>${opts(k).map(o=>`<option>${esc(o)}</option>`).join('')}</select></label>`;
  $('#agentFilters').innerHTML = mk('split','Split') + mk('family','Model family') + mk('cond','Condition') + mk('tissue','Tissue');
  $('#agentFilters').querySelectorAll('select').forEach(s => s.onchange = () => { AF[s.dataset.k] = s.value; agentChart(); agentTable(); });
}
const famColor = f => f === 'gemini' ? css('--series-1') : css('--series-2');   // colour follows the entity
function filtAgents(){ return D.agents.filter(a => Object.entries(AF).every(([k,v]) => !v || a[k] === v)); }
function agentChart(){
  const host = $('#tokChart'); host.innerHTML = '';
  const fams = [...new Set(D.agents.map(a=>a.family))].sort();
  $('#tokLegend').innerHTML = fams.map(f => `<span><i class="sw" style="background:${famColor(f)}"></i>${esc(f)}</span>`).join('') + '<span>bar = total tokens per run</span>';
  const rows = filtAgents(); if (!rows.length) { host.innerHTML = '<p class="note">No runs match.</p>'; return; }
  const W = Math.max(640, host.clientWidth || 900), L = 56, R = 10, top = 10, H = 200, bot = 22;
  const svg = el('svg',{width:W,height:H,role:'img','aria-label':'Tokens per agent run'},host);
  const raw = Math.max(...rows.map(a => a.prompt + a.completion));
  const nice = s => { const e = Math.pow(10, Math.floor(Math.log10(s))); const f = s/e; return (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * e; };
  const step = nice(raw / 4), mx = Math.ceil(raw / step) * step;
  const y = v => H - bot - v/mx*(H-bot-top);
  const ticks = []; for (let t = 0; t <= mx + 1e-9; t += step) ticks.push(t);
  for (const t of ticks) { el('line',{x1:L,x2:W-R,y1:y(t),y2:y(t),stroke:css('--grid'),'stroke-width':1},svg);
    const tx = el('text',{x:L-6,y:y(t)+4,'text-anchor':'end'},svg); tx.textContent = (t/1000).toFixed(0)+'k'; }
  const bw = (W-L-R)/rows.length; const gap = Math.min(2, bw*0.2);
  rows.forEach((a,i) => { const v = a.prompt + a.completion; const x0 = L + i*bw + gap/2, w = Math.max(1, bw-gap), h = H-bot-y(v);
    const r = Math.min(4, w/2, h);
    el('path',{d:`M${x0},${H-bot} v${-(h-r)} a${r},${r} 0 0 1 ${r},${-r} h${w-2*r} a${r},${r} 0 0 1 ${r},${r} v${h-r} z`,fill:famColor(a.family)},svg);
    const hit = el('rect',{x:L+i*bw,y:top,width:bw,height:H-bot-top,fill:'transparent'},svg);
    bindTip(hit, `<b>${esc(a.id)}</b><br>${esc(a.model)} · ${esc(a.cond)} · ${esc(a.status)}<br>tokens ${fint(v)} (prompt ${fint(a.prompt)}, completion ${fint(a.completion)})
      <br>tool calls ${a.tool_calls ?? '—'} · runtime ${a.runtime ?? '—'} s<br>phantom@10 ${f2(a.phantom)} · gate passes ${a.n_gate_pass}/${a.n_claims}`); });
  el('line',{x1:L,x2:W-R,y1:H-bot,y2:H-bot,stroke:css('--axis'),'stroke-width':1},svg);
  const xl = el('text',{x:L,y:H-6},svg); xl.textContent = `${rows.length} runs (hover a bar)`;
}
function agentTable(){
  const rows = filtAgents();
  let h = '<table><thead><tr><th>run</th><th>split</th><th>tissue</th><th>seed</th><th>model</th><th>condition</th><th>status</th><th>tool calls</th><th>tokens</th><th>phantom@10</th><th>gate passes</th><th>transcript</th></tr></thead><tbody>';
  rows.forEach((a,i) => { h += `<tr class="click" data-i="${i}"><td>${esc(a.id)}</td><td>${esc(a.split)}</td><td>${esc(a.tissue)}</td><td class="n">${a.seed}</td>
    <td>${esc(a.model)}</td><td>${esc(a.cond)}</td><td>${badge(a.status==='submitted'?'SUPPORTED':'FALSIFIED').replace(/SUPPORTED|FALSIFIED/, esc(a.status))}</td>
    <td class="n">${a.tool_calls ?? ''}</td><td class="n">${fint(a.prompt+a.completion)}</td><td class="n">${f2(a.phantom)}</td>
    <td class="n">${a.cond==='unguided' ? `${a.n_gate_pass}/${a.n_claims}` : '—'}</td>
    <td><a href="${PRE}results/agent_logs/${encodeURIComponent(a.id)}.json" target="_blank" rel="noopener" onclick="event.stopPropagation()">json</a></td></tr>
    <tr class="det" data-d="${i}" hidden><td colspan="12"></td></tr>`; });
  $('#agentTable').innerHTML = h + '</tbody></table>';
  $('#agentTable').querySelectorAll('tr.click').forEach(tr => tr.onclick = () => {
    const i = +tr.dataset.i, a = rows[i], d = $('#agentTable').querySelector(`tr[data-d="${i}"]`);
    if (!d.hidden) { d.hidden = true; return; }
    const seq = a.seq.reduce((acc, n) => { const last = acc[acc.length-1]; if (last && last[0]===n) last[1]++; else acc.push([n,1]); return acc; }, [])
      .map(([n,c]) => `<span class="pill">${esc(n)}${c>1?' ×'+c:''}</span>`).join(' → ');
    let c = '<table><thead><tr><th>#</th><th>sender → receiver</th><th>ligand–receptor</th><th>truth</th><th>gate</th><th>agent evidence</th></tr></thead><tbody>';
    for (const x of a.claims) c += `<tr><td class="n">${x.rank}</td><td>${esc(x.sender)} → ${esc(x.receiver)}</td><td>${esc(x.ligand)}–${esc(x.receptor)}</td>
      <td>${x.label===1 ? badge('SUPPORTED').replace('SUPPORTED','true') : badge('FALSIFIED').replace('FALSIFIED', esc(x.neg_type||'negative'))}</td>
      <td>${x.gate==='' ? '—' : (x.gate==='True' ? 'pass' : 'fail')}</td><td>${esc(x.evidence)}</td></tr>`;
    d.firstElementChild.innerHTML = `<div style="padding:6px 0"><b>Tool sequence</b><div style="margin:4px 0 10px;line-height:2">${seq || '—'}</div><b>Claims</b>${c}</tbody></table></div>`;
    d.hidden = false; });
}

// ---------- experiments table
function expTable(){
  const q = ($('#expSearch').value || '').toLowerCase(), sp = $('#expSplit').value;
  const rows = D.results.filter(r => (!sp || r.split === sp) && (!q || JSON.stringify([r.id,r.tissue,r.condition,r.split,r.ablation,r.family]).toLowerCase().includes(q))).slice().reverse();
  const key = ['phantom_rate','phantom_rate_after','retention','rel_reduction','n'];
  let h = '<table><thead><tr><th>id</th><th>phase</th><th>split</th><th>tissue</th><th>condition</th><th>K</th>' + key.map(k=>`<th>${k}</th>`).join('') + '<th>exploratory</th><th>time</th><th>git</th><th>versions</th></tr></thead><tbody>';
  for (const r of rows) h += `<tr class="click" title="click to show the raw record"><td>${esc(r.id)}</td><td class="n">${r.phase ?? ''}</td><td>${esc(r.split)}</td><td>${esc(r.tissue)}</td><td>${esc(r.condition)}</td><td class="n">${r.K ?? ''}</td>`
    + key.map(k => `<td class="n">${r.metrics[k] === undefined ? '' : (k==='n' ? r.metrics[k] : f2(r.metrics[k]))}</td>`).join('')
    + `<td>${r.exploratory ? 'yes' : ''}</td><td class="n">${esc((r.timestamp||'').replace('T',' ').slice(0,16))}</td><td><code>${esc(r.git)}</code></td><td class="n">${r.versions}</td></tr>
       <tr hidden><td colspan="${key.length+10}"><pre style="white-space:pre-wrap;font-size:11px;margin:0">${esc(r.raw)}</pre></td></tr>`;
  $('#expTable').innerHTML = h + '</tbody></table>';
  $('#expTable').querySelectorAll('tr.click').forEach(tr => tr.onclick = () => { tr.nextElementSibling.hidden = !tr.nextElementSibling.hidden; });
  $('#expNote').textContent = `${rows.length} of ${D.results.length} unique experiment ids (latest version of each id; ${D.n_result_lines} lines in total). Click a row for the raw JSON record.`;
}
[...new Set(D.results.map(r => r.split).filter(Boolean))].sort().forEach(s => { const o = document.createElement('option'); o.textContent = s; $('#expSplit').appendChild(o); });
$('#expSearch').oninput = expTable; $('#expSplit').onchange = expTable;

// ---------- timeline, reports, figures
$('#progress').innerHTML = D.progress.slice().reverse().map(e => `<li><span class="ts">${esc(e.ts)}</span>${esc(e.text)}</li>`).join('');
$('#commits').innerHTML = D.commits.map(c => `<li><span class="ts">${esc(c.date)} · <code>${esc(c.hash)}</code></span>${esc(c.msg)}</li>`).join('');
$('#morning').innerHTML = D.morning; $('#state').innerHTML = D.state; $('#amend').innerHTML = D.amendments;
$('#gallery').innerHTML = D.figures.map(f => `<figure><a href="${PRE}figures/${esc(f)}" target="_blank" rel="noopener"><img loading="lazy" src="${PRE}figures/${esc(f)}" alt="${esc(f)}"></a><figcaption>${esc(f)}</figcaption></figure>`).join('');

function renderCharts(){ dumbbell(); hardBars(); agentChart(); }
tiles(); verdictTable(); plantedTable(); agentFilters(); agentTable(); expTable(); renderCharts();
let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(renderCharts, 150); });
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', renderCharts);
</script></div></body></html>
"""


def _clean(o):
    """NaN/Infinity are not valid JSON (browser JSON.parse rejects them) -> null."""
    if isinstance(o, float):
        return None if (o != o or o in (float("inf"), float("-inf"))) else o
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    return o


def build(out_html, prefix):
    data = _clean(collect())
    payload = json.dumps(data, default=str, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    page = PAGE.replace("__DATA__", payload).replace("__PREFIX__", prefix)
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(page)
    return data


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", help="write a self-contained site (index.html + copied figures and agent logs) to this dir")
    a = ap.parse_args()
    if a.site:
        site = os.path.abspath(a.site)
        d = build(os.path.join(site, "index.html"), "")
        shutil.copytree(p("figures"), os.path.join(site, "figures"), dirs_exist_ok=True)
        shutil.copytree(p("results", "agent_logs"), os.path.join(site, "results", "agent_logs"), dirs_exist_ok=True)
        open(os.path.join(site, ".nojekyll"), "w").close()
    else:
        d = build(p("dashboard", "index.html"), "../")
    print(f"dashboard built: {len(d['results'])} experiments, {len(d['agents'])} agent runs, "
          f"{len(d['progress'])} progress entries, {len(d['commits'])} commits")
