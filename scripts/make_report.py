"""Build verdict tables + figures from results/results.jsonl ONLY (latest record per id), and write
results/tables/verdicts.json + paper tables. Every number shown comes from results.jsonl."""
import json, math, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = {}
for line in open("results/results.jsonl"):
    line = line.strip()
    if line:
        r = json.loads(line); R[r["id"]] = r


def get(i):
    return R.get(i)


def fmt(x, nd=2):
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}f}"


def h1_verdict(r):
    if r is None:
        return "NOT YET RUN"
    p, lo = round(r["metrics"]["phantom_rate"], 9), round(r["ci"]["phantom_rate"][0], 9)   # avoid 6/60 -> 0.0999.. artefact
    if p >= 0.30 and lo >= 0.20:
        return "SUPPORTED"
    if p < 0.10:
        return "FALSIFIED"
    return "INCONCLUSIVE"


def h2_verdict(r):
    if r is None:
        return "NOT YET RUN"
    m = r["metrics"]; red, ret = m.get("rel_reduction"), m.get("retention")
    nan = lambda v: v is None or (isinstance(v, float) and math.isnan(v))
    if not nan(red) and round(red, 9) < 0.25:
        return "FALSIFIED"                       # PREREG: reduction < 0.25 falsifies regardless of retention
    if not nan(ret) and round(ret, 9) < 0.70:
        return "FALSIFIED"
    if nan(red) or nan(ret):
        return "INCONCLUSIVE"
    red, ret = round(red, 9), round(ret, 9)
    if red >= 0.5 and ret >= 0.8:
        return "SUPPORTED"
    if ret < 0.7 or red < 0.25:
        return "FALSIFIED"
    return "INCONCLUSIVE"


def h1_id(split, tissue, m, K=10):
    a = f"{split}_{'default' if split == 'dev' else 'final'}_{tissue}_H1_expr_{m}_K{K}"
    b = f"{split}_default_H1_expr_{m}_K{K}"
    return a if a in R else (b if split == "dev" and b in R else a)


def h2_id(split, tissue, m, K=10, abl="none"):
    a = f"{split}_{'default' if split == 'dev' else 'final'}_{tissue}_H2_expr_{m}_K{K}_{abl}"
    b = f"{split}_default_H2_expr_{m}_K{K}_{abl}"
    return a if a in R else (b if split == "dev" and b in R else a)


rows = []
sets = [("dev", "breast_rep1"), ("test", "breast_rep1"), ("test", "breast_rep2"), ("test", "pancreas"), ("test", "lung")]
for split, tissue in sets:
    for m in ("cpdb", "prod"):
        r1, r2 = get(h1_id(split, tissue, m)), get(h2_id(split, tissue, m))
        rows.append(dict(split=split, tissue=tissue, condition=f"expr_{m}",
                         pr=r1 and r1["metrics"]["phantom_rate"], pr_ci=r1 and r1["ci"]["phantom_rate"], H1=h1_verdict(r1),
                         pr_after=r2 and r2["metrics"].get("phantom_rate_after"),
                         red=r2 and r2["metrics"].get("rel_reduction"), red_ci=r2 and r2["ci"].get("rel_reduction"),
                         ret=r2 and r2["metrics"].get("retention"), ret_ci=r2 and r2["ci"].get("retention"),
                         H2=h2_verdict(r2), ids=[x for x in (h1_id(split, tissue, m), h2_id(split, tissue, m)) if x in R]))
    for fam in ("gemini", "deepseek", "both"):
        a1 = get(f"agent_{split}_{tissue}_{fam}_unguided_K10"); a2 = get(f"agent_{split}_{tissue}_{fam}_unguided_gated_K10")
        a3 = get(f"agent_{split}_{tissue}_{fam}_gate_tool_K10")
        if a1 is None:
            continue
        rows.append(dict(split=split, tissue=tissue, condition=f"agent_{fam}",
                         pr=a1["metrics"]["phantom_rate"], pr_ci=a1["ci"]["phantom_rate"], H1=h1_verdict(a1),
                         pr_after=a2 and a2["metrics"].get("phantom_rate_after"),
                         red=a2 and a2["metrics"].get("rel_reduction"), red_ci=a2 and a2["ci"].get("rel_reduction"),
                         ret=a2 and a2["metrics"].get("retention"), ret_ci=a2 and a2["ci"].get("retention"),
                         H2=h2_verdict(a2), gate_tool_pr=a3 and a3["metrics"]["phantom_rate"],
                         gate_tool_ci=a3 and a3["ci"]["phantom_rate"], n_runs=a1.get("n_runs"),
                         ids=[x["id"] for x in (a1, a2, a3) if x]))
T = pd.DataFrame(rows)
os.makedirs("paper/tables", exist_ok=True)
T.to_csv("results/tables/verdict_table.csv", index=False)

# H3: frozen gate on held-out breast seeds AND >= 2 other tissues (both expression tools meet H2 criteria)
def h2_ok(tissue):
    v = [h2_verdict(get(h2_id("test", tissue, m))) for m in ("cpdb", "prod")]
    return None if "NOT YET RUN" in v else all(x == "SUPPORTED" for x in v)
others = {t: h2_ok(t) for t in ("pancreas", "lung")}
breast = h2_ok("breast_rep1")
if breast is None or any(v is None for v in others.values()):
    H3 = "NOT YET RUN"
else:
    H3 = "SUPPORTED" if breast and sum(others.values()) >= 2 else "NOT SUPPORTED"
verd = dict(H3=H3, H3_detail=dict(breast_rep1_heldout=breast, rep2=h2_ok("breast_rep2"), **others))
json.dump(dict(table=T.to_dict("records"), **verd), open("results/tables/verdicts.json", "w"), indent=1, default=str)

# ---- Figure 1: phantom rate before/after gate, K=10
fig, ax = plt.subplots(figsize=(10, 3.6))
d = T[T.pr.notna()].reset_index(drop=True)
x = np.arange(len(d))
lab = [f"{r.split}:{r.tissue.replace('breast_', 'br_')}\n{r.condition.replace('expr_', '').replace('agent_', 'LLM-')}" for r in d.itertuples()]
pr = d.pr.astype(float).values
pa = d.pr_after.astype(float).values
err = np.array([[p - c[0], c[1] - p] if isinstance(c, list) else [0, 0] for p, c in zip(pr, d.pr_ci)]).T
ax.bar(x - 0.2, pr, 0.4, yerr=err, capsize=2, color="#c44e52", label="reported (no gate)")
ok_ = ~np.isnan(pa)                                   # review: do not draw missing values as 0
ax.bar((x + 0.2)[ok_], pa[ok_], 0.4, color="#4c72b0", label="after claim gate")
ax.scatter((x + 0.2)[ok_], pa[ok_] + 0.012, marker="v", s=14, color="#4c72b0", zorder=3)
for xi, v in zip(x, pa):
    if not np.isnan(v):
        ax.text(xi + 0.2, v + 0.03, f"{v:.2f}", ha="center", fontsize=5, color="#4c72b0")
ax.axhline(0.30, ls="--", lw=0.8, color="k"); ax.text(len(d) - 0.5, 0.31, "H1 threshold 0.30", ha="right", fontsize=7)
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=6, rotation=90); ax.set_ylabel("phantom rate @10"); ax.set_ylim(0, 1)
ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig("figures/fig_phantom_rates.png", dpi=200); plt.close(fig)

# ---- Figure 2: planted-claim pass rates by category (dev ablations + test tissues)
pp = []
dv = get("dev_gate_planted_pass_rates")
if dv:
    for v, cats in dv["metrics"].items():
        pp.append(dict(setting=f"dev:{v}", **cats))
for t in ("breast_rep1", "breast_rep2", "pancreas", "lung"):
    r = get(f"test_final_planted_pass_rates_{t}")
    if r:
        pp.append(dict(setting=f"test:{t}", **r["metrics"]))
if pp:
    P = pd.DataFrame(pp).set_index("setting")[["true", "phantom", "niche", "random"]]
    P.to_csv("results/tables/planted_pass_rates.csv")
    fig, ax = plt.subplots(figsize=(8, 3.2))
    P.plot.bar(ax=ax, color=["#55a868", "#c44e52", "#dd8452", "#8172b3"], width=0.8)
    ax.set_ylabel("gate pass rate (planted claims)"); ax.set_ylim(0, 1.05)
    ax.legend(["true", "phantom (all 0)" if (P.phantom == 0).all() else "phantom", "niche",
               "random (all 0)" if (P.random == 0).all() else "random"],
              fontsize=7, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right", fontsize=7); fig.tight_layout()
    fig.savefig("figures/fig_planted_pass_rates.png", dpi=200); plt.close(fig)

# ---- Figure 3: exploratory knob sweep on dev
sw = []
for tag in ("default", "effect0.5", "effect0.25", "noise0.5", "r50"):
    for m in ("cpdb", "prod"):
        r1 = get(f"dev_{tag}_breast_rep1_H1_expr_{m}_K10") or get(f"dev_{tag}_H1_expr_{m}_K10")
        r2 = get(f"dev_{tag}_breast_rep1_H2_expr_{m}_K10_none") or get(f"dev_{tag}_H2_expr_{m}_K10_none")
        if r1 and r2:
            sw.append(dict(bench=tag, method=m, pr=r1["metrics"]["phantom_rate"], pr_after=r2["metrics"]["phantom_rate_after"],
                           retention=r2["metrics"]["retention"]))
if sw:
    S = pd.DataFrame(sw); S.to_csv("results/tables/dev_sweep.csv", index=False)
print(T[["split", "tissue", "condition", "pr", "H1", "pr_after", "red", "ret", "H2"]].to_string())
print("H3:", verd)

# ---- LaTeX verdict table for the paper (numbers from results.jsonl via T)
def f2(x):
    return "--" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.2f}"
def fci(c):
    return "" if not isinstance(c, list) else f" [{c[0]:.2f}, {c[1]:.2f}]"
BS = chr(92)  # backslash, avoids escape mangling
def esc(s):
    return s.replace("_", BS + "_")
lines = [BS + "begin{table*}[t]" + BS + "centering" + BS + "small",
         BS + "begin{tabular}{llllllll}" + BS + "toprule",
         "split & tissue & condition & phantom@10 & H1 & after gate & retention & H2" + BS + BS + BS + "midrule"]
for r in T.itertuples():
    if r.pr is None or (isinstance(r.pr, float) and math.isnan(r.pr)):
        continue
    lines.append(f"{r.split} & {esc(r.tissue)} & {esc(r.condition)} & {f2(r.pr)}{fci(r.pr_ci)} & "
                 f"{r.H1} & {f2(r.pr_after)} & {f2(r.ret)}{fci(r.ret_ci)} & {r.H2}" + BS + BS)
lines += [BS + "bottomrule" + BS + "end{tabular}",
          BS + "caption{Pre-registered verdicts (K=10; 95" + BS + "% hierarchical bootstrap CIs). H3: " + verd["H3"]
          + ".}" + BS + "label{tab:verdicts}",
          BS + "end{table*}"]
open("paper/tables/verdicts.tex", "w").write("\n".join(lines) + "\n")
