"""Score agent runs: label claims (first K unique), post-hoc gate unguided claims with the frozen gate,
aggregate with hierarchical bootstrap, log to results.jsonl. Usage: python scripts/agent_score.py <split>"""
import glob, json, os, re, sys
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.generate import get_or_make, label_claim
from bench.pipeline import make_gate, summarize, log_result
from gate.gate import GateConfig, write_audit

K = 10
split = sys.argv[1]
sel = json.load(open("results/tables/gate_selected_dev.json"))
cfg = GateConfig(**sel)
pat = re.compile(rf"^{split}_(.+)_s(\d+)_(gemini|deepseek)_(unguided|gate_tool)$")
runs = []
for f in sorted(glob.glob(f"results/agent_logs/{split}_*.json")):
    rid = os.path.basename(f)[:-5]
    m = pat.match(rid)
    if not m:
        continue
    l = json.load(open(f))
    runs.append(dict(run_id=rid, tissue=m.group(1), seed=int(m.group(2)), family=m.group(3), cond=m.group(4), log=l))
tissues = {}
rows = []
for r in runs:
    if r["tissue"] not in tissues:
        tissues[r["tissue"]] = load_tissue(r["tissue"])
    t = tissues[r["tissue"]]
    inst = get_or_make(t, r["seed"])
    claims, seen = [], set()
    for c in (r["log"].get("claims") or []):
        if not isinstance(c, dict):
            continue
        k = (str(c.get("sender")), str(c.get("receiver")), str(c.get("ligand")), str(c.get("receptor")))
        if k in seen:
            continue
        seen.add(k)
        claims.append(dict(sender=k[0], receiver=k[1], ligand=k[2], receptor=k[3], evidence=str(c.get("evidence", ""))))
    claims = claims[:K]
    labs = [label_claim(inst, c) for c in claims]
    passed = [None] * len(claims)
    if r["cond"] == "unguided" and claims:
        g = make_gate(t, inst, cfg)
        recs = g.check([{k: c[k] for k in ("sender", "receiver", "ligand", "receptor")} for c in claims])
        for rr, (y, sb) in zip(recs, labs):
            rr["truth"] = dict(label=y, neg_type=sb, seed=r["seed"])
        write_audit(recs, f"results/audit/agent_{r['run_id']}_posthoc.jsonl")
        passed = [rr["pass_"] for rr in recs]
    for i, (c, (y, sb), p) in enumerate(zip(claims, labs, passed)):
        rows.append(dict(run_id=r["run_id"], tissue=r["tissue"], seed=r["seed"], family=r["family"], cond=r["cond"],
                         rank=i + 1, **c, label=y, neg_type=sb, gate_pass=p, status=r["log"].get("status"),
                         n_tool_calls=r["log"].get("n_tool_calls"), tokens=sum(r["log"].get("usage", {}).values())))
    if not claims:
        rows.append(dict(run_id=r["run_id"], tissue=r["tissue"], seed=r["seed"], family=r["family"], cond=r["cond"],
                         rank=0, label=None, status=r["log"].get("status")))
df = pd.DataFrame(rows)
df.to_csv(f"results/claims/agent_{split}_claims.csv", index=False)
d = df[df["rank"] > 0]
# tissue groups: each tissue separately, plus pooled 'all'
for tissue in sorted(d.tissue.unique()) + (["all"] if d.tissue.nunique() > 1 else []):
    dd = d if tissue == "all" else d[d.tissue == tissue]
    for fam in sorted(dd.family.unique()) + ["both"]:
        df_ = dd if fam == "both" else dd[dd.family == fam]
        for cond in ("unguided", "gate_tool"):
            x = df_[df_.cond == cond]
            if x.empty:
                continue
            groups = [(g.label.astype(int).values, None) for _, g in x.groupby(["tissue", "seed", "family"])]
            point, ci = summarize(groups)
            cats = x.neg_type.value_counts().to_dict()
            log_result(dict(id=f"agent_{split}_{tissue}_{fam}_{cond}_K{K}", phase=5 if split == "dev" else 6, split=split,
                            tissue=tissue, family=fam, condition=f"agent_{cond}", K=K, hypothesis="H1",
                            n_runs=len(groups), metrics=point, ci=ci, claim_breakdown=cats))
            print(f"{split} {tissue} {fam} {cond}: n_runs={len(groups)} phantom={point['phantom_rate']:.3f} CI={np.round(ci['phantom_rate'], 3)} {cats}")
            if cond == "unguided":
                groups = [(g.label.astype(int).values, g.gate_pass.astype(bool).values) for _, g in x.groupby(["tissue", "seed", "family"])]
                point, ci = summarize(groups)
                log_result(dict(id=f"agent_{split}_{tissue}_{fam}_unguided_gated_K{K}", phase=5 if split == "dev" else 6,
                                split=split, tissue=tissue, family=fam, condition="agent_unguided_gated", K=K,
                                hypothesis="H2", n_runs=len(groups), gate_config=sel, metrics=point, ci=ci))
                print(f"   gated: PR {point['phantom_rate']:.3f}->{point['phantom_rate_after']:.3f} red={point['rel_reduction']:.3f} "
                      f"CI={np.round(ci['rel_reduction'], 3)} ret={point['retention']:.3f} CI={np.round(ci['retention'], 3)}")

# ---- gate_tool compliance: submitted claims vs. the agent's own verify_claim results in its transcript
comp = dict(runs=0, submitted=0, verified_pass=0, verified_fail_submitted=0, unverified=0, n_verify_calls=0,
            verified_fail_submitted_negatives=0)
for r in runs:
    if r["cond"] != "gate_tool":
        continue
    ver = {}
    for s in r["log"]["steps"]:
        for tc in s.get("tool_calls", []):
            if tc["name"] == "verify_claim":
                a = tc["args"]; k = (a.get("sender"), a.get("receiver"), a.get("ligand"), a.get("receptor"))
                try:
                    ver[k] = json.loads(tc["result"]).get("pass_")
                except Exception:
                    pass
    comp["runs"] += 1; comp["n_verify_calls"] += len(ver)
    sub = d[(d.run_id == r["run_id"])]
    for c in sub.itertuples():
        k = (c.sender, c.receiver, c.ligand, c.receptor)
        comp["submitted"] += 1
        if ver.get(k) is True:
            comp["verified_pass"] += 1
        elif ver.get(k) is False:
            comp["verified_fail_submitted"] += 1
            comp["verified_fail_submitted_negatives"] += int(c.label == 0)
        else:
            comp["unverified"] += 1
log_result(dict(id=f"agent_{split}_gate_tool_compliance", phase=5 if split == "dev" else 6, split=split, metrics=comp))
print("gate_tool compliance", comp)
