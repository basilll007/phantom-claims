"""Composition of negatives in the claim lists that the gate was evaluated on (tools: dev + test; agents: dev + test).
Answers: how often do hard negatives (niche, other) actually reach the gate?"""
import glob
import pandas as pd
from bench.pipeline import log_result

rows = []
for f in glob.glob("results/claims/dev_default_tools_s*.csv") + glob.glob("results/claims/test_*_tools_s*.csv"):
    d = pd.read_csv(f)
    d["split"] = "dev" if "/dev_" in f.replace("\\", "/") else "test"
    rows.append(d)
tools = pd.concat(rows)
out = {}
for K in (10, 20):
    x = tools[(tools.K == K) & (tools.label == 0)]
    out[f"tools_K{K}_negatives"] = x.neg_type.value_counts().to_dict()
    out[f"tools_K{K}_n_claims"] = int((tools.K == K).sum())
for split in ("dev", "test"):
    a = pd.read_csv(f"results/claims/agent_{split}_claims.csv")
    a = a[(a["rank"] > 0) & (a.label == 0)]
    for cond in ("unguided", "gate_tool"):
        out[f"agents_{split}_{cond}_negatives"] = a[a.cond == cond].neg_type.value_counts().to_dict()
print(out)
log_result(dict(id="list_negative_composition", phase=6, split="dev+test", metrics=out,
                note="which negative categories actually appear in gated lists (from results/claims/*.csv)"))
