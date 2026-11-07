"""Exploratory harder benchmark (dev only; review finding #1): NICHE expression cap 3 -> 10 so niche-confounded
pairs reach top-K lists. Reports H1/H2-style metrics with the frozen gate and its no-tiles / no-split ablations,
plus the negative composition of the lists."""
import json
import pandas as pd
from gate.gate import GateConfig
from scripts.dev_eval import load_dev, h1_tools, h2
from bench.pipeline import log_result

sel = json.load(open("results/tables/gate_selected_dev.json"))
tag = "nichecap10"
t, data = load_dev(dict(niche_cap=10.0), tag)
h1_tools(data, tag)
h2(t, data, GateConfig(**sel), tag, write=False)
h2(t, data, GateConfig(**sel, use_tiles=False), tag, ablation="no_tiles", write=False)
h2(t, data, GateConfig(**sel, use_split=False), tag, ablation="no_split", write=False)
comp = {}
for s, inst, lists in data:
    d = pd.read_csv(f"results/claims/dev_{tag}_tools_s{s}.csv")
    for K in (10, 20):
        for k, v in d[(d.K == K) & (d.label == 0)].neg_type.value_counts().items():
            comp[f"K{K}_{k}"] = comp.get(f"K{K}_{k}", 0) + int(v)
print("composition", comp)
log_result(dict(id=f"dev_{tag}_list_composition", phase=4, split="dev", tissue="breast_rep1", exploratory=True, metrics=comp))
