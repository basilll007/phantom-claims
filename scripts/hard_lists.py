"""Exploratory 'hard list' analysis: per seed, a list of all 10 TRUE + all 10 NICHE planted claims (PR_before = 0.5),
gated by the full gate and by ablations, from the per-claim decisions in results/tables/gate_planted_claims_dev.csv.
H2-style metrics with hierarchical bootstrap."""
import numpy as np, pandas as pd
from bench.metrics import hier_bootstrap
from bench.pipeline import log_result

d = pd.read_csv("results/tables/gate_planted_claims_dev.csv")
d = d[d.category.isin(["true", "niche"])]
out = {}
for v, g in d.groupby("variant"):
    groups = [((x.category == "true").astype(int).values, x.passed.astype(bool).values) for _, x in g.groupby("seed")]
    point, ci = hier_bootstrap(groups, keys=("phantom_rate_after", "retention", "rel_reduction"))
    out[v] = dict(phantom_rate_before=point["phantom_rate"], phantom_rate_after=point["phantom_rate_after"],
                  retention=point["retention"], rel_reduction=point["rel_reduction"], ci=ci)
    print(v, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in out[v].items() if k != "ci"},
          {k: np.round(c, 3).tolist() for k, c in ci.items()})
log_result(dict(id="dev_hard_list_true_vs_niche", phase=4, split="dev", tissue="breast_rep1", exploratory=True,
                design="per seed: 10 TRUE + 10 NICHE planted claims at their planted type pairs", metrics=out))
