"""Exploratory stress test of gate calibration on REAL structured data (review finding #4).
Breast Rep1, no planting. 200 random claims: random distinct sender/receiver types (>=200 cells), ligand, receptor and
3 'targets' = random real panel genes. Real tissue keeps regional structure and segmentation spillover, so this is a
harder null than within-type-shuffled channels. Reports p<0.05 and full-gate pass rates (frozen gate vs. global null),
stratified by whether the targets look like sender-type markers (spillover risk)."""
import json
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.pipeline import log_result
from gate.gate import Gate, GateConfig, ExprSource

sel = json.load(open("results/tables/gate_selected_dev.json"))
t = load_tissue("breast_rep1")
rng = np.random.default_rng(2026)
types = t.type_list(min_cells=200)
genes = list(t.genes)
X = t.counts.tocsc()
tmean = {k: np.asarray(X[t.types == k].mean(0)).ravel() for k in types}
claims, targets, marker = [], {}, []
for i in range(200):
    a, b = rng.choice(len(types), 2, replace=False)
    g = rng.choice(len(genes), 5, replace=False)
    L, R, T = genes[g[0]], genes[g[1]], [genes[x] for x in g[2:]]
    claims.append(dict(sender=types[a], receiver=types[b], ligand=L, receptor=R))
    targets[(L, R)] = T
    mA = np.mean([tmean[types[a]][genes.index(x)] for x in T]); mB = np.mean([tmean[types[b]][genes.index(x)] for x in T])
    marker.append(bool(mA > 2 * mB + 0.05))
res = {}
for name, kw in [("frozen", {}), ("no_tiles", dict(use_tiles=False))]:
    g = Gate(t, ExprSource(t, t.counts, genes, targets), GateConfig(**{**sel, **kw}))
    recs = g.check(claims)
    p = np.array([r.get("p_perm", 1.0) if r.get("p_perm") is not None else 1.0 for r in recs])
    ps = np.array([r["pass_"] for r in recs]); mk = np.array(marker)
    res[name] = dict(n=len(p), frac_p_lt_0_05=float((p < 0.05).mean()), frac_pass=float(ps.mean()),
                     n_marker_targets=int(mk.sum()),
                     frac_p_lt_0_05_marker=float((p[mk] < 0.05).mean()) if mk.any() else None,
                     frac_p_lt_0_05_nonmarker=float((p[~mk] < 0.05).mean()),
                     frac_pass_marker=float(ps[mk].mean()) if mk.any() else None,
                     frac_pass_nonmarker=float(ps[~mk].mean()))
    print(name, res[name], flush=True)
log_result(dict(id="real_gene_null_breast_rep1", phase=7, split="real", tissue="breast_rep1", exploratory=True,
                gate_config=sel, metrics=res, design="200 random real-gene claims, random type pairs, no planting"))
