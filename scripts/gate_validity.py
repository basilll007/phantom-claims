"""Validity checks on dev (breast Rep1 seeds 0-4): (1) null calibration of the gate's permutation p-values on
claims with no planted signal; (2) gate pass rates on all planted claims by category, full gate vs ablations;
(3) which negatives pass in tool top-20 lists."""
import json, sys
import numpy as np, pandas as pd
from dataclasses import asdict
from bench.tissue import load_tissue
from bench.generate import get_or_make, label_claim
from bench.pipeline import make_gate, log_result
from gate.gate import GateConfig

sel = json.load(open("results/tables/gate_selected_dev.json"))
t = load_tissue("breast_rep1")
rng = np.random.default_rng(123)
calib, planted = [], []
for s in [0, 1, 2, 3, 4]:
    inst = get_or_make(t, s)
    types = t.type_list(min_cells=200)
    # (1) null claims: random-category pairs (no planting at all) at random type pairs
    rnd = [p for p in inst.truth if p["category"] == "random"]
    claims = []
    for p in rnd:
        for _ in range(4):
            a, b = rng.choice(len(types), 2, replace=False)
            claims.append(dict(sender=types[a], receiver=types[b], ligand=p["ligand"], receptor=p["receptor"]))
    for name, kw in [("full", {}), ("no_tiles", dict(use_tiles=False))]:
        g = make_gate(t, inst, GateConfig(**{**sel, **kw}))
        for c in claims:
            r = g.evaluate(c)
            calib.append(dict(seed=s, variant=name, p=r.get("p_perm"), rho=r.get("rho_within_tile")))
    # (2) every planted claim at its planted type pair
    pc = [dict(sender=p["sender"], receiver=p["receiver"], ligand=p["ligand"], receptor=p["receptor"]) for p in inst.truth]
    cats = [p["category"] for p in inst.truth]
    for name, kw in [("full", {}), ("no_tiles", dict(use_tiles=False)), ("no_split", dict(use_split=False)),
                     ("no_min_evidence", dict(use_min_evidence=False)), ("no_perm", dict(use_perm=False))]:
        g = make_gate(t, inst, GateConfig(**{**sel, **kw}))
        recs = g.check(pc)
        for c, r, cat in zip(pc, recs, cats):
            planted.append(dict(seed=s, variant=name, category=cat, passed=r["pass_"], p=r.get("p_perm"),
                                rho=r.get("rho_within_tile"), n_pairs=r.get("n_expressing_pairs_within_r"),
                                reasons=";".join(r.get("reasons", []))))
cal = pd.DataFrame(calib); pl = pd.DataFrame(planted)
cal.to_csv("results/tables/gate_null_calibration_dev.csv", index=False)
pl.to_csv("results/tables/gate_planted_claims_dev.csv", index=False)
cs = cal.dropna().groupby("variant").p.agg(lambda v: dict(n=len(v), frac_p_lt_0_05=float((v < 0.05).mean()),
                                                         frac_p_lt_0_01=float((v < 0.01).mean()),
                                                         mean_p=float(v.mean()))).to_dict()
ps = pl.groupby(["variant", "category"]).passed.mean().unstack().round(3)
print("NULL CALIBRATION", json.dumps(cs, indent=1))
print("PASS RATE BY CATEGORY\n", ps.to_string())
log_result(dict(id="dev_gate_null_calibration", phase=4, split="dev", tissue="breast_rep1", seeds=[0, 1, 2, 3, 4],
                gate_config=sel, metrics=cs))
log_result(dict(id="dev_gate_planted_pass_rates", phase=4, split="dev", tissue="breast_rep1", seeds=[0, 1, 2, 3, 4],
                gate_config=sel, metrics={v: ps.loc[v].to_dict() for v in ps.index}))
