"""Lane A exploratory (DEV seeds only; N2 is seen and must not be used for this): which component of the frozen gate
stops, or fails to stop, the H4 failure modes?

For variants A0 (reference), A4 (no readout) and A5 (mixed + hard negatives) on dev seeds, we compare the full frozen
gate with one-component ablations, using the gate's own switches (no new code paths):
no_tiles, no_split, no_min_evidence, no_perm.
Metrics:
- K=10 lists, cpdb+prod pooled: retention and post-gate phantom rate;
- planted-claim pass rate per category (true / niche / niche_hard / decoy_sender / phantom / random).
This is a descriptive failure map, not tuning: nothing here feeds back into the gate, the variants or the thresholds.

  python -m scripts.laneA_ablate <tissue> <seeds_csv>
"""
import os, sys
from dataclasses import replace
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.pipeline import tool_lists, make_gate, evaluate_records, decide_list
from bench.variants import make_variant, label_claim_v, planted_claims
from bench.metrics import hier_bootstrap
from gate.gate import GateConfig
from scripts.laneA_run import SEEDS, load_base, frozen_cfg
from scripts.laneA_paths import DATA, TABLES

ABL = {"full": {}, "no_tiles": dict(use_tiles=False), "no_split": dict(use_split=False),
       "no_min_evidence": dict(use_min_evidence=False), "no_perm": dict(use_perm=False)}
VARS = ["A0", "A4", "A5"]


def main(tissue, seeds):
    assert set(seeds) <= set(SEEDS["dev"]), "dev seeds only"
    base_cfg = frozen_cfg()
    t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
    rows = []
    for s in seeds:
        base = load_base(t, s)
        lists0 = None
        for v in VARS:
            inst = make_variant(t, base, v, s)
            if v == "A5":
                lists, _ = tool_lists(t, inst, Ks=(10,), seed=s)
            else:
                if lists0 is None:
                    lists0, _ = tool_lists(t, base, Ks=(10,), seed=s)
                lists = lists0
            planted = planted_claims(inst)
            for a, kw in ABL.items():
                cfg = replace(base_cfg, **kw)
                g = make_gate(t, inst, cfg)
                allc = [c for cl in lists.values() for c in cl] + [c for c, _, _ in planted]
                recs = evaluate_records(g, allc)
                for (m, K), cl in lists.items():
                    for c, d in zip(cl, decide_list(g, recs, cl)):
                        lab, sub = label_claim_v(inst, c)
                        rows.append(dict(tissue=tissue, seed=s, variant=v, ablation=a, list=f"{m}_K{K}", label=lab,
                                         category="true" if lab else sub, passed=bool(d["pass_"])))
                for (c, cat, _), d in zip(planted, decide_list(g, recs, [c for c, _, _ in planted])):
                    rows.append(dict(tissue=tissue, seed=s, variant=v, ablation=a, list="planted",
                                     label=int(cat == "true"), category=cat, passed=bool(d["pass_"])))
            print(tissue, s, v, "done", flush=True)
    d = pd.DataFrame(rows)
    os.makedirs(TABLES, exist_ok=True)
    d.to_csv(f"{TABLES}/laneA_dev_ablation_claims_{tissue}.csv", index=False)
    out = []
    for (v, a), g in d.groupby(["variant", "ablation"]):
        L = g[g.list != "planted"]
        groups = [(x.label.values, x.passed.values) for _, x in L.groupby("seed")]
        point, ci = hier_bootstrap(groups, keys=("phantom_rate_after", "retention"), n_boot=1000)
        rec = dict(variant=v, ablation=a, K10_retention=point["retention"], K10_phantom_after=point["phantom_rate_after"],
                   K10_ret_lo=ci["retention"][0], K10_ret_hi=ci["retention"][1],
                   K10_pr_lo=ci["phantom_rate_after"][0], K10_pr_hi=ci["phantom_rate_after"][1])
        P = g[g.list == "planted"]
        for cat, x in P.groupby("category"):
            rec[f"pass_{cat}"] = float(x.passed.mean())
        out.append(rec)
    out = pd.DataFrame(out)
    out.to_csv(f"{TABLES}/laneA_dev_ablation_summary_{tissue}.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(3).to_string(index=False))
    return out


if __name__ == "__main__":
    main(sys.argv[1], [int(x) for x in sys.argv[2].split(",")])
