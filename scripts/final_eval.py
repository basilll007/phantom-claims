"""Phase 6 — final evaluation with the FROZEN gate (PREREG Amendment 2) on held-out test sets. Run once.
Also per-tissue planted-claim pass rates and null calibration (robustness, same frozen gate)."""
import json, os, sys, time
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.generate import get_or_make
from bench.pipeline import tool_lists, make_gate, labels_for, log_result
from gate.gate import GateConfig
from scripts.dev_eval import h1_tools, h2

TEST = {"breast_rep1": [100, 101, 102, 103, 104], "breast_rep2": [100, 101, 102, 103, 104],
        "pancreas": [100, 101, 102, 103, 104], "lung": [100, 101, 102, 103, 104]}
sel = json.load(open("results/tables/gate_selected_dev.json"))
cfg = GateConfig(**sel)
only = sys.argv[1:] or list(TEST)
marker = "results/tables/final_eval_done.json"
done = json.load(open(marker)) if os.path.exists(marker) else {}
for tissue in only:
    if tissue in done:
        print("already evaluated once:", tissue, "- refusing to re-run"); continue
    t0 = time.time()
    t = load_tissue(tissue)
    data = []
    for s in TEST[tissue]:
        inst = get_or_make(t, s)
        lists, _ = tool_lists(t, inst, Ks=(5, 10, 20), seed=s)
        data.append((s, inst, lists))
        rows = []
        for (m, K), cl in lists.items():
            labs, subs = labels_for(inst, cl)
            rows += [dict(method=m, K=K, rank=i + 1, **c, label=y, neg_type=sb) for i, (c, y, sb) in enumerate(zip(cl, labs, subs))]
        pd.DataFrame(rows).to_csv(f"results/claims/test_{tissue}_tools_s{s}.csv", index=False)
    h1_tools(data, "final", split="test", tissue=tissue)
    h2(t, data, cfg, "final", split="test", tissue=tissue)
    # planted-claim pass rates by category (frozen gate)
    pl = []
    for s, inst, _ in data:
        pc = [dict(sender=p["sender"], receiver=p["receiver"], ligand=p["ligand"], receptor=p["receptor"]) for p in inst.truth]
        recs = make_gate(t, inst, cfg).check(pc)
        pl += [dict(seed=s, category=p["category"], passed=r["pass_"]) for p, r in zip(inst.truth, recs)]
    pl = pd.DataFrame(pl)
    rates = pl.groupby("category").passed.mean().round(3).to_dict()
    print(tissue, "planted pass rates", rates)
    log_result(dict(id=f"test_final_planted_pass_rates_{tissue}", phase=6, split="test", tissue=tissue,
                    seeds=TEST[tissue], gate_config=sel, metrics=rates, runtime_s=round(time.time() - t0, 1)))
    done[tissue] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    json.dump(done, open(marker, "w"), indent=1)
