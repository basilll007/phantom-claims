"""Dev-split evaluation (breast Rep1, seeds 0-4): H1 for expression tools, gate tuning grid (PREREG s.4),
H2 on dev with the selected config, ablations, and exploratory knob sweeps.
Usage: python scripts/dev_eval.py [tune|ablate|sweep|all]"""
import copy, itertools, json, os, sys, time
from dataclasses import asdict
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.generate import get_or_make
from bench.pipeline import (tool_lists, make_gate, evaluate_records, decide_list, labels_for, summarize,
                            log_result, METHODS)
from gate.gate import GateConfig, write_audit

DEV_SEEDS = [0, 1, 2, 3, 4]
GRID = dict(r=[20.0, 30.0, 50.0], tile=[150.0, 250.0, 400.0], n_min=[20, 50, 100], rho_min=[0.0, 0.02, 0.05])
DEFAULT = GateConfig()
os.makedirs("results/claims", exist_ok=True); os.makedirs("results/audit", exist_ok=True)


def load_dev(bench_cfg=None, tag="default"):
    t = load_tissue("breast_rep1")
    data = []
    for s in DEV_SEEDS:
        inst = get_or_make(t, s, bench_cfg)
        lists, df = tool_lists(t, inst, Ks=(5, 10, 20), seed=s)
        data.append((s, inst, lists))
        rows = []
        for (m, K), cl in lists.items():
            labs, subs = labels_for(inst, cl)
            for rank, (c, y, sb) in enumerate(zip(cl, labs, subs)):
                rows.append(dict(method=m, K=K, rank=rank + 1, **c, label=y, neg_type=sb))
        pd.DataFrame(rows).to_csv(f"results/claims/dev_{tag}_tools_s{s}.csv", index=False)
    return t, data


def h1_tools(data, tag, split="dev", tissue="breast_rep1"):
    out = {}
    for m in METHODS:
        for K in (5, 10, 20):
            groups, cats = [], {}
            for s, inst, lists in data:
                labs, subs = labels_for(inst, lists[(m, K)])
                groups.append((labs, None))
                for sb in subs:
                    cats[sb] = cats.get(sb, 0) + 1
            point, ci = summarize(groups)
            rec = log_result(dict(id=f"{split}_{tag}_{tissue}_H1_expr_{m}_K{K}", phase=3, split=split, tissue=tissue,
                                  seeds=[d[0] for d in data], condition=f"expr_{m}", K=K, hypothesis="H1",
                                  bench_tag=tag, metrics=point, ci=ci, claim_breakdown=cats))
            out[(m, K)] = rec
            print(f"H1 {m} K={K}: phantom={point['phantom_rate']:.3f} CI={np.round(ci['phantom_rate'], 3)} {cats}")
    return out


def tune(t, data):
    rows = []
    t0 = time.time()
    for r, tile in itertools.product(GRID["r"], GRID["tile"]):
        base = GateConfig(r=r, tile=tile)
        per_seed = []
        for s, inst, lists in data:
            gate = make_gate(t, inst, base)
            union = lists[("cpdb", 20)] + lists[("prod", 20)]
            per_seed.append((s, inst, lists, gate, evaluate_records(gate, union)))
        for n_min, rho_min in itertools.product(GRID["n_min"], GRID["rho_min"]):
            cfg = GateConfig(r=r, tile=tile, n_min=n_min, rho_min=rho_min)
            for s, inst, lists, gate, recs in per_seed:
                for m in METHODS:
                    for K in (10, 20):
                        cl = lists[(m, K)]
                        dec = decide_list(gate, recs, cl, cfg)
                        labs, _ = labels_for(inst, cl)
                        rows.append(dict(r=r, tile=tile, n_min=n_min, rho_min=rho_min, seed=s, method=m, K=K,
                                         labels=labs, passed=[d["pass_"] for d in dec]))
        print(f"grid r={r} tile={tile} done {time.time() - t0:.0f}s", flush=True)
    tab = []
    for key, grp in pd.DataFrame(rows).groupby(["r", "tile", "n_min", "rho_min", "K"]):
        y = np.concatenate([np.array(v) for v in grp.labels]); p = np.concatenate([np.array(v) for v in grp.passed])
        pr0 = 1 - y.mean(); pr1 = 1 - y[p].mean() if p.sum() else 0.0
        tab.append(dict(zip(["r", "tile", "n_min", "rho_min", "K"], key), n=len(y), pr_before=pr0, pr_after=pr1,
                        retention=y[p].sum() / max(y.sum(), 1), n_pass=int(p.sum())))
    tab = pd.DataFrame(tab)
    tab.to_csv("results/tables/gate_tuning_dev.csv", index=False)
    k10 = tab[tab.K == 10].copy()
    ok = k10[k10.retention >= 0.90].copy()
    d = DEFAULT
    ok["dist"] = ((ok.r != d.r).astype(int) + (ok.tile != d.tile).astype(int) + (ok.n_min != d.n_min).astype(int)
                  + (ok.rho_min != d.rho_min).astype(int))
    best = ok.sort_values(["pr_after", "dist"]).iloc[0] if len(ok) else None
    sel = None
    if best is not None:
        sel = dict(r=float(best.r), tile=float(best.tile), n_min=int(best.n_min), rho_min=float(best.rho_min))
        print("SELECTED", sel, best.to_dict())
    log_result(dict(id="dev_gate_tuning", phase=4, split="dev", tissue="breast_rep1", seeds=DEV_SEEDS,
                    grid=GRID, selection_rule="min pr_after (K=10, pooled cpdb+prod) s.t. retention>=0.90; ties->defaults",
                    selected=sel, n_configs=len(k10), selected_row=None if best is None else best.to_dict(),
                    runtime_s=round(time.time() - t0, 1)))
    json.dump(sel, open("results/tables/gate_selected_dev.json", "w"), indent=1)
    return sel


def h2(t, data, cfg, tag, split="dev", tissue="breast_rep1", ablation="none", write=True):
    res = {}
    for m in METHODS:
        for K in (5, 10, 20):
            groups, audits = [], []
            for s, inst, lists in data:
                gate = make_gate(t, inst, cfg)
                cl = lists[(m, K)]
                recs = gate.check(cl)
                labs, subs = labels_for(inst, cl)
                for rr, y, sb in zip(recs, labs, subs):
                    rr["truth"] = dict(label=y, neg_type=sb, seed=s)
                audits += recs
                groups.append((labs, [r_["pass_"] for r_ in recs]))
            point, ci = summarize(groups)
            if write and K == 10:
                write_audit(audits, f"results/audit/{split}_{tag}_{tissue}_{m}_K{K}_{ablation}.jsonl")
            rec = log_result(dict(id=f"{split}_{tag}_{tissue}_H2_expr_{m}_K{K}_{ablation}", phase=4, split=split,
                                  tissue=tissue, seeds=[d[0] for d in data], condition=f"expr_{m}_gated", K=K,
                                  hypothesis="H2", ablation=ablation, bench_tag=tag, gate_config=asdict(cfg),
                                  metrics=point, ci=ci))
            res[(m, K)] = rec
            print(f"H2[{ablation}] {m} K={K}: PR {point['phantom_rate']:.3f}->{point['phantom_rate_after']:.3f} "
                  f"red={point['rel_reduction']:.3f} CI={np.round(ci['rel_reduction'], 3)} "
                  f"ret={point['retention']:.3f} CI={np.round(ci['retention'], 3)}")
    return res


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "all"
    t, data = load_dev()
    if stage in ("tune", "all"):
        h1_tools(data, "default")
        sel = tune(t, data)
    sel = json.load(open("results/tables/gate_selected_dev.json"))
    cfg = GateConfig(**sel)
    if stage in ("tune", "all"):
        h2(t, data, cfg, "default")
    if stage in ("ablate", "all"):
        for name, kw in [("no_tiles", dict(use_tiles=False)), ("no_split", dict(use_split=False)),
                         ("no_min_evidence", dict(use_min_evidence=False)), ("no_perm", dict(use_perm=False))]:
            h2(t, data, GateConfig(**{**sel, **kw}), "default", ablation=name, write=False)
    if stage in ("sweep", "all"):
        for tag, bc in [("effect0.5", dict(effect=0.5)), ("effect0.25", dict(effect=0.25)),
                        ("noise0.5", dict(noise=0.5)), ("r50", dict(r_plant=50.0))]:
            t2, d2 = load_dev(bc, tag)
            h1_tools(d2, tag)
            h2(t2, d2, cfg, tag, write=False)
