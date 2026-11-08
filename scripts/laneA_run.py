"""Lane A (Phase 2, H4) runner: misspecified variants, frozen gate, expression-only baselines.

  python -m scripts.laneA_run run   <tissue> <seed_set> [seeds_csv] [variants_csv]
  python -m scripts.laneA_run pilot <tissue> <seeds_csv>         # dev seeds only: A5 list composition, NO gate

Per (tissue, seed): the default base instance (shared cache, standard name, reusable by other lanes) ->
expression-tool top-K lists (identical for A0-A4 because L/R channels are identical) -> for each variant the
frozen gate on the union of tool lists + all planted claims -> one CSV per (tissue, seed, variant) under
results/claims/laneA/. Aggregation, CIs and logging are in scripts/laneA_aggregate.py.

Env: PC_DATA = data root holding processed/ and bench/ (default 'data').
"""
import json, os, sys, time
from dataclasses import asdict
import numpy as np, pandas as pd
import torch
from bench.tissue import load_tissue
from bench.generate import get_or_make, generate
from bench.pipeline import tool_lists, make_gate, evaluate_records, decide_list, claim_key
from bench.variants import VARIANTS, make_variant, label_claim_v, planted_claims, save_truth
from gate.gate import GateConfig

SEEDS = {"dev": [0, 1, 2, 3, 4], "N2": [200, 201, 202, 203, 204, 205, 206, 207, 208, 209]}
FORBIDDEN = {100, 101, 102, 103, 104}                       # Phase-1 test1 seeds: never used here
FROZEN = dict(r=30.0, tile=250.0, n_perm=1000, alpha_half=0.10, n_min=50, rho_min=0.02, q=0.05)
DATA = os.environ.get("PC_DATA", "data")
OUT = "results/claims/laneA"
KS = (5, 10, 20)


def frozen_cfg():
    cfg = GateConfig()
    assert {k: getattr(cfg, k) for k in FROZEN} == FROZEN, "gate config drifted from frozen defaults"
    sel = json.load(open("results/tables/gate_selected_dev.json"))
    assert all(float(getattr(cfg, k)) == float(v) for k, v in sel.items()), "frozen != Phase-1 selected config"
    return cfg


def load_base(t, seed):
    try:
        return get_or_make(t, seed, None, root=f"{DATA}/bench")
    except Exception as e:                                   # concurrent writer from another lane: rebuild in memory
        print(f"[warn] base cache read failed ({e!r}); regenerating in memory", flush=True)
        return generate(t, seed, None)


def mech_of(inst):
    return {p["ligand"]: p.get("mechanism", "") for p in inst.truth if p["category"] == "true"}


def eval_variant(t, inst, lists, cfg, shared_W):
    gate = make_gate(t, inst, cfg)
    gate._W = shared_W                                       # type-pair distance weights depend on tissue only
    planted = planted_claims(inst)
    allc = [c for cl in lists.values() for c in cl] + [c for c, _, _ in planted]
    recs = evaluate_records(gate, allc)
    mech = mech_of(inst)
    rows = []

    def add_rows(listname, method, K, cl, cats=None):
        dec = decide_list(gate, recs, cl)
        for rank, (c, d) in enumerate(zip(cl, dec)):
            lab, sub = label_claim_v(inst, c)
            rows.append(dict(list=listname, method=method, K=K, rank=rank + 1, **c, label=lab,
                             neg_type=sub if not lab else "", category=cats[rank] if cats else "",
                             mechanism=mech.get(c["ligand"], "") if (lab or (cats and cats[rank] == "decoy_sender")) else "",
                             passed=bool(d["pass_"]), q_bh=d.get("q_bh"), rho=d.get("rho_within_tile"),
                             p_perm=d.get("p_perm"), n_expr_pairs=d.get("n_expressing_pairs_within_r"),
                             effect_d=d.get("effect_d"), reasons=";".join(d.get("reasons", []))))

    for (m, K), cl in lists.items():
        add_rows(f"{m}_K{K}", m, K, cl)
    add_rows("planted", "", 0, [c for c, _, _ in planted], cats=[k for _, k, _ in planted])
    return pd.DataFrame(rows)


def run(tissue, seed_set, seeds=None, variants=None):
    seeds = seeds or SEEDS[seed_set]
    assert set(seeds) <= set(SEEDS[seed_set]) and not (set(seeds) & FORBIDDEN), "seed outside the declared set"
    variants = variants or list(VARIANTS)
    cfg = frozen_cfg()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
    print(f"[{tissue}] loaded n={t.n} in {time.time() - t0:.0f}s", flush=True)
    shared_W = {}
    for s in seeds:
        todo = [v for v in variants if not os.path.exists(f"{OUT}/{seed_set}_{tissue}_s{s}_{v}.csv")]
        if not todo:
            print(f"[{tissue} s{s}] all variants done, skip", flush=True); continue
        ts = time.time()
        base = load_base(t, s)
        lists = None
        for v in todo:
            tv = time.time()
            inst = make_variant(t, base, v, s)
            if v == "A5":
                lv, _ = tool_lists(t, inst, Ks=KS, seed=s)
            else:
                if lists is None:
                    lists, _ = tool_lists(t, base, Ks=KS, seed=s)
                lv = lists
                lr = [inst.gene_index()[g] for p in base.truth for g in (p["ligand"], p["receptor"])]
                assert (inst.added[:, lr] != base.added[:, lr]).nnz == 0, "L/R channels must be identical to base"
            df = eval_variant(t, inst, lv, cfg, shared_W)
            df.insert(0, "variant", v); df.insert(0, "seed", s); df.insert(0, "tissue", tissue)
            df.insert(0, "seed_set", seed_set)
            save_truth(inst, f"{DATA}/bench/laneA/{seed_set}")
            tmp = f"{OUT}/{seed_set}_{tissue}_s{s}_{v}.csv.tmp"
            df.to_csv(tmp, index=False)
            os.replace(tmp, f"{OUT}/{seed_set}_{tissue}_s{s}_{v}.csv")
            k10 = df[df.K == 10]
            vram = torch.cuda.max_memory_allocated() / 1e6 if torch.cuda.is_available() else 0
            print(f"[{tissue} s{s} {v}] {time.time() - tv:.0f}s  K10 rows={len(k10)}  vram_peak={vram:.0f}MB",
                  flush=True)
        print(f"[{tissue} s{s}] seed done in {time.time() - ts:.0f}s (total {time.time() - t0:.0f}s)", flush=True)


def pilot(tissue, seeds, configs):
    """Dev-only benchmark pilot: does each A5 hard-negative config put niche / decoy items into K=10 lists?
    Looks at list composition only; the gate is never run here."""
    assert set(seeds) <= set(SEEDS["dev"])
    t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
    out = []
    for s in seeds:
        base = load_base(t, s)
        for name, nh in configs:
            inst = make_variant(t, base, "A5", s, overrides={"niche_hard": nh})
            lists, _ = tool_lists(t, inst, Ks=(10, 20), seed=s)
            for (m, K), cl in lists.items():
                subs = [label_claim_v(inst, c) for c in cl]
                comp = pd.Series([("true" if y else n) for y, n in subs]).value_counts().to_dict()
                out.append(dict(seed=s, config=name, method=m, K=K, **comp))
                print(s, name, m, K, comp, flush=True)
    return pd.DataFrame(out).fillna(0)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "run":
        seeds = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 and sys.argv[4] else None
        variants = sys.argv[5].split(",") if len(sys.argv) > 5 else None
        run(sys.argv[2], sys.argv[3], seeds, variants)
    elif mode == "pilot":
        PILOT = [("c1_s400_b2_x1", dict(sigma=400.0, blobs=2, boost=1.0, cap=20.0)),
                 ("c2_s400_b2_x2", dict(sigma=400.0, blobs=2, boost=2.0, cap=20.0)),
                 ("c3_s600_b3_x1", dict(sigma=600.0, blobs=3, boost=1.0, cap=20.0)),
                 ("c4_s600_b3_x2", dict(sigma=600.0, blobs=3, boost=2.0, cap=20.0))]
        df = pilot(sys.argv[2], [int(x) for x in sys.argv[3].split(",")], PILOT)
        os.makedirs("results/tables", exist_ok=True)
        df.to_csv("results/tables/laneA_pilot_dev_composition.csv", index=False)
        print(df.groupby(["config", "K"]).mean(numeric_only=True).round(2).to_string())
