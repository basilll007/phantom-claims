"""Lane B (Phase 2, H5) runner: dose-response of expression-only phantom rate and frozen-gate phantom rate.

  python -m scripts.laneB_run <tissue> <seed_set> [seeds_csv] [levels_csv] [list_types_csv]

Per (tissue, seed, effect level), with the effect applied by the unchanged default generator (`effect` x default
planted signal, every planted category):
  default : the default instance itself -> expression-tool top-K lists -> frozen gate     (H5 confirmatory)
  A5m     : gate-matched TRUE response + A5 hard negatives (decoy senders, hard niche)    (exploratory, addendum)
  A5      : the exact Lane-A A5 construction (mixed mechanisms + hard negatives)          (exploratory, continuity)
One CSV per (tissue, seed, level, list type) under $PC_LANEB_CLAIMS (default results/claims/laneB). Prints timings and
row counts only, never metrics (run-once rule, reports/laneB_prespec.md).
"""
import os, sys, time
import pandas as pd
import torch
from bench.tissue import load_tissue
from bench.generate import get_or_make, generate
from bench.pipeline import tool_lists
from bench.variants import make_variant, save_truth
from scripts.laneA_run import SEEDS, FORBIDDEN, frozen_cfg, eval_variant, GIT
from scripts.laneA_paths import DATA, BENCH

OUT = os.environ.get("PC_LANEB_CLAIMS", "results/claims/laneB")
LEVELS = [1.0, 0.5, 0.25, 0.125]
LIST_TYPES = ["default", "A5m", "A5"]
KS = (5, 10, 20)


def base_for(t, seed, e):
    cfg = None if e == 1.0 else {"effect": e}           # effect=1 -> the exact default N2 instance (shared cache)
    try:
        return get_or_make(t, seed, cfg, root=BENCH)
    except Exception as ex:                              # concurrent writer: rebuild in memory
        print(f"[warn] cache read failed ({ex!r}); regenerating in memory", flush=True)
        return generate(t, seed, cfg)


def run(tissue, seed_set, seeds=None, levels=None, list_types=None):
    seeds = seeds or SEEDS[seed_set]
    assert set(seeds) <= set(SEEDS[seed_set]) and not (set(seeds) & FORBIDDEN), "seed outside the declared set"
    levels = levels or LEVELS
    assert set(levels) <= set(LEVELS)
    list_types = list_types or LIST_TYPES
    cfg = frozen_cfg()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
    W = {}
    for s in seeds:
        for e in levels:
            todo = [lt for lt in list_types if not os.path.exists(f"{OUT}/{seed_set}_{tissue}_s{s}_e{e:g}_{lt}.csv")]
            if not todo:
                continue
            ts = time.time()
            base = base_for(t, s, e)
            assert abs(base.cfg["effect"] - e) < 1e-12
            for lt in todo:
                inst = base if lt == "default" else make_variant(t, base, lt, s)
                lists, _ = tool_lists(t, inst, Ks=KS, seed=s)
                df = eval_variant(t, inst, lists, cfg, W)
                for col, val in (("list_type", lt), ("effect", e), ("seed", s), ("tissue", tissue), ("seed_set", seed_set)):
                    df.insert(0, col, val)
                df["git_hash"] = GIT
                if lt != "default":
                    save_truth(inst, f"{BENCH}/laneB/{seed_set}")
                tmp = f"{OUT}/{seed_set}_{tissue}_s{s}_e{e:g}_{lt}.csv"
                df.to_csv(tmp + ".tmp", index=False); os.replace(tmp + ".tmp", tmp)
            vram = torch.cuda.max_memory_allocated() / 1e6 if torch.cuda.is_available() else 0
            print(f"[{tissue} s{s} e{e:g}] {len(todo)} list types in {time.time() - ts:.0f}s "
                  f"(total {time.time() - t0:.0f}s) vram_peak={vram:.0f}MB", flush=True)


if __name__ == "__main__":
    a = sys.argv
    run(a[1], a[2], [int(x) for x in a[3].split(",")] if len(a) > 3 and a[3] else None,
        [float(x) for x in a[4].split(",")] if len(a) > 4 and a[4] else None,
        a[5].split(",") if len(a) > 5 else None)
