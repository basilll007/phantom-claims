"""Lane R (review response), EXPLORATORY diagnostics only: how sensitive are the frozen gate's decisions to its fixed
hyperparameters, its split scheme, its null, border tiles and the response normalisation?

DEV seeds 0-4 of breast Rep1 only (the Phase-1 tuning seeds; N2 is not touched here). Variants A0, A1, A2a, A2b, A3 and
A5 (bench/variants.py). One setting is changed at a time around the frozen gate (Amendment 2 defaults):

  frozen          GateConfig() as frozen
  tile125/tile500 tile size 125 / 500 um                      (GateConfig.tile)
  r15/r60         exposure range 15 / 60 um                   (GateConfig.r; kernel exp(-d/(r/2)), d <= r)
  nmin25/nmin100  minimum-evidence contacts 25 / 100          (GateConfig.n_min)
  split_checker   split halves = tile checkerboard ((tx+ty) mod 2) instead of left/right of the median x
  null_global     global permutation null                     (GateConfig.use_tiles=False; the gate's own switch)
  null_density    permutation (and centring) within density strata inside tiles: groups = tile x tercile of
                  local density (number of cells within 30 um, terciles over all cells of the section)
  drop_edge       receivers in border tiles removed (a tile is a border tile if any of its 8 neighbours holds no
                  cell, i.e. it touches the tissue edge or a hole)
  resp_raw        response z-scores computed on log1p(raw counts) instead of log1p(count x median_lib / lib);
                  the ligand exposure keeps the frozen normalisation

Nothing here changes gate/gate.py (hash-frozen): settings that the GateConfig cannot express are subclasses defined in
this file. Nothing feeds back into the gate, the variants, the thresholds or any verdict.

  python -m scripts.laneR_sensitivity <seed>          # one dev seed -> results/tables/laneR_sens_claims_breast_rep1_s<seed>.csv
"""
import os, sys, time
from dataclasses import replace
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from bench.tissue import load_tissue
from bench.pipeline import tool_lists, evaluate_records, decide_list
from bench.variants import make_variant, label_claim_v, planted_claims
from gate.gate import Gate, GateConfig, ExprSource
from scripts.laneA_run import SEEDS, load_base, frozen_cfg
from scripts.laneA_paths import DATA, TABLES

TISSUE = "breast_rep1"
VARS = ["A0", "A1", "A2a", "A2b", "A3", "A5"]
CFG = {"frozen": {}, "tile125": dict(tile=125.0), "tile500": dict(tile=500.0), "r15": dict(r=15.0),
       "r60": dict(r=60.0), "nmin25": dict(n_min=25), "nmin100": dict(n_min=100), "null_global": dict(use_tiles=False)}
SPECIAL = ["split_checker", "null_density", "drop_edge", "resp_raw"]
SETTINGS = list(CFG) + SPECIAL


class RawZ(ExprSource):
    """Response z-scores on log1p(raw counts); exposure (log/log_multi) keeps the frozen library-size normalisation."""
    def z(self, g):
        if g not in self._z:
            v = np.log1p(self.raw(g)); self._z[g] = (v - v.mean()) / (v.std() + 1e-9)
        return self._z[g]


class CheckerGate(Gate):
    def __init__(self, t, x, cfg):
        super().__init__(t, x, cfg)
        tx = np.floor(t.xy[:, 0] / cfg.tile).astype(np.int64); ty = np.floor(t.xy[:, 1] / cfg.tile).astype(np.int64)
        self.half = ((tx + ty) % 2).astype(np.int64)


class DensityGate(Gate):
    def __init__(self, t, x, cfg, n30):
        super().__init__(t, x, cfg)
        cut = np.quantile(n30, [1 / 3, 2 / 3])
        strata = np.searchsorted(cut, n30, side="right").astype(np.int64)          # 0, 1, 2
        self.tile_id = self.tile_id * 3 + strata


class EdgeGate(Gate):
    def __init__(self, t, x, cfg, edge):
        super().__init__(t, x, cfg)
        self.edge = edge

    def _weights(self, A, B):
        key = ("noedge", A, B)
        if key not in self._W:
            iA, iB, W = super()._weights(A, B)
            if W is None:
                self._W[key] = (iA, iB, None)
            else:
                keep = ~self.edge[iB]
                self._W[key] = (iA, iB[keep], W[keep]) if keep.sum() else (iA, iB[keep], None)
        return self._W[key]


def tissue_features(t, tile):
    n30 = cKDTree(t.xy).query_ball_point(t.xy, 30.0, return_length=True) - 1
    tx = np.floor(t.xy[:, 0] / tile).astype(np.int64); ty = np.floor(t.xy[:, 1] / tile).astype(np.int64)
    occ = set(zip(tx.tolist(), ty.tolist()))
    edge_tile = {}
    for (a, b) in occ:
        edge_tile[(a, b)] = any((a + i, b + j) not in occ for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0))
    edge = np.array([edge_tile[(a, b)] for a, b in zip(tx.tolist(), ty.tolist())])
    return n30, edge


def make(t, inst, setting, feats):
    tg = {(d["ligand"], d["receptor"]): d["targets"] for d in inst.database}
    base = frozen_cfg()
    if setting in CFG:
        return Gate(t, ExprSource(t, inst.counts, inst.syn_genes, tg), replace(base, **CFG[setting]))
    if setting == "resp_raw":
        return Gate(t, RawZ(t, inst.counts, inst.syn_genes, tg), base)
    x = ExprSource(t, inst.counts, inst.syn_genes, tg)
    if setting == "split_checker":
        return CheckerGate(t, x, base)
    if setting == "null_density":
        return DensityGate(t, x, base, feats[0])
    if setting == "drop_edge":
        return EdgeGate(t, x, base, feats[1])
    raise KeyError(setting)


def main(seed):
    assert seed in SEEDS["dev"], "dev seeds only (N2 is not used for this diagnostic)"
    t = load_tissue(TISSUE, f"{DATA}/processed/{TISSUE}.h5ad")
    feats = tissue_features(t, frozen_cfg().tile)
    base = load_base(t, seed)
    lists0, _ = tool_lists(t, base, Ks=(10,), seed=seed)
    rows, t0 = [], time.time()
    for v in VARS:
        inst = make_variant(t, base, v, seed)
        lists = tool_lists(t, inst, Ks=(10,), seed=seed)[0] if v == "A5" else lists0
        planted = planted_claims(inst)
        allc = [c for cl in lists.values() for c in cl] + [c for c, _, _ in planted]
        for s in SETTINGS:
            g = make(t, inst, s, feats)
            recs = evaluate_records(g, allc)
            def add(listname, cl, cats=None):
                for i, (c, d) in enumerate(zip(cl, decide_list(g, recs, cl))):
                    lab, sub = label_claim_v(inst, c)
                    cat = cats[i] if cats else ("true" if lab else sub)
                    iB = np.flatnonzero(t.types == c["receiver"])
                    rows.append(dict(tissue=TISSUE, seed=seed, variant=v, setting=s, list=listname,
                                     sender=c["sender"], receiver=c["receiver"], ligand=c["ligand"], receptor=c["receptor"],
                                     label=int(lab) if not cats else int(cat == "true"), category=cat,
                                     passed=bool(d["pass_"]), rho=d.get("rho_within_tile"), q=d.get("q_bh"),
                                     n_expr=d.get("n_expressing_pairs_within_r"),
                                     frac_receivers_edge=float(feats[1][iB].mean()) if len(iB) else np.nan,
                                     reasons=";".join(d.get("reasons", []))))
            for (m, K), cl in lists.items():
                add(f"{m}_K{K}", cl)
            add("planted", [c for c, _, _ in planted], [cat for _, cat, _ in planted])
            print(f"{TISSUE} s{seed} {v} {s} done ({time.time() - t0:.0f}s)", flush=True)
    os.makedirs(TABLES, exist_ok=True)
    out = f"{TABLES}/laneR_sens_claims_{TISSUE}_s{seed}.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(out, "edge receivers (all cells):", round(float(feats[1].mean()), 4))


if __name__ == "__main__":
    main(int(sys.argv[1]))
