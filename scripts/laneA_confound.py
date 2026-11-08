"""Lane A exploratory: is a gate pass driven by the planted response or by the library-size x density confound?

Background (found in the dev smoke run): Xenium library size falls with local cell density. The gate's response
readout is log1p(count * median_lib / lib), so the normalisation factor is higher in crowded regions, and exposure is
also higher there because there are more neighbours. A response with no planted signal can then correlate weakly with
exposure.

For each claim examined: (i) A4 TRUE planted claims (no response planted), (ii) every planted negative that passed in any
variant, (iii) A0 TRUE planted claims (gate-matched reference). We recompute, with the frozen gate's own exposure and
within-tile permutation test:
  rho_gate  : (E, Y_lognorm)   = the gate statistic
  rho_raw   : (E, Y_raw)       = same response without library-size normalisation
  rho_scale : (E, log(median_lib/lib) of receivers) = exposure vs normalisation factor
Nothing here changes the gate; this is a diagnostic.

  python -m scripts.laneA_confound <seed_set> [tissues_csv]
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.generate import get_or_make
from bench.variants import make_variant
from bench.pipeline import make_gate
from gate.gate import GateConfig
from scipy.stats import spearmanr

from scripts.laneA_paths import DATA, BENCH, CLAIMS, TABLES


def diag(g, t, claim, targets):
    A, B, L, R = claim["sender"], claim["receiver"], claim["ligand"], claim["receptor"]
    iA, iB, W = g._weights(A, B)
    E = W @ g.x.log_multi(L)[iA]
    Y = np.mean([g.x.z(x)[iB] for x in targets], 0)
    Yr = np.mean([g.x.raw(x)[iB] for x in targets], 0)
    sc = np.log(np.median(t.lib) / t.lib)[iB]
    tile = g.tile_id[iB]
    out = {}
    for k, v in (("gate", Y), ("raw", Yr), ("scale", sc)):
        r, p = g._perm_test(E, v, tile, 12345)
        out[f"rho_{k}"], out[f"p_{k}"] = r, p
    return out


def main(seed_set, tissues=None):
    fs = sorted(glob.glob(f"{CLAIMS}/{seed_set}_*.csv"))
    d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    d["passed"] = d.passed.astype(str).str.lower().eq("true")
    pl = d[d.list == "planted"]
    sel = pl[((pl.variant == "A4") & (pl.category == "true")) | ((pl.category != "true") & pl.passed) |
             ((pl.variant == "A0") & (pl.category == "true"))]
    if tissues:
        sel = sel[sel.tissue.isin(tissues)]
    rows, dens = [], []
    for tissue, gt in sel.groupby("tissue"):
        t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
        cnt = np.array([len(x) for x in t.tree.query_ball_point(t.xy[::25], 30.0)])
        rs = spearmanr(t.lib[::25], cnt)
        dens.append(dict(tissue=tissue, spearman_lib_vs_n30=float(rs.statistic), n_cells_sampled=int(len(cnt))))
        W = {}
        for s, gs in gt.groupby("seed"):
            base = get_or_make(t, int(s), None, root=BENCH)
            for v, gv in gs.groupby("variant"):
                inst = make_variant(t, base, v, int(s))
                g = make_gate(t, inst, GateConfig()); g._W = W
                tg = {(x["ligand"], x["receptor"]): x["targets"] for x in inst.database}
                for _, r in gv.iterrows():
                    c = dict(sender=r.sender, receiver=r.receiver, ligand=r.ligand, receptor=r.receptor)
                    rows.append(dict(tissue=tissue, seed=int(s), variant=v, category=r.category, mechanism=r.mechanism,
                                     **c, passed=bool(r.passed), **diag(g, t, c, tg[(r.ligand, r.receptor)])))
            print(tissue, s, len(rows), flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(f"{TABLES}/laneA_{seed_set}_confound.csv", index=False)
    dens = pd.DataFrame(dens)
    dens.to_csv(f"{TABLES}/laneA_{seed_set}_lib_density.csv", index=False)
    out["raw_supports"] = (out.p_raw < 0.05) & (out.rho_raw >= 0.02)
    summ = out.groupby(["variant", "category", "passed"]).agg(
        n=("rho_gate", "size"), rho_gate_med=("rho_gate", "median"), rho_raw_med=("rho_raw", "median"),
        rho_scale_med=("rho_scale", "median"), frac_raw_supports=("raw_supports", "mean")).reset_index()
    print(dens.to_string(index=False)); print(summ.round(3).to_string(index=False))
    summ.to_csv(f"{TABLES}/laneA_{seed_set}_confound_summary.csv", index=False)
    return out, summ, dens


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2].split(",") if len(sys.argv) > 2 else None)
