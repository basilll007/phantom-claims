"""Lane A exploratory: do A5 decoy-sender claims pass the gate when the decoy type co-localises with the true sender?

For each decoy claim (decoy D, receiver B, L, R) the co-localisation index is the Pearson correlation, over receivers
B, between the gate-kernel exposure from D cells and from true-sender A cells. The kernel is exp(-d/15), r=30, and
expression is unweighted, so the index is purely spatial. The index is joined to the planted-claim gate decision in
results/claims/laneA/<seed_set>_*_A5.csv.

  python -m scripts.laneA_decoy_coloc <seed_set>
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from bench.variants import exposure, MECHS

from scripts.laneA_paths import DATA, BENCH, CLAIMS, TABLES


def main(seed_set):
    rows = []
    fs = sorted(glob.glob(f"{CLAIMS}/{seed_set}_*_A5.csv"))
    by_tissue = {}
    for f in fs:
        d = pd.read_csv(f)
        by_tissue.setdefault(d.tissue.iloc[0], []).append(d)
    for tissue, ds in by_tissue.items():
        t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
        for d in ds:
            seed = int(d.seed.iloc[0])
            truth = json.load(open(glob.glob(f"{BENCH}/laneA/{seed_set}/{tissue}_s{seed}_A5_*_truth.json")[0]))
            pl = d[(d.list == "planted") & (d.category == "decoy_sender")]
            for p in [q for q in truth["truth"] if q.get("decoy_sender")]:
                iA = np.flatnonzero(t.types == p["sender"]); iD = np.flatnonzero(t.types == p["decoy_sender"])
                iB = np.flatnonzero(t.types == p["receiver"])
                EA = exposure(t.xy[iA], t.xy[iB], MECHS["exp30"]); ED = exposure(t.xy[iD], t.xy[iB], MECHS["exp30"])
                coloc = float(np.corrcoef(EA, ED)[0, 1]) if EA.std() > 0 and ED.std() > 0 else 0.0
                r = pl[(pl.sender == p["decoy_sender"]) & (pl.ligand == p["ligand"])].iloc[0]
                rows.append(dict(tissue=tissue, seed=seed, ligand=p["ligand"], sender=p["sender"],
                                 decoy=p["decoy_sender"], receiver=p["receiver"], mechanism=p["mechanism"],
                                 coloc=coloc, frac_B_exposed_decoy=float((ED > 0).mean()),
                                 passed=bool(r.passed), rho=r.rho, reasons=r.reasons))
    out = pd.DataFrame(rows)
    out.to_csv(f"{TABLES}/laneA_{seed_set}_decoy_coloc.csv", index=False)
    if len(out):
        out["coloc_bin"] = pd.cut(out.coloc, [-1, 0.1, 0.3, 0.5, 1.0])
        print(out.groupby("coloc_bin", observed=True).passed.agg(["mean", "size"]))
        print(out.groupby("mechanism").passed.agg(["mean", "size"]))
    return out


if __name__ == "__main__":
    main(sys.argv[1])
