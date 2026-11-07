"""Exploratory (real data): are ligand+ sender / receptor+ receiver contacts within 30 um DEPLETED relative to chance?
Null: permute ligand positivity among sender-type cells and receptor positivity among receiver-type cells within
the same 250 um tile (1000 permutations). Two-sided view: observed/expected ratio and lower-tail p."""
import json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from bench.tissue import load_tissue
from bench.pipeline import log_result

claims = pd.read_csv("results/tables/real_breast_gate.csv")
claims = claims[claims.source.isin(["PD-L1/PD-1 case", "positive_control_contact"])][["sender", "receiver", "ligand", "receptor", "source"]].drop_duplicates()
rng = np.random.default_rng(0)
out = []
for rep in ("rep1", "rep2"):
    t = load_tissue(f"breast_{rep}")
    g = {x: i for i, x in enumerate(t.genes)}
    X = t.counts.tocsc()
    tile = (np.floor(t.xy[:, 0] / 250) * 100000 + np.floor(t.xy[:, 1] / 250)).astype(np.int64)
    for c in claims.itertuples():
        iA = np.flatnonzero(t.types == c.sender); iB = np.flatnonzero(t.types == c.receiver)
        D = cKDTree(t.xy[iB]).sparse_distance_matrix(cKDTree(t.xy[iA]), 30.0, output_type="coo_matrix")
        Lp = np.asarray(X[iA, g[c.ligand]].todense()).ravel() > 0
        Rp = np.asarray(X[iB, g[c.receptor]].todense()).ravel() > 0
        obs = int((Lp[D.col] & Rp[D.row]).sum())
        def perm_within(v, tl):
            v = v.copy()
            for u in np.unique(tl):
                m = np.flatnonzero(tl == u); v[m] = v[rng.permutation(m)]
            return v
        null = np.array([int((perm_within(Lp, tile[iA])[D.col] & perm_within(Rp, tile[iB])[D.row]).sum()) for _ in range(1000)])
        out.append(dict(replicate=rep, source=c.source, sender=c.sender, receiver=c.receiver, ligand=c.ligand, receptor=c.receptor,
                        n_ligand_pos=int(Lp.sum()), n_receptor_pos=int(Rp.sum()), observed=obs, expected=float(null.mean()),
                        ratio=float(obs / max(null.mean(), 1e-9)), p_depleted=float((1 + (null <= obs).sum()) / 1001),
                        p_enriched=float((1 + (null >= obs).sum()) / 1001)))
df = pd.DataFrame(out); df.to_csv("results/tables/real_contacts_vs_chance.csv", index=False)
print(df.round(3).to_string())
log_result(dict(id="real_breast_contacts_vs_chance", phase=7, split="real", tissue="breast_rep1+rep2", exploratory=True,
                null="tile-restricted (250 um) within-type permutation of positivity, 1000 perms", records=df.to_dict("records")))
