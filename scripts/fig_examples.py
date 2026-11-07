"""Figure: (a) planted TRUE pair (dev seed 0), (b) planted PHANTOM pair, (c) real CD274+ / PDCD1+ cells (Rep1)."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bench.tissue import load_tissue
from bench.generate import get_or_make

t = load_tissue("breast_rep1")
inst = get_or_make(t, 0)
gi = inst.gene_index()
C = inst.counts.tocsc()
col = lambda g: np.asarray(C[:, gi[g]].todense()).ravel()
tp = next(p for p in inst.truth if p["category"] == "true" and p["n_receivers_exposed"] > 300)
ph = next(p for p in inst.truth if p["category"] == "phantom" and p["n_senders_far"] > 50)
fig, axs = plt.subplots(1, 3, figsize=(13, 4.3))
# (a) TRUE: crop around the densest sender-receiver area
iA = np.flatnonzero(t.types == tp["sender"]); iB = np.flatnonzero(t.types == tp["receiver"])
H, xe, ye = np.histogram2d(t.xy[iB, 0], t.xy[iB, 1], bins=[np.arange(0, 7600, 500), np.arange(0, 5600, 500)])
bx, by = np.unravel_index(H.argmax(), H.shape)
x0, y0 = xe[bx] - 250, ye[by] - 250
win = lambda idx: idx[(t.xy[idx, 0] > x0) & (t.xy[idx, 0] < x0 + 1000) & (t.xy[idx, 1] > y0) & (t.xy[idx, 1] < y0 + 1000)]
ax = axs[0]
allw = win(np.arange(t.n))
ax.scatter(t.xy[allw, 0], t.xy[allw, 1], s=1, c="#dddddd", rasterized=True)
a = win(iA); a = a[col(tp["ligand"])[a] > 0]
b = win(iB); tgt = sum(col(g) for g in tp["responders"])
ax.scatter(t.xy[a, 0], t.xy[a, 1], s=6, c="#dd8452", label=f"{tp['sender']} ligand+", rasterized=True)
sc = ax.scatter(t.xy[b, 0], t.xy[b, 1], s=6, c=tgt[b], cmap="viridis", vmin=0, vmax=np.percentile(tgt[b], 98) + 1e-9,
                label=f"{tp['receiver']} (colour = target counts)", rasterized=True)
ax.set_title(f"(a) planted TRUE: targets rise near senders", fontsize=9); ax.legend(fontsize=6, loc="lower left", markerscale=2)
# (b) PHANTOM
ax = axs[1]
ax.scatter(t.xy[::7, 0], t.xy[::7, 1], s=0.2, c="#dddddd", rasterized=True)
pa = np.flatnonzero((t.types == ph["sender"]) & (inst.added.tocsc()[:, gi[ph["ligand"]]].toarray().ravel() > 0))
pb = np.flatnonzero((t.types == ph["receiver"]) & (inst.added.tocsc()[:, gi[ph["receptor"]]].toarray().ravel() > 0))
ax.scatter(t.xy[pa, 0], t.xy[pa, 1], s=3, c="#c44e52", label=f"{ph['sender']} ligand-planted", rasterized=True)
ax.scatter(t.xy[pb, 0], t.xy[pb, 1], s=3, c="#4c72b0", label=f"{ph['receiver']} receptor-planted", rasterized=True)
ax.set_title("(b) planted PHANTOM: expressing cells never meet", fontsize=9); ax.legend(fontsize=6, loc="lower left", markerscale=3)
# (c) real CD274+ LAMP3+ DCs vs PDCD1+ T cells (the expression tool's top PD-L1/PD-1 type pairs)
ax = axs[2]
g = {x: i for i, x in enumerate(t.genes)}
X = t.counts.tocsc()
pos = lambda gene: np.asarray(X[:, g[gene]].todense()).ravel() > 0
dc = np.flatnonzero((t.types == "LAMP3+_DCs") & pos("CD274"))
tc = np.flatnonzero(np.isin(t.types, ["CD4+_T_Cells", "CD8+_T_Cells"]) & pos("PDCD1"))
ax.scatter(t.xy[::7, 0], t.xy[::7, 1], s=0.2, c="#dddddd", rasterized=True)
ax.scatter(t.xy[dc, 0], t.xy[dc, 1], s=6, c="#c44e52", label=f"CD274+ LAMP3+ DCs (n={len(dc)})", rasterized=True)
ax.scatter(t.xy[tc, 0], t.xy[tc, 1], s=6, c="#4c72b0", label=f"PDCD1+ CD4/CD8 T cells (n={len(tc)})", rasterized=True)
ax.set_title("(c) real Rep1: few contacts, at chance level", fontsize=9); ax.legend(fontsize=6, loc="lower left", markerscale=2)
for ax in axs:
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig("figures/fig_examples.png", dpi=200)
print("pair TRUE", tp["pair_id"], tp["sender"], tp["receiver"], "PHANTOM", ph["pair_id"], ph["sender"], ph["receiver"])
