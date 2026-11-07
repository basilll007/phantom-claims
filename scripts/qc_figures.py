"""Counts table + spatial plots for processed tissues. Usage: python scripts/qc_figures.py tissue1 tissue2 ..."""
import sys, json
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = []
for t in sys.argv[1:]:
    a = ad.read_h5ad(f"data/processed/{t}.h5ad")
    vc = a.obs["cell_type"].value_counts()
    for k, v in vc.items():
        rows.append(dict(tissue=t, cell_type=k, n_cells=int(v)))
    qc = json.loads(a.uns["qc"])
    rows.append(dict(tissue=t, cell_type="__TOTAL__", n_cells=qc["cells_qc"]))
    types = vc.index.tolist()
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(9, 6.5))
    for i, k in enumerate(types):
        m = (a.obs["cell_type"] == k).values
        ax.scatter(a.obs.x[m], a.obs.y[m], s=0.05, color="lightgrey" if k == "Unlabeled" else cmap(i % 20),
                   label=f"{k} ({m.sum()})", rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    ax.set_title(f"{t}: {qc['cells_qc']} cells, {qc['genes']} genes")
    ax.legend(markerscale=40, fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    fig.tight_layout(); fig.savefig(f"figures/spatial_{t}.png", dpi=150); plt.close(fig)
pd.DataFrame(rows).to_csv("results/tables/cell_counts.csv", index=False)
print(pd.DataFrame(rows).pivot_table(index="cell_type", columns="tissue", values="n_cells").to_string())
