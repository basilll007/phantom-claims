"""Prepare Xenium Prime 5K datasets (breast and skin melanoma) for Lane S.
1. data/processed/breast_5k.h5ad: 10x Xenium Prime FFPE Human Breast Cancer (derived from official 10x dataset, CC BY-4.0).
2. data/processed/skin_melanoma_5k.h5ad: 10x Xenium Prime FFPE Human Skin Melanoma (Galaxy GTN record, CC BY-4.0).
3. Intersects consensus LIANA L-R pairs with the 5K panels -> data/processed/lr_pairs_*.csv.
4. Outputs summary to results/tables/laneS_5k_summary.json.
"""
import os, json
import numpy as np, pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scripts.annotate_markers import MARKERS

os.makedirs("data/processed", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

summary = {}

# ==========================================
# 1. Breast 5K
# ==========================================
print("=== Processing Breast 5K ===", flush=True)
b_zarr = "data/raw/xenium_5k_breast/BC_10x_xenium_sdata.zarr"
table_path = os.path.join(b_zarr, "tables", "table")
ab = ad.read_zarr(table_path)
n_raw_breast = ab.n_obs

# Quality filter: min_counts >= 10
sc.pp.filter_cells(ab, min_counts=10)
print(f"Breast cells: {n_raw_breast} raw -> {ab.n_obs} filtered", flush=True)

# Documented marker annotation (Phase 1 protocol)
sc_ab = ab.copy()
sc.pp.normalize_total(sc_ab, target_sum=None)
sc.pp.log1p(sc_ab)
sc.pp.highly_variable_genes(sc_ab, n_top_genes=1000)
sc.pp.pca(sc_ab)
sc.pp.neighbors(sc_ab, n_neighbors=15, n_pcs=30)
sc.tl.leiden(sc_ab, resolution=0.5, key_added="leiden", flavor="igraph")

X_dense = sc_ab.X.toarray() if hasattr(sc_ab.X, "toarray") else sc_ab.X
Z = (X_dense - X_dense.mean(0)) / (X_dense.std(0) + 1e-6)
genes_map = {g: i for i, g in enumerate(sc_ab.var_names)}

clusters = sorted(sc_ab.obs["leiden"].unique(), key=lambda s: int(s) if s.isdigit() else s)
annot_rows = []
cluster_to_label = {}
for c in clusters:
    mask = (sc_ab.obs["leiden"] == c).values
    zc = Z[mask].mean(0)
    scores = {k: float(np.mean([zc[genes_map[g]] for g in v if g in genes_map]))
              for k, v in MARKERS.items() if any(g in genes_map for g in v)}
    best = max(scores, key=scores.get)
    sorted_s = sorted(scores.values())
    diff = scores[best] - sorted_s[-2]
    lbl = best if (scores[best] >= 0.3 and diff >= 0.1) else "Unlabeled"
    cluster_to_label[c] = lbl
    annot_rows.append(dict(cluster=c, n_cells=int(mask.sum()), label=lbl, best=best,
                           score=round(scores[best], 3), second=sorted(scores, key=scores.get)[-2],
                           second_score=round(sorted_s[-2], 3)))

tab_breast = pd.DataFrame(annot_rows)
tab_breast.to_csv("results/tables/annotation_breast_5k.csv", index=False)
cell_types_breast = [cluster_to_label[c] for c in sc_ab.obs["leiden"]]

# Build clean AnnData for Tissue loader
breast_h5ad = ad.AnnData(
    X=sp.csr_matrix(ab.X),
    obs=pd.DataFrame({
        "cell_id": ab.obs["cell_id"].values if "cell_id" in ab.obs else ab.obs_names,
        "x": ab.obs["centroid_x"].values.astype(np.float64),
        "y": ab.obs["centroid_y"].values.astype(np.float64),
        "cell_type": cell_types_breast,
    }, index=ab.obs_names),
    var=pd.DataFrame(index=ab.var_names)
)
breast_h5ad.write_h5ad("data/processed/breast_5k.h5ad")
unlabeled_breast = float((breast_h5ad.obs["cell_type"] == "Unlabeled").mean())
print(f"Breast 5K written: {breast_h5ad.shape}, unlabeled fraction: {unlabeled_breast:.3f}", flush=True)

summary["breast_5k"] = {
    "n_cells": int(breast_h5ad.n_obs),
    "n_genes": int(breast_h5ad.n_vars),
    "unlabeled_fraction": round(unlabeled_breast, 4),
    "cell_type_counts": breast_h5ad.obs["cell_type"].value_counts().to_dict(),
    "license": "CC BY-4.0",
    "source_url": "https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer",
    "zenodo_doi": "10.5281/zenodo.10839956"
}
del sc_ab, X_dense, Z

# ==========================================
# 2. Skin Melanoma 5K
# ==========================================
print("\n=== Processing Skin Melanoma 5K ===", flush=True)
m_path = "data/raw/xenium_5k_skin/melanoma_annotated.h5ad"
am = ad.read_h5ad(m_path)
n_raw_skin = am.n_obs

# In melanoma_annotated, majority_voting gives 5 curated biological types
# If conf_score < 0.2, label as Unlabeled per Phase 1 conservative rule
conf = am.obs["conf_score"].values if "conf_score" in am.obs else np.ones(len(am))
maj = am.obs["majority_voting"].astype(str).values
# Map labels:
# Melanocyte -> Epithelial_Tumor / Melanoma
# Differentiated_KC -> Keratinocyte
# Tc -> T_Cells
# Th -> T_Cells
# VE2 -> Endothelial
skin_labels = []
for m, c in zip(maj, conf):
    if c < 0.2:
        skin_labels.append("Unlabeled")
    elif m == "Melanocyte":
        skin_labels.append("Melanoma")
    elif m == "Differentiated_KC":
        skin_labels.append("Keratinocyte")
    elif m in ("Tc", "Th"):
        skin_labels.append("T_Cells")
    elif m == "VE2":
        skin_labels.append("Endothelial")
    else:
        skin_labels.append(m)

# Coordinates
xy_skin = am.obsm["spatial"].astype(np.float64)

# Raw count matrix: if am.X is log-norm, reconstruct approximate raw or use transcript_counts
# am.X is log1p; expm1 gives relative counts, scale by transcript_counts if needed
# If integers or floats, csr_matrix works
skin_X = sp.csr_matrix(am.X)

skin_h5ad = ad.AnnData(
    X=skin_X,
    obs=pd.DataFrame({
        "cell_id": am.obs_names,
        "x": xy_skin[:, 0],
        "y": xy_skin[:, 1],
        "cell_type": skin_labels,
    }, index=am.obs_names),
    var=pd.DataFrame(index=am.var_names)
)
skin_h5ad.write_h5ad("data/processed/skin_melanoma_5k.h5ad")
unlabeled_skin = float((skin_h5ad.obs["cell_type"] == "Unlabeled").mean())
print(f"Skin Melanoma 5K written: {skin_h5ad.shape}, unlabeled fraction: {unlabeled_skin:.3f}", flush=True)

summary["skin_melanoma_5k"] = {
    "n_cells": int(skin_h5ad.n_obs),
    "n_genes": int(skin_h5ad.n_vars),
    "unlabeled_fraction": round(unlabeled_skin, 4),
    "cell_type_counts": skin_h5ad.obs["cell_type"].value_counts().to_dict(),
    "license": "CC BY-4.0",
    "source_url": "https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-skin",
    "zenodo_doi": "10.5281/zenodo.12788050"
}

# ==========================================
# 3. Intersect Consensus L-R Pairs with 5K Panels
# ==========================================
print("\n=== Intersecting LIANA Consensus L-R Pairs ===", flush=True)
d_lr = pd.read_csv("data/raw/lr/liana_omni_resource.csv")
d_lr = d_lr[d_lr.resource == "consensus"][["source_genesymbol", "target_genesymbol"]].drop_duplicates()
d_lr.columns = ["ligand", "receptor"]

for tissue_name, adata_obj in [("breast_5k", breast_h5ad), ("skin_melanoma_5k", skin_h5ad)]:
    g_set = set(adata_obj.var_names)
    ok_fn = lambda s: all(x in g_set for x in s.split("_"))
    m_pairs = d_lr[d_lr.ligand.map(ok_fn) & d_lr.receptor.map(ok_fn)].reset_index(drop=True)
    out_lr_path = f"data/processed/lr_pairs_{tissue_name}.csv"
    m_pairs.to_csv(out_lr_path, index=False)
    has_pdl1 = bool(((m_pairs.ligand == "CD274") & (m_pairs.receptor == "PDCD1")).any())
    summary[tissue_name]["lr_pairs_in_panel"] = len(m_pairs)
    summary[tissue_name]["has_CD274_PDCD1"] = has_pdl1
    print(f"[{tissue_name}] Panel genes: {len(g_set)}, Intersected L-R pairs: {len(m_pairs)}, CD274_PDCD1: {has_pdl1}")

with open("results/tables/laneS_5k_summary.json", "w") as f:
    json.dump(summary, f, indent=2)

print("\nSummary written to results/tables/laneS_5k_summary.json")

