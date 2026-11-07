"""Process a Xenium tissue into data/processed/<tissue>.h5ad (gene-expression counts, um coordinates, cell types).
Usage: python scripts/process_xenium.py <tissue> <raw_dir> [annotation_csv]"""
import json, sys, time, os
import numpy as np, pandas as pd, scanpy as sc, anndata as ad, scipy.sparse as sp

def load_cells(raw):
    for f in ["cells.csv.gz", "cells.parquet"]:
        p = os.path.join(raw, f)
        if os.path.exists(p):
            return pd.read_csv(p) if f.endswith("gz") else pd.read_parquet(p)
    raise FileNotFoundError(raw)

def main(tissue, raw, annot=None, min_counts=10):
    a = sc.read_10x_h5(os.path.join(raw, "cell_feature_matrix.h5"), gex_only=False)
    a = a[:, a.var["feature_types"] == "Gene Expression"].copy()
    a.var_names_make_unique()
    cells = load_cells(raw)
    cells["cell_id"] = cells["cell_id"].astype(str)
    a.obs_names = a.obs_names.astype(str)
    cells = cells.set_index("cell_id").loc[a.obs_names]
    a.obs["x"] = cells["x_centroid"].values.astype(np.float32)
    a.obs["y"] = cells["y_centroid"].values.astype(np.float32)
    a.obs["cell_area"] = cells["cell_area"].values
    n0 = a.n_obs
    if annot:
        an = pd.read_csv(annot)
        an["Barcode"] = an["Barcode"].astype(float).astype(np.int64).astype(str)
        an = an.set_index("Barcode")["Cluster"]
        a.obs["cell_type"] = an.reindex(a.obs_names).fillna("Unlabeled").values
        a.uns["annotation_source"] = os.path.basename(annot)
    a.X = sp.csr_matrix(a.X, dtype=np.float32)
    a.obs["total_counts"] = np.asarray(a.X.sum(1)).ravel()
    a = a[a.obs["total_counts"] >= min_counts].copy()
    qc = dict(tissue=tissue, cells_raw=int(n0), cells_qc=int(a.n_obs), genes=int(a.n_vars),
              median_counts=float(np.median(a.obs.total_counts)),
              n_types=int(a.obs["cell_type"].nunique()) if "cell_type" in a.obs else None,
              x_range=[float(a.obs.x.min()), float(a.obs.x.max())], y_range=[float(a.obs.y.min()), float(a.obs.y.max())])
    a.uns["qc"] = json.dumps(qc)
    a.write_h5ad(f"data/processed/{tissue}.h5ad", compression="gzip")
    print(json.dumps(qc))
    return a

if __name__ == "__main__":
    main(*sys.argv[1:])
