"""Documented marker-based annotation of 10x XOA graph-based clusters (no published labels available).
Each XOA graphclust cluster gets the lineage whose marker-set score (mean z-scored log-normalised expression,
cluster mean) is highest; score < 0.3 or margin to runner-up < 0.1 -> 'Unlabeled'. Clusters with the same lineage are merged.
Writes data/processed/annot_<tissue>.csv (Barcode, Cluster, xoa_cluster, score) and a marker table.
Usage: python scripts/annotate_markers.py <tissue>
"""
import io, json, sys, tarfile
import numpy as np, pandas as pd, scanpy as sc

MARKERS = {
    "Epithelial_Tumor": ["EPCAM", "KRT7", "KRT20", "CFTR", "EHF", "TMC5", "SPDEF", "CAPN8", "GPRC5A", "AGR3", "MET"],
    "Proliferating": ["MKI67", "TOP2A", "CDK1", "CCNB2", "UBE2C", "CENPF", "PCNA"],
    "Fibroblast": ["PDGFRA", "DPT", "FBLN1", "MFAP5", "SFRP2", "SFRP4", "OGN", "THBS2", "COL5A2", "VCAN", "ASPN", "PCOLCE"],
    "SMC_Pericyte": ["ACTA2", "ACTG2", "MYH11", "DES", "CNN1", "MYLK", "RERGL", "HIGD1B", "PDGFRB"],
    "Endothelial": ["PECAM1", "VWF", "CLEC14A", "EGFL7", "ERG", "ADGRL4", "RAMP2", "SOX17", "SOX18", "ECSCR"],
    "Lymphatic_Endo": ["LYVE1", "PROX1", "MMRN1", "TFPI"],
    "T_Cells": ["CD3D", "CD3E", "CD2", "TRAC", "CD247", "IL7R", "CD8A", "GZMK", "CD28"],
    "NK_Cells": ["NKG7", "GNLY", "KLRD1", "PRF1", "FGFBP2"],
    "B_Cells": ["MS4A1", "CD19", "CD79A", "BANK1", "TCL1A"],
    "Plasma_Cells": ["MZB1", "TNFRSF17", "DERL3", "FKBP11", "PRDM1", "TENT5C"],
    "Macrophages": ["CD68", "CD163", "MRC1", "MARCO", "VSIG4", "FCGR1A", "MS4A4A", "MS4A6A", "TREM2", "AIF1", "MPEG1"],
    "Monocytes_Neutrophils": ["FCN1", "S100A12", "CD300E", "LILRA5", "AQP9", "CD14"],
    "Dendritic_Cells": ["CD1C", "CLEC10A", "FCER1A", "LAMP3", "IRF8", "CD1A", "CD1E"],
    "pDC": ["LILRA4", "IL3RA", "PLD4"],
    "Mast_Cells": ["CPA3", "KIT", "MS4A2", "HPGDS"],
}
TISSUE_MARKERS = {
    "pancreas": {"Acinar": ["AMY2A", "CFB"], "Endocrine": ["INS", "GCG", "SST", "PPY", "CHGA", "PCSK2", "SCGN", "GHRL"]},
    "lung": {"Alveolar_Epi": ["AGER", "SFTA2", "PEBP4", "CA4"],
             "Ciliated": ["FOXJ1", "CCDC39", "C20orf85", "SNTN", "CFAP53", "CCDC78", "C1orf194", "DNAAF1"]},
}


def main(tissue):
    raw = f"data/raw/{tissue}"
    a = sc.read_10x_h5(f"{raw}/cell_feature_matrix.h5")
    a.var_names_make_unique()
    with tarfile.open(f"{raw}/analysis.tar.gz") as tf:
        m = [x for x in tf.getmembers() if x.name.endswith("gene_expression_graphclust/clusters.csv")][0]
        cl = pd.read_csv(io.BytesIO(tf.extractfile(m).read()))
    cl["Barcode"] = cl["Barcode"].astype(str)
    a = a[a.obs_names.isin(cl["Barcode"])].copy()
    a.obs["xoa"] = cl.set_index("Barcode").loc[a.obs_names, "Cluster"].astype(str).values
    sc.pp.normalize_total(a, target_sum=None)
    sc.pp.log1p(a)
    X = a.X.toarray() if hasattr(a.X, "toarray") else a.X
    Z = (X - X.mean(0)) / (X.std(0) + 1e-6)
    genes = {g: i for i, g in enumerate(a.var_names)}
    sets = dict(MARKERS); sets.update(TISSUE_MARKERS.get(tissue, {}))
    clusters = sorted(a.obs["xoa"].unique(), key=lambda s: int(s) if s.isdigit() else s)
    rows = []
    for c in clusters:
        mask = (a.obs["xoa"] == c).values
        zc = Z[mask].mean(0)
        sc_ = {k: float(np.mean([zc[genes[g]] for g in v if g in genes])) for k, v in sets.items()
               if any(g in genes for g in v)}
        best = max(sc_, key=sc_.get)
        rows.append(dict(xoa_cluster=c, n_cells=int(mask.sum()), label=best if (sc_[best] >= 0.3 and sc_[best] - sorted(sc_.values())[-2] >= 0.1) else "Unlabeled",
                         best=best, score=round(sc_[best], 3),
                         second=sorted(sc_, key=sc_.get)[-2], second_score=round(sorted(sc_.values())[-2], 3)))
    tab = pd.DataFrame(rows)
    tab.to_csv(f"results/tables/annotation_{tissue}.csv", index=False)
    lab = dict(zip(tab.xoa_cluster, tab.label))
    out = pd.DataFrame(dict(Barcode=a.obs_names, Cluster=[lab[c] for c in a.obs["xoa"]], xoa_cluster=a.obs["xoa"].values))
    out.to_csv(f"data/processed/annot_{tissue}.csv", index=False)
    print(tab.to_string())
    print(out.Cluster.value_counts().to_string())


if __name__ == "__main__":
    main(sys.argv[1])
