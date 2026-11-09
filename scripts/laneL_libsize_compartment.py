"""Library size against local density, by broad compartment (descriptive; real sections, no planting, all cells).
Supports the appendix text next to Figure figA5 (paper_figures.figA5); compartments as in paper_figures.COMPARTMENT.

    python scripts/laneL_libsize_compartment.py                                    # writes the table
    PC_RESULTS=results/results.jsonl python scripts/laneL_libsize_compartment.py --log
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from bench.tissue import load_tissue

COMPARTMENT = {  # keep in sync with scripts/paper_figures.py
    "epithelial_tumour": ["Invasive_Tumor", "Prolif_Invasive_Tumor", "DCIS_1", "DCIS_2", "Myoepi_ACTA2+",
                          "Myoepi_KRT15+", "Epithelial_Tumor", "Acinar", "Endocrine", "Ciliated"],
    "immune": ["CD4+_T_Cells", "CD8+_T_Cells", "B_Cells", "Macrophages_1", "Macrophages_2", "Macrophages", "T_Cells",
               "IRF7+_DCs", "LAMP3+_DCs", "Dendritic_Cells", "Mast_Cells", "Plasma_Cells"],
    "stroma_vessels": ["Stromal", "Fibroblast", "Endothelial", "Perivascular-Like", "SMC_Pericyte"],
}


def table():
    rows = []
    for n in ("breast_rep1", "breast_rep2", "lung", "pancreas"):
        t = load_tissue(n)
        d = cKDTree(t.xy).query_ball_point(t.xy, 30.0, return_length=True) - 1
        rl = t.lib / np.median(t.lib)
        comp = np.full(t.n, "other", dtype=object)
        for k, ty in COMPARTMENT.items():
            comp[np.isin(t.types, ty)] = k
        for c in ["all"] + list(COMPARTMENT) + ["other"]:
            m = np.ones(t.n, bool) if c == "all" else comp == c
            rows.append(dict(tissue=n, compartment=c, n_cells=int(m.sum()), median_rel_lib=float(np.median(rl[m])),
                             median_n30=float(np.median(d[m])),
                             spearman_lib_vs_n30=float(spearmanr(t.lib[m], d[m]).statistic)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    T = table()
    T.to_csv("results/tables/laneL_libsize_density_by_compartment.csv", index=False)
    print(T.round(3).to_string(index=False))
    if "--log" in sys.argv:
        assert os.environ.get("PC_RESULTS"), "set PC_RESULTS to the shared results.jsonl"
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import log_result as lr
        lr.RESULTS_FILE = os.environ["PC_RESULTS"]
        rid = "laneL_libsize_density_by_compartment"
        if any(f'"id": "{rid}"' in l for l in open(lr.RESULTS_FILE, encoding="utf-8")):
            print(f"{rid} already logged; not logging again")
        else:
            lr.log_result(dict(id=rid, phase="P2", hypothesis="descriptive (failure map: normalisation)",
                               tissues=["breast_rep1", "breast_rep2", "lung", "pancreas"], seeds=[],
                               table="results/tables/laneL_libsize_density_by_compartment.csv",
                               note="Real sections, no planting, all cells. median_rel_lib = median library size / "
                                    "section median; n30 = other cells within 30 um; compartments are display "
                                    "groupings of the published / marker-based cell types (paper_figures.COMPARTMENT).",
                               rows=T.to_dict("records")),
                          lane="L", agent_id="claude-opus-5.5", seed_set="dev", labelled="exploratory")
            print(f"logged {rid}")
