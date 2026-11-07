"""Intersect LIANA+ consensus resource (OmniPath-derived; complexes '_'-joined) with each panel.
A pair is kept iff every subunit of ligand and receptor is in the panel. -> data/processed/lr_pairs_<tissue>.csv"""
import sys, json
import pandas as pd, anndata as ad
d = pd.read_csv("data/raw/lr/liana_omni_resource.csv")
d = d[d.resource == "consensus"][["source_genesymbol", "target_genesymbol"]].drop_duplicates()
d.columns = ["ligand", "receptor"]
summary = {}
for t in sys.argv[1:]:
    g = set(ad.read_h5ad(f"data/processed/{t}.h5ad", backed="r").var_names)
    ok = lambda s: all(x in g for x in s.split("_"))
    m = d[d.ligand.map(ok) & d.receptor.map(ok)].reset_index(drop=True)
    m.to_csv(f"data/processed/lr_pairs_{t}.csv", index=False)
    summary[t] = dict(panel_genes=len(g), consensus_pairs=len(d), pairs_in_panel=len(m),
                      has_CD274_PDCD1=bool(((m.ligand == "CD274") & (m.receptor == "PDCD1")).any()))
print(json.dumps(summary, indent=1))
json.dump(summary, open("results/tables/lr_pool.json", "w"), indent=1)
