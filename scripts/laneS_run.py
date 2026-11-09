"""Lane S runner: Phase 2 Strengthening Experiments.
Evaluates:
1. Frozen gate on fresh seed set N3 (replication of H4 across A0-A5, pooled K=10).
2. Spatial baselines on the same candidate lists (SpatialDM, Squidpy ligrec, Squidpy ligrec+nhood_enrichment).
3. Standalone top-K lists for each spatial baseline.
4. Per-claim class pass rates for the base-rate law.

Outputs per (tissue, seed, variant) CSV under results/claims/laneS/.
"""
import copy, gc, json, os, subprocess, sys, time
from dataclasses import asdict
import numpy as np, pandas as pd
from scipy import sparse
import torch
import anndata as ad
import squidpy as sq
import spatialdm as sdm

from bench.tissue import load_tissue
from bench.generate import get_or_make, generate
from bench.pipeline import tool_lists, make_gate, evaluate_records, decide_list, claim_key
from bench.variants import VARIANTS, make_variant, label_claim_v, planted_claims, save_truth
from gate.gate import GateConfig
from gate.status import status_of, readout_of
from scripts.laneS_paths import DATA, BENCH, CLAIMS as OUT, TABLES

SEEDS = {
    "N3": [300, 301, 302, 303, 304, 305, 306, 307, 308, 309],
    "dev": [0, 1, 2, 3, 4],
}
FORBIDDEN = {0, 1, 2, 3, 4, 100, 101, 102, 103, 104, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209}
FROZEN = dict(r=30.0, tile=250.0, n_perm=1000, alpha_half=0.10, n_min=50, rho_min=0.02, q=0.05)
KS = (5, 10, 20)
TISSUES = ["breast_rep1", "breast_rep2", "pancreas", "lung"]

GIT = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()


def frozen_cfg():
    cfg = GateConfig()
    assert {k: getattr(cfg, k) for k in FROZEN} == FROZEN, "gate config drifted from frozen defaults"
    sel = json.load(open("results/tables/gate_selected_dev.json"))
    assert all(float(getattr(cfg, k)) == float(v) for k, v in sel.items()), "frozen != Phase-1 selected config"
    return cfg


def load_base(t, seed):
    try:
        return get_or_make(t, seed, None, root=BENCH)
    except Exception as e:
        print(f"[warn] base cache read failed ({e!r}); regenerating in memory", flush=True)
        return generate(t, seed, None)


def mech_of(inst):
    return {p["ligand"]: p.get("mechanism", "") for p in inst.truth if p["category"] == "true"}


def precompute_spatial_structures(t):
    """Precompute tissue-level spatial matrices (constant across seeds and variants)."""
    t0 = time.time()
    obs = pd.DataFrame({"cell_type": pd.Categorical(t.types)}, index=[f"cell_{i}" for i in range(t.n)])
    adata_base = ad.AnnData(obs=obs)
    adata_base.obsm["spatial"] = t.xy.astype(np.float64)
    
    # Filter out unlabeled
    valid_mask = (adata_base.obs["cell_type"] != "Unlabeled").values
    adata_sub = adata_base[valid_mask].copy()
    adata_sub.obs["cell_type"] = adata_sub.obs["cell_type"].cat.remove_unused_categories()
    
    # 1. SpatialDM weight matrix
    sdm.weight_matrix(adata_sub, l=30.0, cutoff=0.1, n_nearest_neighbors=6, single_cell=True)
    spatial_W = adata_sub.obsp["weight"]
    spatial_knn = adata_sub.obsp["nearest_neighbors"]
    
    # 2. Squidpy spatial neighborhood enrichment
    sq.gr.spatial_neighbors(adata_sub, coord_type="generic", radius=30.0)
    nhood_res = sq.gr.nhood_enrichment(adata_sub, cluster_key="cell_type", n_perms=1000, show_progress_bar=False, copy=True)
    cluster_names = list(adata_sub.obs["cell_type"].cat.categories)
    cluster_to_idx = {c: i for i, c in enumerate(cluster_names)}
    
    print(f"[{t.name}] spatial precomputation done in {time.time()-t0:.1f}s", flush=True)
    return valid_mask, spatial_W, spatial_knn, nhood_res.zscore, cluster_to_idx


def run_spatial_baselines(t, inst, valid_mask, spatial_W, spatial_knn, nhood_zscore, cluster_to_idx):
    """Run SpatialDM and Squidpy ligrec on this instance."""
    pairs = [(d["ligand"], d["receptor"]) for d in inst.database]
    obs = pd.DataFrame({"cell_type": pd.Categorical(t.types[valid_mask])}, index=[f"cell_{i}" for i in range(valid_mask.sum())])
    obs["cell_type"] = obs["cell_type"].cat.remove_unused_categories()
    
    # --- SpatialDM ---
    t_sdm0 = time.time()
    adata_sdm = ad.AnnData(X=sparse.csr_matrix(inst.counts)[valid_mask], obs=obs)
    adata_sdm.var_names = list(inst.syn_genes)
    adata_sdm.obsm["spatial"] = t.xy[valid_mask].astype(np.float64)
    adata_sdm.obsp["weight"] = spatial_W
    adata_sdm.obsp["nearest_neighbors"] = spatial_knn
    adata_sdm.uns["mean"] = "algebra"
    adata_sdm.uns["single_cell"] = True
    
    lig_df = pd.DataFrame([[p[0]] for p in pairs], columns=["Ligand1"], index=[f"pair_{i}" for i in range(len(pairs))])
    rec_df = pd.DataFrame([[p[1]] for p in pairs], columns=["Receptor1"], index=[f"pair_{i}" for i in range(len(pairs))])
    inter_df = pd.DataFrame({
        "interaction_name": [f"{p[0]}_{p[1]}" for p in pairs],
        "Ligand": [p[0] for p in pairs],
        "Receptor": [p[1] for p in pairs],
        "annotation": ["Secreted Signaling"] * len(pairs)
    }, index=[f"pair_{i}" for i in range(len(pairs))])
    adata_sdm.uns["geneInter"] = inter_df
    adata_sdm.uns["ligand"] = lig_df
    adata_sdm.uns["receptor"] = rec_df
    sdm.spatialdm_global(adata_sdm, method="z-score", n_perm=1000)
    sdm_time = time.time() - t_sdm0
    sdm_res = adata_sdm.uns["global_res"].set_index(["Ligand1", "Receptor1"])
    
    # --- Squidpy ligrec ---
    t_sq0 = time.time()
    adata_sq = adata_sdm.copy()
    adata_sq.var_names = [g.replace("_", "") for g in inst.syn_genes]
    sq_pairs = [(p[0].replace("_", ""), p[1].replace("_", "")) for p in pairs]
    interactions = pd.DataFrame(sq_pairs, columns=["source", "target"])
    res_ligrec = sq.gr.ligrec(
        adata_sq,
        cluster_key="cell_type",
        interactions=interactions,
        threshold=0.01,
        n_perms=1000,
        use_raw=False,
        show_progress_bar=False,
        copy=True
    )
    sq_time = time.time() - t_sq0
    sq_means = res_ligrec["means"]
    sq_pvals = res_ligrec["pvalues"]
    
    return {
        "sdm_res": sdm_res,
        "sdm_time": sdm_time,
        "sq_means": sq_means,
        "sq_pvals": sq_pvals,
        "sq_time": sq_time,
        "nhood_zscore": nhood_zscore,
        "cluster_to_idx": cluster_to_idx,
    }


def make_spatial_top_k(t, inst, spatial_data, expr_df, Ks=(5, 10, 20)):
    """Generate standalone top-K lists for each spatial baseline."""
    sdm_res = spatial_data["sdm_res"]
    sq_means = spatial_data["sq_means"]
    sq_pvals = spatial_data["sq_pvals"]
    nhood_z = spatial_data["nhood_zscore"]
    c2i = spatial_data["cluster_to_idx"]
    
    # expr_df has columns: sender, receiver, ligand, receptor, fracL, fracR, score_cpdb, p_cpdb
    # Candidate pool: all (A, B, L, R) where A != B
    df = expr_df.copy()
    
    # Map SpatialDM stats
    sdm_z_list, sdm_p_list = [], []
    for _, row in df.iterrows():
        try:
            r = sdm_res.loc[(row.ligand, row.receptor)]
            sdm_z_list.append(float(r["z"]))
            sdm_p_list.append(float(r["z_pval"]))
        except KeyError:
            sdm_z_list.append(np.nan)
            sdm_p_list.append(1.0)
    df["sdm_z"] = sdm_z_list
    df["sdm_pval"] = sdm_p_list
    
    # Map Squidpy stats
    sq_m_list, sq_p_list, nh_z_list = [], [], []
    for _, row in df.iterrows():
        L_clean, R_clean = row.ligand.replace("_", ""), row.receptor.replace("_", "")
        A, B = row.sender, row.receiver
        try:
            m = float(sq_means.loc[(L_clean, R_clean), (A, B)])
            p = float(sq_pvals.loc[(L_clean, R_clean), (A, B)])
        except KeyError:
            m, p = 0.0, 1.0
        try:
            z = float(nhood_z[c2i[A], c2i[B]])
        except KeyError:
            z = 0.0
        sq_m_list.append(m)
        sq_p_list.append(p)
        nh_z_list.append(z)
    df["sq_mean"] = sq_m_list
    df["sq_pval"] = sq_p_list
    df["nhood_z"] = nh_z_list
    
    out_lists = {}
    
    # 1. SpatialDM: valid expression (frac >= 0.1) & sdm_pval < 0.05, ranked by sdm_z
    v_sdm = df[(df.fracL >= 0.1) & (df.fracR >= 0.1) & (df.sdm_pval < 0.05) & (df.sdm_z > 0)]
    v_sdm = v_sdm.sort_values(["sdm_z", "score_cpdb", "sender", "receiver", "ligand"], ascending=[False, False, True, True, True])
    for K in Ks:
        out_lists[("spatialdm", K)] = v_sdm.head(K)[["sender", "receiver", "ligand", "receptor"]].to_dict("records")
        
    # 2. Squidpy ligrec: sq_pval < 0.05 & frac >= 0.01, ranked by sq_mean
    v_sq = df[(df.fracL >= 0.01) & (df.fracR >= 0.01) & (df.sq_pval < 0.05)]
    v_sq = v_sq.sort_values(["sq_mean", "sender", "receiver", "ligand"], ascending=[False, True, True, True])
    for K in Ks:
        out_lists[("squidpy_ligrec", K)] = v_sq.head(K)[["sender", "receiver", "ligand", "receptor"]].to_dict("records")
        
    # 3. Squidpy spatial: sq_pval < 0.05 & nhood_z > 0, ranked by sq_mean
    v_sq_sp = df[(df.fracL >= 0.01) & (df.fracR >= 0.01) & (df.sq_pval < 0.05) & (df.nhood_z > 0)]
    v_sq_sp = v_sq_sp.sort_values(["sq_mean", "sender", "receiver", "ligand"], ascending=[False, True, True, True])
    for K in Ks:
        out_lists[("squidpy_spatial", K)] = v_sq_sp.head(K)[["sender", "receiver", "ligand", "receptor"]].to_dict("records")
        
    return out_lists, df


def eval_variant_full(t, inst, lists, cfg, shared_W, spatial_data):
    gate = make_gate(t, inst, cfg)
    gate._W = shared_W
    planted = planted_claims(inst)
    allc = [c for cl in lists.values() for c in cl] + [c for c, _, _ in planted]
    recs = evaluate_records(gate, allc)
    mech = mech_of(inst)
    
    sdm_res = spatial_data["sdm_res"]
    sq_means = spatial_data["sq_means"]
    sq_pvals = spatial_data["sq_pvals"]
    nhood_z = spatial_data["nhood_zscore"]
    c2i = spatial_data["cluster_to_idx"]
    
    rows = []
    def add_rows(listname, method, K, cl, cats=None):
        dec = decide_list(gate, recs, cl)
        for rank, (c, d) in enumerate(zip(cl, dec)):
            lab, sub = label_claim_v(inst, c)
            L, R = c["ligand"], c["receptor"]
            A, B = c["sender"], c["receiver"]
            
            # SpatialDM stats
            try:
                s_row = sdm_res.loc[(L, R)]
                sdm_z = float(s_row["z"])
                sdm_p = float(s_row["z_pval"])
                sdm_pass = bool(sdm_p < 0.05 and sdm_z > 0)
            except KeyError:
                sdm_z, sdm_p, sdm_pass = np.nan, 1.0, False
                
            # Squidpy stats
            L_clean, R_clean = L.replace("_", ""), R.replace("_", "")
            try:
                sq_m = float(sq_means.loc[(L_clean, R_clean), (A, B)])
                sq_p = float(sq_pvals.loc[(L_clean, R_clean), (A, B)])
                sq_pass = bool(sq_p < 0.05)
            except KeyError:
                sq_m, sq_p, sq_pass = 0.0, 1.0, False
                
            # Squidpy nhood stats
            try:
                nh_val = float(nhood_z[c2i[A], c2i[B]])
            except KeyError:
                nh_val = 0.0
            sq_sp_pass = bool(sq_pass and nh_val > 0)
            
            rows.append(dict(
                list=listname, method=method, K=K, rank=rank + 1, **c, label=lab,
                neg_type=sub if not lab else "", category=cats[rank] if cats else "",
                mechanism=mech.get(L, "") if (lab or (cats and cats[rank] == "decoy_sender")) else "",
                passed=bool(d["pass_"]), q_bh=d.get("q_bh"), rho=d.get("rho_within_tile"),
                p_perm=d.get("p_perm"), n_expr_pairs=d.get("n_expressing_pairs_within_r"),
                effect_d=d.get("effect_d"), reasons=";".join(d.get("reasons", [])),
                status=status_of(d), readout=readout_of(d),
                # Spatial baseline fields
                spatialdm_z=sdm_z, spatialdm_pval=sdm_p, spatialdm_pass=sdm_pass,
                squidpy_mean=sq_m, squidpy_pval=sq_p, squidpy_ligrec_pass=sq_pass,
                squidpy_nhood_z=nh_val, squidpy_spatial_pass=sq_sp_pass
            ))

    for (m, K), cl in lists.items():
        add_rows(f"{m}_K{K}", m, K, cl)
    add_rows("planted", "", 0, [c for c, _, _ in planted], cats=[k for _, k, _ in planted])
    return pd.DataFrame(rows)


def run(tissue, seed_set="N3", seeds=None, variants=None):
    seeds = seeds or SEEDS[seed_set]
    assert set(seeds) <= set(SEEDS[seed_set])
    if seed_set == "N3":
        assert not (set(seeds) & FORBIDDEN), "seed overlap with forbidden sets"
    variants = variants or list(VARIANTS)
    cfg = frozen_cfg()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    t = load_tissue(tissue, f"{DATA}/processed/{tissue}.h5ad")
    print(f"[{tissue}] loaded n={t.n} in {time.time() - t0:.0f}s", flush=True)
    
    # Precompute spatial structures once per tissue
    valid_mask, spatial_W, spatial_knn, nhood_zscore, cluster_to_idx = precompute_spatial_structures(t)
    shared_W = {}
    
    for s in seeds:
        todo = [v for v in variants if not os.path.exists(f"{OUT}/{seed_set}_{tissue}_s{s}_{v}.csv")]
        if not todo:
            print(f"[{tissue} s{s}] all variants done, skip", flush=True)
            continue
        ts = time.time()
        base = load_base(t, s)
        lists = None
        base_expr_df = None
        
        for v in todo:
            tv = time.time()
            inst = make_variant(t, base, v, s)
            
            # Expression tool lists
            if v == "A5":
                lv, expr_df = tool_lists(t, inst, Ks=KS, seed=s)
            else:
                if lists is None:
                    lists, base_expr_df = tool_lists(t, base, Ks=KS, seed=s)
                lv = lists
                expr_df = base_expr_df
            
            # Spatial baselines evaluation
            spatial_data = run_spatial_baselines(
                t, inst, valid_mask, spatial_W, spatial_knn, nhood_zscore, cluster_to_idx
            )
            
            # Standalone spatial tool lists
            sp_lists, _ = make_spatial_top_k(t, inst, spatial_data, expr_df, Ks=KS)
            
            # Combine all candidate lists
            all_lists = {**lv, **sp_lists}
            
            df = eval_variant_full(t, inst, all_lists, cfg, shared_W, spatial_data)
            df.insert(0, "variant", v)
            df.insert(0, "seed", s)
            df.insert(0, "tissue", tissue)
            df.insert(0, "seed_set", seed_set)
            df["git_hash"] = GIT
            df["runtime_sdm_s"] = spatial_data["sdm_time"]
            df["runtime_sq_s"] = spatial_data["sq_time"]
            
            save_truth(inst, f"{BENCH}/laneS/{seed_set}")
            tmp = f"{OUT}/{seed_set}_{tissue}_s{s}_{v}.csv.tmp"
            df.to_csv(tmp, index=False)
            os.replace(tmp, f"{OUT}/{seed_set}_{tissue}_s{s}_{v}.csv")
            
            k10 = df[(df.K == 10) & (df.method == "cpdb")]
            vram = torch.cuda.max_memory_allocated() / 1e6 if torch.cuda.is_available() else 0
            print(f"[{tissue} s{s} {v}] {time.time() - tv:.1f}s (SDM {spatial_data['sdm_time']:.1f}s, SQ {spatial_data['sq_time']:.1f}s) K10 rows={len(k10)} vram_peak={vram:.0f}MB", flush=True)
            del df, spatial_data, all_lists, inst
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            
        del base, lists, base_expr_df
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        print(f"[{tissue} s{s}] seed done in {time.time() - ts:.0f}s (total {time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        tissue_arg = sys.argv[2]
        seed_set_arg = sys.argv[3] if len(sys.argv) > 3 else "N3"
        seeds_arg = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 and sys.argv[4] else None
        vars_arg = sys.argv[5].split(",") if len(sys.argv) > 5 else None
        run(tissue_arg, seed_set_arg, seeds_arg, vars_arg)
