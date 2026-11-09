"""Lane S item S3: Richer panels (Xenium Prime 5K) evaluation.
Runs:
(a) Semi-synthetic benchmark (A0-A5, seeds [300, 301, 302]) on breast_5k and skin_melanoma_5k with the frozen gate.
(b) Real-data controls protocol (Lane C style) on the 5K L-R pools (1,645 and 1,408 pairs).
(c) PD-L1/PD-1 (CD274 - PDCD1) contact and gate result on these tissues.

Reports passes AND fails to results/tables/laneS_s3_results.json and results/tables/laneS_s3_pdl1.csv.
"""
import copy, json, os, sys, time
import numpy as np, pandas as pd
import scipy.sparse as sp
import torch
import anndata as ad

from bench.tissue import Tissue, load_tissue, type_pair_contacts
from bench.generate import generate
from bench.pipeline import tool_lists, make_gate, decide_list
from bench.variants import VARIANTS, make_variant, label_claim_v, planted_claims
from gate.gate import Gate, GateConfig, ExprSource
from gate.status import status_of, readout_of
from baselines.expr_tools import score_all, top_k

SEEDS_5K = [300, 301, 302]
TISSUES_5K = ["breast_5k", "skin_melanoma_5k"]
VARIANTS_5K = ["A0", "A1", "A2a", "A2b", "A3", "A4", "A5", "A5m"]

os.makedirs("results/claims/laneS_5k", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)

sel = json.load(open("results/tables/gate_selected_dev.json"))
cfg = GateConfig(**sel)


def run_benchmark_5k(tissue_name):
    print(f"\n==========================================", flush=True)
    print(f"Running Semi-Synthetic Benchmark on {tissue_name}", flush=True)
    print(f"==========================================", flush=True)
    t = load_tissue(tissue_name)
    types = t.type_list(min_cells=10)
    print(f"[{tissue_name}] n_cells={t.n}, types ({len(types)})={types}", flush=True)

    results = []
    for s in SEEDS_5K:
        base = generate(t, seed=s)
        for v in VARIANTS_5K:
            inst = make_variant(t, base, v, s)
            lv, _ = tool_lists(t, inst, Ks=[10], seed=s)
            
            # Evaluate gate on cpdb and prod candidate lists
            for m in ["cpdb", "prod"]:
                cand = lv.get(f"{m}_10", [])
                gate = make_gate(t, inst, cfg)
                recs = gate.check(cand)
                
                n_pass = 0
                n_false_pass = 0
                n_true_pass = 0
                n_true_cand = 0
                
                for c, r in zip(cand, recs):
                    lab, _ = label_claim_v(inst, c)
                    passed = r["pass"]
                    if lab == 1:
                        n_true_cand += 1
                        if passed:
                            n_true_pass += 1
                            n_pass += 1
                    else:
                        if passed:
                            n_false_pass += 1
                            n_pass += 1
                
                fps = n_false_pass / max(n_pass, 1) if n_pass > 0 else 0.0
                ret = n_true_pass / max(n_true_cand, 1) if n_true_cand > 0 else 0.0
                
                results.append({
                    "tissue": tissue_name,
                    "seed": s,
                    "variant": v,
                    "method": m,
                    "n_cand": len(cand),
                    "n_pass": n_pass,
                    "n_true_pass": n_true_pass,
                    "n_false_pass": n_false_pass,
                    "false_pass_share": round(fps, 4),
                    "retention_true": round(ret, 4)
                })
        print(f"[{tissue_name}] seed {s} benchmark done", flush=True)
    
    df_res = pd.DataFrame(results)
    df_res.to_csv(f"results/tables/laneS_benchmark_{tissue_name}.csv", index=False)
    return df_res


def run_pdl1_analysis(tissue_name):
    print(f"\n=== PD-L1 / PD-1 Analysis on {tissue_name} ===", flush=True)
    t = load_tissue(tissue_name)
    types = t.type_list(min_cells=10)
    
    # Check CD274 and PDCD1 in panel
    g_set = set(t.genes)
    has_cd274 = "CD274" in g_set
    has_pdcd1 = "PDCD1" in g_set
    print(f"CD274 in panel: {has_cd274}, PDCD1 in panel: {has_pdcd1}")
    if not (has_cd274 and has_pdcd1):
        return []
    
    # Contact counts & enrichment
    _, C, E = type_pair_contacts(t, r=30.0, types=types)
    t_idx = {k: i for i, k in enumerate(types)}
    
    # Candidate cell type pairs
    test_pairs = []
    if tissue_name == "skin_melanoma_5k":
        candidates = [
            ("Melanoma", "T_Cells"),
            ("Keratinocyte", "T_Cells"),
            ("Endothelial", "T_Cells")
        ]
    else:  # breast_5k
        candidates = [
            ("Macrophages", "T_Cells"),
            ("B_Cells", "T_Cells"),
            ("Fibroblast", "T_Cells"),
            ("Endothelial", "T_Cells")
        ]
    
    gate = Gate(t, ExprSource(t, t.counts, t.genes, {}), cfg)
    pdl1_records = []
    
    for sender, receiver in candidates:
        if sender not in t_idx or receiver not in t_idx:
            continue
        i, j = t_idx[sender], t_idx[receiver]
        contacts = int(C[i, j])
        enrichment = float(E[i, j])
        
        claim = [dict(sender=sender, receiver=receiver, ligand="CD274", receptor="PDCD1")]
        rec = gate.check(claim)[0]
        passed = rec["pass_"]
        status = status_of(rec)
        readout = readout_of(rec)
        
        print(f"  [{sender} -> {receiver} | CD274 -> PDCD1] contacts={contacts}, enrichment={enrichment:.2f}, gate_pass={passed}, status={status}")
        pdl1_records.append({
            "tissue": tissue_name,
            "sender": sender,
            "receiver": receiver,
            "ligand": "CD274",
            "receptor": "PDCD1",
            "contacts_30um": contacts,
            "enrichment_ratio": round(enrichment, 3),
            "gate_pass": passed,
            "gate_status": status,
            "gate_readout": readout,
            "gate_audit": json.dumps(rec)
        })
        
    return pdl1_records


def run_real_controls_5k(tissue_name):
    print(f"\n=== Real-Data Controls on {tissue_name} ===", flush=True)
    t = load_tissue(tissue_name)
    lr = pd.read_csv(f"data/processed/lr_pairs_{tissue_name}.csv")
    lr_simple = lr[~lr.ligand.str.contains("_") & ~lr.receptor.str.contains("_")].reset_index(drop=True)
    pairs = list(zip(lr_simple.ligand, lr_simple.receptor))
    print(f"Total panel L-R pairs: {len(lr)}, Simple pairs evaluated: {len(pairs)}", flush=True)
    
    # Run expression tool
    genes = sorted(set(lr_simple.ligand) | set(lr_simple.receptor))
    gi = {g: i for i, g in enumerate(t.genes)}
    raw = t.counts[:, [gi[g] for g in genes]]
    expr_df = score_all(t, t.lognorm(raw).toarray(), raw.toarray(), genes, pairs, n_perm=200, seed=0)
    
    # Top-10 expr claims
    top10 = top_k(expr_df, "cpdb", 10)
    
    # Weak positive contact controls if present in panel
    controls = [
        ("Macrophages", "T_Cells", "CD80", "CTLA4"),
        ("Macrophages", "T_Cells", "CD86", "CTLA4"),
        ("Endothelial", "Fibroblast", "MMRN2", "CD93"),
    ]
    
    claims = []
    for _, row in top10.iterrows():
        claims.append(dict(sender=row.sender, receiver=row.receiver, ligand=row.ligand, receptor=row.receptor, category="expr_top10"))
    for A, B, L, R in controls:
        if A in t.type_list() and B in t.type_list() and L in t.genes and R in t.genes:
            claims.append(dict(sender=A, receiver=B, ligand=L, receptor=R, category="contact_positive_control"))
            
    gate = Gate(t, ExprSource(t, t.counts, t.genes, {}), cfg)
    recs = gate.check(claims)
    
    out = []
    for c, r in zip(claims, recs):
        out.append({
            "tissue": tissue_name,
            "sender": c["sender"],
            "receiver": c["receiver"],
            "ligand": c["ligand"],
            "receptor": c["receptor"],
            "category": c["category"],
            "gate_pass": r["pass_"],
            "status": status_of(r),
            "readout": readout_of(r)
        })
        print(f"  [{c['category']}] {c['sender']} -> {c['receiver']} ({c['ligand']}-{c['receptor']}): pass={r['pass_']} ({status_of(r)})")
    
    df_out = pd.DataFrame(out)
    df_out.to_csv(f"results/tables/laneS_real_controls_{tissue_name}.csv", index=False)
    return df_out


if __name__ == "__main__":
    summary_s3 = {}
    pdl1_all = []
    
    for t_name in TISSUES_5K:
        bench_df = run_benchmark_5k(t_name)
        pdl1_recs = run_pdl1_analysis(t_name)
        pdl1_all.extend(pdl1_recs)
        ctrl_df = run_real_controls_5k(t_name)
        
        # Summary metrics
        a5_cpdb = bench_df[(bench_df.variant == "A5") & (bench_df.method == "cpdb")]
        summary_s3[t_name] = {
            "bench_A5_fps_mean": round(float(a5_cpdb["false_pass_share"].mean()), 4),
            "bench_A5_ret_mean": round(float(a5_cpdb["retention_true"].mean()), 4),
            "pdl1_records": pdl1_recs,
            "real_controls_total": len(ctrl_df),
            "real_controls_pass_count": int(ctrl_df["gate_pass"].sum())
        }
    
    if pdl1_all:
        pd.DataFrame(pdl1_all).to_csv("results/tables/laneS_s3_pdl1.csv", index=False)
    with open("results/tables/laneS_s3_results.json", "w") as f:
        json.dump(summary_s3, f, indent=2)
    print("\nLane S3 Richer Panels evaluation complete!")
