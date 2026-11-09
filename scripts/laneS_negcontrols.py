"""Anatomical Negative Controls on Real Data (Lane S, S4).
Implements the pre-specified protocol from reports/negcontrols_prespec.md.

Selects L-R pairs between cell types occupying mutually exclusive anatomical compartments
where both L and R are genuinely expressed (fraction >= 0.05 in respective types).
Runs the frozen gate and audits any false passes against library size and local density confounds.

Outputs:
    results/tables/laneS_anatomical_negcontrols.csv
    results/audit/laneS_anatomical_negcontrols.jsonl
"""
import json, os, sys, time
import numpy as np, pandas as pd
from scipy.stats import beta, spearmanr
from scipy.spatial import cKDTree

from bench.tissue import load_tissue
from gate.gate import Gate, GateConfig, ExprSource, write_audit
from gate.status import status_of, readout_of

FROZEN = dict(r=30.0, tile=250.0, n_perm=1000, alpha_half=0.10, n_min=50, rho_min=0.02, q=0.05)
AGENT = "gemini-3.8-flash"

SEGREGATED_PAIRS = {
    "breast_rep1": [
        ("DCIS_1", "B_Cells"), ("DCIS_2", "B_Cells"),
        ("DCIS_1", "Mast_Cells"), ("DCIS_2", "Mast_Cells"),
        ("Myoepi_ACTA2+", "B_Cells"), ("Myoepi_KRT15+", "B_Cells"),
    ],
    "breast_rep2": [
        ("DCIS_1", "B_Cells"), ("DCIS_2", "B_Cells"),
        ("DCIS_1", "Mast_Cells"), ("DCIS_2", "Mast_Cells"),
        ("Myoepi_ACTA2+", "B_Cells"), ("Myoepi_KRT15+", "B_Cells"),
    ],
    "pancreas": [
        ("Acinar", "Epithelial_Tumor"), ("Epithelial_Tumor", "Acinar"),
        ("Endocrine", "Epithelial_Tumor"), ("Epithelial_Tumor", "Endocrine"),
        ("Acinar", "Plasma_Cells"), ("Plasma_Cells", "Acinar"),
    ],
    "lung": [
        ("Ciliated", "Epithelial_Tumor"), ("Epithelial_Tumor", "Ciliated"),
        ("Ciliated", "Plasma_Cells"), ("Plasma_Cells", "Ciliated"),
        ("Ciliated", "B_Cells"), ("B_Cells", "Ciliated"),
    ],
}


def clopper_pearson_ci(k, n, alpha=0.05):
    """Exact Clopper-Pearson binomial confidence interval."""
    if n == 0:
        return [0.0, 0.0]
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return [lo, hi]


def select_claims(t, lr_df, min_frac=0.05):
    tname = t.name
    type_pairs = SEGREGATED_PAIRS.get(tname, [])
    tix = {g: i for i, g in enumerate(t.genes)}
    raw = t.counts
    
    valid_claims = []
    for sender, receiver in type_pairs:
        s_mask = (t.types == sender)
        r_mask = (t.types == receiver)
        if s_mask.sum() == 0 or r_mask.sum() == 0:
            continue
        for _, row in lr_df.iterrows():
            L, R = row.ligand, row.receptor
            if L not in tix or R not in tix:
                continue
            frac_L = float((raw[s_mask, tix[L]] > 0).mean())
            frac_R = float((raw[r_mask, tix[R]] > 0).mean())
            mean_L = float(raw[s_mask, tix[L]].mean())
            mean_R = float(raw[r_mask, tix[R]].mean())
            if frac_L >= min_frac and frac_R >= min_frac:
                valid_claims.append({
                    "tissue": tname,
                    "sender": sender,
                    "receiver": receiver,
                    "ligand": L,
                    "receptor": R,
                    "fracL": frac_L,
                    "fracR": frac_R,
                    "meanL": mean_L,
                    "meanR": mean_R,
                })
    return valid_claims


def audit_confounds(t, passed_claims):
    """Inspect passed claims for library size and local density confounds."""
    tree = t.tree
    # Local density = number of cells within 50 um
    print("  Computing local cell densities for confound audit...")
    densities = np.array([len(tree.query_ball_point(pt, 50.0)) for pt in t.xy])
    lib_sizes = np.array(t.counts.sum(axis=1)).flatten()
    
    confound_records = []
    tix = {g: i for i, g in enumerate(t.genes)}
    for c in passed_claims:
        sender, receiver = c["sender"], c["receiver"]
        L, R = c["ligand"], c["receptor"]
        r_idx = np.where(t.types == receiver)[0]
        if len(r_idx) == 0:
            continue
        # Receptor counts in receivers
        r_counts = np.array(t.counts[r_idx, tix[R]].toarray()).flatten()
        r_libs = lib_sizes[r_idx]
        r_dens = densities[r_idx]
        
        rho_lib, _ = spearmanr(r_counts, r_libs)
        rho_den, _ = spearmanr(r_counts, r_dens)
        
        confound_records.append({
            "claim": f"{sender}->{receiver}:{L}-{R}",
            "tissue": c["tissue"],
            "spearman_receptor_vs_libsize": float(rho_lib) if not np.isnan(rho_lib) else 0.0,
            "spearman_receptor_vs_density": float(rho_den) if not np.isnan(rho_den) else 0.0,
        })
    return confound_records


def run_negcontrols(do_log=False):
    cfg = GateConfig(**FROZEN)
    os.makedirs("results/tables", exist_ok=True)
    os.makedirs("results/audit", exist_ok=True)
    
    all_claims = []
    tissues = ["breast_rep1", "breast_rep2", "pancreas", "lung"]
    
    print("=== Selecting Anatomical Negative Controls across 4 Tissues ===")
    for tis in tissues:
        t = load_tissue(tis, f"data/processed/{tis}.h5ad")
        lr_file = f"data/processed/lr_pairs_{tis}.csv"
        lr_df = pd.read_csv(lr_file)
        claims = select_claims(t, lr_df, min_frac=0.05)
        print(f"[{tis}] {len(claims)} candidate pairs met expression threshold (frac >= 0.05)")
        all_claims.extend(claims)
        
    print(f"\nTotal pre-specified anatomical negative claims: {len(all_claims)}")
    
    # Run the frozen gate per tissue
    all_results = []
    all_audits = []
    
    for tis in tissues:
        tis_claims = [c for c in all_claims if c["tissue"] == tis]
        if not tis_claims:
            continue
        t = load_tissue(tis, f"data/processed/{tis}.h5ad")
        gate = Gate(t, ExprSource(t, t.counts, t.genes, {}), cfg)
        
        gate_cl = [{"sender": c["sender"], "receiver": c["receiver"],
                    "ligand": c["ligand"], "receptor": c["receptor"]} for c in tis_claims]
        
        recs = gate.check(gate_cl)
        for c, rec in zip(tis_claims, recs):
            rec["tissue"] = tis
            all_audits.append(rec)
            all_results.append({
                **c,
                "passed": bool(rec["pass_"]),
                "q_bh": rec.get("q_bh"),
                "rho": rec.get("rho_within_tile"),
                "p_perm": rec.get("p_perm"),
                "n_expr_pairs": rec.get("n_expressing_pairs_within_r"),
                "reasons": ";".join(rec.get("reasons", [])),
                "status": status_of(rec),
                "readout": readout_of(rec),
            })
            
    df_res = pd.DataFrame(all_results)
    df_res.to_csv("results/tables/laneS_anatomical_negcontrols.csv", index=False)
    write_audit(all_audits, "results/audit/laneS_anatomical_negcontrols.jsonl")
    
    # Calculate pass rate and CI
    n_total = len(df_res)
    n_pass = int(df_res.passed.sum())
    pass_rate = float(n_pass / n_total) if n_total > 0 else 0.0
    ci = clopper_pearson_ci(n_pass, n_total)
    
    print("\n=== Anatomical Negative Controls Results ===")
    print(f"Total Claims Evaluated: {n_total}")
    print(f"Claims Passed:           {n_pass}")
    print(f"Pass Rate:               {pass_rate:.4f} ({pass_rate*100:.2f}%)")
    print(f"95% CI (Clopper-Pearson): [{ci[0]:.4f}, {ci[1]:.4f}]")
    
    passed_df = df_res[df_res.passed]
    confound_analysis = []
    if len(passed_df) > 0:
        print(f"\n[WARNING] {len(passed_df)} claims passed. Inspecting confounds...")
        for tis in tissues:
            t_pass = passed_df[passed_df.tissue == tis]
            if len(t_pass) > 0:
                t = load_tissue(tis, f"data/processed/{tis}.h5ad")
                confounds = audit_confounds(t, t_pass.to_dict("records"))
                confound_analysis.extend(confounds)
        print(pd.DataFrame(confound_analysis).to_string(index=False))
    else:
        print("\n[PERFECT PASS] Exactly 0 claims passed the frozen gate! Anatomical segregation holds with 0 false passes.")
        
    summary = {
        "n_total": n_total,
        "n_passed": n_pass,
        "pass_rate": pass_rate,
        "ci_95": ci,
        "passed_claims": passed_df[["tissue", "sender", "receiver", "ligand", "receptor", "n_expr_pairs", "rho"]].to_dict("records"),
        "confound_analysis": confound_analysis,
    }
    json.dump(summary, open("results/tables/laneS_anatomical_negcontrols_summary.json", "w"), indent=2, default=float)
    
    if do_log:
        import scripts.log_result as lr
        lr.log_result({
            "id": "laneS_anatomical_negcontrols_real_data",
            "lane": "S",
            "agent_id": AGENT,
            "phase": "P2",
            "hypothesis": "S4_anatomical_negative_controls",
            "verdict": "SUPPORTED" if pass_rate <= 0.05 else "FALSIFIED",
            "metrics": summary,
            "gate": "frozen (Amendment 2 defaults)",
            "design": "reports/negcontrols_prespec.md",
            "table": "results/tables/laneS_anatomical_negcontrols.csv",
            "note": "Real sections, 4 cancer tissues, mutually exclusive anatomical compartments; seed_set 'dev' is logger placeholder",
        }, lane="S", agent_id=AGENT, seed_set="dev", labelled="confirmatory")
        print("Logged anatomical negative controls to results.jsonl successfully.")


if __name__ == "__main__":
    do_log = "--log" in sys.argv
    run_negcontrols(do_log)
