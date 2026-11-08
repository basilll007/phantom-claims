"""Compare a Lane A rebuild (scripts/laneA_reproduce.sh) with the committed Lane A results.

  python -m scripts.laneA_compare <rebuild_dir> [seed_set]

<rebuild_dir> holds claims/, tables/, figures/ written by the rebuild. Reference = committed results/claims/laneA,
results/tables/laneA_<seed_set>_*, figures/fig_laneA_h4_<seed_set>.png.
Decisions (labels, gate pass/fail, list membership and order, verdict) must match exactly. Continuous columns report the max
abs difference against a tolerance. Exit code 1 if any decision differs or any tolerance is exceeded.
"""
import glob, hashlib, json, os, sys
import numpy as np, pandas as pd

KEY = ["tissue", "seed", "variant", "list", "rank", "sender", "receiver", "ligand", "receptor"]
DECISION = ["label", "neg_type", "category", "mechanism", "passed"]
CONT = ["rho", "p_perm", "q_bh", "n_expr_pairs", "effect_d"]
TOL = 1e-9          # continuous columns: float summation-order noise only


def claims(d, seed_set):
    fs = sorted(glob.glob(f"{d}/{seed_set}_*.csv"))
    x = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in fs], ignore_index=True)
    return x, len(fs)


def txt(s):
    """Text view with missing values made explicit: with pandas' string dtype NaN survives astype(str), and NaN != NaN."""
    return s.astype(object).where(s.notna(), "<NA>").astype(str).values


def num_diff(a, b, cols):
    out = {}
    for c in cols:
        if c in a and c in b and pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
            v = np.abs(a[c].astype(float).values - b[c].astype(float).values)
            out[c] = float(np.nanmax(v)) if len(v) else 0.0
    return out


def main(rb, seed_set="N2"):
    ok = True
    report = {}
    ref_c, n_ref = claims("results/claims/laneA", seed_set)
    new_c, n_new = claims(f"{rb}/claims", seed_set)
    report["claim_files"] = dict(reference=n_ref, rebuild=n_new)
    ok &= n_ref == n_new
    a = ref_c.sort_values(KEY).reset_index(drop=True); b = new_c.sort_values(KEY).reset_index(drop=True)
    same_keys = len(a) == len(b) and all((txt(a[k]) == txt(b[k])).all() for k in KEY)
    report["claim_rows"] = dict(reference=len(a), rebuild=len(b), same_lists_and_order=bool(same_keys))
    ok &= bool(same_keys)
    if same_keys:
        mism = {c: int((txt(a[c]) != txt(b[c])).sum()) for c in DECISION}
        report["decision_mismatches"] = mism
        ok &= sum(mism.values()) == 0
        nd = num_diff(a, b, CONT)
        report["continuous_max_abs_diff"] = nd
        ok &= all(v <= TOL for v in nd.values())

    for name in ["summary", "list_composition", "planted_pass_rates", "true_fail_reasons", "decoy_coloc", "confound",
                 "confound_summary", "lib_density"]:
        fr, fn = f"results/tables/laneA_{seed_set}_{name}.csv", f"{rb}/tables/laneA_{seed_set}_{name}.csv"
        if not (os.path.exists(fr) and os.path.exists(fn)):
            report[name] = "missing: " + ", ".join(p for p in (fr, fn) if not os.path.exists(p)); ok = False; continue
        r, n = pd.read_csv(fr), pd.read_csv(fn)
        if r.shape != n.shape or list(r.columns) != list(n.columns):
            report[name] = f"shape/columns differ {r.shape} vs {n.shape}"; ok = False; continue
        tcols = [c for c in r.columns if not pd.api.types.is_numeric_dtype(r[c])]
        tdiff = int(sum((txt(r[c]) != txt(n[c])).sum() for c in tcols))
        nd = num_diff(r, n, [c for c in r.columns if c not in tcols])
        mx = max(nd.values()) if nd else 0.0
        report[name] = dict(text_cells_differing=tdiff, max_abs_numeric_diff=mx)
        ok &= tdiff == 0 and mx <= TOL

    vr = json.load(open(f"results/tables/laneA_{seed_set}_h4_verdict.json"))
    vn = json.load(open(f"{rb}/tables/laneA_{seed_set}_h4_verdict.json"))
    same_v = vr["verdict"]["verdict"] == vn["verdict"]["verdict"] and \
        {u: x["met"] for u, x in vr["verdict"]["units"].items()} == {u: x["met"] for u, x in vn["verdict"]["units"].items()}
    report["verdict"] = dict(reference=vr["verdict"]["verdict"], rebuild=vn["verdict"]["verdict"], same=bool(same_v))
    ok &= bool(same_v)

    h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()[:16] if os.path.exists(p) else None
    report["figure_sha256"] = dict(reference=h(f"figures/fig_laneA_h4_{seed_set}.png"),
                                   rebuild=h(f"{rb}/figures/fig_laneA_h4_{seed_set}.png"))   # informative only
    report["RESULT"] = "MATCH" if ok else "MISMATCH"
    print(json.dumps(report, indent=1, default=str))
    json.dump(report, open(f"{rb}/compare_{seed_set}.json", "w"), indent=1, default=str)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "N2")
