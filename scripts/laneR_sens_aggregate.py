"""Lane R, EXPLORATORY: aggregate the sensitivity diagnostics (scripts/laneR_sensitivity.py, dev seeds 0-4, breast Rep1)
and the per-tissue normalisation diagnostic on the N2 A5 hard-niche false passes (N2 is SEEN: descriptive only).

Outputs
  results/tables/laneR_sens_summary.csv        setting x group (A0-A3 pooled | A5): K=10 cpdb+prod retention,
                                               post-gate phantom rate (false share of passes), n_pass, 95% CIs
                                               (hierarchical bootstrap over seeds, then claims; 2000 resamples)
  results/tables/laneR_sens_planted.csv        setting x variant group x planted category: pass rate, n
  results/tables/laneR_edge.json               border-tile dependence of decisions (frozen vs drop_edge)
  results/tables/laneR_N2_A5_niche_norm_by_tissue.csv
  --log  appends exploratory records to results/results.jsonl via scripts/log_result.py (lane R).

  python -m scripts.laneR_sens_aggregate [--log]
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from bench.metrics import hier_bootstrap

TABLES = "results/tables"
ORDER = ["frozen", "tile125", "tile500", "r15", "r60", "nmin25", "nmin100", "split_checker", "null_global",
         "null_density", "drop_edge", "resp_raw"]
GROUPS = {"A0-A3": ["A0", "A1", "A2a", "A2b", "A3"], "A5": ["A5"]}


def load():
    fs = sorted(glob.glob(f"{TABLES}/laneR_sens_claims_breast_rep1_s*.csv"))
    d = pd.concat([pd.read_csv(f) for f in fs])
    d["passed"] = d.passed.astype(str).str.lower().eq("true")
    return d, sorted(d.seed.unique().tolist())


def summary(d):
    L = d[d.list != "planted"]
    out = []
    for s in ORDER:
        for gname, vs in GROUPS.items():
            x = L[(L.setting == s) & L.variant.isin(vs)]
            if not len(x):
                continue
            groups = [(g.label.values.astype(int), g.passed.values) for _, g in x.groupby(["seed", "variant"])]
            point, ci = hier_bootstrap(groups, keys=("phantom_rate_after", "retention"), n_boot=2000)
            out.append(dict(setting=s, group=gname, n_claims=len(x), n_true=int(x.label.sum()), n_pass=int(x.passed.sum()),
                            n_false_pass=int((x.passed & (x.label == 0)).sum()),
                            retention=point["retention"], retention_lo=ci["retention"][0], retention_hi=ci["retention"][1],
                            phantom_after=point["phantom_rate_after"], phantom_after_lo=ci["phantom_rate_after"][0],
                            phantom_after_hi=ci["phantom_rate_after"][1]))
    return pd.DataFrame(out)


def planted(d):
    P = d[d.list == "planted"].copy()
    P["group"] = np.where(P.variant == "A5", "A5", "A0-A3")
    t = P.groupby(["setting", "group", "category"]).passed.agg(["mean", "size"]).reset_index()
    return t.rename(columns={"mean": "pass_rate", "size": "n"})


def edge(d):
    key = ["seed", "variant", "list", "sender", "receiver", "ligand", "receptor"]
    a = d[d.setting == "frozen"].drop_duplicates(key); b = d[d.setting == "drop_edge"].drop_duplicates(key)
    m = a.merge(b, on=key, suffixes=("_f", "_e"))
    flip = m[m.passed_f != m.passed_e]
    return dict(n_claims=int(len(m)), n_flipped=int(len(flip)), n_pass_to_fail=int((flip.passed_f & ~flip.passed_e).sum()),
                n_fail_to_pass=int((~flip.passed_f & flip.passed_e).sum()),
                flipped_by_category=flip.assign(cat=np.where(flip.label_f == 1, "true", flip.category_f))
                .groupby("cat").size().to_dict(),
                median_frac_receivers_in_edge_tiles=float(m.frac_receivers_edge_f.median()),
                note="unique claims per (seed, variant, list); decisions with the frozen gate vs the same gate with "
                     "receivers in border tiles removed (border tile: any of 8 neighbours empty)")


def norm_by_tissue():
    c = pd.read_csv(f"{TABLES}/laneA_N2_confound.csv")
    c["passed"] = c.passed.astype(str).str.lower().eq("true")
    c["raw_supports"] = (c.p_raw < 0.05) & (c.rho_raw >= 0.02)
    rows = []
    for (cat, var), sel in (("niche_hard", "A5"), "hard-niche false passes"), (("true", "A0"), "true passes (reference)"):
        x = c[(c.category == cat) & (c.variant == var) & c.passed]
        for t, g in x.groupby("tissue"):
            rows.append(dict(set=sel, tissue=t, n=len(g), frac_raw_supported=float(g.raw_supports.mean()),
                             median_rho_gate=float(g.rho_gate.median()), median_rho_raw=float(g.rho_raw.median()),
                             median_rho_scale=float(g.rho_scale.median())))
    return pd.DataFrame(rows)


def main(log=False):
    d, seeds = load()
    S, P, E, N = summary(d), planted(d), edge(d), norm_by_tissue()
    S.to_csv(f"{TABLES}/laneR_sens_summary.csv", index=False)
    P.to_csv(f"{TABLES}/laneR_sens_planted.csv", index=False)
    json.dump(E, open(f"{TABLES}/laneR_edge.json", "w"), indent=1)
    N.to_csv(f"{TABLES}/laneR_N2_A5_niche_norm_by_tissue.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(3).to_string(index=False)); print(json.dumps(E, indent=1)); print(N.round(3).to_string(index=False))
    pv = P.pivot_table(index=["setting", "group"], columns="category", values="pass_rate").round(3)
    print(pv.to_string())
    if not log:
        return
    import scripts.log_result as lr
    assert os.environ.get("PC_RESULTS"), "set PC_RESULTS to the shared results.jsonl (main tree)"
    lr.RESULTS_FILE = os.environ["PC_RESULTS"]
    have = {json.loads(l)["id"] for l in open(lr.RESULTS_FILE, encoding="utf-8") if l.strip()}

    def log_result(rec, **kw):                                   # once per id
        if rec["id"] in have:
            print("already logged, skip:", rec["id"]); return
        lr.log_result(rec, **kw)
    base = dict(lane="R", agent_id="claude-opus-5.5-laneR", phase="P2")
    log_result(dict(base, id="laneR_dev_sensitivity", hypothesis="review response (exploratory diagnostics)",
                    tissue="breast_rep1", seeds=seeds, settings=ORDER, groups=GROUPS,
                    note="One setting changed at a time around the frozen gate (scripts/laneR_sensitivity.py); K=10 "
                         "cpdb+prod lists; 95% hierarchical bootstrap. Diagnostics only: nothing feeds back into the "
                         "gate, the variants, the thresholds or any verdict.",
                    rows=S.to_dict("records"), planted=P.to_dict("records"),
                    table="results/tables/laneR_sens_summary.csv"),
               lane="R", seed_set="dev", labelled="exploratory")
    log_result(dict(base, id="laneR_dev_edge_handling", hypothesis="review response (exploratory diagnostics)",
                    tissue="breast_rep1", seeds=seeds, metrics=E, table="results/tables/laneR_edge.json"),
               lane="R", seed_set="dev", labelled="exploratory")
    log_result(dict(base, id="laneR_N2_A5_niche_normalisation_by_tissue", hypothesis="review response (exploratory)",
                    note="N2 is SEEN (descriptive re-analysis of results/tables/laneA_N2_confound.csv; no new gate "
                         "runs). raw_supported = p_raw<0.05 and rho_raw>=0.02 (exposure vs raw target counts).",
                    rows=N.to_dict("records"), table="results/tables/laneR_N2_A5_niche_norm_by_tissue.csv"),
               lane="R", seed_set="N2", labelled="exploratory")
    print("logged 3 records")


if __name__ == "__main__":
    main("--log" in sys.argv)
