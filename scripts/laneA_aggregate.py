"""Lane A aggregation (H4): pooled K=10 metrics with hierarchical bootstrap CIs, verdict per PREREG Amendment 4 as
operationalised in reports/laneA_prespec.md, list composition, planted pass rates, gate failure reasons.

  python -m scripts.laneA_aggregate <seed_set> [--log]

Reads results/claims/laneA/<seed_set>_*.csv (written by scripts/laneA_run.py). Writes results/tables/laneA_*.
--log appends records to the shared results.jsonl (env PC_RESULTS) via scripts/log_result.py, once per id, and
only when the design is complete (4 tissues x 10 seeds x 7 variants for N2). Confirmatory = N2 primary only.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from bench.metrics import hier_bootstrap, list_metrics
from scripts.laneA_paths import CLAIMS, TABLES

TISSUES = ["breast_rep1", "breast_rep2", "pancreas", "lung"]
VARS = ["A0", "A1", "A2a", "A2b", "A3", "A4", "A5"]
H4_UNITS = {"A1": ["A1"], "A2": ["A2a", "A2b"], "A3": ["A3"], "A4": ["A4"], "A5": ["A5"]}
THR = dict(ret_min=0.60, pr_max=0.10, n_needed=4, fals_ret=0.50, fals_pr=0.15)
KEYS = ("phantom_rate", "phantom_rate_after", "retention", "rel_reduction")
AGENT = "claude-opus-5.5-laneA"
NSEEDS = {"N2": 10, "dev": 5}


def load(seed_set):
    fs = sorted(glob.glob(f"{CLAIMS}/{seed_set}_*.csv"))
    if not fs:
        sys.exit(f"no claim files for {seed_set}")
    d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in fs], ignore_index=True)
    d["passed"] = d.passed.astype(str).str.lower().eq("true")
    d["label"] = d.label.astype(int)
    return d


def groups_for(d, v, K, methods=("cpdb", "prod"), tissues=None):
    x = d[(d.variant == v) & (d.K == K) & (d.method.isin(methods))]
    if tissues:
        x = x[x.tissue.isin(tissues)]
    return [(g.label.values, g.passed.values) for _, g in x.groupby(["tissue", "seed"], sort=True)]


def metric_block(groups, n_boot=2000):
    point, ci = hier_bootstrap(groups, keys=KEYS, n_boot=n_boot)
    return {k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in point.items()}, ci


def rate_ci(groups_bool, n_boot=2000, seed=0):
    """Pass rate over claims, CI by resampling (tissue, seed) groups then claims."""
    rng = np.random.default_rng(seed)
    allv = np.concatenate(groups_bool) if groups_bool else np.array([])
    if len(allv) == 0:
        return None, [None, None], 0
    S = len(groups_bool); bs = []
    for _ in range(n_boot):
        xs = [groups_bool[i][rng.integers(0, len(groups_bool[i]), len(groups_bool[i]))]
              for i in rng.integers(0, S, S) if len(groups_bool[i])]
        bs.append(np.concatenate(xs).mean() if xs else np.nan)
    return float(allv.mean()), [float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))], int(len(allv))


def verdict(prim):
    units = {}
    for u, vs in H4_UNITS.items():
        ret = min(prim[v]["metrics"]["retention"] for v in vs)
        pr = max(prim[v]["metrics"]["phantom_rate_after"] for v in vs)
        units[u] = dict(variants=vs, retention=ret, phantom_rate_after=pr,
                        met=bool(ret >= THR["ret_min"] and pr <= THR["pr_max"]))
    n_met = sum(u["met"] for u in units.values())
    n_low_ret = sum(u["retention"] < THR["fals_ret"] for u in units.values())
    n_high_pr = sum(u["phantom_rate_after"] > THR["fals_pr"] for u in units.values())
    if n_met >= THR["n_needed"]:
        v = "SUPPORTED"
    elif n_low_ret >= 2 or n_high_pr >= 2:
        v = "FALSIFIED"
    else:
        v = "INCONCLUSIVE"
    return dict(verdict=v, n_units_met=n_met, n_units_retention_below_050=n_low_ret,
                n_units_phantom_above_015=n_high_pr, units=units, thresholds=THR)


def main(seed_set, do_log):
    d = load(seed_set)
    os.makedirs(TABLES, exist_ok=True)
    have = d.groupby("variant")[["tissue", "seed"]].apply(lambda x: len(x.drop_duplicates())).to_dict()
    complete = all(have.get(v, 0) == len(TISSUES) * NSEEDS[seed_set] for v in VARS) and \
        set(d.tissue) == set(TISSUES)
    print("design cells per variant:", have, "complete:", complete)

    prim, rows = {}, []
    for v in VARS:
        for K in (10, 20):
            m, ci = metric_block(groups_for(d, v, K))
            rec = dict(variant=v, K=K, scope="pooled", metrics=m, ci=ci)
            if K == 10:
                prim[v] = rec
            rows.append(dict(variant=v, K=K, tissue="ALL", method="cpdb+prod", **{k: m[k] for k in KEYS},
                             **{f"{k}_lo": ci[k][0] for k in KEYS}, **{f"{k}_hi": ci[k][1] for k in KEYS},
                             n_claims=m["n"], n_true=m["n_true"], n_pass=m["n_pass"]))
            for tis in sorted(set(d.tissue)):
                mt, cit = metric_block(groups_for(d, v, K, tissues=[tis]), n_boot=1000)
                rows.append(dict(variant=v, K=K, tissue=tis, method="cpdb+prod", **{k: mt[k] for k in KEYS},
                                 **{f"{k}_lo": cit[k][0] for k in KEYS}, **{f"{k}_hi": cit[k][1] for k in KEYS},
                                 n_claims=mt["n"], n_true=mt["n_true"], n_pass=mt["n_pass"]))
            for meth in ("cpdb", "prod"):
                mm, cim = metric_block(groups_for(d, v, K, methods=(meth,)), n_boot=1000)
                rows.append(dict(variant=v, K=K, tissue="ALL", method=meth, **{k: mm[k] for k in KEYS},
                                 **{f"{k}_lo": cim[k][0] for k in KEYS}, **{f"{k}_hi": cim[k][1] for k in KEYS},
                                 n_claims=mm["n"], n_true=mm["n_true"], n_pass=mm["n_pass"]))
    summ = pd.DataFrame(rows)
    summ.to_csv(f"{TABLES}/laneA_{seed_set}_summary.csv", index=False)
    ver = verdict(prim)
    print(json.dumps({k: v for k, v in ver.items() if k != "units"}, indent=1))
    for u, x in ver["units"].items():
        print(f"  {u}: ret={x['retention']:.3f} PR_after={x['phantom_rate_after']:.3f} met={x['met']}")

    # list composition (all items in tool lists, by neg_type; TRUE counted as 'true')
    comp = []
    lists = d[d.list != "planted"]
    for (v, K), g in lists.groupby(["variant", "K"]):
        vc = g.apply(lambda r: "true" if r.label == 1 else r.neg_type, axis=1).value_counts()
        comp.append(dict(variant=v, K=K, n_lists=g.groupby(["tissue", "seed", "method"]).ngroups, **vc.to_dict()))
    comp = pd.DataFrame(comp).fillna(0)
    comp.to_csv(f"{TABLES}/laneA_{seed_set}_list_composition.csv", index=False)
    print(comp.to_string(index=False))

    # planted-claim pass rates by category, and by mechanism for TRUE pairs
    pl = d[d.list == "planted"]
    prow = []
    for (v, cat), g in pl.groupby(["variant", "category"]):
        gb = [x.passed.values for _, x in g.groupby(["tissue", "seed"])]
        r, ci, n = rate_ci(gb)
        prow.append(dict(variant=v, category=cat, mechanism="", pass_rate=r, lo=ci[0], hi=ci[1], n=n))
    for (v, mech), g in pl[(pl.category.isin(["true", "decoy_sender"])) & (pl.variant == "A5")].groupby(
            ["variant", "mechanism"]):
        for cat, gg in g.groupby("category"):
            gb = [x.passed.values for _, x in gg.groupby(["tissue", "seed"])]
            r, ci, n = rate_ci(gb)
            prow.append(dict(variant=v, category=cat, mechanism=mech, pass_rate=r, lo=ci[0], hi=ci[1], n=n))
    prow = pd.DataFrame(prow)
    prow.to_csv(f"{TABLES}/laneA_{seed_set}_planted_pass_rates.csv", index=False)
    print(prow.round(3).to_string(index=False))

    # why TRUE claims fail (K=10 lists, pooled tools) — first failing reason class
    rs = []
    tr = d[(d.list != "planted") & (d.K == 10) & (d.label == 1) & (~d.passed)]
    for v, g in tr.groupby("variant"):
        cnt = {}
        for r in g.reasons:
            for c in str(r).split(";"):
                c = c.split("(")[0]
                if c:
                    cnt[c] = cnt.get(c, 0) + 1
        rs.append(dict(variant=v, n_true_failed=len(g), **cnt))
    rs = pd.DataFrame(rs).fillna(0)
    rs.to_csv(f"{TABLES}/laneA_{seed_set}_true_fail_reasons.csv", index=False)
    print(rs.to_string(index=False))
    json.dump(dict(verdict=ver, complete=complete, design=have), open(f"{TABLES}/laneA_{seed_set}_h4_verdict.json", "w"),
              indent=1, default=float)

    if do_log:
        if not complete:
            print("design incomplete: nothing logged"); return
        log(seed_set, d, prim, ver, summ, comp, prow, rs)


def log(seed_set, d, prim, ver, summ, comp, prow, rs):
    import scripts.log_result as lr
    lr.RESULTS_FILE = os.environ.get("PC_RESULTS", lr.RESULTS_FILE)
    have = set()
    if os.path.exists(lr.RESULTS_FILE):
        for line in open(lr.RESULTS_FILE, encoding="utf-8"):
            try:
                have.add(json.loads(line).get("id"))
            except Exception:
                pass
    conf = "confirmatory" if seed_set == "N2" else "exploratory"
    seeds = sorted(set(int(s) for s in d.seed))
    common = dict(lane="A", agent_id=AGENT, seed_set=seed_set, phase="P2", hypothesis="H4", tissues=TISSUES,
                  seeds=seeds, gate="frozen (Amendment 2 defaults)", design="reports/laneA_prespec.md",
                  claims_dir="results/claims/laneA")

    def put(rec, labelled):
        if rec["id"] in have:
            print("exists, skip:", rec["id"]); return
        lr.log_result({**common, **rec}, lane="A", agent_id=AGENT, seed_set=seed_set, labelled=labelled)
        print("logged", rec["id"], labelled)

    for v in VARS:
        lab = conf if v != "A0" else "exploratory"          # A0 = gate-matched reference, not part of H4
        put(dict(id=f"laneA_{seed_set}_H4_{v}_K10_pooled", variant=v, K=10, methods=["cpdb", "prod"],
                 metrics=prim[v]["metrics"], ci=prim[v]["ci"],
                 note="phantom_rate = expression-only list before gate; phantom_rate_after/retention = frozen gate"),
            lab)
    put(dict(id=f"laneA_{seed_set}_H4_verdict", metrics=ver, verdict=ver["verdict"]), conf)
    put(dict(id=f"laneA_{seed_set}_H4_secondary_summary", table=f"results/tables/laneA_{seed_set}_summary.csv",
             rows=summ.round(4).to_dict("records")), "exploratory")
    put(dict(id=f"laneA_{seed_set}_list_composition", table=f"results/tables/laneA_{seed_set}_list_composition.csv",
             rows=comp.to_dict("records")), "exploratory")
    put(dict(id=f"laneA_{seed_set}_planted_pass_rates", table=f"results/tables/laneA_{seed_set}_planted_pass_rates.csv",
             rows=prow.round(4).to_dict("records")), "exploratory")
    put(dict(id=f"laneA_{seed_set}_true_fail_reasons", table=f"results/tables/laneA_{seed_set}_true_fail_reasons.csv",
             rows=rs.to_dict("records")), "exploratory")


if __name__ == "__main__":
    main(sys.argv[1], "--log" in sys.argv)
