"""Lane B aggregation (H5, PREREG Amendment 4; operationalised in reports/laneB_prespec.md).

  python -m scripts.laneB_aggregate <seed_set> [--log]

H5a (expression-only baseline): for each (tissue, seed) unit, Spearman rho between effect level and the phantom rate of
the unit's pooled cpdb+prod K=10 default lists across the 4 levels. A unit whose phantom rate is constant across levels
gets rho = 0, i.e. no monotone increase (conservative). Criterion: median unit rho <= -0.80.
H5b (frozen gate): pooled post-gate phantom rate (K=10, cpdb+prod, all units) <= 0.15 at every level.
Verdict: SUPPORTED if H5a and H5b; FALSIFIED if median rho > -0.50 or post-gate phantom rate > 0.20 at any level;
otherwise INCONCLUSIVE.
A5m / A5 list types are exploratory and reported separately (addendum).
--log appends to $PC_RESULTS once per id, only when the design is complete.
"""
import glob, json, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from bench.metrics import hier_bootstrap

TISSUES = ["breast_rep1", "breast_rep2", "pancreas", "lung"]
LEVELS = [1.0, 0.5, 0.25, 0.125]
LTYPES = ["default", "A5m", "A5"]
KEYS = ("phantom_rate", "precision", "phantom_rate_after", "retention", "rel_reduction")
THR = dict(rho_max=-0.80, pr_max=0.15, fals_rho=-0.50, fals_pr=0.20)
AGENT = "claude-opus-5.5-laneB"
NSEEDS = {"N2": 10, "dev": 5}
CLAIMS = os.environ.get("PC_LANEB_CLAIMS", "results/claims/laneB")
TABLES = os.environ.get("PC_TABLES", "results/tables")


def load(seed_set):
    fs = sorted(glob.glob(f"{CLAIMS}/{seed_set}_*.csv"))
    if not fs:
        sys.exit("no claim files")
    d = pd.concat([pd.read_csv(f, keep_default_na=False, na_values=[""]) for f in fs], ignore_index=True)
    d["passed"] = d.passed.astype(str).str.lower().eq("true")
    d["label"] = d.label.astype(int)
    return d


def groups(d, lt, e, K, tissues=None, methods=("cpdb", "prod")):
    x = d[(d.list_type == lt) & (d.effect == e) & (d.K == K) & d.method.isin(methods)]
    if tissues:
        x = x[x.tissue.isin(tissues)]
    return [(g.label.values, g.passed.values) for _, g in x.groupby(["tissue", "seed"], sort=True)]


def block(gs, n_boot=2000):
    point, ci = hier_bootstrap(gs, keys=KEYS, n_boot=n_boot)
    return {k: (None if isinstance(v, float) and np.isnan(v) else v) for k, v in point.items()}, ci


def unit_rhos(d, lt):
    x = d[(d.list_type == lt) & (d.K == 10) & d.method.isin(["cpdb", "prod"])]
    pr = x.groupby(["tissue", "seed", "effect"]).label.apply(lambda y: 1 - y.mean()).reset_index(name="pr")
    rows = []
    for (tis, s), g in pr.groupby(["tissue", "seed"]):
        g = g.sort_values("effect")
        if len(g) < len(LEVELS):
            continue
        if g.pr.nunique() == 1:
            rho = 0.0                                           # constant -> no monotone increase (conservative)
        else:
            # 12-decimal rounding: Spearman over 4 levels is discrete (-1, -0.8, ...) and -0.8 must not become
            # -0.7999999999999999 and fail "<= -0.80" by float error
            rho = round(float(spearmanr(g.effect, g.pr).statistic), 12)
        rows.append(dict(tissue=tis, seed=int(s), rho=rho, **{f"pr_e{e:g}": float(v) for e, v in zip(g.effect, g.pr)}))
    return pd.DataFrame(rows)


def main(seed_set, do_log):
    d = load(seed_set)
    cells = d.groupby(["list_type", "effect"])[["tissue", "seed"]].apply(lambda x: len(x.drop_duplicates()))
    complete = all(cells.get((lt, e), 0) == len(TISSUES) * NSEEDS[seed_set] for lt in LTYPES for e in LEVELS) and \
        set(d.tissue) == set(TISSUES)
    print("design complete:", complete)
    os.makedirs(TABLES, exist_ok=True)
    rows, prim = [], {}
    for lt in LTYPES:
        for e in LEVELS:
            for K in (5, 10, 20):
                m, ci = block(groups(d, lt, e, K))
                if K == 10:
                    prim[(lt, e)] = (m, ci)
                rows.append(dict(list_type=lt, effect=e, K=K, tissue="ALL", **{k: m.get(k) for k in KEYS},
                                 **{f"{k}_lo": ci[k][0] for k in KEYS}, **{f"{k}_hi": ci[k][1] for k in KEYS},
                                 n_claims=m["n"], n_true=m["n_true"], n_pass=m["n_pass"]))
            for tis in sorted(set(d.tissue)):
                m, ci = block(groups(d, lt, e, 10, tissues=[tis]), n_boot=1000)
                rows.append(dict(list_type=lt, effect=e, K=10, tissue=tis, **{k: m.get(k) for k in KEYS},
                                 **{f"{k}_lo": ci[k][0] for k in KEYS}, **{f"{k}_hi": ci[k][1] for k in KEYS},
                                 n_claims=m["n"], n_true=m["n_true"], n_pass=m["n_pass"]))
    summ = pd.DataFrame(rows)
    summ.to_csv(f"{TABLES}/laneB_{seed_set}_summary.csv", index=False)

    rho = {lt: unit_rhos(d, lt) for lt in LTYPES}
    pd.concat([r.assign(list_type=lt) for lt, r in rho.items()]).to_csv(f"{TABLES}/laneB_{seed_set}_unit_spearman.csv",
                                                                       index=False)
    r0 = rho["default"]
    med = float(r0.rho.median())
    pr_after = {e: round(prim[("default", e)][0]["phantom_rate_after"], 12) for e in LEVELS}
    h5a = med <= THR["rho_max"]
    h5b = all(v <= THR["pr_max"] for v in pr_after.values())
    if h5a and h5b:
        v = "SUPPORTED"
    elif med > THR["fals_rho"] or any(x > THR["fals_pr"] for x in pr_after.values()):
        v = "FALSIFIED"
    else:
        v = "INCONCLUSIVE"
    ver = dict(verdict=v, H5a_met=bool(h5a), H5b_met=bool(h5b), median_unit_rho=med, n_units=int(len(r0)),
               frac_units_rho_le_m080=float((r0.rho <= THR["rho_max"]).mean()), min_unit_rho=float(r0.rho.min()),
               max_unit_rho=float(r0.rho.max()), post_gate_phantom_rate_by_level={f"{e:g}": pr_after[e] for e in LEVELS},
               n_pass_by_level={f"{e:g}": prim[("default", e)][0]["n_pass"] for e in LEVELS},
               retention_by_level={f"{e:g}": prim[("default", e)][0]["retention"] for e in LEVELS}, thresholds=THR)
    print(json.dumps(ver, indent=1, default=float))

    comp = []
    for (lt, e, K), g in d[d.list != "planted"].groupby(["list_type", "effect", "K"]):
        vc = g.apply(lambda r: "true" if r.label == 1 else r.neg_type, axis=1).value_counts()
        comp.append(dict(list_type=lt, effect=e, K=K, n_lists=g.groupby(["tissue", "seed", "method"]).ngroups, **vc.to_dict()))
    comp = pd.DataFrame(comp).fillna(0)
    comp.to_csv(f"{TABLES}/laneB_{seed_set}_list_composition.csv", index=False)
    pl = d[d.list == "planted"].groupby(["list_type", "effect", "category"]).passed.agg(["mean", "size"]).reset_index()
    pl.columns = ["list_type", "effect", "category", "pass_rate", "n"]
    pl.to_csv(f"{TABLES}/laneB_{seed_set}_planted_pass_rates.csv", index=False)
    fp = []
    for (lt, e), g in d[(d.list != "planted") & (d.K == 10) & d.passed & (d.label == 0)].groupby(["list_type", "effect"]):
        fp.append(dict(list_type=lt, effect=e, **g.neg_type.value_counts().to_dict()))
    fp = pd.DataFrame(fp).fillna(0)
    fp.to_csv(f"{TABLES}/laneB_{seed_set}_false_pass_composition.csv", index=False)
    json.dump(dict(verdict=ver, complete=complete), open(f"{TABLES}/laneB_{seed_set}_h5_verdict.json", "w"), indent=1,
              default=float)
    pd.set_option("display.width", 250)
    print(summ[(summ.K == 10) & (summ.tissue == "ALL")][["list_type", "effect", "phantom_rate", "phantom_rate_after",
                                                         "retention", "n_pass"]].round(3).to_string(index=False))
    print(comp[comp.K == 10].to_string(index=False)); print(fp.to_string(index=False))
    if do_log:
        if not complete:
            print("design incomplete: nothing logged"); return
        log(seed_set, d, prim, ver, summ, rho, comp, pl, fp)


def log(seed_set, d, prim, ver, summ, rho, comp, pl, fp):
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
    common = dict(lane="B", agent_id=AGENT, seed_set=seed_set, phase="P2", hypothesis="H5", tissues=TISSUES,
                  seeds=sorted(int(s) for s in set(d.seed)), gate="frozen (Amendment 2 defaults)",
                  design="reports/laneB_prespec.md", claims_dir=CLAIMS)

    def put(rec, labelled):
        if rec["id"] in have:
            print("exists, skip:", rec["id"]); return
        lr.log_result({**common, **rec}, lane="B", agent_id=AGENT, seed_set=seed_set, labelled=labelled)
        print("logged", rec["id"], labelled)

    for lt in LTYPES:
        for e in LEVELS:
            m, ci = prim[(lt, e)]
            put(dict(id=f"laneB_{seed_set}_H5_{lt}_e{e:g}_K10_pooled", list_type=lt, effect=e, K=10, methods=["cpdb", "prod"],
                     metrics=m, ci=ci, note="phantom_rate/precision = expression-only list before gate; "
                                            "phantom_rate_after/retention = frozen gate"),
                conf if lt == "default" else "exploratory")
    put(dict(id=f"laneB_{seed_set}_H5_unit_spearman", list_type="default",
             metrics=dict(median_unit_rho=ver["median_unit_rho"], frac_units_rho_le_m080=ver["frac_units_rho_le_m080"],
                          min_unit_rho=ver["min_unit_rho"], max_unit_rho=ver["max_unit_rho"], n_units=ver["n_units"]),
             units=rho["default"].round(4).to_dict("records"), table=f"{TABLES}/laneB_{seed_set}_unit_spearman.csv"), conf)
    put(dict(id=f"laneB_{seed_set}_H5_verdict", metrics=ver, verdict=ver["verdict"]), conf)
    for lt in ("A5m", "A5"):
        r = rho[lt]
        put(dict(id=f"laneB_{seed_set}_{lt}_unit_spearman", list_type=lt,
                 metrics=dict(median_unit_rho=float(r.rho.median()), frac_units_rho_le_m080=float((r.rho <= -0.8).mean()))),
            "exploratory")
    put(dict(id=f"laneB_{seed_set}_secondary_summary", table=f"{TABLES}/laneB_{seed_set}_summary.csv",
             rows=summ.round(4).to_dict("records")), "exploratory")
    put(dict(id=f"laneB_{seed_set}_list_composition", rows=comp.to_dict("records")), "exploratory")
    put(dict(id=f"laneB_{seed_set}_planted_pass_rates", rows=pl.round(4).to_dict("records")), "exploratory")
    put(dict(id=f"laneB_{seed_set}_false_pass_composition", rows=fp.to_dict("records")), "exploratory")


if __name__ == "__main__":
    main(sys.argv[1], "--log" in sys.argv)
