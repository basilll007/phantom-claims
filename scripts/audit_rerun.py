"""Lane L audit: re-run a random 10% of a lane's logged N2 design cells on a DIFFERENT seed, to catch silent bugs
(results that hinge on one seed, cached state, or a mis-wired variant).

  python -m scripts.audit_rerun laneA [frac=0.10] [rng_seed=20261108]

laneA: cells = (tissue, N2 seed, variant). Each sampled cell is re-run with the same tissue and variant on one dev seed
(0-4, drawn at random; dev is the only non-N2, non-test1 seed set; N3 stays untouched for a future amendment). Output goes
to a separate claims directory (results/claims/audit_laneA). Check: the re-run's K=10 retention and post-gate phantom rate
(cpdb+prod pooled) must fall inside the range spanned by the 10 N2 seeds of the same tissue x variant, widened by 0.10.
A miss is a FLAG to investigate, not proof of a bug: dev seeds are a different sample.
Writes results/tables/audit_laneA_rerun.csv. Logging is done separately by lane L.
"""
import glob, os, sys
import numpy as np, pandas as pd

MARGIN = 0.10


def k10_metrics(d):
    x = d[(d.K == 10) & d.method.isin(["cpdb", "prod"])]
    y = x.label.astype(int).values; p = x.passed.astype(str).str.lower().eq("true").values
    ret = float(y[p].sum() / y.sum()) if y.sum() else np.nan
    pr = float(1 - y[p].mean()) if p.sum() else 0.0
    return ret, pr


def main(lane="laneA", frac=0.10, rng_seed=20261108):
    assert lane == "laneA", "only lane A has logged N2 results so far"
    out_claims = "results/claims/audit_laneA"
    os.environ["PC_LANEA_CLAIMS"] = out_claims          # read by scripts.laneA_paths at import
    from scripts import laneA_run
    ref = {}
    for f in glob.glob("results/claims/laneA/N2_*.csv"):
        d = pd.read_csv(f)
        ref[(d.tissue.iloc[0], int(d.seed.iloc[0]), d.variant.iloc[0])] = k10_metrics(d)
    cells = sorted(ref)
    rng = np.random.default_rng(rng_seed)
    pick = [cells[i] for i in sorted(rng.choice(len(cells), size=int(round(frac * len(cells))), replace=False))]
    plan = [(t, s, v, int(rng.integers(0, 5))) for t, s, v in pick]
    for t, s, v, ds in plan:
        laneA_run.run(t, "dev", [ds], [v])
    rows = []
    for t, s, v, ds in plan:
        d = pd.read_csv(f"{out_claims}/dev_{t}_s{ds}_{v}.csv")
        ret, pr = k10_metrics(d)
        n2 = [ref[(t, ss, v)] for ss in range(200, 210)]
        rets = [a for a, _ in n2 if not np.isnan(a)]; prs = [b for _, b in n2]
        ok_r = np.isnan(ret) or (min(rets) - MARGIN <= ret <= max(rets) + MARGIN)
        ok_p = min(prs) - MARGIN <= pr <= max(prs) + MARGIN
        rows.append(dict(tissue=t, n2_seed=s, variant=v, dev_seed=ds, retention=ret, n2_ret_min=min(rets), n2_ret_max=max(rets),
                         n2_ret_this_seed=ref[(t, s, v)][0], phantom_after=pr, n2_pr_min=min(prs), n2_pr_max=max(prs),
                         n2_pr_this_seed=ref[(t, s, v)][1], flag=not (ok_r and ok_p)))
    res = pd.DataFrame(rows)
    res.to_csv("results/tables/audit_laneA_rerun.csv", index=False)
    pd.set_option("display.width", 250)
    print(res.round(3).to_string(index=False))
    print(f"cells re-run: {len(res)} of {len(cells)}; flagged: {int(res.flag.sum())}")
    return res


if __name__ == "__main__":
    a = sys.argv
    main(a[1] if len(a) > 1 else "laneA", float(a[2]) if len(a) > 2 else 0.10, int(a[3]) if len(a) > 3 else 20261108)
