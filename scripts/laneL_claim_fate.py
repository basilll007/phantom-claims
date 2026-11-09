"""Where each planted claim stops: a descriptive breakdown of the frozen gate's committed N2 decisions on the Lane A
instances (4 tissues x seeds 200-209 x variants A0-A5, planted-claim list). Nothing is re-run: the per-claim rows in
results/claims/laneA/N2_*.csv are read and their logged reasons are tallied.

    python scripts/laneL_claim_fate.py            # writes results/tables/laneL_N2_claim_fate.csv
    PC_RESULTS=results/results.jsonl python scripts/laneL_claim_fate.py --log   # also logs id laneL_N2_claim_fate

Two views per claim group:
  * any_<check>:   share of claims that fail that check (a claim can fail several; the gate logs all of them);
  * first_<check>: share whose FIRST failed check, in the order evidence -> effect size -> tile null -> split halves,
                   is that check (a presentational order, stated in the figure caption; the gate itself applies all
                   checks and passes a claim only if none fails).
Labelled exploratory: a post-hoc description of confirmatory decisions, not a pre-registered analysis.
"""
import glob, os, sys
import pandas as pd

CHECKS = [("contacts", "too_few_expressing_pairs"), ("effect", "effect_below_min"),
          ("tile_null", "perm_null_not_rejected"), ("halves", "split_inconsistent")]
OTHER = ("gene_not_in_state", "cell_type_absent")
GROUPS = [  # (group id, label used in the figure, variants, planted category)
    ("true_A0", "true, matched kernel (A0)", ["A0"], "true"),
    ("true_A1A3", "true, other kernels (A1–A3)", ["A1", "A2a", "A2b", "A3"], "true"),
    ("true_A4", "true, no target readout (A4)", ["A4"], "true"),
    ("true_A5", "true, in A5 lists", ["A5"], "true"),
    ("phantom", "phantom (cells never meet)", None, "phantom"),
    ("random", "random pair", None, "random"),
    ("niche", "niche (A0–A4)", ["A0", "A1", "A2a", "A2b", "A3", "A4"], "niche"),
    ("niche_hard", "hard niche (A5)", ["A5"], "niche_hard"),
    ("decoy_sender", "decoy sender (A5)", ["A5"], "decoy_sender"),
]


def load(src="results/claims/laneA"):
    fs = sorted(glob.glob(f"{src}/N2_*.csv"))
    assert len(fs) == 280, f"expected 280 N2 files (4 tissues x 10 seeds x 7 variants), found {len(fs)}"
    D = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
    return D[D.list == "planted"].copy()


def tally(P):
    rs = P.reasons.fillna("")
    for k, pref in CHECKS:
        P[f"f_{k}"] = rs.str.contains(pref, regex=False)
    P["f_other"] = rs.apply(lambda s: any(o in s for o in OTHER))
    first = pd.Series("passed", index=P.index)
    for k, _ in reversed(CHECKS):                       # earliest check in the order wins
        first[P[f"f_{k}"]] = k
    first[P.f_other] = "other"
    bad = (first == "passed") & ~P.passed.astype(bool)
    assert not bad.any(), "a failed claim with no recognised reason"
    assert not ((first != "passed") & P.passed.astype(bool)).any(), "a passed claim with a failure reason"
    P["first"] = first
    rows = []
    for gid, lab, vs, cat in GROUPS:
        g = P[(P.category == cat) & (P.variant.isin(vs) if vs else True)]
        r = dict(group=gid, label=lab, category=cat, variants="all" if vs is None else "+".join(vs), n=int(len(g)),
                 pass_rate=float(g.passed.astype(bool).mean()))
        for k, _ in CHECKS + [("other", None)]:
            r[f"first_{k}"] = float((g["first"] == k).mean())
            if k != "other":
                r[f"any_{k}"] = float(g[f"f_{k}"].mean())
                r[f"only_{k}"] = float((g[f"f_{k}"] & (g[[f"f_{c}" for c, _ in CHECKS]].sum(axis=1) == 1)).mean())
        rows.append(r)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    P = load()
    T = tally(P)
    os.makedirs("results/tables", exist_ok=True)
    T.to_csv("results/tables/laneL_N2_claim_fate.csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(T.round(3).to_string(index=False))
    if "--log" in sys.argv:
        assert os.environ.get("PC_RESULTS"), "set PC_RESULTS to the shared results.jsonl"
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import log_result as lr
        lr.RESULTS_FILE = os.environ["PC_RESULTS"]
        rid = "laneL_N2_claim_fate"
        have = [l for l in open(lr.RESULTS_FILE, encoding="utf-8") if f'"id": "{rid}"' in l]
        if have:
            print(f"{rid} already logged; not logging again")
        else:
            lr.log_result(dict(id=rid, phase="P2", hypothesis="descriptive (failure map)",
                               tissues=["breast_rep1", "breast_rep2", "pancreas", "lung"],
                               seeds=list(range(200, 210)), source="results/claims/laneA/N2_*.csv (list == planted)",
                               table="results/tables/laneL_N2_claim_fate.csv",
                               note="Post-hoc description of the frozen gate's committed N2 decisions; no re-run. "
                                    "any_* = share failing that check (non-exclusive); first_* = first failed check "
                                    "in the order contacts -> effect -> tile_null -> halves; only_* = that check is "
                                    "the claim's sole failure.",
                               rows=T.to_dict("records")),
                          lane="L", agent_id="claude-opus-5.5", seed_set="N2", labelled="exploratory")
            print(f"logged {rid}")
