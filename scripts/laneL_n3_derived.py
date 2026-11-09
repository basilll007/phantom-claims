"""Lane L: two descriptive summaries of lane S's N3 / real-data results, for the paper.

1. Base-rate law, odds form (derived from the logged H9 record; no new gate run).
   With pass rate r_T on true claims and f on a class of negatives, the false share of a filter's passes on a list
   whose true share is pi_T obeys
        FPS / (1 - FPS) = ((1 - pi_T) / pi_T) * (f / r_T),
   so a filter multiplies the odds that a listed claim is true by Lambda = r_T / f. For a target false share phi the
   list must be at least
        pi_min(phi) = f (1 - phi) / (phi r_T + f (1 - phi))
   true. Negative mixes use the compositions of the pre-registered H9 grid: hard = niche : decoy_sender = 0.6 : 0.4,
   easy = phantom : random = 0.10 : 0.05. Class rates are those logged in laneS_N3_H9_base_rate_law_verdict (pooled
   over variants A0-A5m, four tissues, seeds 300-309; niche includes hard niche).

2. Real-data anatomical negative controls (S4): which checks the 138 pre-specified claims failed
   (from results/tables/laneS_anatomical_negcontrols.csv; reasons are non-exclusive).

    PC_RESULTS=results/results.jsonl python scripts/laneL_n3_derived.py --log
Both records are labelled exploratory (post-hoc summaries of logged results).
"""
import json, os, sys
import pandas as pd

RES = os.environ.get("PC_RESULTS", "results/results.jsonl")
R = {}
for line in open(RES, encoding="utf-8"):
    if line.strip():
        r = json.loads(line); R[r["id"]] = r

h9 = R["laneS_N3_H9_base_rate_law_verdict"]["metrics"]
METHODS = [("gate", "class_rates"), ("spatialdm", "spatialdm_class_rates"), ("squidpy_ligrec", "squidpy_class_rates"),
           ("squidpy_ligrec_nhood", "squidpy_spatial_class_rates")]
HARD = {"niche": 0.6, "decoy_sender": 0.4}
EASY = {"phantom": 0.10 / 0.15, "random": 0.05 / 0.15}
PHI = 0.10
rows = []
for name, key in METHODS:
    c = h9[key]
    rT = c["true"]
    f_hard = sum(w * c[k] for k, w in HARD.items())
    f_easy = sum(w * c[k] for k, w in EASY.items())
    pi_min = lambda f: f * (1 - PHI) / (PHI * rT + f * (1 - PHI))
    rows.append(dict(method=name, r_true=rT, f_niche=c["niche"], f_decoy=c["decoy_sender"], f_phantom=c["phantom"],
                     f_random=c["random"], f_hard=f_hard, f_easy=f_easy,
                     lr_hard=rT / f_hard if f_hard > 0 else None, lr_easy=rT / f_easy if f_easy > 0 else None,
                     pi_min_hard_phi10=pi_min(f_hard), pi_min_easy_phi10=pi_min(f_easy)))
T = pd.DataFrame(rows)
os.makedirs("results/tables", exist_ok=True)
T.to_csv("results/tables/laneL_N3_base_rate_odds.csv", index=False)
print(T.round(4).to_string(index=False))

d = pd.read_csv("results/tables/laneS_anatomical_negcontrols.csv")
assert len(d) == 138, len(d)
rs = d.reasons.fillna("")
neg = dict(n_total=int(len(d)), n_passed=int(d.passed.astype(str).str.lower().eq("true").sum()),
           n_fail_contacts=int(rs.str.contains("too_few_expressing_pairs").sum()),
           n_fail_effect=int(rs.str.contains("effect_below_min").sum()),
           n_fail_tile_null=int(rs.str.contains("perm_null_not_rejected").sum()),
           n_fail_halves=int(rs.str.contains("split_inconsistent").sum()),
           n_contacts_ge50=int((d.n_expr_pairs >= 50).sum()),
           n_untestable=int((d.status == "untestable").sum()),
           n_receptor_only=int((d.readout == "receptor_only").sum()),
           n_by_tissue={k: int(v) for k, v in d.tissue.value_counts().items()},
           n_type_pairs=int(d[["tissue", "sender", "receiver"]].drop_duplicates().shape[0]))
print(neg)

# 3. A5 top-10 lists on N3 (the H8 comparison lists): which true claims the frozen gate lost, by planted mechanism
import glob
a5 = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("results/claims/laneS/N3_*_A5.csv"))], ignore_index=True)
assert a5.variant.eq("A5").all() and a5.seed_set.eq("N3").all()
a5 = a5[a5["list"].isin(["cpdb_K10", "prod_K10"]) & (a5.label == 1)].copy()
a5["passed"] = a5.passed.astype(str).str.lower().eq("true")
loss = (a5.groupby("mechanism").passed.agg(n_true="size", n_passed="sum").reset_index())
loss["n_lost"] = loss.n_true - loss.n_passed
print(loss.to_string(index=False))
loss_metrics = dict(n_true=int(loss.n_true.sum()), n_lost=int(loss.n_lost.sum()),
                    by_mechanism={r.mechanism: dict(n_true=int(r.n_true), n_passed=int(r.n_passed),
                                                    n_lost=int(r.n_lost)) for r in loss.itertuples()})

if "--log" in sys.argv:
    assert os.environ.get("PC_RESULTS"), "set PC_RESULTS to the shared results.jsonl"
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import log_result as lr
    lr.RESULTS_FILE = os.environ["PC_RESULTS"]
    recs = [
        (dict(id="laneL_N3_base_rate_odds", phase="P2", hypothesis="descriptive (H9, odds form)",
              source="laneS_N3_H9_base_rate_law_verdict (class rates)", table="results/tables/laneL_N3_base_rate_odds.csv",
              note="Derived from the logged H9 class rates; no gate run. lr_* = r_true / f (odds multiplier); "
                   "pi_min_*_phi10 = minimum true share of a list for a false share of passes <= 0.10. hard = niche "
                   ": decoy_sender 0.6 : 0.4, easy = phantom : random 2 : 1 (the H9 grid compositions).",
              rows=T.to_dict("records")), "N3"),
        (dict(id="laneL_S4_negcontrol_reasons", phase="P2", hypothesis="descriptive (S4 real-data negative controls)",
              source="results/tables/laneS_anatomical_negcontrols.csv", metrics=neg,
              note="Fail reasons are non-exclusive. Real sections, no planting: seed_set is not applicable and 'dev' "
                   "is the logger placeholder."), "dev"),
        (dict(id="laneL_N3_A5_true_losses_by_mechanism", phase="P2", hypothesis="descriptive (H8 lists)",
              source="results/claims/laneS/N3_*_A5.csv (lists cpdb_K10, prod_K10; label == 1)", metrics=loss_metrics,
              note="True claims in the A5 top-10 lists used for H8, by planted mechanism (none = A4-type, no responding "
                   "target; para80 = A3-type; exp15/exp60 = A2; ring = A1); n_lost = true claims the frozen gate "
                   "did not pass."), "N3"),
    ]
    for rec, ss in recs:
        if any(f'"id": "{rec["id"]}"' in l for l in open(lr.RESULTS_FILE, encoding="utf-8")):
            print(f'{rec["id"]} already logged; not logging again')
            continue
        lr.log_result(rec, lane="L", agent_id="claude-opus-5.5", seed_set=ss, labelled="exploratory")
        print("logged", rec["id"])
