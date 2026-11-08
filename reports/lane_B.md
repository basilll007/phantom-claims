# Lane B report: H5 dose-response
Agent `claude-opus-5.5-laneB`. This is the same agent as lane A and lane L, taking lane B as a second task at the human's
request. Branch `lane-B`. Every number below is from `results/results.jsonl` (ids in brackets). Status: **N2 run complete
and logged.**

## 1. Verdict
**H5: FALSIFIED (confirmatory, N2), and reported as a miss.** [`laneB_N2_H5_verdict`]
- **H5a, expression-only tools: met, strongly.** Their phantom rate climbs as the effect falls. The median per-unit
  Spearman rho is -1.0, and all 40 (tissue, seed) units are <= -0.80 (31 units at -1.0, the rest -0.95 or -0.89)
  [`laneB_N2_H5_unit_spearman`].
- **H5b, frozen gate <= 0.15 at every level: not met.** At effect 1/8 the post-gate phantom rate is 0.447, above both the
  0.15 target and the 0.20 falsification line.

Design (pre-spec `reports/laneB_prespec.md`, commit `e330eb7`, 12:13:02, before any effect < 1 N2 instance existed):
4 levels x 40 units (4 tissues x N2 seeds 200-209); frozen gate; K=10 lists of the two expression-only tools pooled;
95% hierarchical bootstrap CI. [`laneB_N2_H5_default_e<level>_K10_pooled`]

| effect | expr-only phantom rate [CI] | post-gate phantom rate [CI] | retention [CI] | TRUE in lists / passed |
|---|---|---|---|---|
| 1 | 0.190 [0.152, 0.230] | 0.000 [0, 0] | 0.994 [0.982, 1.000] | 648 / 644 |
| 1/2 | 0.304 [0.257, 0.349] | 0.016 [0.002, 0.034] | 0.982 [0.957, 0.998] | 557 / 556 |
| 1/4 | 0.758 [0.711, 0.801] | 0.123 [0.049, 0.211] | 0.809 [0.714, 0.897] | 194 / 179 |
| 1/8 | 0.965 [0.944, 0.982] | **0.447 [0.154, 0.706]** | 0.929 [0.706, 1.000] | 28 / 47 |

Figure: `figures/fig_dose_response.png`.

## 2. Findings, negatives first
1. **At low effect, the gate's output is mostly false. The cause is base rate, not a per-claim weakness.**
   [`laneB_N2_negative_pass_rate_by_level`]
   - The gate's pass rate per negative claim in the default lists is about flat: 0.037 at 1/2, 0.036 at 1/4, 0.027 at
     1/8 (0 at effect 1, where the negatives are phantoms).
   - What collapses is the share of TRUE claims in the lists: 0.81 at effect 1, 0.035 at 1/8.
   - With 28 TRUE items among 800, a ~3% negative pass rate yields 21 false vs 26 true passes.
   - The false passes at 1/8 are 18 "other" (wrong type pair for a planted ligand) and 3 random
     [`laneB_N2_false_pass_composition`].
   - Practical reading: a per-claim filter cannot rescue a candidate list that is ~97% wrong. The gate's output is only
     as clean as the list's prevalence allows.
2. **Hard-negative lists (A5m: gate-matched TRUE + decoy senders + hard niche; exploratory).**
   [`laneB_N2_H5_A5m_e<level>_K10_pooled`]
   - Post-gate phantom rate: 0.164 at effect 1, 0.203 at 1/2, 0.465 at 1/4, 0.652 at 1/8.
   - So even with a correctly specified response, hard negatives keep the false-pass share above 0.10 at full effect.
   - This agrees with lane A's A5 miss.
   - Exact-A5 lists: 0.204, 0.364, 0.625, 0.606 [`laneB_N2_H5_A5_e<level>_K10_pooled`].
3. **Per tissue, the 1/8 level is unstable** (14 and 12 passes in the two breast tissues). Post-gate phantom rate there: breast Rep1
   0.786, Rep2 0.667, lung 0.167, pancreas 0.067 [`laneB_N2_secondary_summary`]. The pooled CI [0.154, 0.706] reflects this.
4. **Retention stays high** (0.81-0.99), but the TRUE count behind it shrinks to 28 at 1/8. The expression tools barely
   surface TRUE pairs at low effect, which is H5a.
5. **Cross-lane consistency.** Lane B's exact-A5 lists at effect 1 match lane A's A5 claim for claim: 4,600 rows, 0 label
   or pass/fail mismatches [`laneB_N2_A5_e1_consistency_with_laneA`].

## 3. Caveats
- **The effect knob scales every planted category** (TRUE, phantom, niche) together, as in Phase 1. It is not "TRUE signal
  only".
- **H5b's falsification rests on one level (1/8)** with 47 passed claims and a wide CI. At 1/4 the point estimate (0.123)
  meets the 0.15 target, but its CI reaches 0.211. The rule says "at every level", so the verdict stands.
- **The A5-style lists are high-prevalence-of-negatives by construction.** Their expression-only phantom rates (0.65-0.98)
  are not H1/H5 baseline results.
- **Effect 1 = the default N2 lists already seen in lane A** (disclosed in the pre-spec). This does not affect the
  pre-registered H5 rule.

## 4. Reproduce
```
PC_DATA=<data root> python -m scripts.laneB_run <tissue> N2     # each tissue
python -m scripts.laneB_aggregate N2                             # tables + verdict (add --log only for the first logging)
python -m scripts.laneB_figure N2
```
Per-claim decisions: `results/claims/laneB/`. Tables: `results/tables/laneB_N2_*`.
Instances are cached in `data/bench/` (effect-level default instances) and `data/bench/laneB/N2/` (A5m/A5 truth).

## 5. For the paper (lane L)
- **H5 row:** FALSIFIED. H5a met (median rho -1.0); H5b missed at effect 1/8 (0.447).
- **Mechanism sentence:** the base rate (per-negative pass rate about 3% while TRUE prevalence falls to 3.5%).
- **Failure map entry:** "low-prevalence lists". The gate filters claims; it cannot create signal the list lacks.
