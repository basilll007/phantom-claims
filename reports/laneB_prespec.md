# Lane B pre-specification: H5 dose-response (fixed before any effect < 1 N2 data)

- Written 2026-11-08 by lane B (agent `claude-opus-5.5-laneB`; the same agent that ran lane A and is lane L), branch
  `lane-B`, parent `563d129`. The commit of this file precedes the first N2 instance with effect < 1. The effect = 1
  default N2 lists were already seen through lane A's A0 reference: expression-only phantom rate 0.190, post-gate 0.000
  (`laneA_N2_H4_A0_K10_pooled`). That is disclosed here; it does not affect the H5 rule, which was fixed in Amendment 4.
- Governs PREREG Amendment 4, H5. Levels, seeds, tissues and thresholds are copied from the amendment; only the
  operational details below are added. `STATE.md`: `THRESHOLDS: default`.
- Gate: frozen `GateConfig()`, asserted by the runner. No gate or generator parameter is changed or tuned. N2 is
  pre-registered for H5; nothing is tuned on it.

## 1. Design
- Effect levels e in {1, 0.5, 0.25, 0.125}. Each level is the unchanged default generator with `effect = e`, which scales
  the planted Poisson means of every planted category (TRUE ligand/receptor/targets, phantom, niche). This is the same
  knob Phase 1 used for its exploratory effect sweeps. e = 1 is the cached default N2 instance.
- Seeds N2 = 200..209, tissues breast Rep1, breast Rep2, pancreas, lung: 40 (tissue, seed) units x 4 levels.
- Expression-only tools: `cpdb` and `prod` (Phase-1 baselines, unchanged), top-K with K in {5, 10, 20}; tool permutation
  seed = instance seed. The gate decides each list separately (BH within list), as in H2/H4.

## 2. H5 (confirmatory; "default" list type)
- **Unit statistic.** For each unit, the expression-only phantom rate at each level is the non-TRUE share of the unit's
  pooled cpdb+prod K=10 lists (20 items). rho_unit = Spearman(effect level, phantom rate) over the 4 levels, with average
  ranks for ties. If the phantom rate is constant over the 4 levels, rho_unit := 0 (no monotone increase; conservative).
  rho_unit and the post-gate rates are rounded to 12 decimals before any threshold comparison. Spearman over 4 levels is
  discrete, and the dev smoke test showed an exact -0.8 computed as -0.7999999999999999 (fixed before any N2 data;
  `tests/test_laneB.py`).
- **H5a** ("increases monotonically as effect size falls, Spearman rho <= -0.80 over levels, per seed"): met if the
  **median of rho_unit over the 40 units <= -0.80**. Also reported, without changing the verdict: the fraction of units
  with rho_unit <= -0.80, min/max rho_unit, and every unit's rho.
- **H5b** ("gated phantom rate stays <= 0.15 at every level"): pooled post-gate phantom rate (K=10, cpdb+prod, 40 units)
  <= 0.15 at each of the 4 levels. The post-gate phantom rate is the non-TRUE share of passed claims, and 0 when nothing passes
  (`bench/metrics.list_metrics`). Because a gate that passes nothing trivially "meets" H5b, the number of passed claims and
  the retention at each level are always reported next to it.
- **Verdict (Amendment 4, unchanged):** SUPPORTED if H5a and H5b; FALSIFIED if median rho_unit > -0.50 or post-gate phantom
  rate > 0.20 at any level; otherwise INCONCLUSIVE. Point estimates decide; 95% hierarchical bootstrap CIs
  ((tissue, seed) groups, then claims, 2000 resamples) are reported with every rate.
- Also reported for every level: precision@K (= 1 - expression-only phantom rate), retention and relative reduction with CIs,
  per-tissue values, K = 5 and K = 20, list composition, planted-claim pass rates per category.

## 3. A5-style lists (addendum; exploratory, reported separately, never part of the H5 verdict)
- **A5m** (`bench/variants.py`): the level-e default instance, with gate-matched TRUE response (exp kernel, r = 30)
  plus the lane-A A5 hard negatives (5 decoy senders, 10 hard-niche pairs; same construction and parameters as lane A),
  all scaled by the same effect e. This separates effect size from mechanism misspecification.
- **A5** (exact lane-A construction: mixed mechanisms + hard negatives) at every level. At e = 1 this is identical to
  lane A's A5, so it is also a consistency check against `laneA_N2_H4_A5_K10_pooled`.
- Same metrics as above, including false-pass composition by negative type.
- The A5-style expression-only phantom rates are high by construction (hard negatives strengthened to reach the lists).

## 4. Run-once rule
- `scripts/laneB_run.py` prints timings only. N2 is run once; after a crash, only missing
  (tissue, seed, level, list type) files are re-run.
- `scripts/laneB_aggregate.py --log` logs only for the complete design: 40 units x 4 levels x 3 list types.
  Confirmatory records = the 4 default-list level records, the unit-Spearman record and the H5 verdict.
- Dev seed 0 (breast Rep1) is used only for a code smoke test. It changes nothing here.

## 5. Written-down expectations (before the data)
- **Phase-1 exploratory dev sweeps** (expression-only K=10 phantom rate rising steeply as effect falls; retention
  falling) make H5a likely.
- **H5b's main risk is the low-effect end.** At e = 0.125 the gate may pass almost nothing; then H5b is "met" trivially.
  This will be stated next to the verdict, with retention and n_pass.
- **A5-style lists:** lane A found 0.204 post-gate phantom rate for A5 at e = 1 on N2. A5m at e = 1 is expected in the
  same range. No threshold applies; this is a description.
