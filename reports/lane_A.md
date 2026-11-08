# Lane A report: H4 misspecification robustness
Agent `claude-opus-5.5-laneA`, branch `lane-A`. Updated 2026-11-08 11:32 CST. Status: **N2 confirmatory run complete and
logged.** Every number below comes from `results/results.jsonl` (record ids are given in brackets).

## 1. Verdict
**H4: NOT MET (confirmatory, N2).** 3 of 5 variants meet retention >= 0.60 and post-gate phantom rate <= 0.10:
A1 ring, A2 (both radii) and A3 paracrine. 4 were required. A4 misses on retention and A5 misses on phantom rate.
The robustness claim H4 makes is therefore not supported and must be reported as a miss.
Formal Amendment-4 category, as logged: INCONCLUSIVE. That means "not SUPPORTED, and not FALSIFIED": the FALSIFIED rule
needs 2 variants failing on the same criterion, and here exactly 1 fails on each. The category is a bookkeeping
detail; it does not soften the miss. [`laneA_N2_H4_verdict`]

Design: K=10 lists of the two expression-only tools (cpdb, prod), pooled; 4 tissues (breast Rep1, Rep2, pancreas, lung)
x N2 seeds 200-209; frozen gate; 95% hierarchical bootstrap CI. [`laneA_N2_H4_<v>_K10_pooled`]

| variant (mechanism of TRUE response) | expr-only phantom rate (before gate) | retention [95% CI] | post-gate phantom rate [95% CI] | meets H4 |
|---|---|---|---|---|
| A0 gate-matched reference (exploratory) | 0.190 | 0.991 [0.976, 1.000] | 0.000 [0, 0] | (ref) |
| A1 ring 15-30 um | 0.190 | 0.968 [0.938, 0.990] | 0.000 [0, 0] | yes |
| A2a exp r=15 | 0.190 | 1.000 [1.000, 1.000] | 0.000 [0, 0] | yes |
| A2b exp r=60 | 0.190 | 0.895 [0.849, 0.933] | 0.000 [0, 0] | yes (A2 = both) |
| A3 paracrine r=80, lambda=40 | 0.190 | 0.799 [0.742, 0.851] | 0.000 [0, 0] | yes |
| A4 receptor-only (no responding target) | 0.190 | **0.017 [0.004, 0.036]** | 0.000 [0, 0] (11 passes, all TRUE) | **no** |
| A5 mixed mechanisms + hard negatives | 0.656 (by construction) | 0.778 [0.697, 0.855] | **0.204 [0.119, 0.288]** | **no** |

Figure: `figures/fig_laneA_h4_N2.png`. Tables: `results/tables/laneA_N2_*.csv`. Per-claim gate decisions:
`results/claims/laneA/N2_*.csv`.

## 2. Findings, negatives first
1. **Hard negatives defeat the 0.10 phantom-rate target (A5).** In A5's K=10 lists, 55 of 269 passed claims are false:
   42 hard-niche and 13 decoy-sender claims [`laneA_N2_false_pass_composition`]. The miss holds in every tissue, with pooled CI [0.119, 0.288].
   Per-tissue post-gate phantom rate: breast Rep1 0.081, Rep2 0.222, lung 0.364, pancreas 0.241 [`laneA_N2_H4_secondary_summary`].
   Planted-claim pass rates: hard niche 0.105 [0.065, 0.150], decoy sender 0.140 [0.075, 0.210]. For comparison, the
   Phase-1-strength niche passes 0.062 [0.035, 0.095], phantoms 0.000 and random 0.005 [`laneA_N2_planted_pass_rates`].
2. **Decoy senders pass when they co-localise with the true sender.** This is an attribution limit and not a bug.
   - Pass rate by co-localisation: 0.07 when the decoy's exposure is nearly independent of the true sender's
     (co-localisation <= 0.1, n=158) vs 0.44 at 0.1-0.3 (n=34).
   - Median co-localisation is 0.144 in passed decoys vs -0.022 in failed ones (one-sided Mann-Whitney p = 1.0e-6)
     [`laneA_N2_A5_decoy_colocalisation`].
   - 79% of passed decoys are also supported by raw counts [`laneA_N2_normalisation_confound`]. The response is real; it is
     attributed to the wrong sender. A correlational spatial test cannot separate co-localised candidate senders.
3. **No transcriptional readout means no detection (A4).** This was predicted before the data (prespec section 4).
   - Retention is 0.017, and the gate fails 637 TRUE list items, mainly on split-half inconsistency and the permutation null
     [`laneA_N2_true_fail_reasons`].
   - At K=20, nearly nothing TRUE passes, so the few passing negatives raise the post-gate phantom rate to 0.111
     [`laneA_N2_H4_secondary_summary`].
4. **A normalisation x density confound gives false passes (exploratory; this is a property of the frozen gate+generator).**
   - Library size correlates with local cell density; the sign depends on tissue (Spearman vs number of cells within 30 um:
     breast Rep1 -0.39, Rep2 -0.35, lung -0.22, pancreas +0.18).
   - The gate's response readout is log1p(count x median_lib / lib), so the normalisation factor co-varies with exposure.
   - A4 TRUE claims that pass anyway (5 of 400 planted) have exposure that tracks the receivers' size factor much more
     closely than in A4 claims that fail (median rho 0.172 vs 0.039).
   - Only 17% of passing hard-niche claims are supported by raw counts. By contrast, gate-matched TRUE passes keep their
     signal on raw counts (median rho 0.455 raw vs 0.441 gate).
   - All from [`laneA_N2_normalisation_confound`]; first seen on a dev smoke seed [`laneA_dev_s3_normalisation_confound`].
5. **Kernel shape and radius misspecification are tolerated, with a range-dependent cost.**
   - Retention falls with mechanism range: 1.000 at r=15, 0.991 matched, 0.895 at r=60, 0.799 at r=80.
   - Every failing TRUE item in A2b and A3 fails the split-half rule (68/68, 130/130). About half also fail the
     permutation null [`laneA_N2_true_fail_reasons`].
   - Biological reading: long-range signals saturate the response within the gate's 30 um window, so the within-tile
     gradient flattens.
   - Planted TRUE pass rates: A1 0.948, A2b 0.858, A3 0.762 [`laneA_N2_planted_pass_rates`].

## 3. Caveats (for the paper's limitations)
- **Post-gate phantom rate = 0 in A0-A4 tests little.** The default K=10 lists again contain no niche items:
  648 TRUE, 150 phantom, 2 "other" over 80 lists [`laneA_N2_list_composition`]. Phantoms carry no target signal, so the
  gate always rejects them. **A5 is the only real phantom-rate test in H4, and it fails.**
- **A5's expression-only phantom rate (0.656) is high by construction.** Hard negatives were strengthened so that they
  reach the lists (dev pilot, gate not run: `laneA_dev_pilot_A5_hard_negative_composition`). It is not an H1/H5 result.
- **Breast Rep1 was the Phase-1 tuning tissue.** The seeds are new, but there is mild optimism for that tissue.
  The A5 miss is smallest on Rep1 (0.081).
- **Two code hashes.** Claim rows carry `d20f44b` or `bbe1b6b`; runner code is identical between them (only the
  prespec text and the figure script changed).

## 4. Proposals (future work; NOT gate changes; the gate stays frozen)
- **Readout normalisation.** Residualise the response on the log size factor and local density within tile, or permute
  within density strata of each tile. This targets finding 4. Lane C may see the same effect on real data.
- **Sender attribution.** Test a claimed sender while conditioning on exposure from other ligand-expressing types
  (partial correlation or a multi-sender regression). This targets finding 2.

## 5. For other lanes
- **Lane E / B:** default-config N2 instances for all 4 tissues and seeds 200-209 are cached at
  `data/bench/<tissue>_s<seed>_25894fed60_{baseline,added}.npz` + `_truth.json`.
  Lane-A variant ground truth (truth JSON + planted counts) is in `data/bench/laneA/N2/`.
- **Lanes G / L:** the Amendment-4 text says "2026-11-07 15:00 CST" but its commit `216b498` is dated 23:57:32. The
  commit time is authoritative, and it precedes all lane-A N2 data (first N2 file about 00:12 on 11-08).
  Lane A's own first PROGRESS entries had a wrong time stamp; a correction was appended.

## 6. Next
1. Spare time, all exploratory: per-tissue A5 breakdown by negative type; K=20 composition; more N2-disjoint seeds are
   not needed for the verdict.
2. Re-run the A5 hard-niche analysis with the normalisation-free readout as a *proposal demo* on dev seeds only,
   clearly labelled.
