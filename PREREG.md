# PREREG — Phantom Claims: claim-level governance for agentic spatial biology

**Frozen:** 2026-11-06 23:58 CST (local). Frozen before any evaluation code was written or run.
The git hash of the commit that adds this file is recorded in `STATE.md` and `PROGRESS.md`
(a file cannot contain its own commit hash). **This file is never edited after that commit; only the
`## Amendments` section at the bottom may be appended to, with timestamp and reason.**

## 0. Framing
Tissue = world state (cells with positions, types, expression). A claim
`(sender type A, receiver type B, ligand L, receptor R)` = an action proposed on that state.
A state-grounded verification gate tests every claim against the state before it is reported and
writes a JSON audit record. The gate is plain statistics (no LLM inside), deterministic given a seed.

## 1. Hypotheses
- **H1 (problem is real):** On the semi-synthetic benchmark, expression-only tools and an unguided
  LLM agent report phantom (= not planted-true) interactions in >= 30% of their top-K claims.
- **H2 (gate works):** Applied post hoc to a condition's top-K claims, the gate reduces the phantom
  rate by >= 50% (relative) while retaining >= 80% of the planted-true claims in that list.
- **H3 (transfers):** The gate frozen after dev tuning meets the H2 criteria on held-out seeds of breast
  Rep1, and on >= 2 other public Xenium cancer tissues, with no retuning.
- **Falsification / negative-result triggers:** baseline phantom rate < 10%; or gate retention < 70%
  (loses > 30% of true claims); or H2 criteria fail on held-out data. Any of these is reported as-is.

## 2. Benchmark (semi-synthetic, `bench/`)
Real cells, positions, cell types and real expression are kept. Planted signal lives on **added,
anonymized gene channels** (`LIG_k`, `REC_k`, `TGT_k_m`) so that (i) every claim in the evaluation
universe has a known label and (ii) LLMs cannot use memorized priors. Each channel's baseline is a real
panel gene's counts permuted **within cell type** (realistic per-type level/sparsity, no spatial signal).
Per seed (defaults): 10 TRUE, 10 PHANTOM, 10 NICHE, 10 RANDOM pairs; each pair has 1 ligand, 1 receptor,
3 listed targets.
- **TRUE:** type pair (A,B) with >= 500 A-B cell pairs within r_plant. +Poisson(mu_L) ligand counts in
  A cells, +Poisson(mu_R) receptor counts in B cells; for each B cell, exposure
  E_i = sum_{j in A, d_ij <= r_plant} exp(-d_ij/lambda); 2 of 3 listed targets get
  +Poisson(mu_T * (1 - exp(-E_i))) (dose-response by distance); the 3rd listed target is a non-responder
  (imperfect prior).
- **PHANTOM ("never meet"):** (A,B) preferring segregated pairs. Ligand raised only in A cells with no B
  cell within r_seg = 3 r_plant, receptor only in B cells with no A cell within r_seg; per-cell means
  scaled so the type-mean increase matches TRUE (multiplier capped at 3). No target effect.
- **NICHE-CONFOUNDED:** region field rho(x) in [0,1] = sum of 2 Gaussian blobs (sigma 400 um), ligand in
  A and receptor in B scaled by rho, all 3 listed targets raised in all cells by rho (shared regional
  program, no signaling). Label negative (hard).
- **RANDOM:** no planting (baseline only).
- Knobs: effect multiplier, r_plant, noise (binomial thinning of planted counts), #true vs #decoy, panel
  donor subset. Ground truth saved with each benchmark instance.
- **Claim universe:** all 40 planted pairs x all ordered pairs of distinct cell types (excluding
  `Unlabeled`). Positive iff (A,B,L,R) equals a TRUE plant. Everything else is negative; negatives are
  sub-labelled phantom / niche / random / other (wrong type pair).

## 3. Conditions
1. `expr_cpdb`: CellPhoneDB-style. Log-normalised means; claim valid if >= 10% of A express L and >= 10%
   of B express R; score = (mean_L(A) + mean_R(B))/2; p-value from 1000 cell-label permutations; rank
   claims with p < 0.05 by score.
2. `expr_prod`: LIANA/CellChat-style. score = mean_L(A) * mean_R(B), same validity/p-value filter.
3. `agent_unguided` (LLM, tool-using): asked "Report the top-K cell-cell interactions in this tissue with
   evidence." Models: Gemini (`gemini-2.5-flash`) and DeepSeek (`deepseek-flash`); temperature 0.2;
   3 runs per dataset per model.
4. `*_gated`: the same claim lists filtered post hoc by the gate.
5. (optional) `agent_gate_tool`: the gate is available to the agent as a tool.
- **K = 10 primary**, K = 5 and 20 secondary.

## 4. The gate (`gate/`), per claim (A,B,L,R)
- (a) **Spatial statistic + tile-restricted permutation null.** For receiver cells i in B:
  E_i = sum_{j in A, d_ij <= r} L_j exp(-d_ij/(r/2)); response Y_i = mean z-scored log expression of the
  pair's listed targets (receptor if no targets are listed, e.g. real data). Statistic = within-tile
  Pearson correlation of (E, Y) (both centred within square tiles). Null: permute B cells' responses among
  B cells **within the same tile** (= permuting positions within type within tile), n_perm = 1000,
  one-sided p (enrichment).
- (b) **Split consistency:** same statistic in the two spatial halves (split at median x; for real breast
  data also Rep1 vs Rep2): positive sign and p_half < alpha_half in both.
- (c) **Minimum evidence:** >= n_min (A,B) cell pairs within r with L > 0 in A and R > 0 in B; within-tile
  correlation >= rho_min.
- BH-FDR at q over the claims being gated. Pass = BH-significant (a) AND (b) AND (c).
- Output: one JSON audit record per claim (claim, statistics, p-values, effect sizes, pass/fail, reasons).
- **Pre-registered defaults:** r = 30 um, tile = 250 um, n_perm = 1000, alpha_half = 0.10, n_min = 50,
  rho_min = 0.02, q = 0.05.
- **Dev tuning grid (dev only):** r in {20,30,50}, tile in {150,250,400}, n_min in {20,50,100},
  rho_min in {0,0.02,0.05}. Selection rule: lowest post-gate phantom rate pooled over the two expression
  tools at K = 10 among configs with retention >= 0.90; ties -> closest to defaults. The selected config is
  frozen and recorded as an amendment before any test data are touched.

## 5. Metrics
Phantom rate@K = fraction of reported claims that are not planted-TRUE. Precision@K = 1 - phantom rate.
Recall = TRUE plants reported / TRUE plants. Retention = TRUE claims passing gate / TRUE claims in the
list. Relative reduction = 1 - PR_after / PR_before (PR_after over claims that pass; if nothing passes,
PR_after = 0 and retention = 0). CIs: 95% hierarchical bootstrap (resample seeds, then claims within
seed), 2000 resamples. Breakdown of phantoms by negative category.

## 6. Splits (dev vs test)
- **Dev (tuning allowed):** breast Rep1, benchmark seeds 0,1,2,3,4. Agent dev runs on seeds 0,1,2.
- **Test (touched once, final evaluation only):** breast Rep1 seeds 100-104 (held-out seeds); breast Rep2
  seeds 100-104 (held-out section, same block); >= 2 other public 10x Xenium cancer tissues, seeds 100-104
  each. Agent test runs on seeds 100-102 of each test dataset.
- Test benchmark instances are not generated until the gate config is frozen.

## 7. Verdict rules
- H1 per condition: SUPPORTED iff phantom rate@10 >= 0.30 and its 95% CI lower bound >= 0.20;
  FALSIFIED iff point < 0.10; else INCONCLUSIVE.
- H2 per condition: SUPPORTED iff relative reduction >= 0.50 and retention >= 0.80 (point estimates),
  CIs reported; FALSIFIED iff retention < 0.70 or reduction < 0.25; else INCONCLUSIVE.
- H3: SUPPORTED iff H2 criteria hold (both expression tools) on held-out breast seeds AND on >= 2 other
  tissues; else NOT SUPPORTED (report which).
- Known limitation, stated up front: phantoms are matched to TRUE plants in expression by design, so H1 on
  the benchmark measures whether standard expression-only scoring/filters separate them; the real-data
  PD-L1/PD-1 check (Phase 7) is the non-synthetic evidence.

## 8. Agent budget
Per model family: hard cap 4,000,000 total tokens and USD 5 estimated at list price; the agent condition
stops when either is hit and this is reported.

## Amendments
(append-only; timestamp + reason)

### Amendment 1 — 2026-11-07 00:07 CST
**What:** benchmark donor pool restricted to real panel genes with overall mean 0.02–0.30 counts/cell and
max per-type mean <= 0.5 counts/cell (generator knobs donor_max_mean, donor_max_type_mean); previously
any gene with mean >= 0.05. Gate: complex-subunit splitting on '_' only applies when the full name is not a
gene (bug fix: anonymized names contain '_').
**Why:** first dev smoke run (seed 0, before any gate tuning or gate metric was inspected) showed the
expression-only top-10 contained zero planted pairs: donor genes with strong type-specific baselines (e.g.
DST in myoepithelium) dominated type means, so top-K was all wrong-type-pair claims and retention would be
undefined. The restricted range matches real L-R genes on the panel (e.g. CD274 0.037, CTLA4 0.13, EGFR 0.30).
**Effect on hypotheses/verdict rules:** none. Seeds/splits unchanged.

### Amendment 2 — 2026-11-07 00:13 CST
**Gate config frozen (dev tuning complete):** r = 30 um, tile = 250 um, n_min = 50, rho_min = 0.02,
n_perm = 1000, alpha_half = 0.10, q = 0.05 — identical to the pre-registered defaults (selection rule s.4:
many grid configs tied at post-gate phantom rate 0 with retention >= 0.90; tie broken toward defaults).
Recorded in results/tables/gate_selected_dev.json and results.jsonl id=dev_gate_tuning. No test benchmark
instance had been generated at this point. Ablation bug fix (no effect on the frozen gate): with use_perm=False
the split check uses only the sign of rho in both halves.

### Amendment 3 — 2026-11-07 00:53 CST
**Documentation of post-freeze code edits (no change to the frozen gate's behaviour).** After Amendment 2,
gate/gate.py gained (i) an exploratory flag use_targets (default True = frozen behaviour) used only for the
receptor-only control, and (ii) a change to the split rule that applies only when use_perm=False (ablation).
An independent read-only code review (subagent, this session) confirmed via git diff that both are behaviour-preserving
for the frozen configuration, so final-eval results (including the lung re-run after an OOM crash) are unaffected.
Also post hoc, after all agent runs: tool-name whitelist in agent/harness.py (security hardening; no logged run called an
undeclared tool) and report-script fixes (float rounding in verdict thresholds; PREREG-conformant FALSIFIED rule when
retention is NaN; missing values no longer drawn as 0). None of these changes any reported verdict.
