# Lane A pre-specification: H4 misspecification robustness (fixed before any N2 data)

- Written 2026-11-08 (first committed in 95cf788 at 00:08:52 CST; an earlier draft line said '~00:30', which was wrong) by lane A (agent `claude-opus-5.5-laneA`), branch `lane-A`, parent commit 216b498
  (Amendment 4 commit). This file is committed **before** any N2 instance is generated or any gate is run on a variant.
- Governs: PREREG Amendment 4, H4. Here I only make the amendment's variant descriptions concrete; thresholds and verdict
  rules are copied from the amendment unchanged. `STATE.md` has `THRESHOLDS: default`.
- Gate: frozen `GateConfig()` (r=30, tile=250, n_perm=1000, alpha_half=0.10, n_min=50, rho_min=0.02, q=0.05).
  `scripts/laneA_run.py` asserts this and that it equals `results/tables/gate_selected_dev.json`. Gate code untouched.

## 1. Paired design
For each (tissue, seed), a **base** instance comes from the unchanged default generator (`bench/generate.py`, `DEFAULTS`).
This is the gate-matched mechanism: exponential kernel, lambda = r/2, r_plant = 30 um, 2 of 3 targets respond.
Variants (`bench/variants.py`) keep real genes, ligand and receptor channels, and every PHANTOM / NICHE / RANDOM pair
**byte-identical**. They re-plant only the TRUE pairs' target response, using the same responder genes, the same response law
`Poisson(mu_T * effect * (1 - exp(-E)))`, and a different exposure kernel `E = sum_j w(d_ij)`. Consequences:
(i) expression-only tool lists are identical across A0-A4 (unit-tested); (ii) the variants differ only in the mechanism.
Variant noise is drawn from `default_rng([seed, variant_code, 20261108])`.

| variant | mechanism of TRUE response | biological stand-in |
|---|---|---|
| A0 (reference, not part of H4) | exp kernel, r=30, lambda=15 (= default generator) | gate-matched contact/short-range signal |
| A1 kernel shape | ring/annulus: w=1 for 15 <= d <= 30 um, else 0 | influence not monotone from contact (e.g. ligand consumed by immediate neighbours) |
| A2a radius 0.5x | exp, r=15, lambda=7.5 | strictly juxtacrine (about one cell diameter) |
| A2b radius 2x | exp, r=60, lambda=30 | short-range paracrine beyond the gate's radius |
| A3 diffusive paracrine | exp, r=80, lambda=40 (amendment values); response saturates -> weak gradient | long-range diffusible signal |
| A4 receptor-only | n_responders = 0: no target responds; receptor unchanged | signalling with no transcriptional readout on the panel |
| A5 mixed + hard negatives | 10 TRUE pairs: 2 each of {ring, exp15, exp60, para80, none} (seeded permutation; no gate-matched pair) + hard negatives below | heterogeneous mechanisms in one tissue |

**A5 hard negatives.**
- (a) **Decoy senders (wrong-type-pair):** in 5 of the 10 TRUE pairs, chosen at random, the ligand is also expressed
  at the true-sender level (`Poisson(mu_L)`) by one decoy cell type. The decoy type is drawn uniformly from the eligible
  types, excluding the pair's sender and receiver; it is not filtered for distance. Receivers respond only to the true
  sender. The claim (decoy, receiver, L, R) is a negative of sub-type `decoy_sender`.
- (b) **Hard niche:** the 10 NICHE pairs are re-planted with the Phase-1 niche mechanism (L, R and all targets driven by
  a shared regional field). The field is broader and stronger: sigma=600 um, 3 blobs, sender/receiver scale =
  min(2.0 / mean field, 20).
- **Pilot on dev seeds, done before this file.** Seeds 0-2, breast Rep1, tool lists only; the gate was **never** run in the pilot.
  The listed configs, in order, were: c1 (sigma 400, 2 blobs, x1), c2 (400, 2, x2), c3 (600, 3, x1), c4 (600, 3, x2).
  The aim was the first config that puts niche negatives into K=10 lists. c1-c3 put 0 niche items into K=10 lists;
  c4 put 3.3 niche + 1.7 decoy items per K=10 list (table `results/tables/laneA_pilot_dev_composition.csv`). So c4 is fixed.
  Phantom items mostly drop out of A5 K=10 lists, displaced by the hard negatives. Phantom rejection is therefore tested in A0-A4.
- Note: in A5, the expression-only phantom rate is high **by construction**, because hard negatives were made to reach the
  lists. It must not be read as an H1/H5-type baseline result.

## 2. Evaluation
- Seeds: N2 = 200..209 (all 10). Tissues: breast Rep1, breast Rep2, pancreas, lung (all 4, pooled). The runner refuses
  test1 seeds and seeds outside the declared set.
- Tool lists: expression-only `cpdb` and `prod` (Phase-1 baselines, unchanged), top-K, K in {5, 10, 20}. Tool
  permutation seed = instance seed.
- The gate is evaluated on every unique claim. Each list is decided separately: BH over that list's claims (as in Phase-1 H2).
- **Primary metric (H4):** K=10, both tools pooled, i.e. each (tissue, seed) group = the cpdb and prod K=10 lists.
  - `retention` = fraction of TRUE claims in the lists that pass the gate.
  - `phantom rate` (post-gate) = fraction of passed claims that are not TRUE (any negative sub-type; 0 if nothing passes).
    The definitions are those of `bench/metrics.list_metrics`.
  - Point estimates are pooled over 4 tissues x 10 seeds. 95% CIs come from a hierarchical bootstrap: (tissue, seed)
    groups, then claims, 2000 resamples (`bench/metrics.hier_bootstrap`).
- **Verdict (Amendment 4, unchanged):** a variant *meets* H4 if retention >= 0.60 and post-gate phantom rate <= 0.10.
  - SUPPORTED if >= 4 of 5 variants (A1-A5) meet it.
  - FALSIFIED if retention < 0.50 in >= 2 variants or phantom rate > 0.15 in >= 2 variants.
  - Otherwise INCONCLUSIVE.
- **Operationalisation of A2** (two radii, one variant): A2 meets H4 only if **both** A2a and A2b meet it. For the
  falsification counts, A2 contributes retention = min(A2a, A2b) and phantom rate = max(A2a, A2b). This is the conservative choice.
- Verdicts use point estimates; CIs are always reported next to them. No p-values, so no multiplicity correction is
  needed for the verdict. All variants and both A2 radii are reported regardless of outcome.
- **Secondary / exploratory (reported, never relabelled confirmatory):**
  - per tissue, per tool, and K=20;
  - A0 reference;
  - planted-claim pass rates per category (true / phantom / niche / niche_hard / random / decoy_sender), with BH over
    all planted claims of an instance;
  - A5 TRUE and decoy pass rates per mechanism;
  - list composition;
  - the failing-reason breakdown for TRUE claims.

## 3. Run-once rule
- `scripts/laneA_run.py` prints only timings and row counts, never metrics.
- N2 is evaluated once. After a crash, only missing (tissue, seed, variant) CSVs are re-run.
- No variant, threshold, kernel or list parameter is changed after N2 data exist. Bug fixes are logged in PROGRESS.md
  with the reason. If a fix changes any N2 output, the run is redone in full and the change is reported.
- Aggregation (`scripts/laneA_aggregate.py`) logs to `results/results.jsonl` only when the design is complete
  (4 tissues x 10 seeds x 7 variants). Confirmatory records = A1-A5 pooled K=10 records + the H4 verdict.

## 4. Expectations written before the data (so they cannot be adjusted afterwards)
- **A4 is expected to fail retention by construction.** The gate tests an exposure-dependent response, and A4 has none.
  Phase-1 exploratory dev already showed retention 0.00-0.03 with 0/3 responding targets. So H4 can be SUPPORTED only
  if A1, A2 (both radii), A3 and A5 all meet the thresholds.
- **A5 retention ceiling is about 0.8,** because 2 of its 10 TRUE pairs are receptor-only.
  If A4 and A5 both fall below 0.50 retention, H4 is FALSIFIED by the amendment rule. That outcome will be reported as such.
- **Decoy senders can pass the gate.** This can happen when the decoy type co-localises with the true sender near the
  receivers. A correlation-based spatial test cannot separate co-localised senders. This is reported as a limitation,
  not engineered away.
