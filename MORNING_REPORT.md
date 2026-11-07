# MORNING REPORT — Phantom Claims (rewritten every cycle)

**Status:** All phases (0–8) done; the paper draft needs a human LaTeX compile and citation and licence checks. ~95% of plan done. Last update 2026-11-07 00:59 CST.

## 1. Verdicts (pre-registered rules; K = 10; 95% hierarchical bootstrap CIs; all numbers from results/results.jsonl)
| Hyp | Verdict | Key numbers (result ids) |
|---|---|---|
| **H1** expression tools | **INCONCLUSIVE** (everywhere) | dev 0.30 [0.16, 0.46] cpdb / 0.24 prod; test: breast Rep1 0.20, Rep2 0.18, pancreas 0.12, lung 0.14 (`*_H1_expr_*_K10`). None meets "≥0.30 and CI-low ≥0.20". |
| **H1** LLM agents | **INCONCLUSIVE** (pooled) | dev 0.20 [0.07, 0.35]; test, 24 runs: 0.13 [0.08, 0.18] (Gemini 0.15, DeepSeek 0.11) (`agent_test_all_both_unguided_K10`). Per-tissue pancreas (0.08) is FALSIFIED (<0.10). |
| **H2** gate works | **SUPPORTED** (every split and condition) | tools: dev → 0.00, retention 0.97 [0.88, 1.0]; test → 0.00, retention 1.00. Agents post hoc on test: → 0.00, retention 0.995 [0.98, 1.0] |
| **H3** transfer | **SUPPORTED** | frozen gate meets H2 on held-out breast seeds, breast Rep2, **pancreatic cancer** and **lung cancer**, with no retuning (`results/tables/verdicts.json`). K=20: 0.51–0.57 → 0.00–0.02, retention 0.98–1.00 |

## 2. Top findings (surprises and negatives included)
1. **H1 failed at the pre-registered defaults; the benchmark is easy at K=10.** Expression tools rank TRUE plants above the matched never-meet decoys, and agents with spatial tools avoid most phantoms. Exploratory: at ¼ effect size the expression-only phantom rate is 0.72, which the gate reduces to 0.08 (retention 0.86), so H1 depends on effect size.
2. **The gate is calibrated even on real structured tissue, and its tile-restricted null matters.** On 200 random real-gene claims on real breast tissue, 6% reach p<0.05 (8% when the targets are sender markers) and 0% pass the full gate. With a global null the rate is 11%, and 19% for marker targets (`real_gene_null_breast_rep1`). NICHE false passes are 0.12 dev and 0.04–0.06 test (0.32 with a global null). TRUE pass rates are 0.96–1.00, and PHANTOM/RANDOM 0.
3. **Agents don't fully enforce their own verifier:** on test, 5 of 240 submitted claims (4 of them false) had been rejected by the agent's own `verify_claim`, and 13 were never verified. The phantom rate is 0.07 with the gate as a tool vs 0.00 post hoc (`agent_test_gate_tool_compliance`).
4. **Negative: the gate's power is entirely the response readout.** With 0/3 responding targets retention is 0.00–0.03; in receptor-only mode (the real-data readout) it is 0.20. On real breast tissue: 0/6 PD-L1/PD-1 claims, 0/10 expression-only top claims and **0/4 positive controls** pass. The PD-L1/PD-1 contacts are at **chance level** (obs/exp 0.53–1.58), *not depleted*, so "never meet" is **not** supported at these type pairs (`real_breast_contacts_vs_chance`).

## 3. Figures
`figures/fig_examples.png` (TRUE / PHANTOM / real PD-L1–PD-1), `figures/fig_phantom_rates.png`, `figures/fig_planted_pass_rates.png`, `figures/spatial_{breast_rep1,breast_rep2,pancreas,lung}.png`. Tables: `results/tables/verdict_table.csv`, `paper/tables/verdicts.tex`.

## 4. Blockers / what the human must do
1. **Compile `paper/main.tex` on Overleaf** with the official AAAI-27 author kit (no LaTeX here). Swap the `\documentclass` line as marked. Numbers come from `paper/numbers.tex` (auto-generated; provenance in `paper/numbers_provenance.json`).
2. **Verify citations** (`paper/refs.bib`, written from memory) and the 3 `\todo{}` items: the medical-world-model and agent citations, and the parent-thesis PD-L1/PD-1 finding, **which our analysis does not reproduce at the type-pair level**.
3. **Check dataset licences** in `paper/DATA_SOURCES.md` (marked TODO-VERIFY).
4. Optional third non-breast cancer (colorectal or ovarian, IO panel): the 10x pages returned HTTP 429, so I didn't guess URLs. Paste an official `outs.zip` URL and run `scripts/zip_remote_extract.py`.
5. Submission deadline is **Nov 15 AoE** (OpenReview link not yet announced).

## 5. Risks / possible bugs (harshest-reviewer view)
- **Independent code review (subagent, read-only):** no leakage and no bugs that inflate results. The permutation test, BH, metrics, bootstrap, verdict rules and dev/test separation were all verified. Fixed: the agent tool whitelist, a verdict edge case and a figure NaN. Post-freeze edits are documented in PREREG Amendment 3 and are behavior-preserving.
- **What H2 actually tests:** the K=10 tool lists contain 88 never-meet phantoms, 3 wrong-type-pair negatives and **0 niche** negatives (`list_negative_composition`), so H2 at K=10 mostly tests the contact and response rules. Harder negatives appear only at K=20 (124 wrong-type-pair, 9 niche) and in the planted-claim analysis. Phantoms carry no target signal, so any target check rejects them, and the gate's exposure kernel equals the generator's. Together, H2/H3 "supported" is partly by construction; the paper states all of this.
- 3 agent runs per family per dataset is small; agent H1 per tissue is noisy (CIs ±0.15).
- Pancreas and lung cell types are marker-based (Unlabeled after QC: pancreas 11.1%, lung 15.7%; `results/tables/cell_counts.csv`).
- Two amendments were made (the donor-gene range after a degenerate first smoke run, and freezing the gate config). Both are logged before any test data were generated.
- OOM crash during lung final eval: nothing had been logged, so it was rerun. Documented in PROGRESS.md.

**Added this cycle (exploratory, dev): hard-list ablation.** Per seed, a list of 10 TRUE + 10 NICHE planted claims (phantom rate 0.50). The full gate brings it to 0.11 [0.02, 0.22], retention 0.98. Without tiles it is 0.24, without the permutation test 0.22, and without split consistency 0.15 (`dev_hard_list_true_vs_niche`). This is the evidence that the spatial components matter. Lifting the niche expression cap did *not* get niche pairs into tool lists (`dev_nichecap10_*`).

## 6. Next 3 actions
1. Done: paper read-through (abstract, amendments, phantom definition, Figure 1 moved into the main text).
2. Human: compile the paper (expected ~6 pages plus appendix), verify citations and licences, and decide whether to add a third non-breast cancer.
3. Optional: more agent runs (budget used: Gemini ~1.8M and DeepSeek ~1.4M of 4M tokens each) to tighten the per-tissue CIs.

Done this cycle: the r_plant=50 sweep (retention 0.90), the independent review, the real-gene null, and an appendix with audit records, the agent protocol and the annotation rules.
