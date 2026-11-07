# MORNING REPORT — Phantom Claims (rewritten every cycle)

**Status:** Phase 6 (final eval) running; phases 0–5 and 7 done. ~70% of plan done. Last update 2026-11-07 00:29 CST.

## 1. Verdicts (all numbers from results/results.jsonl; K = 10; 95% hierarchical bootstrap CIs)
| Hyp | Verdict | Key numbers (result id) |
|---|---|---|
| H1 tools (dev) | **INCONCLUSIVE** | cpdb phantom 0.30 [0.16, 0.46]; prod 0.24 [0.12, 0.38] (`dev_default_H1_expr_*_K10`). Rule needs CI-low ≥ 0.20. |
| H1 tools (test breast) | **INCONCLUSIVE** | Rep1 held-out seeds 0.20 [0.08, 0.32]; Rep2 0.18 [0.06, 0.32] (`test_final_breast_rep*_H1_expr_cpdb_K10`) |
| H1 LLM agents (dev) | **INCONCLUSIVE** | both families 0.20 [0.07, 0.35]; Gemini 0.27, DeepSeek 0.13 (`agent_dev_breast_rep1_*_unguided_K10`) |
| H2 (dev) | **SUPPORTED** | tools: phantom → 0.00, retention 0.97 [0.88, 1.0]; agents post-hoc: → 0.00, retention 0.96 [0.87, 1.0] |
| H2 (test breast Rep1/Rep2) | **SUPPORTED** | phantom → 0.00, retention 1.00 (`test_final_breast_rep*_H2_expr_*_K10_none`) |
| H3 | **NOT YET RUN** | breast held-out passes; pancreas + lung final eval in progress |

## 2. Top findings (including the surprises and the negative ones)
1. **The benchmark is easier than intended at K = 10.** Expression-only tools rank planted TRUE pairs above the matched never-meet decoys (with log-normalised means, concentrating a matched count increase into fewer cells lowers the score). So H1 is *not* supported at K = 10. At K = 20 the phantom rate is 0.57 [0.48, 0.67] (test Rep1), but that is partly mechanical, since only 10 TRUE claims exist.
2. **The gate removes phantoms while keeping the true claims, and is calibrated:** null p < 0.05 for 5.0% of unplanted claims (`dev_gate_null_calibration`). Planted pass rates: TRUE 0.98, PHANTOM 0.00, NICHE 0.12. With a global instead of a tile-restricted null, NICHE rises to 0.32 and the null p < 0.05 rate to 9%, so tiles are the component that matters for hard negatives.
3. **Agents ignore their own verifier:** with `verify_claim` available, agents still submitted 5 of 60 claims the tool had rejected (3 of them were negatives) and 5 they never verified (`agent_dev_gate_tool_compliance`). The phantom rate was 0.15 with the gate as a tool vs 0.00 when the gate is applied post hoc. Enforcement has to be outside the agent.
4. **Real data, a negative result:** on real breast tissue the frozen gate passes 0/6 PD-L1/PD-1 claims (0–12 expressing contacts within 30 µm, i.e. the never-meet pattern), but it also passes **0/4 contact-dependent positive controls** (`real_breast_gate_check`). Without target priors the readout is too weak to validate positives.

## 3. Figures
`figures/spatial_{breast_rep1,breast_rep2,pancreas,lung}.png`, `figures/fig_phantom_rates.png`, `figures/fig_planted_pass_rates.png`

## 4. Blockers / what the human must do
- Colorectal and ovarian (Immuno-Oncology panel) data: the 10x dataset pages return HTTP 429 to automated fetches, and I did not guess URLs. To add a third non-breast cancer: open the 10x dataset page, copy the `outs.zip` URL into `configs/`, then run `scripts/zip_remote_extract.py`.
- LaTeX isn't installed locally: compile `paper/` on Overleaf.
- If you have the official Janesick annotations for pancreas or lung (there are none): ours are marker-based (`scripts/annotate_markers.py`).

## 5. Risks / possible bugs (harshest-reviewer view)
- H2 is near-perfect partly by construction: never-meet decoys have ~0 expressing contacts, and TRUE plants have a strong default effect (within-tile ρ ≈ 0.2–0.5). Exploratory sweeps over effect 0.5 / 0.25, noise and plant radius are running to find where it breaks.
- The gate uses the database's target lists. TRUE targets are 2/3 correct by design, which is optimistic compared with real priors.
- Gate false positives come from niche pairs (6–12%) and wrong-sender claims next to a TRUE sender.
- Pancreas and lung annotations are approximate (marker-scored XOA clusters; 17–20% Unlabeled).
- One Gemini smoke run (249k tokens) used the uncapped harness; it is excluded from metrics (`results/agent_logs/smoke/`).

## 6. Next 3 actions
1. Finish the pancreas + lung final eval, then compute the H3 verdict.
2. Run test agents (breast running; pancreas/lung next), then `scripts/agent_score.py test`.
3. Finish the knob sweeps, then a paper draft (`paper/`) generated from the results.
