# TASKS (claims table)
| task | owner | status | timestamp |
|---|---|---|---|
| P0 setup env/repo/prereg | claude-main | done | 2026-11-06 23:58 |
| P1 breast data | claude-main | done | 2026-11-07 00:00 |
| P1 pancreas + lung cancer | claude-main | done | 2026-11-07 00:05 |
| P1 colorectal/ovarian IO panel | claude-main | blocked (10x 429, no guessed URLs) | 2026-11-07 00:28 CST |
| P2 bench generator + tests | claude-main | done | 2026-11-07 00:02 |
| P3 baselines + H1 dev | claude-main | done | 2026-11-07 00:12 |
| P4 gate + tuning + validity + ablations | claude-main | done | 2026-11-07 00:20 |
| P5 agent harness + dev runs | claude-main | done | 2026-11-07 00:28 CST |
| P6 final eval tools (all 4 test sets) | claude-main | done | 2026-11-07 00:39 CST (commit e91b18b; row first stamped 00:28) |
| P6 final eval agents (24+24 runs) | claude-main | done | 2026-11-07 00:45 CST (commit fcdf63f; row first stamped 00:28) |
| P7 real-data check | claude-main | done (negative on positive controls) | 2026-11-07 00:22 |
| P8 paper draft | claude-main | done (needs human compile/cite check) | 2026-11-07 00:45 CST |
| exploratory sweeps, real-gene null, hard lists, code review | claude-main | done | 2026-11-07 00:58 CST (commit 9e5786b; row first stamped 00:28) |
| Phase 2 Lane L setup (Amendment 4, check_frozen, log_result) | lane-L (gemini-lead) | done | 2026-11-07 23:57 CST (commit 216b498; row first stamped 15:00) |
| Phase 2 Lane A Misspecified benchmark (H4, variants A1-A5, N2) | lane-A (claude-opus-5.5-laneA) | done: H4 NOT MET (3/5; logged category INCONCLUSIVE); merged into main 62cd982; one-command rebuild verified (MATCH); report reports/lane_A.md | 2026-11-08 12:05 CST |
| Phase 2 Lane B Dose-response sweep (H5, levels 1..1/8, N2) | lane-B (claude-opus-5.5; second task of the lane-L agent, at the human's request) | done: H5 FALSIFIED (H5a met; H5b missed at 1/8, 0.447); merged into main 9286401; report reports/lane_B.md | 2026-11-08 13:08 CST |
| Phase 2 Lane C Real positive controls (H7, pre-spec, frozen gate) | lane-C | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane D PD-L1/PD-1 reconciliation (resolutions & radii) | lane-D | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane E Agent enforced harness (H6, models, N2) | lane-E | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane F Reproduce script, citation audit, licences | lane-F | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane G Red team audit & adversarial tests | lane-G | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane R Review response (text/definitions, output-only status field, exploratory sensitivity, future work, release plan) | lane-R (claude-opus-5.5; human pasted the lane R brief) | done: merged d13c6a4 (report reports/lane_R.md; status equivalence IDENTICAL; lane F citation sign-off and the review text itself still pending) | 2026-11-08 17:01 CST |
| Phase 2 Lane L Lead integration, 10% re-run audit, morning report & paper | lane-L (claude-opus-5.5; took over from gemini-lead at the human's request) | in progress: lane-A + lane-B merged (62cd982, 9286401); timestamp audit; DEADLINE verified; 10% re-run audit of lane A (0/28 flagged). Paper rewritten on the AAAI-27 kit with H1-H5, failure map, caveats, Amendment-4 disclosure (human asked to start before C/E); needs a real compile; H6/H7 + MORNING_REPORT wait for lanes E/C. 23:47: Figure 5 (claim fate) + appendix figures A4-A6 added | 2026-11-08 23:47 CST |
| Phase 2 Lane S Strengthening experiments (H8 spatial baselines, H9 base-rate law, 5K panel, neg controls) | lane-S (gemini-3.8-flash) | done: S0 prereg Amendment 5 committed (5f37cd3); S4 neg controls logged (0/138 passed); S1 H8 SUPPORTED (gate FPS 0.213 vs baselines 0.447-0.643, margin >= 0.10); S2 H9 SUPPORTED (Spearman rho 1.000, error 0.0058 <= 0.05, fig_base_rate.png & calculator shipped); H4 N3 replication INCONCLUSIVE (3/5 units met, replicates Lane A); S3 5K richer panels evaluated & logged (1,645 & 1,408 pairs; PD-L1 audited); report reports/lane_S.md complete | 2026-11-09 01:28 CST |
