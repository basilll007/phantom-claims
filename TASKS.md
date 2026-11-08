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
| Phase 2 Lane B Dose-response sweep (H5, levels 1..1/8, N2) | lane-B | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane C Real positive controls (H7, pre-spec, frozen gate) | lane-C | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane D PD-L1/PD-1 reconciliation (resolutions & radii) | lane-D | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane E Agent enforced harness (H6, models, N2) | lane-E | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane F Reproduce script, citation audit, licences | lane-F | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane G Red team audit & adversarial tests | lane-G | pending | 2026-11-07 23:57 CST (commit 216b498; first stamped 15:00) |
| Phase 2 Lane L Lead integration, 10% re-run audit, morning report & paper | lane-L (claude-opus-5.5; took over from gemini-lead at the human's request; was idle since setup) | in progress: merge lane-A, commit shared files, fix timestamps; paper edits wait for B/C/E | 2026-11-08 12:03 CST |
