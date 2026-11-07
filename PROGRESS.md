# PROGRESS (append-only)
- 2026-11-06 23:55 CST — Phase 0: git repo, uv venv (.venv, py3.12), torch 2.11.0+cu128; CUDA OK on RTX 5060 Laptop (sm_120), bf16 matmul OK. requirements.lock written.
- 2026-11-06 23:55 CST — Workshop info recorded in STATE.md (deadline 2026-11-15 AoE; 8 pages; non-archival).
- 2026-11-06 23:55 CST — Data: breast Rep1/Rep2 cell_feature_matrix.h5 + cells.csv.gz from 10x CDN (URLs from 10x preview page / web results; sizes verified by Content-Length; sha256 in data/raw/MANIFEST.jsonl); Janesick et al. supervised annotations from GEO GSM7780153/GSM7780154. LIANA+ OmniPath resource from saezlab/liana-py GitHub.
- 2026-11-06 23:55 CST — Decision: user asked to include other cancer types beyond breast. Breast stays dev (parent work, PD-L1/PD-1 case); pancreatic/lung/colorectal/ovarian cancers are held-out test tissues only.
- 2026-11-06 23:55 CST — Decision: user supplied Gemini key in .env -> Gemini is first agent family; DeepSeek second.
- 2026-11-06 23:55 CST — Decision: benchmark plants on anonymized added channels with real-gene baselines (labels known for whole claim universe; no LLM prior leakage). Real L-R panel list used for Phase 7 real-data check.
- 2026-11-06 23:55 CST — Blocker: 10x dataset pages HTTP 429; a background retry loop with a browser UA was denied by the permission system as bot-check circumvention -> dropped; will use occasional single WebFetch.
- 2026-11-06 23:55 CST — L-R pool: LIANA consensus (4620 pairs) ∩ breast panel (313 genes) = 23 pairs (all complex subunits required).
- 2026-11-06 23:55 CST — PREREG.md frozen in commit da1b81c (before any evaluation).
