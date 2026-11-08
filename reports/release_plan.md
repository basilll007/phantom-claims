# Release plan (lane R, R5). PLAN ONLY: nothing is published, pushed to a public place, or uploaded.

Written 2026-11-08 by lane R (claude-opus-5.5). Decisions about licences and publication belong to the repository owner;
this file lists what would be released, how, and what is blocked.

## 1. What would be released (anonymised folder `release/phantom-claims-anon/`, built from tracked files only)
| component | source in repo | notes |
|---|---|---|
| Benchmark generator | `bench/` (`tissue.py`, `generate.py`, `variants.py`, `metrics.py`, `pipeline.py`), `configs/` | deterministic given (tissue, seed, config hash); default config hash `25894fed60` |
| Seeds and pre-registration | `PREREG.md` (base frozen at `da1b81c`, Amendments 1-4, Amendment-4 erratum) | dev 0-4, test1 100-104, N2 200-209 |
| Frozen gate | `gate/gate.py` (SHA-256 `f61003c6...7066`, checked by `scripts/check_frozen.py`), `GateConfig` defaults | r=30 um, tile=250 um, n_perm=1000, alpha_half=0.10, n_min=50, rho_min=0.02, q=0.05 |
| Audit status (output-only) | `gate/status.py` | adds `status` (pass/fail/untestable), `readout`, `status_note`; equivalence test: `results/tables/laneR_status_equivalence.json` |
| Audit-record schema | section 3 below | one JSON object per claim |
| Baselines and agent harness | `baselines/`, `agent/` | tool definitions and prompts; no keys (keys only ever came from environment variables) |
| Results | `results/results.jsonl`, `results/tables/`, `results/claims/`, `results/audit/` | every paper number is a macro from `results.jsonl` (`paper/numbers_provenance.json`) |
| Paper sources | `paper/main.tex`, `paper/numbers.tex`, `paper/refs.bib`, `paper/figures/`, `paper/tables/` | `scripts/flatten_paper.py` makes the single-file source |
| Tests and checks | `tests/`, `scripts/check_frozen.py`, `scripts/check_paper.py`, `scripts/citation_audit.py` | |

## 2. `make reproduce` (targets to be written as a Makefile; they wrap scripts that exist today)
| target | runs | needs |
|---|---|---|
| `make data` | `scripts/download.py`, `scripts/process_xenium.py`; verifies SHA-256 / CRC32 against `data/raw/MANIFEST.jsonl` | network; licence-verified sources (section 4) |
| `make check` | `pytest`, `scripts/check_frozen.py`, `scripts/check_paper.py`, `scripts/citation_audit.py` | network for the citation check |
| `make phase1` | `scripts/dev_eval.py`, `scripts/final_eval.py`, `scripts/make_report.py` | GPU optional |
| `make laneA` | `scripts/laneA_reproduce.sh` (clean rebuild + diff; last result MATCH, id `laneA_N2_rebuild_check`) | about 20 min with 2 tissues in parallel |
| `make laneB` | `scripts/laneB_run.py` per tissue, `scripts/laneB_aggregate.py`, `scripts/laneB_figure.py` | |
| `make laneR` | `scripts/laneR_sensitivity.py` (dev seeds), `scripts/laneR_sens_aggregate.py`, `scripts/laneR_status_equivalence.py` | |
| `make paper` | `scripts/build_paper.sh` (numbers -> figures -> checks -> single-file source -> pdflatex if present) | TeX Live for the PDF |
| `make reproduce` | `data check phase1 laneA laneB laneR paper` | |
LLM agent runs are not part of `make reproduce` (they cost API money and are not deterministic); their transcripts and
scores are released as recorded (section 4).

## 3. Audit-record schema (fields written by the frozen gate, plus the output-only status)
`claim` {sender, receiver, ligand, receptor}; `tissue`; `gate_config`; `n_senders`; `n_receivers`; `response`
({targets: [...]} or {receptor: R}); `n_pairs_within_r`; `n_expressing_pairs_within_r`; `rho_within_tile`; `p_perm`;
`effect_d`; `halves` [{half, n, rho, p}]; `runtime_s`; `q_bh`; `pass_`; `reasons` [...]; and from `gate/status.py`:
`status` (pass | fail | untestable), `readout` (targets | receptor_only | none), `status_note`.

## 4. Licence notes and what is BLOCKED
| item | status | blocking question |
|---|---|---|
| Code licence | **BLOCKED** | there is no LICENSE file; the owner must choose one |
| Breast Xenium Rep1/Rep2 (10x preview page; GEO GSE243280) | **BLOCKED** | licence marked TODO-VERIFY in `paper/DATA_SOURCES.md`; until verified, release download scripts + checksums only, not data |
| Pancreas and lung Xenium (10x datasets pages) | **BLOCKED** | "CC BY 4.0" so far only from secondary sources (TODO-VERIFY on the 10x pages); same rule: scripts + checksums only |
| Processed `.h5ad` and benchmark instances | **BLOCKED** | derived from the data above; same licence question |
| LIANA/OmniPath consensus resource file | **BLOCKED** | repository licence TODO-VERIFY (GPL-3.0?); plan: download at build time instead of shipping the file |
| LLM agent transcripts (`results/agent_logs/`, 1.9 MB) | **BLOCKED** | provider output terms (Gemini, DeepSeek) TODO-VERIFY before redistributing raw transcripts; scores can be released |
| AAAI-27 kit files (`aaai2027.sty/.bst`) | **BLOCKED** | redistribution terms not checked; the release can point to https://aaai.org/authorkit27/ instead |
| Everything else in section 1 | ready once the code licence is chosen | |

## 5. Anonymisation checklist (before any anonymous release)
- New history: export tracked files with `git archive` into a fresh repository with a neutral author; do not publish
  the current history (it carries the owner's account name and e-mail).
- Remove or neutralise strings found in tracked files: the GitHub owner/repository name (in `STATE.md`,
  `PROGRESS.md`, `dashboard/index.html`); local absolute paths in logs and reports.
- Exclude `.env`, `logs/`, `tools/`, `data/`, `.worktrees/`, `dashboard/` (it embeds the repository URL).
- Run `scripts/check_paper.py` and a key scan (last scan: 0 key matches in tracked files, PROGRESS 2026-11-07 14:38).
- The workshop review is single-blind (submissions need not be anonymised), so this applies only if an anonymous
  artefact is wanted.
