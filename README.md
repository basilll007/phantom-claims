# Phantom Claims — claim-level governance for agentic spatial biology

Target venue: First Workshop on Medical World Models @ AAAI-27 (submission deadline 2026-11-15 AoE).

A single-cell spatial transcriptomics section (10x Xenium) is treated as a **world state**; a biological claim
*"cell type A signals to B via ligand L / receptor R"* is an **action** on that state. A deterministic, LLM-free
**claim gate** checks each claim against the state (tile-restricted spatial permutation null, split consistency,
minimum-evidence rule, BH-FDR) and writes a JSON audit record per claim.

## Dashboard

Open **`dashboard/index.html`** in a browser (self-contained; no server or internet needed). It shows the
pre-registered verdicts, the phantom rate before/after the gate for every split and condition, gate ablations, planted
pass rates, every LLM agent run (tokens, tool calls, claims, link to the full transcript), all experiments in
`results/results.jsonl`, the progress log, commits, the morning report and all figures. Append `?theme=light` or
`?theme=dark` to the URL to force a theme.

Rebuild after any new results:

```bash
python scripts/build_dashboard.py            # -> dashboard/index.html (stdlib only)
```

The GitHub Actions workflow `.github/workflows/dashboard.yml` rebuilds it on every push (artifact `dashboard-site`)
and deploys it to GitHub Pages when the repository variable `PAGES_ENABLED` is `true`.

## Current results (details: `MORNING_REPORT.md`, numbers: `results/results.jsonl`)

| Hypothesis (pre-registered, `PREREG.md`) | Verdict |
|---|---|
| H1: expression-only tools / LLM agents report >= 30% phantom claims in their top-10 | INCONCLUSIVE |
| H2: the gate cuts the phantom rate >= 50% while keeping >= 80% of true claims | SUPPORTED |
| H3: the frozen gate transfers to held-out seeds, a held-out breast section, pancreatic and lung cancer | SUPPORTED |

Read the caveats in `MORNING_REPORT.md` (benchmark difficulty, dependence on target readouts, real-data negative
results) before citing any of these.

## Repository layout

| path | contents |
|---|---|
| `PREREG.md` | frozen pre-registration + append-only amendments |
| `STATE.md`, `TASKS.md`, `PROGRESS.md`, `MORNING_REPORT.md` | loop state, task claims, append-only log, current report |
| `bench/` | tissue loader, semi-synthetic benchmark generator, metrics, pipeline |
| `gate/` | the claim gate |
| `baselines/` | expression-only CCC baselines (CellPhoneDB-style, LIANA/CellChat-style) |
| `agent/` | tool-using LLM agent harness (Gemini, DeepSeek via OpenAI-compatible APIs) |
| `scripts/` | data download / processing, evaluation, report, paper numbers, dashboard |
| `tests/` | unit tests (`pytest tests/`) |
| `results/` | `results.jsonl`, tables, claim lists, audit records, agent transcripts |
| `figures/`, `paper/` | figures; LaTeX draft whose numbers are generated from `results/` |
| `dashboard/` | the built dashboard |

## Reproduce

```bash
uv venv .venv --python 3.12
uv pip install -r requirements.lock      # torch from the CUDA 12.8 index: --index-url https://download.pytorch.org/whl/cu128
python scripts/download.py configs/data_breast.json
python scripts/download.py configs/data_lr.json
# pancreas / lung: scripts/zip_remote_extract.py <official outs.zip URL> data/raw/<tissue> cell_feature_matrix.h5 cells.csv.gz analysis.tar.gz
python scripts/process_xenium.py breast_rep1 data/raw/breast_rep1 data/raw/breast_rep1/Xenium_R1_Fig1-5_supervised.csv.gz
PYTHONPATH=. python -m pytest tests/
PYTHONPATH=. python scripts/dev_eval.py tune
```

Raw data are **not** in the repository (`data/raw/` is git-ignored); sources, URLs and checksums are in
`paper/DATA_SOURCES.md` and `data/raw/MANIFEST.jsonl` (local). API keys are read from the environment / a local
`.env` that is git-ignored; no key is stored in this repository.

## Data and licences

10x Genomics Xenium datasets (Janesick et al. 2023, Nat Commun; GEO GSE243280; 10x pancreas and lung cancer
datasets, CC BY 4.0 per the 10x pages) and the LIANA+ / OmniPath ligand-receptor resource. Verify each licence in
`paper/DATA_SOURCES.md` before redistribution. No licence has yet been chosen for this code.
