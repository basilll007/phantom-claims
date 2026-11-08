#!/usr/bin/env bash
# One-command rebuild of the Lane A (H4) results: per-claim gate decisions, tables, confound/decoy diagnostics, figure.
# It then compares the rebuild with the committed results (scripts/laneA_compare.py; exit 1 on any mismatch).
#
#   scripts/laneA_reproduce.sh [OUT_DIR]        (default OUT_DIR=repro/laneA)
#
# Env:
#   PYTHON   interpreter (default: python)
#   PC_DATA  data root holding processed/<tissue>.h5ad (default: data)
#   PC_BENCH benchmark cache (default: $OUT_DIR/bench, i.e. EVERY instance is regenerated from processed data)
#   PAR      tissues run in parallel (default 2)
# Input: data/processed/{breast_rep1,breast_rep2,pancreas,lung}.h5ad. Building those from the raw 10x downloads is
# covered by the project-wide reproduce target (lane F): scripts/download.py -> scripts/process_xenium.py /
# scripts/annotate_markers.py. This script never logs to results/results.jsonl and never touches committed outputs.
set -euo pipefail
OUT="${1:-repro/laneA}"
PY="${PYTHON:-python}"
PAR="${PAR:-2}"
export PC_DATA="${PC_DATA:-data}"
export PC_BENCH="${PC_BENCH:-$OUT/bench}"
export PC_LANEA_CLAIMS="$OUT/claims" PC_LANEA_TABLES="$OUT/tables" PC_LANEA_FIGS="$OUT/figures"
mkdir -p "$PC_BENCH" "$PC_LANEA_CLAIMS" "$PC_LANEA_TABLES" "$PC_LANEA_FIGS" "$OUT/logs"
for t in breast_rep1 breast_rep2 pancreas lung; do
  [ -f "$PC_DATA/processed/$t.h5ad" ] || { echo "missing $PC_DATA/processed/$t.h5ad (build it first, see header)"; exit 2; }
done
echo "[laneA-repro] $(date +'%Y-%m-%d %H:%M:%S') out=$OUT bench=$PC_BENCH git=$(git rev-parse --short HEAD)"
"$PY" scripts/check_frozen.py
"$PY" -m pytest -q tests/test_bench.py tests/test_variants.py
printf '%s\n' pancreas lung breast_rep1 breast_rep2 | \
  xargs -P "$PAR" -I{} sh -c "\"$PY\" -m scripts.laneA_run run {} N2 > \"$OUT/logs/run_{}.log\" 2>&1 || { echo 'run {} FAILED'; exit 255; }"
"$PY" -m scripts.laneA_aggregate N2            > "$OUT/logs/aggregate.log" 2>&1
"$PY" -m scripts.laneA_decoy_coloc N2          > "$OUT/logs/decoy.log" 2>&1
"$PY" -m scripts.laneA_confound N2             > "$OUT/logs/confound.log" 2>&1
"$PY" -m scripts.laneA_figure N2
echo "[laneA-repro] $(date +'%Y-%m-%d %H:%M:%S') rebuild done; comparing with committed results"
"$PY" -m scripts.laneA_compare "$OUT" N2
