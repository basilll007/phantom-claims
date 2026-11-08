"""Lane A I/O locations. Defaults reproduce the committed layout. Override via environment variables so that a
rebuild (scripts/laneA_reproduce.sh) can write to a scratch directory and be compared with the committed results.

  PC_DATA          data root holding processed/<tissue>.h5ad         (default: data)
  PC_BENCH         benchmark-instance cache                          (default: $PC_DATA/bench)
  PC_LANEA_CLAIMS  per-claim gate decisions written by laneA_run     (default: results/claims/laneA)
  PC_LANEA_TABLES  aggregated tables                                 (default: results/tables)
  PC_LANEA_FIGS    figures                                           (default: figures)
"""
import os

DATA = os.environ.get("PC_DATA", "data")
BENCH = os.environ.get("PC_BENCH", f"{DATA}/bench")
CLAIMS = os.environ.get("PC_LANEA_CLAIMS", "results/claims/laneA")
TABLES = os.environ.get("PC_LANEA_TABLES", "results/tables")
FIGS = os.environ.get("PC_LANEA_FIGS", "figures")
