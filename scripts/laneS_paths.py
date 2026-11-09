"""Lane S I/O locations. Defaults reproduce the committed layout."""
import os

DATA = os.environ.get("PC_DATA", "data")
BENCH = os.environ.get("PC_BENCH", f"{DATA}/bench")
CLAIMS = os.environ.get("PC_LANES_CLAIMS", "results/claims/laneS")
TABLES = os.environ.get("PC_LANES_TABLES", "results/tables")
FIGS = os.environ.get("PC_LANES_FIGS", "figures")

