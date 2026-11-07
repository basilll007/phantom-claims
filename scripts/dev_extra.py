"""Exploratory negative controls on dev (breast Rep1 seeds 0-4), frozen gate otherwise:
(1) target-prior misspecification: TRUE plants with 1 or 0 responding listed targets;
(2) receptor-only gate mode (use_targets=False) on the default benchmark (= the real-data readout)."""
import json
from gate.gate import GateConfig
from scripts.dev_eval import load_dev, h1_tools, h2

sel = json.load(open("results/tables/gate_selected_dev.json"))
t, data = load_dev()
h2(t, data, GateConfig(**sel, use_targets=False), "default", ablation="receptor_only", write=False)
for tag, bc in [("resp1", dict(n_responders=1)), ("resp0", dict(n_responders=0))]:
    t2, d2 = load_dev(bc, tag)
    h1_tools(d2, tag)
    h2(t2, d2, GateConfig(**sel), tag, write=False)
