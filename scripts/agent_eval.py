"""Run LLM agent conditions. Resumable (skips runs whose log exists).
Usage: python scripts/agent_eval.py <split> <family> <conditions,comma> <tissue:seed,seed> [<tissue:seeds> ...]
e.g.   python scripts/agent_eval.py dev gemini unguided,gate_tool breast_rep1:0,1,2"""
import json, os, sys, time
from bench.tissue import load_tissue
from bench.generate import get_or_make
from baselines.expr_tools import run_on_bench
from bench.pipeline import make_gate
from gate.gate import GateConfig
from agent.harness import Env, run_agent, over_budget

split, family, conds = sys.argv[1], sys.argv[2], sys.argv[3].split(",")
sel = json.load(open("results/tables/gate_selected_dev.json"))
for spec in sys.argv[4:]:
    tissue, seeds = spec.split(":")
    t = load_tissue(tissue)
    for s in [int(x) for x in seeds.split(",")]:
        inst = get_or_make(t, s)
        df = None
        for cond in conds:
            run_id = f"{split}_{tissue}_s{s}_{family}_{cond}"
            if os.path.exists(f"results/agent_logs/{run_id}.json"):
                print("skip", run_id, flush=True); continue
            if over_budget(family):
                print("BUDGET EXHAUSTED", family, flush=True); sys.exit(0)
            if df is None:
                df = run_on_bench(t, inst, seed=s)
            gate = make_gate(t, inst, GateConfig(**sel)) if cond == "gate_tool" else None
            env = Env(t, inst, gate=gate, expr_df=df)
            log = run_agent(env, family, cond, run_id, temperature=0.2, seed=s)
            print(run_id, log.get("status"), log.get("usage"), log.get("n_tool_calls"), log.get("runtime_s"), flush=True)
