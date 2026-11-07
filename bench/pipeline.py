"""Shared evaluation pipeline: benchmark instance -> tool claim lists -> gate records -> labels -> results.jsonl."""
import copy, json, os, subprocess, time
from dataclasses import asdict
import numpy as np
import torch
from bench.generate import get_or_make, label_claim
from bench.metrics import list_metrics, hier_bootstrap
from baselines.expr_tools import run_on_bench, top_k
from gate.gate import Gate, GateConfig, ExprSource

RESULTS = "results/results.jsonl"
METHODS = ("cpdb", "prod")


def git_hash():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        return "unknown"


def log_result(rec):
    rec = dict(rec)
    rec.setdefault("git_hash", git_hash())
    rec.setdefault("timestamp", time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    if torch.cuda.is_available():
        rec.setdefault("vram_peak_mb", round(torch.cuda.max_memory_allocated() / 1e6, 1))
    with open(RESULTS, "a") as f:
        f.write(json.dumps(rec, default=float) + "\n")
    return rec


def claim_key(c):
    return (c["sender"], c["receiver"], c["ligand"], c["receptor"])


def tool_lists(t, inst, Ks=(5, 10, 20), n_perm=1000, seed=0):
    df = run_on_bench(t, inst, n_perm=n_perm, seed=seed)
    out = {}
    for m in METHODS:
        for K in Ks:
            out[(m, K)] = top_k(df, m, K)[["sender", "receiver", "ligand", "receptor"]].to_dict("records")
    return out, df


def make_gate(t, inst, cfg: GateConfig):
    tg = {(d["ligand"], d["receptor"]): d["targets"] for d in inst.database}
    return Gate(t, ExprSource(t, inst.counts, inst.syn_genes, tg), cfg)


def evaluate_records(gate, claims):
    """Raw (un-decided) records for unique claims; returns dict key -> record."""
    out = {}
    for c in claims:
        k = claim_key(c)
        if k not in out:
            out[k] = gate.evaluate(c)
    return out


def decide_list(gate, recs_by_key, claims, cfg=None):
    """BH + pass rule applied over this list's claims only (a gate run gates one list)."""
    g = copy.copy(gate)
    if cfg is not None:
        g.cfg = cfg
    recs = [copy.deepcopy(recs_by_key[claim_key(c)]) for c in claims]
    return g.decide(recs)


def labels_for(inst, claims):
    out = [label_claim(inst, c) for c in claims]
    return [o[0] for o in out], [o[1] for o in out]


def summarize(groups, n_boot=2000):
    keys = ("phantom_rate", "phantom_rate_after", "retention", "rel_reduction")
    if groups and groups[0][1] is None:
        keys = ("phantom_rate",)
    point, ci = hier_bootstrap(groups, keys=keys, n_boot=n_boot)
    return point, ci
