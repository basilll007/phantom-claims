"""Lane R: prove that the output-only `status` field (gate/status.py) changes no gate decision.

1. prepare  : copy the cached default N2 instances (the Lane A base instances) into a scratch cache
              data/bench/_laneR_equiv/ so that the rerun never rewrites the shared cache.
2. run      : (per tissue) rerun the frozen gate on every Lane A N2 cell (variants A0-A5, seeds 200-209) with the
              current code, which writes `status` and `readout` next to each decision, into
              results/claims/laneR_equiv/ (env PC_BENCH / PC_LANEA_CLAIMS point the Lane A runner there).
3. compare  : join every rerun row to the committed Lane A row (results/claims/laneA/N2_*.csv) on list, rank and claim;
              decisions (label, neg_type, category, mechanism, passed) must be identical; continuous columns are
              reported as max |diff|. Also checks status == "pass" <=> passed, and tabulates status x readout.
              Writes results/tables/laneR_status_equivalence.json.

  python -m scripts.laneR_status_equivalence prepare
  PC_BENCH=data/bench/_laneR_equiv PC_LANEA_CLAIMS=results/claims/laneR_equiv python -m scripts.laneR_status_equivalence run <tissue>
  python -m scripts.laneR_status_equivalence compare
"""
import glob, json, os, shutil, sys
import numpy as np, pandas as pd

ROOT = os.environ.get("PC_DATA", "data")
SCRATCH = f"{ROOT}/bench/_laneR_equiv"
REF, NEW = "results/claims/laneA", "results/claims/laneR_equiv"
VARS = ["A0", "A1", "A2a", "A2b", "A3", "A4", "A5"]
KEY = ["seed_set", "tissue", "seed", "variant", "list", "rank", "sender", "receiver", "ligand", "receptor"]
DECISION = ["label", "neg_type", "category", "mechanism", "passed"]
CONT = ["q_bh", "rho", "p_perm", "n_expr_pairs", "effect_d"]


def prepare():
    os.makedirs(SCRATCH, exist_ok=True)
    fs = [f for f in glob.glob(f"{ROOT}/bench/*_s20[0-9]_25894fed60_*") if os.path.isfile(f)]
    for f in fs:
        shutil.copy2(f, SCRATCH)
    print(f"copied {len(fs)} base-instance files to {SCRATCH}")


def run(tissue):
    from scripts import laneA_run
    assert laneA_run.BENCH.replace("\\", "/").rstrip("/").endswith("_laneR_equiv"), "PC_BENCH must be the scratch cache"
    assert laneA_run.OUT.replace("\\", "/").rstrip("/").endswith("laneR_equiv"), "PC_LANEA_CLAIMS must be the scratch dir"
    laneA_run.run(tissue, "N2", variants=VARS)


def txt(s):
    return s.astype(object).where(s.notna(), "<NA>").astype(str)


def compare():
    fr = sorted(glob.glob(f"{REF}/N2_*.csv")); fn = sorted(glob.glob(f"{NEW}/N2_*.csv"))
    names = lambda fs: sorted(os.path.basename(f) for f in fs)
    assert names(fr) == names(fn), f"file sets differ: ref {len(fr)} new {len(fn)}"
    a = pd.concat([pd.read_csv(f) for f in fr]); b = pd.concat([pd.read_csv(f) for f in fn])
    for d in (a, b):
        d["rank"] = d["rank"].astype(int); d["seed"] = d["seed"].astype(int)
    m = a.merge(b, on=KEY, how="outer", suffixes=("_ref", "_new"), indicator=True)
    unmatched = int((m["_merge"] != "both").sum())
    m = m[m["_merge"] == "both"]
    mism = {c: int((txt(m[f"{c}_ref"]) != txt(m[f"{c}_new"])).sum()) for c in DECISION}
    cont = {c: float(np.nanmax(np.abs(m[f"{c}_ref"].astype(float) - m[f"{c}_new"].astype(float)))) for c in CONT}
    passed = m["passed_new"].astype(str).str.lower().eq("true")
    status_ok = int(((m["status"] == "pass") == passed).all())
    tab = (m.groupby(["variant", "status", "readout"]).size().rename("n").reset_index().to_dict("records"))
    by_list = (m.assign(kind=np.where(m["list"] == "planted", "planted", "lists"))
                .groupby(["kind", "status"]).size().rename("n").reset_index().to_dict("records"))
    out = dict(files_ref=len(fr), files_new=len(fn), rows_ref=len(a), rows_new=len(b), rows_matched=int(len(m)),
               rows_unmatched=unmatched, decision_mismatches=mism, continuous_max_abs_diff=cont,
               status_consistent_with_passed=bool(status_ok), status_by_variant=tab, status_by_kind=by_list,
               result="IDENTICAL" if unmatched == 0 and sum(mism.values()) == 0 and status_ok else "DIFFERENT")
    os.makedirs("results/tables", exist_ok=True)
    json.dump(out, open("results/tables/laneR_status_equivalence.json", "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ("status_by_variant",)}, indent=1, default=float))
    return out


def log():
    """Append the equivalence result (and, as a secondary field, the status of every real-data audit record) once."""
    from collections import Counter
    import scripts.log_result as lr
    from gate.status import annotate
    assert os.environ.get("PC_RESULTS"), "set PC_RESULTS to the shared results.jsonl (main tree)"
    lr.RESULTS_FILE = os.environ["PC_RESULTS"]
    have = {json.loads(l)["id"] for l in open(lr.RESULTS_FILE, encoding="utf-8") if l.strip()}
    rid = "laneR_status_equivalence_laneA_N2"
    if rid in have:
        print("already logged:", rid); return
    eq = json.load(open("results/tables/laneR_status_equivalence.json"))
    real = [json.loads(l) for l in open("results/audit/real_breast_gate.jsonl")]
    rc = Counter(f"{annotate(r)['status']}|{annotate(r)['readout']}" for r in real)
    lr.log_result(dict(
        id=rid, lane="R", agent_id="claude-opus-5.5-laneR", phase="P2", hypothesis="review response R2 (output-only)",
        note="Frozen gate rerun on every Lane A N2 cell (cached instances, variants A0-A5) with gate/status.py writing "
             "status/readout after each decision; joined to the committed Lane A rows on list, rank and claim. "
             "Secondary field real_data_status: status of the 40 Phase-1 real-data audit records "
             "(results/audit/real_breast_gate.jsonl; real tissue, no seeds; re-annotated, gate not rerun).",
        result=eq["result"], files=eq["files_new"], rows_matched=eq["rows_matched"], rows_unmatched=eq["rows_unmatched"],
        decision_mismatches=eq["decision_mismatches"], continuous_max_abs_diff=eq["continuous_max_abs_diff"],
        status_consistent_with_passed=eq["status_consistent_with_passed"], status_by_kind=eq["status_by_kind"],
        status_by_variant=eq["status_by_variant"], real_data_status=dict(rc), n_real=len(real),
        table="results/tables/laneR_status_equivalence.json"), lane="R", seed_set="N2", labelled="exploratory")
    print("logged", rid)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prepare": prepare()
    elif cmd == "run": run(sys.argv[2])
    elif cmd == "compare": compare()
    elif cmd == "log": log()
