"""Phase 7 — real-data check (no planting) on breast Rep1 and Rep2 with the frozen gate.
Response = receptor expression in receivers (no target priors on a 313-gene panel). Claims:
(i) expression-only (CellPhoneDB-like) top-10 over the panel L-R pool, (ii) the PD-L1/PD-1 case (CD274-PDCD1 and
PDCD1LG2-PDCD1) at the 3 sender/receiver type pairs the expression tool scores highest, (iii) weak positive
controls: contact-dependent pairs CD80-CTLA4, CD86-CTLA4 (APC -> T), MMRN2-CLEC14A / MMRN2-CD93 (endothelial).
Pass in both replicates = replicate-consistent."""
import json
import numpy as np, pandas as pd
from bench.tissue import load_tissue
from baselines.expr_tools import score_all, top_k
from bench.pipeline import log_result
from gate.gate import Gate, GateConfig, ExprSource, write_audit

sel = json.load(open("results/tables/gate_selected_dev.json"))
cfg = GateConfig(**sel)
lr = pd.read_csv("data/processed/lr_pairs_breast_rep1.csv")
assert not lr.ligand.str.contains("_").any() and not lr.receptor.str.contains("_").any()
pairs = list(zip(lr.ligand, lr.receptor))
genes = sorted(set(lr.ligand) | set(lr.receptor))
tis = {r: load_tissue(f"breast_{r}") for r in ("rep1", "rep2")}
tool = {}
for r, t in tis.items():
    gi = {g: i for i, g in enumerate(t.genes)}
    raw = t.counts[:, [gi[g] for g in genes]]
    tool[r] = score_all(t, t.lognorm(raw).toarray(), raw.toarray(), genes, pairs, n_perm=1000, seed=0)
    tool[r].to_csv(f"results/tables/real_{r}_expr_scores.csv", index=False)

claims = {}
for _, row in top_k(tool["rep1"], "cpdb", 10).iterrows():
    claims[(row.sender, row.receiver, row.ligand, row.receptor)] = "expr_top10_rep1"
for L, R in [("CD274", "PDCD1"), ("PDCD1LG2", "PDCD1")]:
    d = tool["rep1"][(tool["rep1"].ligand == L) & (tool["rep1"].receptor == R)].sort_values("score_cpdb", ascending=False)
    for _, row in d.head(3).iterrows():
        claims.setdefault((row.sender, row.receiver, L, R), "PD-L1/PD-1 case")
for L, R, A, B in [("CD80", "CTLA4", "Macrophages_1", "CD4+_T_Cells"), ("CD86", "CTLA4", "Macrophages_1", "CD4+_T_Cells"),
                   ("CD86", "CTLA4", "LAMP3+_DCs", "CD4+_T_Cells"), ("MMRN2", "CLEC14A", "Endothelial", "Endothelial"),
                   ("MMRN2", "CD93", "Endothelial", "Perivascular-Like")]:
    if A == B:
        continue      # gate defines sender != receiver types
    claims.setdefault((A, B, L, R), "positive_control_contact")
cl = [dict(sender=k[0], receiver=k[1], ligand=k[2], receptor=k[3]) for k in claims]
rows, audits = [], []
for r, t in tis.items():
    g = Gate(t, ExprSource(t, t.counts, t.genes, {}), cfg)
    recs = g.check(cl)
    for c, rec in zip(cl, recs):
        rec["replicate"] = r; rec["claim_source"] = claims[(c["sender"], c["receiver"], c["ligand"], c["receptor"])]
        audits.append(rec)
        e = tool[r][(tool[r].sender == c["sender"]) & (tool[r].receiver == c["receiver"]) &
                     (tool[r].ligand == c["ligand"]) & (tool[r].receptor == c["receptor"])]
        rows.append(dict(replicate=r, source=rec["claim_source"], **c, expr_score=float(e.score_cpdb.iloc[0]) if len(e) else None,
                         expr_p=float(e.p_cpdb.iloc[0]) if len(e) else None, gate_pass=rec["pass_"],
                         rho=rec.get("rho_within_tile"), p_perm=rec.get("p_perm"), q_bh=rec.get("q_bh"),
                         n_expr_pairs=rec.get("n_expressing_pairs_within_r"), reasons=";".join(rec.get("reasons", []))))
write_audit(audits, "results/audit/real_breast_gate.jsonl")
df = pd.DataFrame(rows)
df.to_csv("results/tables/real_breast_gate.csv", index=False)
wide = df.pivot_table(index=["source", "sender", "receiver", "ligand", "receptor"], columns="replicate",
                      values=["gate_pass", "rho", "n_expr_pairs"], aggfunc="first")
print(wide.to_string())
both = df.groupby(["sender", "receiver", "ligand", "receptor", "source"]).gate_pass.all().reset_index()
summ = both.groupby("source").gate_pass.agg(["sum", "count"]).to_dict("index")
log_result(dict(id="real_breast_gate_check", phase=7, split="real", tissue="breast_rep1+rep2", gate_config=sel,
                metrics={k: dict(pass_both=int(v["sum"]), n=int(v["count"])) for k, v in summ.items()},
                pdl1_claims=df[df.source == "PD-L1/PD-1 case"][["replicate", "sender", "receiver", "ligand", "receptor",
                                                                 "gate_pass", "rho", "n_expr_pairs", "reasons"]].to_dict("records")))
print(json.dumps(summ, indent=1))
