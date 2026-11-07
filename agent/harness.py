"""Tool-using LLM agent over a benchmark tissue (OpenAI-compatible API: Gemini, DeepSeek).
Tools are structured, deterministic functions over the state (no free-form code execution: safer on a shared
workstation; documented in PROGRESS.md). Conditions: 'unguided' and 'gate_tool' (agent may call verify_claim).
Post-hoc gating of unguided claims is done offline. Full transcripts -> results/agent_logs/.
API keys are read at runtime (agent/llm_keys.py) and never logged."""
import json, os, time, traceback
import numpy as np
from scipy.spatial import cKDTree
from openai import OpenAI
from agent.llm_keys import get_key, ENDPOINTS

MODELS = {"gemini": "gemini-2.5-flash", "deepseek": "deepseek-flash"}
# Conservative list-price ESTIMATES (USD per 1M tokens) used only for the budget cap.
PRICE = {"gemini": (0.30, 2.50), "deepseek": (0.50, 2.00)}
BUDGET = {"tokens": 4_000_000, "usd": 5.0}
BUDGET_FILE = "results/agent_logs/budget.json"
MAX_STEPS = 12          # model turns
MAX_TOOL_CALLS = 24
K = 10

SYSTEM = ("You are a computational biology agent analysing a spatial transcriptomics section (10x Xenium, "
          "single-cell resolution, coordinates in micrometres). Cells are annotated with cell types. A ligand-"
          "receptor (L-R) database is provided; gene names are anonymized (LIG_xx / REC_xx / TGT_xx_y) and each "
          "pair lists putative downstream target genes. Use the tools to analyse the tissue. Be efficient: at most "
          f"{MAX_TOOL_CALLS} tool calls in total. Finish by calling submit_claims exactly once.")
TASK = (f"Report the top-{K} cell-cell interactions in this tissue with evidence. Each interaction is "
        "(sender cell type, receiver cell type, ligand, receptor) with the L-R pair taken from the database. "
        f"Call submit_claims with exactly {K} claims ranked from most to least confident, each with a one-line "
        "evidence string.")
GATE_HINT = (" A verification tool, verify_claim, tests a claim against the spatial state (permutation test of "
             "receiver target response vs. ligand exposure, split consistency, minimum evidence). Only submit "
             "claims you have verified if possible.")


def _schema(name, desc, props, req):
    return {"type": "function", "function": {"name": name, "description": desc,
            "parameters": {"type": "object", "properties": props, "required": req}}}


STR = {"type": "string"}
TOOLS = [
    _schema("list_cell_types", "Cell types and their cell counts.", {}, []),
    _schema("list_lr_database", "The ligand-receptor database: pair_id, ligand, receptor, listed targets.", {}, []),
    _schema("expression_summary", "Per cell type mean log-normalised expression and fraction of expressing cells "
            "for up to 12 genes.", {"genes": {"type": "array", "items": STR},
                                     "cell_types": {"type": "array", "items": STR}}, ["genes"]),
    _schema("top_expressed_pairs", "Expression-only ranking (CellPhoneDB-like): for every L-R pair and sender/"
            "receiver type, score = mean(L in sender)+mean(R in receiver) over 2; returns the top n.",
            {"n": {"type": "integer"}}, []),
    _schema("spatial_proximity", "Number of sender-receiver cell pairs within radius (um) and enrichment vs. "
            "random placement.", {"sender": STR, "receiver": STR, "radius_um": {"type": "number"}},
            ["sender", "receiver"]),
    _schema("colocalization", "Counts of ligand-expressing sender cells and receptor-expressing receiver cells "
            "that lie within radius of each other.", {"sender": STR, "receiver": STR, "ligand": STR,
                                                       "receptor": STR, "radius_um": {"type": "number"}},
            ["sender", "receiver", "ligand", "receptor"]),
    _schema("compare_near_far", "Mean expression of a gene in receiver cells within radius of any sender cell vs. "
            "receiver cells farther away.", {"sender": STR, "receiver": STR, "gene": STR,
                                             "radius_um": {"type": "number"}}, ["sender", "receiver", "gene"]),
    _schema("submit_claims", f"Submit the final ranked list of {K} interactions.",
            {"claims": {"type": "array", "items": {"type": "object", "properties": {
                "sender": STR, "receiver": STR, "ligand": STR, "receptor": STR, "evidence": STR},
                "required": ["sender", "receiver", "ligand", "receptor"]}}}, ["claims"]),
]
VERIFY = _schema("verify_claim", "Verify one claim against the spatial state; returns pass/fail and reasons.",
                 {"sender": STR, "receiver": STR, "ligand": STR, "receptor": STR},
                 ["sender", "receiver", "ligand", "receptor"])


class Env:
    def __init__(self, t, inst, gate=None, expr_df=None):
        self.t, self.inst, self.gate = t, inst, gate
        self.types = t.type_list()
        C = inst.counts
        self.raw = C.toarray()
        self.log = t.lognorm(C).toarray()
        self.gi = inst.gene_index()
        self.expr_df = expr_df
        self.by_type = {k: np.flatnonzero(t.types == k) for k in self.types}
        self._trees = {}

    def _tree(self, k):
        if k not in self._trees:
            self._trees[k] = cKDTree(self.t.xy[self.by_type[k]])
        return self._trees[k]

    def _chk_types(self, *ks):
        bad = [k for k in ks if k not in self.by_type]
        if bad:
            raise ValueError(f"unknown cell type(s) {bad}; valid: {self.types}")

    def _chk_genes(self, *gs):
        bad = [g for g in gs if g not in self.gi]
        if bad:
            raise ValueError(f"unknown gene(s) {bad}")

    def list_cell_types(self):
        return {k: int(len(v)) for k, v in self.by_type.items()}

    def list_lr_database(self):
        return self.inst.database

    def expression_summary(self, genes, cell_types=None):
        genes = genes[:12]; self._chk_genes(*genes)
        cts = cell_types or self.types; self._chk_types(*cts)
        out = {}
        for k in cts:
            idx = self.by_type[k]
            out[k] = {g: [round(float(self.log[idx, self.gi[g]].mean()), 3),
                          round(float((self.raw[idx, self.gi[g]] > 0).mean()), 3)] for g in genes}
        return {"format": "gene: [mean_log_expr, frac_expressing]", "by_cell_type": out}

    def top_expressed_pairs(self, n=20):
        n = int(min(max(n, 1), 50))
        d = self.expr_df[(self.expr_df.fracL >= 0.1) & (self.expr_df.fracR >= 0.1)]
        d = d.sort_values("score_cpdb", ascending=False).head(n)
        return [dict(sender=r.sender, receiver=r.receiver, ligand=r.ligand, receptor=r.receptor,
                     score=round(float(r.score_cpdb), 3), perm_p=round(float(r.p_cpdb), 4)) for r in d.itertuples()]

    def spatial_proximity(self, sender, receiver, radius_um=30.0):
        self._chk_types(sender, receiver); r = float(min(max(radius_um, 5), 200))
        iA, iB = self.by_type[sender], self.by_type[receiver]
        n_pairs = int(self._tree(receiver).count_neighbors(self._tree(sender), r))
        dens = len(iB) / (np.ptp(self.t.xy[:, 0]) * np.ptp(self.t.xy[:, 1]))
        exp_ = len(iA) * dens * np.pi * r * r
        frac = float((self._tree(receiver).query(self.t.xy[iA], k=1, distance_upper_bound=r)[0] < np.inf).mean())
        return dict(n_sender_receiver_pairs=n_pairs, expected_if_random=round(float(exp_), 1),
                    enrichment=round(n_pairs / max(exp_, 1e-9), 3), frac_senders_with_receiver_within_r=round(frac, 3))

    def colocalization(self, sender, receiver, ligand, receptor, radius_um=30.0):
        self._chk_types(sender, receiver); self._chk_genes(ligand, receptor); r = float(min(max(radius_um, 5), 200))
        iA, iB = self.by_type[sender], self.by_type[receiver]
        Lp = iA[self.raw[iA, self.gi[ligand]] > 0]; Rp = iB[self.raw[iB, self.gi[receptor]] > 0]
        if len(Lp) == 0 or len(Rp) == 0:
            return dict(n_ligand_pos_senders=int(len(Lp)), n_receptor_pos_receivers=int(len(Rp)), n_pairs_within_r=0)
        tL, tR = cKDTree(self.t.xy[Lp]), cKDTree(self.t.xy[Rp])
        n = int(tR.count_neighbors(tL, r))
        fL = float((tR.query(self.t.xy[Lp], k=1, distance_upper_bound=r)[0] < np.inf).mean())
        fR = float((tL.query(self.t.xy[Rp], k=1, distance_upper_bound=r)[0] < np.inf).mean())
        return dict(n_ligand_pos_senders=int(len(Lp)), n_receptor_pos_receivers=int(len(Rp)), n_pairs_within_r=n,
                    frac_Lpos_senders_near_Rpos_receiver=round(fL, 3), frac_Rpos_receivers_near_Lpos_sender=round(fR, 3))

    def compare_near_far(self, sender, receiver, gene, radius_um=30.0):
        self._chk_types(sender, receiver); self._chk_genes(gene); r = float(min(max(radius_um, 5), 200))
        iB = self.by_type[receiver]
        d = self._tree(sender).query(self.t.xy[iB], k=1, distance_upper_bound=r)[0]
        near = d < np.inf
        v = self.log[iB, self.gi[gene]]
        return dict(n_near=int(near.sum()), n_far=int((~near).sum()),
                    mean_near=round(float(v[near].mean()), 4) if near.any() else None,
                    mean_far=round(float(v[~near].mean()), 4) if (~near).any() else None)

    def verify_claim(self, sender, receiver, ligand, receptor):
        c = dict(sender=sender, receiver=receiver, ligand=ligand, receptor=receptor)
        rec = self.gate.decide([self.gate.evaluate(c)])[0]
        return {k: rec.get(k) for k in ["pass_", "reasons", "rho_within_tile", "p_perm",
                                        "n_expressing_pairs_within_r", "halves"]}

    def call(self, name, args):
        if name == "verify_claim" and self.gate is None:
            raise ValueError("verify_claim not available")
        fn = getattr(self, name)
        return fn(**args)


def _budget():
    if os.path.exists(BUDGET_FILE):
        return json.load(open(BUDGET_FILE))
    return {f: {"prompt": 0, "completion": 0, "usd_est": 0.0, "runs": 0} for f in MODELS}


def _save_budget(b):
    json.dump(b, open(BUDGET_FILE, "w"), indent=1)


def over_budget(family):
    b = _budget()[family]
    return (b["prompt"] + b["completion"]) >= BUDGET["tokens"] or b["usd_est"] >= BUDGET["usd"]


def run_agent(env, family, condition, run_id, temperature=0.2, seed=0):
    os.makedirs("results/agent_logs", exist_ok=True)
    if over_budget(family):
        return dict(status="BUDGET_EXHAUSTED")
    client = OpenAI(api_key=get_key(family), base_url=ENDPOINTS[family], timeout=180, max_retries=3)
    tools = TOOLS + ([VERIFY] if condition == "gate_tool" else [])
    msgs = [{"role": "system", "content": SYSTEM},
            {"role": "user", "content": TASK + (GATE_HINT if condition == "gate_tool" else "")}]
    log = dict(run_id=run_id, family=family, model=MODELS[family], condition=condition, temperature=temperature,
               seed=seed, tissue=env.t.name, bench_seed=env.inst.seed, started=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
               steps=[], usage={"prompt": 0, "completion": 0})
    claims, status = None, "no_submit"
    t0 = time.time()
    forced = False
    for step in range(MAX_STEPS + 3):
        if over_budget(family):
            status = "BUDGET_EXHAUSTED"; break
        kw = dict(model=MODELS[family], messages=msgs, tools=tools, temperature=temperature)
        n_calls = sum(len(s_.get("tool_calls", [])) for s_ in log["steps"])
        if step >= MAX_STEPS or n_calls >= MAX_TOOL_CALLS:
            kw["tool_choice"] = {"type": "function", "function": {"name": "submit_claims"}}
            if not forced:
                msgs.append({"role": "user", "content": "Tool budget reached. Call submit_claims now."})
                forced = True
        try:
            resp = client.chat.completions.create(**kw, seed=seed)
        except Exception as e:
            try:
                if "tool_choice" in str(e) and "tool_choice" in kw:
                    kw.pop("tool_choice")       # e.g. DeepSeek thinking mode: no forced tool_choice
                resp = client.chat.completions.create(**kw)
            except Exception as e2:
                log["steps"].append(dict(step=step, error=f"{type(e2).__name__}: {str(e2)[:500]}"))
                status = "api_error"; break
        u = resp.usage
        pt, ct = (u.prompt_tokens or 0), (u.completion_tokens or 0)
        log["usage"]["prompt"] += pt; log["usage"]["completion"] += ct
        b = _budget(); b[family]["prompt"] += pt; b[family]["completion"] += ct
        b[family]["usd_est"] += pt / 1e6 * PRICE[family][0] + ct / 1e6 * PRICE[family][1]; _save_budget(b)
        m = resp.choices[0].message
        entry = dict(step=step, content=m.content, tool_calls=[], usage=[pt, ct])
        am = {"role": "assistant", "content": m.content or ""}
        if m.tool_calls:
            am["tool_calls"] = [{"id": tc.id, "type": "function",
                                 "function": {"name": tc.function.name, "arguments": tc.function.arguments}}
                                for tc in m.tool_calls]
        msgs.append(am)
        if not m.tool_calls:
            log["steps"].append(entry)
            msgs.append({"role": "user", "content": "Please continue using the tools and finish with submit_claims."})
            continue
        done = False
        for tc in m.tool_calls:
            name = tc.function.name
            try:
                args = json.loads(tc.function.arguments or "{}")
            except Exception:
                args = {}
            if name == "submit_claims":
                claims = args.get("claims", [])
                out = {"status": "received", "n": len(claims)}; done = True
            else:
                try:
                    out = env.call(name, args)
                except Exception as e:
                    out = {"error": str(e)[:300]}
            txt = json.dumps(out, default=float)[:6000]
            entry["tool_calls"].append(dict(name=name, args=args, result=txt))
            msgs.append({"role": "tool", "tool_call_id": tc.id, "content": txt})
        log["steps"].append(entry)
        if done:
            status = "submitted"; break
    log.update(status=status, claims=claims, runtime_s=round(time.time() - t0, 1),
               n_tool_calls=sum(len(s.get("tool_calls", [])) for s in log["steps"]))
    b = _budget(); b[family]["runs"] += 1; _save_budget(b)
    json.dump(log, open(f"results/agent_logs/{run_id}.json", "w"), indent=1, default=str)
    return log
