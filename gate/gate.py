"""State-grounded claim gate (plain statistics, deterministic; no LLM inside).

For a claim (A, B, L, R):
 (a) receivers i in B: exposure E_i = sum_{j in A, d_ij <= r} L_j exp(-d_ij / (r/2)); response Y_i = mean of
     z-scored listed targets (or receptor R if no targets). Statistic rho_w = within-tile Pearson correlation of
     (E, Y). Null: permute Y among B cells within the same tile (= permuting B positions within type within tile).
 (b) split consistency: rho > 0 and p < alpha_half in both spatial halves (median x of the tissue).
 (c) minimum evidence: >= n_min (A,B) pairs within r with L>0 in the sender and R>0 in the receiver;
     rho_w >= rho_min.
 BH-FDR at q over the gated claims. pass = (a: q-value <= q) & (b) & (c).
"""
import hashlib, json, time
from dataclasses import dataclass, asdict
import numpy as np
import torch
from scipy.spatial import cKDTree


@dataclass
class GateConfig:
    r: float = 30.0
    tile: float = 250.0
    n_perm: int = 1000
    alpha_half: float = 0.10
    n_min: int = 50
    rho_min: float = 0.02
    q: float = 0.05
    seed: int = 0
    use_tiles: bool = True          # ablation: False -> global permutation
    use_split: bool = True          # ablation
    use_min_evidence: bool = True   # ablation
    use_perm: bool = True           # ablation: False -> skip (a)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1.0)
    return out


class ExprSource:
    """Gene -> (lognorm vector, raw count vector) over all cells."""
    def __init__(self, tissue, counts, genes, targets=None):
        self.t = tissue
        self.counts = counts.tocsc()
        self.gi = {g: i for i, g in enumerate(genes)}
        self.scale = np.median(tissue.lib) / tissue.lib
        self.targets = targets or {}       # (L, R) -> [target genes]
        self._z = {}

    def raw(self, g):
        return np.asarray(self.counts[:, self.gi[g]].todense()).ravel()

    def log(self, g):
        return np.log1p(self.raw(g) * self.scale)

    def z(self, g):
        if g not in self._z:
            v = self.log(g); self._z[g] = (v - v.mean()) / (v.std() + 1e-9)
        return self._z[g]

    def _parts(self, g):
        return [g] if g in self.gi else g.split("_")

    def has(self, g):
        return all(x in self.gi for x in self._parts(g))

    def log_multi(self, g):
        """Complexes 'A_B': geometric-style min over subunits (CellPhoneDB convention)."""
        return np.min(np.stack([self.log(x) for x in self._parts(g)]), 0)

    def raw_multi(self, g):
        return np.min(np.stack([self.raw(x) for x in self._parts(g)]), 0)


class Gate:
    def __init__(self, tissue, expr: ExprSource, cfg: GateConfig = None, device=None):
        self.t, self.x, self.cfg = tissue, expr, cfg or GateConfig()
        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self._W = {}
        self.xmid = float(np.median(tissue.xy[:, 0]))
        c = self.cfg
        tx = np.floor(tissue.xy[:, 0] / c.tile).astype(np.int64)
        ty = np.floor(tissue.xy[:, 1] / c.tile).astype(np.int64)
        self.tile_id = tx * 100000 + ty if c.use_tiles else np.zeros(tissue.n, np.int64)
        self.half = (tissue.xy[:, 0] >= self.xmid).astype(np.int64)

    def _weights(self, A, B):
        key = (A, B)
        if key not in self._W:
            iA = np.flatnonzero(self.t.types == A); iB = np.flatnonzero(self.t.types == B)
            if len(iA) == 0 or len(iB) == 0:
                self._W[key] = (iA, iB, None)
            else:
                D = cKDTree(self.t.xy[iB]).sparse_distance_matrix(cKDTree(self.t.xy[iA]), self.cfg.r,
                                                                    output_type="coo_matrix")
                W = D.tocsr(); W.data = np.exp(-W.data / (self.cfg.r / 2))
                # sparse_distance_matrix drops exact-zero distances; acceptable (coincident centroids are rare)
                self._W[key] = (iA, iB, W)
        return self._W[key]

    def _perm_test(self, E, Y, groups, seed):
        """Within-group centred correlation and one-sided within-group permutation p-value."""
        dev = self.device
        g = torch.as_tensor(groups, device=dev)
        _, gidx = torch.unique(g, return_inverse=True)
        G = int(gidx.max()) + 1
        E_t = torch.as_tensor(E, dtype=torch.float64, device=dev)
        Y_t = torch.as_tensor(Y, dtype=torch.float64, device=dev)
        cnt = torch.bincount(gidx, minlength=G).double()
        def centre(v):
            s = torch.zeros(G, dtype=torch.float64, device=dev).index_add_(0, gidx, v)
            return v - (s / cnt)[gidx]
        Ec, Yc = centre(E_t), centre(Y_t)
        den = torch.sqrt((Ec ** 2).sum() * (Yc ** 2).sum())
        if float(den) <= 0:
            return 0.0, 1.0
        obs = float((Ec * Yc).sum() / den)
        if not self.cfg.use_perm:
            return obs, np.nan
        order0 = torch.argsort(gidx, stable=True)
        gs = gidx[order0].double()
        gen = torch.Generator(device=dev); gen.manual_seed(int(seed) % (2 ** 62))
        ge, done, chunk = 0, 0, 200
        n = len(E)
        while done < self.cfg.n_perm:
            m = min(chunk, self.cfg.n_perm - done)
            keys = gs.unsqueeze(0) * 2.0 + torch.rand((m, n), generator=gen, device=dev, dtype=torch.float64)
            perm = torch.argsort(keys, dim=1)                         # random order within each group
            src = order0[perm]                                        # (m, n) source index per sorted slot
            Yp = torch.empty((m, n), dtype=torch.float64, device=dev)
            Yp[:, order0] = Yc[src]
            null = (Yp @ Ec) / den
            ge += int((null >= obs - 1e-12).sum())
            done += m
        return obs, (1 + ge) / (1 + self.cfg.n_perm)

    def _claim_seed(self, claim):
        s = f"{claim['sender']}|{claim['receiver']}|{claim['ligand']}|{claim['receptor']}|{self.cfg.seed}"
        return int(hashlib.sha1(s.encode()).hexdigest()[:15], 16)

    def evaluate(self, claim):
        t0 = time.time()
        c = self.cfg
        A, B, L, R = claim["sender"], claim["receiver"], claim["ligand"], claim["receptor"]
        rec = dict(claim=dict(sender=A, receiver=B, ligand=L, receptor=R), tissue=self.t.name,
                   gate_config=asdict(c))
        if not (self.x.has(L) and self.x.has(R)):
            rec.update(pass_=False, reasons=["gene_not_in_state"]); return rec
        iA, iB, W = self._weights(A, B)
        rec["n_senders"], rec["n_receivers"] = int(len(iA)), int(len(iB))
        if W is None:
            rec.update(pass_=False, reasons=["cell_type_absent"]); return rec
        Lv = self.x.log_multi(L)[iA]
        E = W @ Lv
        tg = [g for g in self.x.targets.get((L, R), []) if self.x.has(g)]
        if tg:
            Y = np.mean([self.x.z(g)[iB] for g in tg], 0); rec["response"] = {"targets": tg}
        else:
            Y = self.x.log_multi(R)[iB]; rec["response"] = {"receptor": R}
        # (c) minimum evidence: expressing sender-receiver pairs within r
        Lr = self.x.raw_multi(L)[iA] > 0; Rr = self.x.raw_multi(R)[iB] > 0
        Wc = W.tocoo()
        n_pairs = int((Lr[Wc.col] & Rr[Wc.row]).sum())
        rec["n_pairs_within_r"] = int(Wc.nnz); rec["n_expressing_pairs_within_r"] = n_pairs
        seed = self._claim_seed(claim)
        rho, p = self._perm_test(E, Y, self.tile_id[iB], seed)
        rec["rho_within_tile"], rec["p_perm"] = rho, p
        # effect size: within-tile centred response, exposed vs unexposed receivers (standardised)
        exp_ = E > 0
        if exp_.any() and (~exp_).any():
            rec["effect_d"] = float((Y[exp_].mean() - Y[~exp_].mean()) / (Y.std() + 1e-9))
        else:
            rec["effect_d"] = None
        halves = []
        for h in (0, 1):
            m = self.half[iB] == h
            if m.sum() < 10:
                halves.append(dict(half=h, n=int(m.sum()), rho=None, p=None)); continue
            rh, ph = self._perm_test(E[m], Y[m], self.tile_id[iB][m], seed + h + 1)
            halves.append(dict(half=h, n=int(m.sum()), rho=rh, p=ph))
        rec["halves"] = halves
        rec["runtime_s"] = round(time.time() - t0, 3)
        return rec

    def decide(self, recs):
        """Apply BH over the claims and the pass rule. Mutates and returns recs."""
        c = self.cfg
        ps = [r.get("p_perm") for r in recs]
        idx = [i for i, p in enumerate(ps) if p is not None and not (isinstance(p, float) and np.isnan(p))]
        qv = bh([ps[i] for i in idx])
        for i, qq in zip(idx, qv):
            recs[i]["q_bh"] = float(qq)
        for r in recs:
            if "reasons" in r and r.get("pass_") is False and "rho_within_tile" not in r:
                continue
            reasons = []
            if c.use_perm and not (r.get("q_bh", 1.0) <= c.q):
                reasons.append(f"perm_null_not_rejected(q={r.get('q_bh', 1.0):.3g}>{c.q})")
            if c.use_split:
                ok = all(h["rho"] is not None and h["rho"] > 0 and h["p"] is not None and h["p"] < c.alpha_half
                         for h in r["halves"])
                if not ok:
                    reasons.append("split_inconsistent")
            if c.use_min_evidence:
                if r["n_expressing_pairs_within_r"] < c.n_min:
                    reasons.append(f"too_few_expressing_pairs({r['n_expressing_pairs_within_r']}<{c.n_min})")
                if r["rho_within_tile"] < c.rho_min:
                    reasons.append(f"effect_below_min(rho={r['rho_within_tile']:.3g}<{c.rho_min})")
            r["pass_"] = len(reasons) == 0
            r["reasons"] = reasons
        return recs

    def check(self, claims):
        recs = [self.evaluate(cl) for cl in claims]
        return self.decide(recs)


def write_audit(recs, path):
    with open(path, "w") as f:
        for r in recs:
            f.write(json.dumps(r, default=float) + "\n")
