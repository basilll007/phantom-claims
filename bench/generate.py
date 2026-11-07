"""Semi-synthetic benchmark generator: plant TRUE / PHANTOM / NICHE / RANDOM ligand-receptor pairs on
anonymized added gene channels of a real tissue. Real cells, positions, types and real counts are untouched.

Usage: python -m bench.generate <tissue> <seed> [config.json]
"""
import hashlib, json, os, sys, time
from dataclasses import dataclass, asdict, field
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from bench.tissue import Tissue, load_tissue, type_pair_contacts

DEFAULTS = dict(n_true=10, n_phantom=10, n_niche=10, n_random=10, n_targets=3, n_responders=2,
                r_plant=30.0, lam=None, mu_L=1.0, mu_R=1.0, mu_T=2.0, effect=1.0, noise=0.0,
                r_seg_mult=3.0, phantom_cap=3.0, niche_sigma=400.0, niche_blobs=2, niche_cap=3.0,
                min_type_cells=200, min_true_contacts=500, min_far_frac=0.3, donor_min_mean=0.05)
CATS = ("true", "phantom", "niche", "random")


def cfg_hash(cfg):
    return hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:10]


@dataclass
class BenchInstance:
    tissue: str
    seed: int
    cfg: dict
    syn_genes: list                 # channel names
    baseline: sp.csr_matrix         # (n, c) int counts
    added: sp.csr_matrix            # (n, c) planted counts
    truth: list                     # per-pair dicts incl. category (hidden from tools)
    database: list                  # what tools/agents see: pair_id, ligand, receptor, targets
    meta: dict = field(default_factory=dict)

    @property
    def counts(self):
        return (self.baseline + self.added).tocsr()

    def gene_index(self):
        return {g: i for i, g in enumerate(self.syn_genes)}

    def claim_universe(self, types):
        """All (A,B,L,R) with A != B over database pairs; label 1 iff planted TRUE."""
        true_set = {(p["sender"], p["receiver"], p["ligand"], p["receptor"]) for p in self.truth
                    if p["category"] == "true"}
        cat = {(p["ligand"], p["receptor"]): p for p in self.truth}
        out = []
        for d in self.database:
            p = cat[(d["ligand"], d["receptor"])]
            for A in types:
                for B in types:
                    if A == B:
                        continue
                    key = (A, B, d["ligand"], d["receptor"])
                    lab = int(key in true_set)
                    if lab:
                        sub = "true"
                    elif p["category"] != "true" and (A, B) == (p["sender"], p["receiver"]):
                        sub = p["category"]
                    elif p["category"] == "random":
                        sub = "random"
                    else:
                        sub = "other"
                    out.append(dict(sender=A, receiver=B, ligand=d["ligand"], receptor=d["receptor"],
                                    label=lab, neg_type=None if lab else sub))
        return out


def label_claim(inst, claim):
    """Ground-truth label and negative sub-type for a single claim dict."""
    cat = {(p["ligand"], p["receptor"]): p for p in inst.truth}
    p = cat.get((claim["ligand"], claim["receptor"]))
    if p is None:
        return 0, "not_in_db"
    same_pair = (claim["sender"], claim["receiver"]) == (p["sender"], p["receiver"])
    if p["category"] == "true" and same_pair:
        return 1, "true"
    if p["category"] != "true" and same_pair:
        return 0, p["category"]
    if p["category"] == "random":
        return 0, "random"
    return 0, "other"


def _permute_within_type(col, types, rng):
    out = col.copy()
    for t in np.unique(types):
        idx = np.flatnonzero(types == t)
        out[idx] = col[rng.permutation(idx)]
    return out


def _nearest_dist(src_xy, dst_xy, ub):
    tree = cKDTree(dst_xy)
    d, _ = tree.query(src_xy, k=1, distance_upper_bound=ub)
    return d  # inf if none within ub


def generate(t: Tissue, seed: int, cfg=None) -> BenchInstance:
    c = dict(DEFAULTS)
    c.update(cfg or {})
    lam = c["lam"] or c["r_plant"] / 2
    rng = np.random.default_rng(seed)
    t0 = time.time()
    types = t.type_list(min_cells=c["min_type_cells"])
    T = len(types)
    tl, C, Enr = type_pair_contacts(t, c["r_plant"], types)
    by_type = {k: np.flatnonzero(t.types == k) for k in types}
    r_seg = c["r_seg_mult"] * c["r_plant"]

    # --- candidate type pairs
    true_c = [(i, j) for i in range(T) for j in range(T)
              if i != j and C[i, j] >= c["min_true_contacts"] and Enr[i, j] >= 1.0]
    far = {}
    ph_c = []
    for i in range(T):
        for j in range(T):
            if i == j:
                continue
            fa = np.isinf(_nearest_dist(t.xy[by_type[types[i]]], t.xy[by_type[types[j]]], r_seg)).mean()
            fb = np.isinf(_nearest_dist(t.xy[by_type[types[j]]], t.xy[by_type[types[i]]], r_seg)).mean()
            far[(i, j)] = (fa, fb)
            if fa >= c["min_far_frac"] and fb >= c["min_far_frac"]:
                ph_c.append((i, j))
    ph_c.sort(key=lambda p: Enr[p])                     # prefer segregated pairs
    all_c = [(i, j) for i in range(T) for j in range(T) if i != j]

    def pick(cands, n, prefer_first=None):
        cands = list(cands)
        if prefer_first:
            pool = cands[:max(n * 2, prefer_first)]
        else:
            pool = cands
        if len(pool) == 0:
            pool = all_c
        replace = len(pool) < n
        idx = rng.choice(len(pool), size=n, replace=replace)
        return [pool[k] for k in idx]

    plan = ([("true", p) for p in pick(true_c, c["n_true"])] +
            [("phantom", p) for p in pick(ph_c, c["n_phantom"], prefer_first=c["n_phantom"] * 2)] +
            [("niche", p) for p in pick(all_c, c["n_niche"])] +
            [("random", p) for p in pick(all_c, c["n_random"])])
    order = rng.permutation(len(plan))                    # anonymize: ids carry no category info
    nT = c["n_targets"]
    syn_genes, truth = [], []
    real_means = np.asarray(t.counts.mean(0)).ravel()
    donor_pool = np.flatnonzero(real_means >= c["donor_min_mean"])
    n = t.n
    nch = len(plan) * (2 + nT)
    base = np.zeros((n, nch), dtype=np.int32)
    add = np.zeros((n, nch), dtype=np.int32)
    keep = 1.0 - c["noise"]
    eff = c["effect"]

    for slot, k in enumerate(order):
        cat, (i, j) = plan[k]
        pid = f"P{slot + 1:02d}"
        A, B = types[i], types[j]
        Lg, Rg = f"LIG_{pid[1:]}", f"REC_{pid[1:]}"
        tg = [f"TGT_{pid[1:]}_{x}" for x in "abc"[:nT]]
        col0 = len(syn_genes)
        syn_genes += [Lg, Rg] + tg
        donors = rng.choice(donor_pool, size=2 + nT, replace=False)
        for q, g in enumerate(donors):
            col = np.asarray(t.counts[:, g].todense()).ravel().astype(np.int32)
            base[:, col0 + q] = _permute_within_type(col, t.types, rng)
        iA, iB = by_type[A], by_type[B]
        rec = dict(pair_id=pid, category=cat, sender=A, receiver=B, ligand=Lg, receptor=Rg, targets=tg,
                   responders=[], donor_genes=[str(t.genes[g]) for g in donors],
                   contacts=float(C[i, j]), enrichment=float(Enr[i, j]))
        if cat == "true":
            add[iA, col0] = rng.poisson(c["mu_L"] * eff, size=len(iA))
            add[iB, col0 + 1] = rng.poisson(c["mu_R"] * eff, size=len(iB))
            treeA = cKDTree(t.xy[iA])
            D = cKDTree(t.xy[iB]).sparse_distance_matrix(treeA, c["r_plant"], output_type="coo_matrix")
            E = np.zeros(len(iB))
            np.add.at(E, D.row, np.exp(-D.data / lam))
            resp = sorted(rng.choice(nT, size=c["n_responders"], replace=False).tolist())
            for q in resp:
                add[iB, col0 + 2 + q] = rng.poisson(c["mu_T"] * eff * (1 - np.exp(-E)))
            rec["responders"] = [tg[q] for q in resp]
            rec["n_receivers_exposed"] = int((E > 0).sum())
        elif cat == "phantom":
            dA = _nearest_dist(t.xy[iA], t.xy[iB], r_seg)
            dB = _nearest_dist(t.xy[iB], t.xy[iA], r_seg)
            fA, fB = np.isinf(dA), np.isinf(dB)
            sA = min(len(iA) / max(fA.sum(), 1), c["phantom_cap"])
            sB = min(len(iB) / max(fB.sum(), 1), c["phantom_cap"])
            add[iA[fA], col0] = rng.poisson(c["mu_L"] * eff * sA, size=fA.sum())
            add[iB[fB], col0 + 1] = rng.poisson(c["mu_R"] * eff * sB, size=fB.sum())
            rec.update(n_senders_far=int(fA.sum()), n_receivers_far=int(fB.sum()), r_seg=r_seg)
        elif cat == "niche":
            centers = t.xy[rng.choice(np.r_[iA, iB], size=c["niche_blobs"], replace=False)]
            d2 = ((t.xy[:, None, :] - centers[None]) ** 2).sum(-1)
            rho = np.exp(-d2 / (2 * c["niche_sigma"] ** 2)).max(1)
            sA = min(1 / max(rho[iA].mean(), 1e-6), c["niche_cap"])
            sB = min(1 / max(rho[iB].mean(), 1e-6), c["niche_cap"])
            add[iA, col0] = rng.poisson(c["mu_L"] * eff * sA * rho[iA])
            add[iB, col0 + 1] = rng.poisson(c["mu_R"] * eff * sB * rho[iB])
            for q in range(nT):
                add[:, col0 + 2 + q] = rng.poisson(c["mu_T"] * eff * rho)
            rec.update(niche_centers=centers.tolist(), niche_sigma=c["niche_sigma"])
        if keep < 1.0:
            sl = slice(col0, col0 + 2 + nT)
            add[:, sl] = rng.binomial(add[:, sl], keep)
        truth.append(rec)

    truth.sort(key=lambda p: p["pair_id"])
    database = [dict(pair_id=p["pair_id"], ligand=p["ligand"], receptor=p["receptor"], targets=p["targets"])
                for p in truth]
    meta = dict(cfg_hash=cfg_hash(c), types=types, gen_seconds=round(time.time() - t0, 1),
                n_true_candidates=len(true_c), n_phantom_candidates=len(ph_c))
    return BenchInstance(tissue=t.name, seed=seed, cfg=c, syn_genes=syn_genes,
                         baseline=sp.csr_matrix(base), added=sp.csr_matrix(add), truth=truth,
                         database=database, meta=meta)


def save(inst: BenchInstance, root="data/bench"):
    os.makedirs(root, exist_ok=True)
    stem = f"{root}/{inst.tissue}_s{inst.seed}_{inst.meta['cfg_hash']}"
    sp.save_npz(stem + "_baseline.npz", inst.baseline)
    sp.save_npz(stem + "_added.npz", inst.added)
    json.dump(dict(tissue=inst.tissue, seed=inst.seed, cfg=inst.cfg, syn_genes=inst.syn_genes,
                   truth=inst.truth, database=inst.database, meta=inst.meta),
              open(stem + "_truth.json", "w"), indent=1)
    return stem


def load(stem):
    j = json.load(open(stem + "_truth.json"))
    return BenchInstance(tissue=j["tissue"], seed=j["seed"], cfg=j["cfg"], syn_genes=j["syn_genes"],
                         baseline=sp.load_npz(stem + "_baseline.npz").tocsr(),
                         added=sp.load_npz(stem + "_added.npz").tocsr(),
                         truth=j["truth"], database=j["database"], meta=j["meta"])


def get_or_make(t: Tissue, seed, cfg=None, root="data/bench"):
    c = dict(DEFAULTS); c.update(cfg or {})
    stem = f"{root}/{t.name}_s{seed}_{cfg_hash(c)}"
    if os.path.exists(stem + "_truth.json"):
        return load(stem)
    inst = generate(t, seed, cfg)
    save(inst, root)
    return inst


if __name__ == "__main__":
    tissue, seed = sys.argv[1], int(sys.argv[2])
    cfg = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else None
    inst = generate(load_tissue(tissue), seed, cfg)
    print(save(inst), inst.meta)
