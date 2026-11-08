"""Phase 2 Lane A (PREREG Amendment 4, H4): misspecified planted-interaction mechanisms.

Paired design. A variant is derived from a *base* instance made by the unchanged default generator
(bench/generate.py, gate-matched exponential kernel, r_plant = 30 um). The variant keeps every real gene, every
ligand/receptor channel, and every PHANTOM / NICHE / RANDOM pair byte-identical and re-plants ONLY the target
response of the TRUE pairs with a different mechanism. So, for a given tissue and seed, A0-A4 differ from each
other only in how receivers respond; the expression-only tool lists are identical across A0-A4.

A5 (mixed) additionally (i) assigns the misspecified mechanisms to the TRUE pairs of one tissue (2 pairs each),
(ii) adds DECOY senders (wrong-type-pair negatives): the ligand of some TRUE pairs is also expressed by a second,
non-interacting cell type, and (iii) re-plants the NICHE pairs with the Phase-1 niche mechanism but stronger
expression so that niche-confounded negatives actually reach top-K tool lists.

Nothing in the database (what tools / agents / the gate see) reveals a mechanism, a category or a type.
The frozen gate is not touched. Parameters below are fixed in reports/laneA_prespec.md before any N2 run.
"""
import copy, hashlib, json, os
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from bench.generate import BenchInstance, label_claim

# Response kernels: exposure E_i of receiver i = sum over senders j of w(d_ij); response of each responder target
# ~ Poisson(mu_T * effect * (1 - exp(-E_i))) exactly as in the default generator.
MECHS = {
    "exp30": dict(kernel="exp", r=30.0, lam=15.0),     # = default generator = gate-matched (reference A0)
    "ring": dict(kernel="ring", r_in=15.0, r_out=30.0),  # A1: annulus, w = 1 for 15 <= d <= 30 um, else 0
    "exp15": dict(kernel="exp", r=15.0, lam=7.5),      # A2a: radius 0.5x
    "exp60": dict(kernel="exp", r=60.0, lam=30.0),     # A2b: radius 2x
    "para80": dict(kernel="exp", r=80.0, lam=40.0),    # A3: diffusive paracrine, long range, weak gradient
    "none": dict(kernel=None),                         # A4: receptor-only response, no target responds
}
MIX = ["ring", "exp15", "exp60", "para80", "none"]      # A5: every misspecified mechanism, 2 TRUE pairs each

# A5 hard negatives (fixed after the dev-seed list-composition pilot; see reports/laneA_prespec.md)
DECOY = dict(n_decoy=5, boost=1.0)                     # 5 of the 10 TRUE ligands also expressed by a decoy type
NICHE_HARD = dict(sigma=600.0, blobs=3, boost=2.0, cap=20.0)   # pilot config c4 (only one with niche in K=10)

VARIANTS = {
    "A0": dict(mech="exp30"),
    "A1": dict(mech="ring"),
    "A2a": dict(mech="exp15"),
    "A2b": dict(mech="exp60"),
    "A3": dict(mech="para80"),
    "A4": dict(mech="none"),
    "A5": dict(mech="mix", decoy=DECOY, niche_hard=NICHE_HARD),
    # Lane B (H5 secondary, exploratory): gate-matched TRUE response + the A5 hard negatives, so that the dose-response of
    # the hard-negative false-pass share is not confounded with mechanism misspecification. Appended last: VCODEs of
    # A0-A5 (and therefore their outputs) are unchanged.
    "A5m": dict(mech="exp30", decoy=DECOY, niche_hard=NICHE_HARD),
}
VCODE = {k: i + 1 for i, k in enumerate(VARIANTS)}


def exposure(xy_send, xy_recv, mech):
    """Per-receiver exposure under a mechanism (no ligand weighting, as in the default generator)."""
    k = mech["kernel"]
    if k is None:
        return np.zeros(len(xy_recv))
    r = mech["r"] if k == "exp" else mech["r_out"]
    D = cKDTree(xy_recv).sparse_distance_matrix(cKDTree(xy_send), r, output_type="coo_matrix")
    if k == "exp":
        w = np.exp(-D.data / mech["lam"])
    elif k == "ring":
        w = (D.data >= mech["r_in"]).astype(np.float64)
    else:
        raise ValueError(k)
    E = np.zeros(len(xy_recv))
    np.add.at(E, D.row, w)
    return E


def _cols(syn_genes):
    return {g: i for i, g in enumerate(syn_genes)}


def _niche_field(xy, centers, sigma):
    d2 = ((xy[:, None, :] - centers[None]) ** 2).sum(-1)
    return np.exp(-d2 / (2 * sigma ** 2)).max(1)


def make_variant(t, base: BenchInstance, variant: str, seed: int, overrides=None) -> BenchInstance:
    """Derive variant `variant` from `base` (same tissue `t`). Deterministic in (base, variant, seed).
    `overrides` (dev pilot only) replaces entries of the variant spec, e.g. {'niche_hard': {...}}."""
    spec = copy.deepcopy(VARIANTS[variant])
    spec.update(overrides or {})
    assert base.tissue == t.name and base.seed == seed
    rng = np.random.default_rng([int(seed), VCODE[variant], 20261108])
    c = base.cfg
    eff, mu_T, mu_L, mu_R = c["effect"], c["mu_T"], c["mu_L"], c["mu_R"]
    gi = _cols(base.syn_genes)
    add = base.added.toarray().astype(np.int32)
    truth = copy.deepcopy(base.truth)
    trues = sorted([p for p in truth if p["category"] == "true"], key=lambda p: p["pair_id"])

    if spec["mech"] == "mix":
        order = rng.permutation(len(trues))
        mech_of = {trues[k]["pair_id"]: MIX[pos % len(MIX)] for pos, k in enumerate(order)}
    else:
        mech_of = {p["pair_id"]: spec["mech"] for p in trues}

    for p in trues:
        m = mech_of[p["pair_id"]]
        iA = np.flatnonzero(t.types == p["sender"]); iB = np.flatnonzero(t.types == p["receiver"])
        for tg in p["targets"]:                    # wipe the base response (all targets of TRUE pairs)
            add[:, gi[tg]] = 0
        E = exposure(t.xy[iA], t.xy[iB], MECHS[m])
        resp = [] if m == "none" else list(p["responders"])   # same responder identity as base (paired)
        for tg in resp:
            add[iB, gi[tg]] = rng.poisson(mu_T * eff * (1 - np.exp(-E)))
        p["mechanism"] = m
        p["mechanism_params"] = MECHS[m]
        p["responders"] = resp
        p["n_receivers_exposed"] = int((E > 0).sum())
        p["exposure_mean_exposed"] = float(E[E > 0].mean()) if (E > 0).any() else 0.0

    types = list(base.meta["types"])
    if spec.get("decoy"):
        dc = spec["decoy"]
        pick = rng.choice(len(trues), size=min(dc["n_decoy"], len(trues)), replace=False)
        for k in sorted(pick.tolist()):
            p = trues[k]
            elig = [x for x in types if x not in (p["sender"], p["receiver"])]
            dA = elig[int(rng.integers(len(elig)))]
            iD = np.flatnonzero(t.types == dA)
            add[iD, gi[p["ligand"]]] += rng.poisson(mu_L * eff * dc["boost"], size=len(iD)).astype(np.int32)
            p["decoy_sender"] = dA
    if spec.get("niche_hard"):
        nh = spec["niche_hard"]
        for p in [q for q in truth if q["category"] == "niche"]:
            iA = np.flatnonzero(t.types == p["sender"]); iB = np.flatnonzero(t.types == p["receiver"])
            for g in [p["ligand"], p["receptor"]] + p["targets"]:
                add[:, gi[g]] = 0
            centers = t.xy[rng.choice(np.r_[iA, iB], size=nh["blobs"], replace=False)]
            rho = _niche_field(t.xy, centers, nh["sigma"])
            sA = min(nh["boost"] / max(rho[iA].mean(), 1e-6), nh["cap"])
            sB = min(nh["boost"] / max(rho[iB].mean(), 1e-6), nh["cap"])
            add[iA, gi[p["ligand"]]] = rng.poisson(mu_L * eff * sA * rho[iA])
            add[iB, gi[p["receptor"]]] = rng.poisson(mu_R * eff * sB * rho[iB])
            for g in p["targets"]:
                add[:, gi[g]] = rng.poisson(mu_T * eff * rho)
            p.update(niche_centers=centers.tolist(), niche_sigma=nh["sigma"], niche_hard=True,
                     niche_scale=[float(sA), float(sB)])

    cfg = dict(base.cfg)
    cfg["laneA_variant"] = variant
    if overrides:
        cfg["laneA_overrides"] = overrides
    meta = dict(base.meta)
    meta.update(laneA_variant=variant, base_cfg_hash=base.meta["cfg_hash"],
                cfg_hash=hashlib.sha1(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:10],
                mechanisms={k: v for k, v in sorted(mech_of.items())})
    return BenchInstance(tissue=base.tissue, seed=seed, cfg=cfg, syn_genes=list(base.syn_genes),
                         baseline=base.baseline, added=sp.csr_matrix(add), truth=truth,
                         database=copy.deepcopy(base.database), meta=meta)


def label_claim_v(inst, claim):
    """label_claim plus the A5 decoy-sender negative sub-type (a TRUE ligand claimed from its decoy type)."""
    lab, sub = label_claim(inst, claim)
    if lab == 0 and sub == "other":
        for p in inst.truth:
            if (p["ligand"], p["receptor"]) == (claim["ligand"], claim["receptor"]) and \
                    p.get("decoy_sender") == claim["sender"] and p["receiver"] == claim["receiver"]:
                return 0, "decoy_sender"
    return lab, sub


def planted_claims(inst):
    """All planted claims at their planted type pairs, plus A5 decoy-sender claims. (claim, category, mechanism)"""
    out = []
    for p in inst.truth:
        out.append((dict(sender=p["sender"], receiver=p["receiver"], ligand=p["ligand"], receptor=p["receptor"]),
                    p["category"] + ("_hard" if p.get("niche_hard") else ""), p.get("mechanism", "")))
        if p.get("decoy_sender"):
            out.append((dict(sender=p["decoy_sender"], receiver=p["receiver"], ligand=p["ligand"],
                             receptor=p["receptor"]), "decoy_sender", p.get("mechanism", "")))
    return out


def save_truth(inst, root):
    """Ground truth of a variant: truth/meta JSON + the planted (added) counts of all synthetic channels."""
    os.makedirs(root, exist_ok=True)
    stem = f"{root}/{inst.tissue}_s{inst.seed}_{inst.meta['laneA_variant']}_{inst.meta['cfg_hash']}"
    sp.save_npz(stem + "_added.npz", inst.added)
    json.dump(dict(tissue=inst.tissue, seed=inst.seed, cfg=inst.cfg, syn_genes=inst.syn_genes, truth=inst.truth,
                   database=inst.database, meta=inst.meta), open(stem + "_truth.json", "w"), indent=1)
    return stem
