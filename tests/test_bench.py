"""Phase-2 gate: labels are correct and injection raises exactly the intended genes."""
import numpy as np
import scipy.sparse as sp
import pytest
from scipy.spatial import cKDTree
from bench.tissue import Tissue
from bench.generate import generate, label_claim


def toy_tissue(seed=0, n=6000):
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, 2000, size=(n, 2))
    # 4 types: T0/T1 intermixed on the left, T2 on the right, T3 everywhere
    types = np.where(xy[:, 0] < 900, rng.choice(["T0", "T1"], n), "T2").astype(object)
    types[rng.random(n) < 0.25] = "T3"
    types[rng.random(n) < 0.02] = "Unlabeled"
    counts = sp.csr_matrix(rng.poisson(0.5, size=(n, 30)).astype(np.float32))
    lib = np.asarray(counts.sum(1)).ravel() + 1
    return Tissue("toy", xy, types.astype(str), counts, np.array([f"G{i}" for i in range(30)]), lib)


CFG = dict(n_true=3, n_phantom=3, n_niche=3, n_random=3, min_type_cells=100, min_true_contacts=50,
           niche_sigma=300.0, donor_min_mean=0.0)


@pytest.fixture(scope="module")
def inst():
    return generate(toy_tissue(), seed=7, cfg=CFG)


@pytest.fixture(scope="module")
def tis():
    return toy_tissue()


def cols(inst, names):
    gi = inst.gene_index()
    return [gi[g] for g in names]


def test_determinism(tis):
    a = generate(tis, 3, CFG); b = generate(tis, 3, CFG)
    assert (a.added != b.added).nnz == 0 and (a.baseline != b.baseline).nnz == 0
    assert a.truth == b.truth
    c = generate(tis, 4, CFG)
    assert (a.added != c.added).nnz > 0


def test_categories_and_names_anonymized(inst):
    cats = [p["category"] for p in inst.truth]
    assert sorted(cats) == sorted(["true"] * 3 + ["phantom"] * 3 + ["niche"] * 3 + ["random"] * 3)
    for d in inst.database:                                  # db exposes no category / type info
        assert set(d) == {"pair_id", "ligand", "receptor", "targets"}
    assert len(set(inst.syn_genes)) == len(inst.syn_genes) == 12 * 5


def test_random_pairs_untouched(inst):
    for p in inst.truth:
        if p["category"] == "random":
            c = cols(inst, [p["ligand"], p["receptor"]] + p["targets"])
            assert inst.added[:, c].nnz == 0


def test_true_injection_exact(inst, tis):
    for p in [q for q in inst.truth if q["category"] == "true"]:
        iA = tis.types == p["sender"]; iB = tis.types == p["receiver"]
        L, R = cols(inst, [p["ligand"], p["receptor"]])
        aL = inst.added[:, L].toarray().ravel(); aR = inst.added[:, R].toarray().ravel()
        assert aL[~iA].sum() == 0 and aL[iA].sum() > 0       # ligand only in senders
        assert aR[~iB].sum() == 0 and aR[iB].sum() > 0       # receptor only in receivers
        # targets: only responders change, only in receivers, only if a sender within r_plant
        xyA = tis.xy[iA]; tree = cKDTree(xyA)
        dmin, _ = tree.query(tis.xy, k=1)
        near = dmin <= inst.cfg["r_plant"]
        for tg in p["targets"]:
            a = inst.added[:, cols(inst, [tg])[0]].toarray().ravel()
            if tg in p["responders"]:
                assert a[~iB].sum() == 0
                assert a[iB & ~near].sum() == 0                # no exposure -> no effect
                assert a[iB & near].sum() > 0
                # dose-response: closer receivers respond more on average
                close = iB & (dmin <= inst.cfg["r_plant"] / 3); mid = iB & near & ~close
                if close.sum() > 20 and mid.sum() > 20:
                    assert a[close].mean() > a[mid].mean()
            else:
                assert a.sum() == 0


def test_phantom_never_meet(inst, tis):
    for p in [q for q in inst.truth if q["category"] == "phantom"]:
        L, R = cols(inst, [p["ligand"], p["receptor"]])
        sL = np.flatnonzero(inst.added[:, L].toarray().ravel() > 0)
        sR = np.flatnonzero(inst.added[:, R].toarray().ravel() > 0)
        assert np.all(tis.types[sL] == p["sender"]) and np.all(tis.types[sR] == p["receiver"])
        if len(sL) and len(sR):
            d, _ = cKDTree(tis.xy[sR]).query(tis.xy[sL], k=1)
            assert d.min() > p["r_seg"]                       # planted L and R cells never meet
        for tg in p["targets"]:
            assert inst.added[:, cols(inst, [tg])[0]].nnz == 0


def test_niche_regional(inst, tis):
    for p in [q for q in inst.truth if q["category"] == "niche"]:
        c = np.array(p["niche_centers"])
        d = np.sqrt(((tis.xy[:, None] - c[None]) ** 2).sum(-1)).min(1)
        far = d > 5 * p["niche_sigma"]
        L, R = cols(inst, [p["ligand"], p["receptor"]])
        aL = inst.added[:, L].toarray().ravel()
        assert aL[tis.types != p["sender"]].sum() == 0
        for tg in p["targets"]:
            a = inst.added[:, cols(inst, [tg])[0]].toarray().ravel()
            assert a[far].sum() == 0 and a[d < p["niche_sigma"]].mean() > a[d > 2 * p["niche_sigma"]].mean()


def test_real_genes_untouched(inst, tis):
    before = tis.counts.copy()
    generate(tis, 11, CFG)
    assert (before != tis.counts).nnz == 0


def test_labels(inst, tis):
    U = inst.claim_universe(tis.type_list())
    assert sum(u["label"] for u in U) == 3
    for u in U:
        assert label_claim(inst, u) == (u["label"], "true" if u["label"] else u["neg_type"])
    for p in inst.truth:
        lab, sub = label_claim(inst, dict(sender=p["sender"], receiver=p["receiver"],
                                          ligand=p["ligand"], receptor=p["receptor"]))
        assert lab == (p["category"] == "true") and sub == p["category"]
    # swapped direction of a true pair is a negative
    p = next(q for q in inst.truth if q["category"] == "true")
    assert label_claim(inst, dict(sender=p["receiver"], receiver=p["sender"], ligand=p["ligand"],
                                  receptor=p["receptor"]))[0] == 0
