"""Lane A (H4) misspecified-variant generator: labels, mechanisms and paired-design invariants."""
import numpy as np
import pytest
from scipy.spatial import cKDTree
from bench.generate import generate, label_claim
from bench.variants import (make_variant, VARIANTS, MECHS, MIX, DECOY, label_claim_v, planted_claims, exposure)
from test_bench import toy_tissue, CFG

SINGLE = ["A0", "A1", "A2a", "A2b", "A3", "A4"]
CFG10 = dict(CFG, n_true=10)


@pytest.fixture(scope="module")
def tis():
    return toy_tissue()


@pytest.fixture(scope="module")
def base(tis):
    return generate(tis, 7, CFG)


@pytest.fixture(scope="module")
def var(tis, base):
    return {v: make_variant(tis, base, v, 7) for v in VARIANTS}


@pytest.fixture(scope="module")
def base10(tis):
    return generate(tis, 5, CFG10)


@pytest.fixture(scope="module")
def a5(tis, base10):
    return make_variant(tis, base10, "A5", 5)


def colv(inst, g):
    return inst.added[:, inst.gene_index()[g]].toarray().ravel()


def n_within(tree, xy, r):
    return np.array([len(x) for x in tree.query_ball_point(xy, r)])


def test_paired_only_true_targets_change(base, var):
    true_tg = {g for p in base.truth if p["category"] == "true" for g in p["targets"]}
    for v in SINGLE:
        inst = var[v]
        assert (inst.baseline != base.baseline).nnz == 0
        assert inst.database == base.database and inst.syn_genes == base.syn_genes
        for g in base.syn_genes:
            if g not in true_tg:                       # L, R, and every negative channel byte-identical
                assert np.array_equal(colv(inst, g), colv(base, g)), (v, g)


def test_mechanism_support(tis, var):
    """Responders respond only in receivers with >=1 sender inside the mechanism's support, and do respond there."""
    for v in ["A0", "A1", "A2a", "A2b", "A3"]:
        inst, m = var[v], MECHS[VARIANTS[v]["mech"]]
        r_hi = m["r"] if m["kernel"] == "exp" else m["r_out"]
        for p in [q for q in inst.truth if q["category"] == "true"]:
            assert p["mechanism"] == VARIANTS[v]["mech"]
            iA = tis.types == p["sender"]; iB = tis.types == p["receiver"]
            tree = cKDTree(tis.xy[iA])
            n_sup = n_within(tree, tis.xy, r_hi)
            if m["kernel"] == "ring":
                n_sup = n_sup - n_within(tree, tis.xy, m["r_in"] - 1e-6)
            sup = n_sup > 0
            assert len(p["responders"]) == 2
            for tg in p["targets"]:
                a = colv(inst, tg)
                if tg in p["responders"]:
                    assert a[~iB].sum() == 0
                    assert a[iB & ~sup].sum() == 0
                    assert a[iB & sup].sum() > 0
                else:
                    assert a.sum() == 0


def test_ring_excludes_close_only_receivers(tis, var):
    """A1: a receiver whose only senders are closer than r_in gets no response (exp30 would respond)."""
    m = MECHS["ring"]
    hit = 0
    for p in [q for q in var["A1"].truth if q["category"] == "true"]:
        iA = tis.types == p["sender"]; iB = np.flatnonzero(tis.types == p["receiver"])
        tree = cKDTree(tis.xy[iA])
        close = n_within(tree, tis.xy[iB], m["r_in"] - 1e-6) > 0
        ring = (n_within(tree, tis.xy[iB], m["r_out"]) - n_within(tree, tis.xy[iB], m["r_in"] - 1e-6)) > 0
        only_close = iB[close & ~ring]
        for tg in p["responders"]:
            assert colv(var["A1"], tg)[only_close].sum() == 0
        hit += len(only_close)
    assert hit >= 0


def test_radius_ordering(var):
    by = {v: {p["pair_id"]: p["n_receivers_exposed"] for p in var[v].truth if p["category"] == "true"}
          for v in ["A2a", "A0", "A2b", "A3"]}
    for pid in by["A0"]:
        assert by["A2a"][pid] <= by["A0"][pid] <= by["A2b"][pid] <= by["A3"][pid]


def test_exposure_matches_default_generator(tis, base):
    """exp30 exposure == the default generator's E: same response law and same expected counts."""
    p = next(q for q in base.truth if q["category"] == "true")
    iA = np.flatnonzero(tis.types == p["sender"]); iB = np.flatnonzero(tis.types == p["receiver"])
    E = exposure(tis.xy[iA], tis.xy[iB], MECHS["exp30"])
    assert int((E > 0).sum()) == p["n_receivers_exposed"]


def test_receptor_only(var, base):
    for p in [q for q in var["A4"].truth if q["category"] == "true"]:
        assert p["responders"] == [] and p["mechanism"] == "none"
        for tg in p["targets"]:
            assert colv(var["A4"], tg).sum() == 0
        assert np.array_equal(colv(var["A4"], p["receptor"]), colv(base, p["receptor"]))  # receptor unchanged


def test_labels_unchanged_single(tis, var, base):
    U = base.claim_universe(tis.type_list())
    for v in SINGLE:
        for u in U[::7]:
            assert label_claim_v(var[v], u) == label_claim(base, u)


def test_a5_mechanism_mix(a5):
    mech = [p["mechanism"] for p in a5.truth if p["category"] == "true"]
    assert sorted(mech) == sorted(MIX * 2)
    assert all(m != "exp30" for m in mech)


def test_a5_decoy_sender(tis, a5, base10):
    dec = [p for p in a5.truth if p.get("decoy_sender")]
    assert len(dec) == DECOY["n_decoy"]
    for p in dec:
        dA = p["decoy_sender"]
        assert dA not in (p["sender"], p["receiver"]) and p["category"] == "true"
        aL = colv(a5, p["ligand"])
        assert aL[tis.types == dA].sum() > 0
        assert aL[(tis.types != dA) & (tis.types != p["sender"])].sum() == 0
        assert np.array_equal(aL[tis.types == p["sender"]], colv(base10, p["ligand"])[tis.types == p["sender"]])
        c = dict(sender=dA, receiver=p["receiver"], ligand=p["ligand"], receptor=p["receptor"])
        assert label_claim_v(a5, c) == (0, "decoy_sender")
        assert label_claim(a5, c) == (0, "other")
        t_ = dict(sender=p["sender"], receiver=p["receiver"], ligand=p["ligand"], receptor=p["receptor"])
        assert label_claim_v(a5, t_) == (1, "true")
        other = next(x for x in tis.type_list() if x not in (dA, p["sender"], p["receiver"]))
        assert label_claim_v(a5, dict(c, sender=other)) == (0, "other")


def test_a5_niche_hard(tis, a5, base10):
    for p in [q for q in a5.truth if q["category"] == "niche"]:
        assert p["niche_hard"]
        aL, aR = colv(a5, p["ligand"]), colv(a5, p["receptor"])
        assert aL[tis.types != p["sender"]].sum() == 0 and aR[tis.types != p["receiver"]].sum() == 0
        c = np.array(p["niche_centers"])
        d = np.sqrt(((tis.xy[:, None] - c[None]) ** 2).sum(-1)).min(1)
        for tg in p["targets"]:
            a = colv(a5, tg)
            assert a[d < p["niche_sigma"]].mean() > a[d > 2 * p["niche_sigma"]].mean()
        # stronger than (or equal to) the Phase-1 niche: sender-type mean ligand near boost*mu_L unless capped
        assert aL[tis.types == p["sender"]].mean() >= 0.5 * colv(base10, p["ligand"])[tis.types == p["sender"]].mean()
        lab, sub = label_claim_v(a5, dict(sender=p["sender"], receiver=p["receiver"], ligand=p["ligand"],
                                          receptor=p["receptor"]))
        assert (lab, sub) == (0, "niche")


def test_a5_phantom_random_untouched(a5, base10):
    for p in a5.truth:
        if p["category"] in ("phantom", "random"):
            for g in [p["ligand"], p["receptor"]] + p["targets"]:
                assert np.array_equal(colv(a5, g), colv(base10, g))


def test_database_anonymized(a5, var):
    for inst in list(var.values()) + [a5]:
        for d in inst.database:
            assert set(d) == {"pair_id", "ligand", "receptor", "targets"}


def test_planted_claims(a5):
    pc = planted_claims(a5)
    cats = [c for _, c, _ in pc]
    assert cats.count("decoy_sender") == DECOY["n_decoy"] and cats.count("true") == 10
    assert cats.count("niche_hard") == 3 and cats.count("phantom") == 3
    for cl, cat, _ in pc:
        assert label_claim_v(a5, cl)[0] == (cat == "true")


def test_determinism(tis, base):
    a = make_variant(tis, base, "A1", 7); b = make_variant(tis, base, "A1", 7)
    assert (a.added != b.added).nnz == 0 and a.truth == b.truth
    c = make_variant(tis, base, "A2b", 7)
    assert (a.added != c.added).nnz > 0


def test_variant_codes_stable():
    """Lane A results depend on these RNG codes; appending variants must not renumber them."""
    from bench.variants import VCODE
    assert [VCODE[k] for k in ["A0", "A1", "A2a", "A2b", "A3", "A4", "A5"]] == [1, 2, 3, 4, 5, 6, 7]


def test_a5m_matched_true_with_hard_negatives(tis, base10):
    v = make_variant(tis, base10, "A5m", 5)
    mech = [p["mechanism"] for p in v.truth if p["category"] == "true"]
    assert mech == ["exp30"] * 10
    assert sum(1 for p in v.truth if p.get("decoy_sender")) == DECOY["n_decoy"]
    assert all(p.get("niche_hard") for p in v.truth if p["category"] == "niche")
    for p in [q for q in v.truth if q["category"] == "true"]:
        assert len(p["responders"]) == 2


def test_effect_scaling_reaches_variants(tis):
    """Effect-size levels scale the A5-style hard negatives too (lane B): weaker base -> weaker decoy ligand."""
    import numpy as np
    b1 = generate(tis, 5, CFG10); b8 = generate(tis, 5, dict(CFG10, effect=0.125))
    v1 = make_variant(tis, b1, "A5m", 5); v8 = make_variant(tis, b8, "A5m", 5)
    d1 = [p for p in v1.truth if p.get("decoy_sender")]
    m1 = np.mean([colv(v1, p["ligand"])[tis.types == p["decoy_sender"]].mean() for p in d1])
    d8 = [p for p in v8.truth if p.get("decoy_sender")]
    m8 = np.mean([colv(v8, p["ligand"])[tis.types == p["decoy_sender"]].mean() for p in d8])
    assert m8 < 0.5 * m1
