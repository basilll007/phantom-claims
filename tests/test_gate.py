"""Gate unit tests: determinism, BH, null calibration, detection of a planted contact effect, never-meet failure."""
import numpy as np
import scipy.sparse as sp
from bench.tissue import Tissue
from gate.gate import Gate, GateConfig, ExprSource, bh


def toy(seed=0, n=8000, signal=0.0, never_meet=False):
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, 2000, size=(n, 2))
    types = rng.choice(["A", "B", "C"], n)
    if never_meet:                        # A only left, B only right
        types = np.where(xy[:, 0] < 900, rng.choice(["A", "C"], n), rng.choice(["B", "C"], n))
    L = np.where(types == "A", rng.poisson(1.0, n), 0)
    R = np.where(types == "B", rng.poisson(1.0, n), 0)
    T = rng.poisson(0.3, n)
    if signal > 0:
        iA = np.flatnonzero(types == "A"); iB = np.flatnonzero(types == "B")
        from scipy.spatial import cKDTree
        d, _ = cKDTree(xy[iA]).query(xy[iB], k=1)
        T[iB] += rng.poisson(signal * (d <= 30) * np.exp(-d / 15))
    counts = sp.csr_matrix(np.c_[L, R, T, rng.poisson(1.0, (n, 5))].astype(np.float32))
    lib = np.asarray(counts.sum(1)).ravel() + 5
    t = Tissue("toy", xy, types.astype(str), counts, np.array(["L", "R", "T"] + [f"G{i}" for i in range(5)]), lib)
    return t, ExprSource(t, counts, list(t.genes), {("L", "R"): ["T"]})


CL = dict(sender="A", receiver="B", ligand="L", receptor="R")


def test_bh():
    q = bh([0.01, 0.04, 0.03, 0.5])
    assert np.allclose(q, [0.04, 0.0533333, 0.0533333, 0.5], atol=1e-4)


def test_deterministic():
    t, x = toy(signal=1.5)
    a = Gate(t, x, GateConfig(n_perm=200)).check([CL])[0]
    b = Gate(t, x, GateConfig(n_perm=200)).check([CL])[0]
    assert a["p_perm"] == b["p_perm"] and a["rho_within_tile"] == b["rho_within_tile"]


def test_detects_planted_contact_effect():
    t, x = toy(signal=2.0)
    r = Gate(t, x, GateConfig(n_perm=300)).check([CL])[0]
    assert r["pass_"], r["reasons"]


def test_never_meet_fails():
    t, x = toy(never_meet=True)
    r = Gate(t, x, GateConfig(n_perm=200)).check([CL])[0]
    assert not r["pass_"]
    assert any("too_few_expressing_pairs" in s for s in r["reasons"])


def test_null_calibration():
    ps = []
    for s in range(30):
        t, x = toy(seed=100 + s, n=4000, signal=0.0)
        ps.append(Gate(t, x, GateConfig(n_perm=200)).evaluate(CL)["p_perm"])
    ps = np.array(ps)
    assert (ps < 0.05).mean() <= 0.2 and 0.3 < ps.mean() < 0.7
