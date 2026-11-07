"""Claim-list metrics + hierarchical (seed -> claim) bootstrap CIs."""
import numpy as np


def list_metrics(labels, passed=None, n_true_planted=None):
    """labels: 1 = planted TRUE claim. passed: gate pass flags (optional)."""
    y = np.asarray(labels, dtype=float)
    out = dict(n=len(y), n_true=int(y.sum()))
    out["phantom_rate"] = float(1 - y.mean()) if len(y) else np.nan
    out["precision"] = float(y.mean()) if len(y) else np.nan
    if n_true_planted:
        out["recall"] = float(y.sum() / n_true_planted)
    if passed is not None:
        p = np.asarray(passed, dtype=bool)
        out["n_pass"] = int(p.sum())
        out["phantom_rate_after"] = float(1 - y[p].mean()) if p.sum() else 0.0
        out["retention"] = float(y[p].sum() / y.sum()) if y.sum() else np.nan
        pr0 = out["phantom_rate"]
        out["rel_reduction"] = float(1 - out["phantom_rate_after"] / pr0) if pr0 > 0 else np.nan
        if n_true_planted:
            out["recall_after"] = float(y[p].sum() / n_true_planted)
    return out


def pooled(groups):
    """groups: list of (labels, passed) per seed -> pooled metrics."""
    y = np.concatenate([np.asarray(g[0], float) for g in groups]) if groups else np.array([])
    p = None
    if groups and groups[0][1] is not None:
        p = np.concatenate([np.asarray(g[1], bool) for g in groups])
    return list_metrics(y, p)


def hier_bootstrap(groups, keys=("phantom_rate",), n_boot=2000, seed=0):
    """Resample seeds with replacement, then claims within each resampled seed."""
    rng = np.random.default_rng(seed)
    S = len(groups)
    vals = {k: [] for k in keys}
    for _ in range(n_boot):
        gs = []
        for s in rng.integers(0, S, size=S):
            y, p = groups[s]
            y = np.asarray(y)
            if len(y) == 0:
                continue
            idx = rng.integers(0, len(y), size=len(y))
            gs.append((y[idx], None if p is None else np.asarray(p)[idx]))
        m = pooled(gs) if gs else {}
        for k in keys:
            v = m.get(k, np.nan)
            vals[k].append(np.nan if v is None else v)
    point = pooled(groups)
    ci = {k: [float(np.nanpercentile(vals[k], 2.5)), float(np.nanpercentile(vals[k], 97.5))] for k in keys}
    return point, ci
