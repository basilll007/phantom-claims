"""Lane B (H5) aggregation rules: Spearman unit statistic, tie/constant handling, float-safe thresholds."""
import pandas as pd
from scripts.laneB_aggregate import unit_rhos, LEVELS


def frame(prs, tissue="t", seed=1):
    rows = []
    for e, pr in zip(LEVELS, prs):
        n_neg = int(round(pr * 20))
        for i in range(20):
            rows.append(dict(list_type="default", K=10, method="cpdb" if i < 10 else "prod", tissue=tissue, seed=seed,
                             effect=e, label=int(i >= n_neg)))
    return pd.DataFrame(rows)


def test_monotone_increase_is_minus_one():
    r = unit_rhos(frame([0.1, 0.3, 0.6, 0.9]), "default")
    assert r.rho.iloc[0] == -1.0


def test_one_swap_is_exactly_minus_080():
    r = unit_rhos(frame([0.3, 0.1, 0.6, 0.9]), "default")
    assert r.rho.iloc[0] == -0.8 and r.rho.iloc[0] <= -0.80


def test_constant_is_zero():
    r = unit_rhos(frame([0.2, 0.2, 0.2, 0.2]), "default")
    assert r.rho.iloc[0] == 0.0


def test_ties_use_average_ranks():
    r = unit_rhos(frame([0.1, 0.1, 0.5, 0.9]), "default")
    assert -1.0 < r.rho.iloc[0] < -0.9
