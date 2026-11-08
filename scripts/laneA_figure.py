"""Lane A figure: per-variant retention and post-gate phantom rate (K=10, cpdb+prod pooled, 95% hierarchical
bootstrap CI) with the Amendment-4 H4 thresholds. Reads results/tables/laneA_<seed_set>_summary.csv.
Two panels with their own axes (no dual axis). Palette: dataviz reference slots 1 (gate) and 2 (expression-only).

  python -m scripts.laneA_figure <seed_set>
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GATE, EXPR = "#2a78d6", "#eb6834"
ORDER = ["A0", "A1", "A2a", "A2b", "A3", "A4", "A5"]
NAMES = {"A0": "A0\nmatched\n(ref)", "A1": "A1\nring", "A2a": "A2a\nr 15um", "A2b": "A2b\nr 60um",
         "A3": "A3\nparacrine", "A4": "A4\nreceptor-\nonly", "A5": "A5\nmixed +\nhard neg."}


def main(seed_set):
    s = pd.read_csv(f"results/tables/laneA_{seed_set}_summary.csv")
    s = s[(s.K == 10) & (s.tissue == "ALL") & (s.method == "cpdb+prod")].set_index("variant").loc[ORDER]
    x = range(len(ORDER))
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK})
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6), facecolor=SURF)
    panels = [("retention", "TRUE-claim retention", 0.60, "H4: retention >= 0.60", [("retention", GATE, "frozen gate")]),
              ("phantom", "phantom rate (non-TRUE share)", 0.10, "H4: post-gate <= 0.10",
               [("phantom_rate", EXPR, "expression-only (before gate)"), ("phantom_rate_after", GATE, "frozen gate")])]
    for ax, (_, ylab, thr, thr_lab, series) in zip(axes, panels):
        ax.set_facecolor(SURF)
        ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
        ax.axhline(thr, color=INK2, lw=1.0, ls=(0, (4, 3)))
        ax.text(len(ORDER) - 0.45, thr, thr_lab, ha="right", va="bottom", fontsize=7, color=INK2)
        ax.axvspan(-0.5, 0.5, color=GRID, alpha=0.5, lw=0)          # A0 = reference, not part of H4
        off = [-0.12, 0.12] if len(series) == 2 else [0.0]
        for (k, col, lab), o in zip(series, off):
            y = s[k].values; lo = s[f"{k}_lo"].values; hi = s[f"{k}_hi"].values
            xs = [i + o for i in x]
            ax.vlines(xs, lo, hi, color=col, lw=2)
            ax.scatter(xs, y, s=30, color=col, edgecolor=SURF, linewidth=1.5, zorder=3, label=lab)
        ax.set_xticks(list(x)); ax.set_xticklabels([NAMES[v] for v in ORDER], fontsize=6.5)
        ax.set_ylim(-0.03, 1.03); ax.set_ylabel(ylab)
        if len(series) > 1:
            ax.legend(frameon=False, fontsize=7, loc="upper left")
    fig.suptitle(f"H4 misspecification ({seed_set}, 4 tissues x {s['n_claims'].iloc[0] and ''}seeds; K=10, 2 tools pooled; 95% CI)",
                 fontsize=8, color=INK2, y=0.99)
    fig.tight_layout()
    out = f"figures/fig_laneA_h4_{seed_set}.png"
    fig.savefig(out, dpi=300, facecolor=SURF)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1])
