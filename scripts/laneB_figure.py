"""Lane B figure: figures/fig_dose_response.png. Phantom rate vs effect size (expression-only before gate vs frozen gate)
for the default lists (H5) and the A5m hard-negative lists (exploratory), plus retention, K=10, cpdb+prod pooled,
95% hierarchical bootstrap CI. Three panels with their own axes (no dual axis).
Palette: dataviz reference slots 1 (gate) / 2 (expression-only); retention uses line style to separate list types.

  python -m scripts.laneB_figure <seed_set>
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
GATE, EXPR = "#2a78d6", "#eb6834"
TABLES = os.environ.get("PC_TABLES", "results/tables")
FIGS = os.environ.get("PC_FIGS", "figures")
LEVELS = [1.0, 0.5, 0.25, 0.125]


def main(seed_set):
    s = pd.read_csv(f"{TABLES}/laneB_{seed_set}_summary.csv")
    s = s[(s.K == 10) & (s.tissue == "ALL")]
    plt.rcParams.update({"font.size": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK})
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8), facecolor=SURF)

    def style(ax, title, ylab):
        ax.set_facecolor(SURF); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
        ax.set_xscale("log", base=2); ax.set_xticks(LEVELS); ax.set_xticklabels(["1", "1/2", "1/4", "1/8"])
        ax.invert_xaxis(); ax.set_xlabel("planted effect size (x default)"); ax.set_ylabel(ylab)
        ax.set_ylim(-0.03, 1.03); ax.set_title(title, fontsize=8, color=INK2)

    def series(ax, d, k, col, lab, ls="-"):
        d = d.set_index("effect").loc[LEVELS]
        ax.plot(LEVELS, d[k], color=col, lw=2, ls=ls, label=lab, zorder=2)
        ax.vlines(LEVELS, d[f"{k}_lo"], d[f"{k}_hi"], color=col, lw=2, zorder=2)
        ax.scatter(LEVELS, d[k], s=30, color=col, edgecolor=SURF, linewidth=1.5, zorder=3)

    for ax, lt, title in [(axes[0], "default", "default lists (H5)"), (axes[1], "A5m", "A5m hard-negative lists (exploratory)")]:
        style(ax, title, "phantom rate (non-TRUE share)")
        d = s[s.list_type == lt]
        series(ax, d, "phantom_rate", EXPR, "expression-only")
        series(ax, d, "phantom_rate_after", GATE, "frozen gate")
        if lt == "default":
            ax.axhline(0.15, color=INK2, lw=1.0, ls=(0, (4, 3)), label="H5 threshold: gated <= 0.15")
    style(axes[2], "gate retention of TRUE claims", "retention")
    series(axes[2], s[s.list_type == "default"], "retention", GATE, "default lists")
    series(axes[2], s[s.list_type == "A5m"], "retention", GATE, "A5m lists", ls=(0, (2, 2)))
    axes[2].legend(frameon=False, fontsize=6.5, loc="lower left")
    fig.suptitle(f"H5 dose-response ({seed_set} seeds x 4 tissues, K=10, 2 tools pooled, 95% CI)", fontsize=7, color=INK2)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=7)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    out = f"{FIGS}/fig_dose_response.png" if seed_set == "N2" else f"{FIGS}/fig_dose_response_{seed_set}.png"
    fig.savefig(out, dpi=300, facecolor=SURF)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1])
