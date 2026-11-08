"""All paper figures, in one place: paper/figures/*.png (600 dpi bitmaps, drawn at their printed size).

    PC_DATA=<data root> python -m scripts.paper_figures            # every figure
    PC_DATA=<data root> python -m scripts.paper_figures fig2 fig3  # a subset

AAAI-27 author-kit rules applied here (AnonymousSubmission2027.tex, "Illustrations and Figures"):
  * figures are drawn at their printed size (figure* = 7.0 in = \\textwidth; figure = 3.3 in = \\columnwidth) and
    included at that width, so every label is >= 9 pt in print;
  * 600 dpi PNG (kit: .png/.jpg/.pdf, >= 300 dpi). Not PDF: matplotlib embeds vector text as CID/Identity-H
    fonts, which the kit says must be outlined or removed even inside figures; a bitmap carries no fonts;
  * no cropping in LaTeX; captions are written in main.tex, under the figure.
Numbers come from results/results.jsonl records (latest per id) or results/tables/*.csv derived from them. The
overview figure illustrates the frozen gate on dev seed 0 (breast Rep1); dev seeds are illustration-only and no
confirmatory (N2) data is used or re-run there.
Palette: validated categorical slots 1-3 (blue = frozen gate, orange = expression-only tool, aqua = third series),
checked with the dataviz validator against a white surface (all pairs pass; aqua always direct-labelled).
"""
import ast, dataclasses, json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

OUT = "paper/figures"
W2, W1 = 7.0, 3.3                    # \textwidth, \columnwidth (inches)
GATE, EXPR, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, BASE, CELL = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#d9d8d2"
BLUES = LinearSegmentedColormap.from_list("blues", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6",
                                                    "#1c5cab", "#104281", "#0d366b"])
FS = 9
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "font.size": FS,
    "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
    "axes.edgecolor": BASE, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.major.size": 3,
    "ytick.major.size": 3, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "axes.titlelocation": "left",
    "axes.titlepad": 4, "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white"})
LEVELS = [1.0, 0.5, 0.25, 0.125]
LEVLAB = ["1", "1/2", "1/4", "1/8"]


def records():
    R = {}
    for line in open("results/results.jsonl"):
        if line.strip():
            r = json.loads(line); R[r["id"]] = r
    return R


def tidy(ax, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.6); ax.set_axisbelow(True)


def title(ax, letter, text):
    ax.set_title(f"$\\bf{{{letter}}}$  {text}", loc="left", fontsize=FS, color=INK)


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    fig.savefig(f"{OUT}/{name}.png", dpi=600)
    plt.close(fig)
    print(f"{OUT}/{name}.png")


def ci_err(p, lo, hi):
    return np.array([[p - lo], [hi - p]])


# ------------------------------------------------------------------------------------------------ Figure 1 (overview)
def null_samples(E, Y, groups, n_perm=1000, seed=0):
    """Within-group permutation null of the gate's statistic (centred correlation), as in gate.Gate._perm_test."""
    _, g = np.unique(groups, return_inverse=True)
    def centre(v):
        s = np.bincount(g, weights=v); c = np.bincount(g)
        return v - (s / c)[g]
    Ec, Yc = centre(E.astype(float)), centre(Y.astype(float))
    den = np.sqrt((Ec ** 2).sum() * (Yc ** 2).sum())
    obs = float((Ec * Yc).sum() / den)
    rng = np.random.default_rng(seed)
    order = np.argsort(g, kind="stable"); gs = g[order]
    null = np.empty(n_perm)
    for k in range(n_perm):
        key = gs * 2.0 + rng.random(len(gs))
        src = order[np.argsort(key, kind="stable")]
        Yp = np.empty_like(Yc); Yp[order] = Yc[src]
        null[k] = (Yp @ Ec) / den
    return obs, null


def fig1():
    from bench.tissue import load_tissue
    from bench.generate import get_or_make
    from bench.pipeline import make_gate
    from gate.gate import GateConfig
    t = load_tissue("breast_rep1"); inst = get_or_make(t, 0)
    cfg = GateConfig(); gate = make_gate(t, inst, cfg)
    P = {p["pair_id"]: p for p in inst.truth}
    picks = [("true", P["P34"]), ("phantom", P["P05"]), ("niche", P["P09"])]
    claim = lambda p: dict(sender=p["sender"], receiver=p["receiver"], ligand=p["ligand"], receptor=p["receptor"])
    recs = gate.decide([gate.evaluate(claim(p)) for _, p in picks])        # gated as one list (BH over 3)
    for (lab, p), r in zip(picks, recs):
        print(lab, p["pair_id"], r["pass_"], r.get("rho_within_tile"), r.get("p_perm"), r.get("q_bh"), r["reasons"])
    gi = inst.gene_index(); C = inst.counts.tocsc(); Add = inst.added.tocsc()
    col = lambda g: np.asarray(C[:, gi[g]].todense()).ravel()
    fig = plt.figure(figsize=(W2, 3.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.3, 1.25], height_ratios=[2.2, 1.0], wspace=0.22, hspace=0.78,
                          left=0.01, right=0.99, bottom=0.02, top=0.91)
    TY = 0.958                                                                  # common title baseline

    def ftitle(ax, letter, text):
        fig.canvas.draw(); x = ax.get_position().x0
        fig.text(x, TY, f"$\\bf{{{letter}}}$  {text}", fontsize=FS, va="baseline")

    # (a) TRUE claim, 600 x 600 um window around the densest exposed-receiver neighbourhood
    tp = P["P34"]
    iA = np.flatnonzero(t.types == tp["sender"]); iB = np.flatnonzero(t.types == tp["receiver"])
    Wt = gate._weights(tp["sender"], tp["receiver"])[2]
    E = Wt @ gate.x.log_multi(tp["ligand"])[iA]
    Y = np.mean([gate.x.z(g)[iB] for g in tp["targets"]], 0)
    H, xe, ye = np.histogram2d(t.xy[iB, 0], t.xy[iB, 1], bins=[np.arange(0, 7600, 300), np.arange(0, 5600, 300)],
                               weights=(E > 0).astype(float))
    bx, by = np.unravel_index(H.argmax(), H.shape)
    S = 600.0; x0, y0 = xe[bx] + 150 - S / 2, ye[by] + 150 - S / 2
    inwin = lambda idx: idx[(t.xy[idx, 0] > x0) & (t.xy[idx, 0] < x0 + S) & (t.xy[idx, 1] > y0) & (t.xy[idx, 1] < y0 + S)]
    ax = fig.add_subplot(gs[0, 0])
    allw = inwin(np.arange(t.n))
    ax.scatter(t.xy[allw, 0], t.xy[allw, 1], s=1.2, c=CELL, lw=0, rasterized=True)
    a = inwin(iA); a = a[col(tp["ligand"])[a] > 0]
    bmask = np.isin(iB, inwin(iB)); b = iB[bmask]
    yv = np.clip(Y[bmask], -1, 3); o = np.argsort(yv)
    ax.scatter(t.xy[b[o], 0], t.xy[b[o], 1], s=5, c=yv[o], cmap=BLUES, vmin=-1, vmax=3, lw=0, rasterized=True)
    ax.scatter(t.xy[a, 0], t.xy[a, 1], s=5, c=EXPR, lw=0, rasterized=True)
    ax.add_patch(Rectangle((x0 + S - 185, y0 + S - 105), 185, 105, facecolor="white", edgecolor="none", alpha=0.9,
                           zorder=4))
    ax.plot([x0 + S - 135, x0 + S - 35], [y0 + S - 25, y0 + S - 25], color=INK, lw=1.5, solid_capstyle="butt", zorder=5)
    ax.text(x0 + S - 85, y0 + S - 40, "100 µm", ha="center", va="bottom", fontsize=FS, zorder=5)
    ax.set_xlim(x0, x0 + S); ax.set_ylim(y0 + S, y0); ax.set_aspect("equal"); ax.axis("off")
    ax.legend([Line2D([], [], ls="", marker="o", ms=4, color=EXPR), Line2D([], [], ls="", marker="o", ms=4, color="#104281")],
              ["ligand$^+$ sender", "receiver (dark = responds)"], loc="upper left",
              bbox_to_anchor=(-0.03, -0.02), handletextpad=0.1, borderaxespad=0, labelspacing=0.3)
    ftitle(ax, "a", "A true claim")
    # (b) PHANTOM claim, whole section
    ph = P["P05"]
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(t.xy[::6, 0], t.xy[::6, 1], s=0.15, c=CELL, lw=0, rasterized=True)
    pa = np.flatnonzero((t.types == ph["sender"]) & (Add[:, gi[ph["ligand"]]].toarray().ravel() > 0))
    pb = np.flatnonzero((t.types == ph["receiver"]) & (Add[:, gi[ph["receptor"]]].toarray().ravel() > 0))
    ax.scatter(t.xy[pa, 0], t.xy[pa, 1], s=2.5, c=EXPR, lw=0, rasterized=True)
    ax.scatter(t.xy[pb, 0], t.xy[pb, 1], s=2.5, c=GATE, lw=0, rasterized=True)
    xm, ym = t.xy[:, 0].max(), t.xy[:, 1].max()
    ax.plot([xm - 1100, xm - 100], [ym + 250, ym + 250], color=INK, lw=1.5, solid_capstyle="butt", clip_on=False)
    ax.text(xm - 600, ym + 420, "1 mm", ha="center", va="top", fontsize=FS)
    ax.set_xlim(0, xm); ax.set_ylim(ym + 100, 0); ax.set_aspect("equal"); ax.axis("off")
    ax.legend([Line2D([], [], ls="", marker="o", ms=4, color=EXPR), Line2D([], [], ls="", marker="o", ms=4, color=GATE)],
              ["ligand-planted sender", "receptor-planted receiver"], loc="upper left", bbox_to_anchor=(-0.02, -0.02),
              handletextpad=0.1, borderaxespad=0, labelspacing=0.3)
    ftitle(ax, "b", "A phantom: the cells never meet")
    # (c) why the null is tile-restricted: NICHE claim under a global vs a within-tile permutation null
    nh = P["P09"]
    iA = np.flatnonzero(t.types == nh["sender"]); iB = np.flatnonzero(t.types == nh["receiver"])
    Wn = gate._weights(nh["sender"], nh["receiver"])[2]
    En = Wn @ gate.x.log_multi(nh["ligand"])[iA]; Yn = np.mean([gate.x.z(g)[iB] for g in nh["targets"]], 0)
    og, ng = null_samples(En, Yn, np.zeros(len(iB), int), seed=1)
    ot, nt = null_samples(En, Yn, gate.tile_id[iB], seed=2)
    pg = (1 + (ng >= og - 1e-12).sum()) / (1 + len(ng)); pt = (1 + (nt >= ot - 1e-12).sum()) / (1 + len(nt))
    print("NICHE P09 global obs %.3f p %.3f | tile obs %.3f p %.3f" % (og, pg, ot, pt))
    ax = fig.add_subplot(gs[0, 2]); tidy(ax, grid=None)
    bins = np.linspace(-0.06, 0.09, 46)
    for yoff, null, obs, p, colr in ((1.3, ng, og, pg, MUTED), (0.0, nt, ot, pt, GATE)):
        h, _ = np.histogram(null, bins=bins); h = h / h.max()
        ax.fill_between(bins[:-1], yoff, yoff + 0.9 * h, step="post", color=colr, alpha=0.3, lw=0)
        ax.step(bins[:-1], yoff + 0.9 * h, where="post", color=colr, lw=1.0)
        ax.plot([obs, obs], [yoff, yoff + 1.0], color=INK, lw=1.5)
        right = obs > 0.04
        ax.text(obs + (-0.003 if right else 0.003), yoff + 1.0, f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}",
                fontsize=FS, ha="right" if right else "left", va="top")
    ax.text(-0.066, 1.3 + 0.75, "global\nnull", fontsize=FS, color=INK2, va="center", ha="left", linespacing=1.0)
    ax.text(0.1, 0.45, "within-\ntile\nnull", fontsize=FS, color=GATE, va="center", ha="right", linespacing=1.0)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.066, 0.1); ax.set_ylim(0, 2.4); ax.set_xticks([-0.05, 0, 0.05])
    ax.set_xticklabels(["−0.05", "0", "0.05"])
    ax.set_xlabel("exposure–response ρ (black: observed)", labelpad=2)
    ftitle(ax, "c", "A niche claim under two nulls")
    # (d) the audit record: the gate's checks for the three claims (frozen gate, gated as one list)
    ax = fig.add_subplot(gs[1, :]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.canvas.draw()
    fig.text(ax.get_position().x0, ax.get_position().y1 + 0.03, "$\\bf{d}$  What the gate writes to the audit record "
             "for these three claims (frozen gate; BH over the three)", fontsize=FS, va="baseline")
    hdr = ["claim (sender → receiver)", "L$^+$R$^+$ contacts ≤ 30 µm\n(need ≥ 50)",
           "within-tile ρ\n(need ≥ 0.02)", "tile-null q\n(need ≤ 0.05)", "both halves\n(p < 0.10)", "verdict"]
    xs = [0.005, 0.40, 0.56, 0.69, 0.82, 0.95]
    for x, h in zip(xs, hdr):
        ax.text(x, 0.98, h, ha="left" if x < 0.1 else "center", va="top", fontsize=FS, color=INK2, linespacing=1.0)
    ax.plot([0, 1], [0.6, 0.6], color=BASE, lw=0.8)
    names = {"true": "true", "phantom": "phantom", "niche": "niche"}
    tick = lambda ok: ("✓", GATE) if ok else ("✗", INK)
    for i, ((lab, p), r) in enumerate(zip(picks, recs)):
        yy = 0.44 - i * 0.2
        halves = r.get("halves", [])
        hp = all(h.get("p") is not None and h["p"] < cfg.alpha_half and (h.get("rho") or 0) > 0 for h in halves)
        q = r.get("q_bh") if r.get("q_bh") is not None else 1.0
        rho = r.get("rho_within_tile") or 0.0
        n = r.get("n_expressing_pairs_within_r", 0)
        cells = [(f"{n:,}", n >= cfg.n_min), (f"{rho:.2f}", rho >= cfg.rho_min),
                 (f"{q:.3f}" if q < 0.01 else f"{q:.2f}", q <= cfg.q), ("", hp)]
        ax.text(xs[0], yy, f"{names[lab]}:  {p['sender'].replace('_', ' ')} → {p['receiver'].replace('_', ' ')}",
                ha="left", va="center", fontsize=FS)
        for x, (v, ok) in zip(xs[1:5], cells):
            g, colr = tick(ok)
            if v:
                ax.text(x + 0.005, yy, v, ha="right", va="center", fontsize=FS)
            ax.text(x + (0.015 if v else 0.0), yy, g, ha="left" if v else "center", va="center", fontsize=FS + 1,
                    color=colr, fontfamily="DejaVu Sans", fontweight="bold")
        ax.text(xs[5], yy, "pass" if r["pass_"] else "fail", ha="center", va="center", fontsize=FS, fontweight="bold",
                color=GATE if r["pass_"] else INK)
    save(fig, "fig1_overview")


# --------------------------------------------------------------------------------------- Figure 2 (H4, misspecified)
VAR = [("A0", "A0  matched (ref.)", dict(kernel="exp", r=30, lam=15)),
       ("A1", "A1  ring 15–30 µm", dict(kernel="ring", r_in=15, r_out=30)),
       ("A2a", "A2a  range 15 µm", dict(kernel="exp", r=15, lam=7.5)),
       ("A2b", "A2b  range 60 µm", dict(kernel="exp", r=60, lam=30)),
       ("A3", "A3  paracrine 80 µm", dict(kernel="exp", r=80, lam=40)),
       ("A4", "A4  no target readout", dict(kernel=None)),
       ("A5", "A5  mixed + hard neg.", dict(kernel="mix"))]


def kern(d, m):
    if m["kernel"] == "exp":
        return np.where(d <= m["r"], np.exp(-d / m["lam"]), 0.0)
    if m["kernel"] == "ring":
        return ((d >= m["r_in"]) & (d <= m["r_out"])).astype(float)
    return np.zeros_like(d)


def fig2(R):
    hv = R["laneA_N2_H4_verdict"]["metrics"]; th = hv["thresholds"]
    rows = [(v, lab, m, R[f"laneA_N2_H4_{v}_K10_pooled"]) for v, lab, m in VAR]
    fig = plt.figure(figsize=(W2, 2.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1.0, 1.0], wspace=0.10, left=0.205, right=0.99, bottom=0.16,
                          top=0.90)
    n = len(rows); ys = np.arange(n)[::-1].astype(float)
    # (a) the planted response kernel of each variant vs the kernel the gate assumes
    ax = fig.add_subplot(gs[0, 0]); tidy(ax, grid=None)
    d = np.linspace(0, 100, 401)
    for y, (v, lab, m, r) in zip(ys, rows):
        ax.axvspan(0, 30, ymin=0, ymax=1, color=GATE, alpha=0.06, lw=0) if v == "A0" else None
        if m["kernel"] in ("exp", "ring"):
            k = kern(d, m)
            ax.fill_between(d, y - 0.38, y - 0.38 + 0.76 * k, color=INK2, alpha=0.25, lw=0)
            ax.plot(d, y - 0.38 + 0.76 * k, color=INK2, lw=1.0)
        elif m["kernel"] is None:
            ax.text(2, y, "receptor only; no target responds", fontsize=FS, va="center", color=INK2)
        else:
            ax.text(2, y, "2 TRUE pairs each of A1–A4", fontsize=FS, va="center", color=INK2)
        gk = kern(d, dict(kernel="exp", r=30, lam=15))
        if m["kernel"] in ("exp", "ring"):
            ax.plot(d, y - 0.38 + 0.76 * gk, color=GATE, lw=1.0, ls=(0, (2, 1.5)))
    ax.set_yticks(ys); ax.set_yticklabels([lab for _, lab, _, _ in rows])
    ax.set_ylim(-0.6, n + 0.25); ax.set_xlim(0, 100); ax.set_xticks([0, 30, 60, 90])
    ax.set_xlabel("sender–receiver distance (µm)")
    ax.tick_params(axis="y", length=0)
    title(ax, "a", "Planted response kernel")
    ax.legend([Line2D([], [], color=INK2, lw=1.0), Line2D([], [], color=GATE, lw=1.0, ls=(0, (2, 1.5)))],
              ["planted", "gate assumes"], loc="upper right", bbox_to_anchor=(1.0, 1.03), handlelength=1.6,
              borderaxespad=0.1, ncol=2, columnspacing=0.8)
    # (b) retention, (c) post-gate phantom rate
    for k, (key, thr, letter, ttl, xl) in enumerate((("retention", th["ret_min"], "b", "Retention of true claims", (0, 1.05)),
                                                     ("phantom_rate_after", th["pr_max"], "c", "Phantom rate after gate", (-0.01, 0.52)))):
        ax = fig.add_subplot(gs[0, 1 + k]); tidy(ax, grid="x")
        good_side = (lambda x: x >= thr) if key == "retention" else (lambda x: x <= thr)
        ax.axvspan(thr if key == "retention" else xl[0], xl[1] if key == "retention" else thr, color=GATE,
                   alpha=0.07, lw=0)
        ax.axvline(thr, color=INK2, lw=1.0, ls=(0, (4, 2)))
        for y, (v, lab, m, r) in zip(ys, rows):
            p = r["metrics"][key]; lo, hi = r["ci"][key]
            met = (r["metrics"]["retention"] >= th["ret_min"]) and (r["metrics"]["phantom_rate_after"] <= th["pr_max"])
            ok = good_side(p)
            ax.plot([lo, hi], [y, y], color=GATE, lw=1.5, solid_capstyle="butt")
            ax.scatter([p], [y], s=36, color=GATE if ok else "white", edgecolor=GATE, lw=1.5, zorder=3)
            txt = f"{p:.3f}" if key == "phantom_rate_after" or p < 0.1 else f"{p:.2f}"
            if not ok:
                ax.text(hi + 0.02 * (xl[1] - xl[0]), y, f"{txt}  miss", va="center", fontsize=FS, color=INK)
        ax.set_yticks(ys); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
        ax.set_ylim(-0.6, n + 0.25); ax.set_xlim(*xl)
        ax.text(thr + (-0.01 if key == "retention" else 0.006), n - 0.25,
                f"H4: {'≥' if key == 'retention' else '≤'} {thr:.2f}", fontsize=FS, color=INK2,
                ha="right" if key == "retention" else "left", va="center")
        ax.set_xlabel("retention (95% CI)" if key == "retention" else "phantom rate (95% CI)")
        title(ax, letter, ttl)
    save(fig, "fig2_h4_misspecification")


# ------------------------------------------------------------------------------------------- Figure 3 (H5, dose-response)
def fig3(R):
    get = lambda lt, e: R[f"laneB_N2_H5_{lt}_e{e:g}_K10_pooled"]
    th = R["laneB_N2_H5_verdict"]["metrics"]["thresholds"]
    nb = pd.DataFrame(R["laneB_N2_negative_pass_rate_by_level"]["rows"])
    fig = plt.figure(figsize=(W2, 2.75))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.0], wspace=0.36, left=0.075, right=0.99, bottom=0.2,
                          top=0.79)
    xs = np.arange(4)

    def rates(ax, lt, letter, ttl, show_thr):
        tidy(ax)
        if show_thr:
            ax.axhline(th["pr_max"], color=INK2, lw=1.0, ls=(0, (4, 2)))
            ax.axhline(th["fals_pr"], color=INK2, lw=1.0, ls=(0, (1, 1.5)))
        for key, colr in (("phantom_rate", EXPR), ("phantom_rate_after", GATE)):
            p = np.array([get(lt, e)["metrics"][key] for e in LEVELS])
            lo = np.array([get(lt, e)["ci"][key][0] for e in LEVELS]); hi = np.array([get(lt, e)["ci"][key][1] for e in LEVELS])
            ax.plot(xs, p, color=colr, lw=1.5, zorder=2)
            ax.vlines(xs, lo, hi, color=colr, lw=1.5, zorder=2)
            ax.scatter(xs, p, s=36, color=colr, edgecolor="white", lw=1.2, zorder=3)
        ax.set_xticks(xs); ax.set_xticklabels(LEVLAB); ax.set_xlim(-0.25, 3.25); ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel("planted effect (× default)")
        title(ax, letter, ttl)

    ax = fig.add_subplot(gs[0, 0]); rates(ax, "default", "a", "Default lists (H5)", True)
    ax.set_ylabel("phantom rate (95% CI)")
    # (b) base rate: what the gate passes at each level; tick labels carry the true share of the lists
    ax = fig.add_subplot(gs[0, 1]); tidy(ax)
    d = nb[nb.list_type == "default"].set_index("effect").loc[LEVELS]
    tp, fp = d.true_passed.values.astype(float), d.neg_passed.values.astype(float)
    ax.bar(xs, tp, 0.62, color=GATE)
    ax.bar(xs, fp, 0.62, bottom=tp, color="white", edgecolor=INK2, hatch="////", lw=0.8)
    for x, a, b in zip(xs, tp, fp):
        ax.text(x, a + b + 12, f"{int(b)}/{int(a + b)}", ha="center", va="bottom", fontsize=FS)
    prev = [f"{p:.2f}" if p >= 0.1 else f"{p:.3f}" for p in d.prevalence_true.values]
    ax.set_xticks(xs); ax.set_xticklabels([f"{l}\n{p}" for l, p in zip(LEVLAB, prev)])
    ax.set_ylim(0, 740); ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel("claims passed (false/all)", labelpad=1)
    ax.set_xlabel("effect / true share of list", labelpad=2)
    title(ax, "b", "What the gate lets through")
    ax = fig.add_subplot(gs[0, 2]); rates(ax, "A5m", "c", "Hard-negative lists", False)
    handles = [Line2D([], [], color=EXPR, lw=1.5, marker="o", ms=6, mec="white"),
               Line2D([], [], color=GATE, lw=1.5, marker="o", ms=6, mec="white"),
               Line2D([], [], color=INK2, lw=1.0, ls=(0, (4, 2))), Line2D([], [], color=INK2, lw=1.0, ls=(0, (1, 1.5))),
               Patch(facecolor=GATE), Patch(facecolor="white", edgecolor=INK2, hatch="////", lw=0.8)]
    labels = ["expression-only tool", "frozen gate", f"H5 target ≤ {th['pr_max']:.2f}",
              f"falsified if > {th['fals_pr']:.2f}", "true passes", "false passes"]
    fig.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.53, 1.0), handlelength=1.8,
               columnspacing=1.6, labelspacing=0.3)
    save(fig, "fig3_h5_dose_response")


# ---------------------------------------------------------------------------------------------- Figure 4 (failure map)
def fig4(R):
    fig = plt.figure(figsize=(W2, 3.9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1, 1], wspace=0.36, hspace=0.62,
                          left=0.175, right=0.985, bottom=0.105, top=0.935)
    # (a) which gate component holds back the hard negatives (dev seeds, A5 lists)
    ab = pd.DataFrame(R["laneA_dev_ablation_failure_map"]["metrics"])
    ab = ab[ab.variant == "A5"].set_index("ablation")
    order = [("full", "full gate"), ("no_min_evidence", "− min. evidence"), ("no_perm", "− permutation null"),
             ("no_split", "− split halves"), ("no_tiles", "− tile restriction")]
    ax = fig.add_subplot(gs[0, 0]); tidy(ax, grid="x")
    y = np.arange(len(order))[::-1]
    for yi, (k, lab) in zip(y, order):
        p, lo, hi = ab.loc[k, "K10_phantom_after"], ab.loc[k, "K10_pr_lo"], ab.loc[k, "K10_pr_hi"]
        ax.plot([lo, hi], [yi, yi], color=GATE, lw=1.5, solid_capstyle="butt")
        ax.scatter([p], [yi], s=36, color=GATE if k == "full" else "white", edgecolor=GATE, lw=1.5, zorder=3)
        ax.text(hi + 0.012, yi, f"{p:.2f}", va="center", fontsize=FS, color=INK)
    ax.set_yticks(y); ax.set_yticklabels([lab for _, lab in order]); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 0.7); ax.set_xticks([0, 0.2, 0.4, 0.6]); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("phantom rate after gate (95% CI)")
    title(ax, "a", "Remove one gate component (dev, A5)")
    # (b) decoy senders pass when they co-localise with the true sender
    dc = R["laneA_N2_A5_decoy_colocalisation"]["by_coloc_bin"]
    ax = fig.add_subplot(gs[0, 1]); tidy(ax)
    labs = ["≤ 0.1", "0.1–0.3", "0.3–0.5", "> 0.5"]
    x = np.arange(len(dc))
    ax.bar(x, [b["pass_rate"] for b in dc], 0.6, color=GATE)
    for xi, b in zip(x, dc):
        ax.text(xi, b["pass_rate"] + 0.015, f"{b['pass_rate']:.2f}\n(n = {b['n']})", ha="center", va="bottom",
                fontsize=FS, color=INK, linespacing=1.0)
    ax.set_xticks(x); ax.set_xticklabels(labs)
    ax.set_ylim(0, 0.65); ax.set_yticks([0, 0.2, 0.4, 0.6]); ax.set_ylabel("pass rate")
    ax.set_xlabel("co-localisation of decoy and true sender")
    title(ax, "b", "Decoy senders: attribution")
    # (c) is the passing signal in the raw counts? (normalisation x density confound)
    cf = R["laneA_N2_normalisation_confound"]["metrics"]
    fm = lambda v, c: next(m for m in cf if m["variant"] == v and m["category"] == c and m["passed"])
    bars = [("true claims (A0)", fm("A0", "true")), ("decoy senders (A5)", fm("A5", "decoy_sender")),
            ("niche (A0)", fm("A0", "niche")), ("hard niche (A5)", fm("A5", "niche_hard"))]
    ax = fig.add_subplot(gs[1, 0]); tidy(ax, grid="x")
    y = np.arange(len(bars))[::-1]
    ax.barh(y, [m["frac_raw_supports"] for _, m in bars], 0.6, color=GATE)
    for yi, (lab, m) in zip(y, bars):
        ax.text(m["frac_raw_supports"] + 0.02, yi, f"{100 * m['frac_raw_supports']:.0f}%  (n = {m['n']})",
                va="center", fontsize=FS)
    ax.set_yticks(y); ax.set_yticklabels([lab for lab, _ in bars]); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1.32); ax.set_xticks([0, 0.5, 1.0]); ax.set_xticklabels(["0", "50%", "100%"])
    ax.set_ylim(-0.6, len(bars) - 0.4)
    ax.set_xlabel("gate passes that raw counts also support")
    title(ax, "c", "Is the passed signal in raw counts?")
    # (d) base rate: post-gate phantom rate vs the true share of the list
    nb = pd.DataFrame(R["laneB_N2_negative_pass_rate_by_level"]["rows"])
    ax = fig.add_subplot(gs[1, 1]); tidy(ax, grid="both")
    pi = np.logspace(np.log10(0.006), 0, 200)
    for f in (0.01, 0.03, 0.10):
        pr = (1 - pi) * f / ((1 - pi) * f + pi * 0.95)
        ax.plot(pi, pr, color=MUTED, lw=0.9, zorder=1)
        above = f == 0.10
        ax.text(0.0064, pr[0] + (0.03 if above else -0.04), f"f = {100 * f:.0f}%", fontsize=FS, color=INK2,
                va="bottom" if above else "top", ha="left")
    mk = {"default": ("o", GATE, "default lists"), "A5m": ("s", AQUA, "hard negatives"),
          "A5": ("^", EXPR, "hard neg. + mixed")}
    for lt, (m, colr, lab) in mk.items():
        d = nb[nb.list_type == lt]
        prv = d.prevalence_true.values; pr = d.neg_passed / (d.neg_passed + d.true_passed)
        ax.scatter(prv, pr, marker=m, s=36, color=colr, edgecolor="white", lw=0.8, zorder=3, label=lab)
    ax.set_xscale("log"); ax.set_xlim(0.006, 1); ax.set_ylim(-0.03, 1.1)
    ax.set_xticks([0.01, 0.1, 1]); ax.set_xticklabels(["0.01", "0.1", "1"])
    ax.set_xlabel("true share of the candidate list")
    ax.set_ylabel("phantom rate after gate", labelpad=2)
    ax.legend(loc="upper right", handletextpad=0.2, borderaxespad=0.2, labelspacing=0.25)
    title(ax, "d", "Base rate of the list")
    save(fig, "fig4_failure_map")


# ---------------------------------------------------------------------------------------------- appendix figures
def figA1():
    T = pd.read_csv("results/tables/verdict_table.csv")
    T = T[T.condition.isin(["expr_cpdb", "expr_prod", "agent_both"])].copy()
    lab = {"expr_cpdb": "cpdb", "expr_prod": "prod", "agent_both": "agents"}
    tis = {"breast_rep1": "breast Rep1", "breast_rep2": "breast Rep2", "pancreas": "pancreas", "lung": "lung"}
    fig = plt.figure(figsize=(W2, 3.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.9], wspace=0.5, left=0.235, right=0.985, bottom=0.14, top=0.84)
    ax = fig.add_subplot(gs[0, 0]); tidy(ax, grid="x")
    y = np.arange(len(T))[::-1]
    for yi, r in zip(y, T.itertuples()):
        lo, hi = ast.literal_eval(r.pr_ci)
        ax.plot([lo, hi], [yi, yi], color=EXPR, lw=1.5)
        ax.scatter([r.pr], [yi], s=30, color=EXPR, zorder=3)
        ax.scatter([r.pr_after], [yi], s=30, marker="D", color=GATE, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([f"{r.split} {tis[r.tissue]} · {lab[r.condition]}" for r in T.itertuples()],
                                         fontsize=FS)
    ax.tick_params(axis="y", length=0); ax.set_xlim(-0.02, 0.6)
    ax.axvline(0.30, color=INK2, lw=1.0, ls=(0, (4, 2)))
    ax.text(0.31, len(T) - 0.6, "H1 ≥ 0.30", fontsize=FS, color=INK2, va="bottom")
    ax.set_xlabel("phantom rate @10")
    fig.legend([Line2D([], [], marker="o", color=EXPR, lw=1.5), Line2D([], [], marker="D", color=GATE, lw=0)],
               ["phantom rate before gate (95% CI)", "after gate"], loc="upper left", ncol=2,
               bbox_to_anchor=(0.03, 1.0), handlelength=1.4)
    title(ax, "a", "Phase 1 lists")
    P = pd.read_csv("results/tables/planted_pass_rates.csv").set_index("setting")
    ax = fig.add_subplot(gs[0, 1]); tidy(ax, grid="x")
    y = np.arange(len(P))[::-1]
    for c, colr, m, nm in (("true", GATE, "o", "true"), ("niche", EXPR, "s", "niche"), ("phantom", INK2, "x", "phantom"),
                           ("random", AQUA, "^", "random")):
        ax.scatter(P[c].values, y, marker=m, s=30, color=colr, lw=1.2, zorder=3, label=nm)
    ax.set_yticks(y); ax.set_yticklabels([s.replace("dev:", "dev ").replace("test:", "test ").replace("_", " ")
                                          .replace("no ", "− ").replace("min evidence", "min. evid.")
                                          .replace("breast rep", "breast Rep") for s in P.index])
    ax.tick_params(axis="y", length=0); ax.set_xlim(-0.03, 1.05); ax.set_xlabel("gate pass rate, planted claims")
    ax.legend(loc="center", bbox_to_anchor=(0.66, 0.5), handletextpad=0.2, borderaxespad=0.1)
    title(ax, "b", "Planted claims by category")
    save(fig, "figA1_phase1")


def figA2():
    from bench.tissue import load_tissue
    t = load_tissue("breast_rep1")
    g = {x: i for i, x in enumerate(t.genes)}; X = t.counts.tocsc()
    pos = lambda gene: np.asarray(X[:, g[gene]].todense()).ravel() > 0
    dc = np.flatnonzero((t.types == "LAMP3+_DCs") & pos("CD274"))
    tc = np.flatnonzero(np.isin(t.types, ["CD4+_T_Cells", "CD8+_T_Cells"]) & pos("PDCD1"))
    fig = plt.figure(figsize=(W1, 2.45)); ax = fig.add_axes([0.0, 0.2, 1.0, 0.8])
    ax.scatter(t.xy[::6, 0], t.xy[::6, 1], s=0.15, c=CELL, lw=0, rasterized=True)
    ax.scatter(t.xy[tc, 0], t.xy[tc, 1], s=3, c=GATE, lw=0, rasterized=True)
    ax.scatter(t.xy[dc, 0], t.xy[dc, 1], s=4, c=EXPR, lw=0, rasterized=True)
    xm, ym = t.xy[:, 0].max(), t.xy[:, 1].max()
    ax.plot([200, 1200], [ym + 250, ym + 250], color=INK, lw=1.5, solid_capstyle="butt")
    ax.text(700, ym + 330, "1 mm", ha="center", va="top", fontsize=FS)
    ax.set_xlim(0, xm); ax.set_ylim(ym + 700, 0); ax.set_aspect("equal"); ax.axis("off")
    fig.text(0.02, 0.03, f"orange: CD274$^+$ LAMP3$^+$ DCs (n = {len(dc)})\nblue: PDCD1$^+$ CD4/CD8 T cells (n = {len(tc)})",
             fontsize=FS, color=INK2)
    save(fig, "figA2_real_pdl1")


SENS_LAB = {"frozen": "frozen gate", "tile125": "tiles 125 µm", "tile500": "tiles 500 µm", "r15": "range 15 µm",
            "r60": "range 60 µm", "nmin25": "min. contacts 25", "nmin100": "min. contacts 100",
            "split_checker": "checkerboard halves", "null_global": "global null",
            "null_density": "density-stratified null", "drop_edge": "border tiles dropped",
            "resp_raw": "raw-count response"}


def figA3():
    """Lane R sensitivity diagnostics (exploratory, dev seeds 0-4, breast Rep1)."""
    S = pd.read_csv("results/tables/laneR_sens_summary.csv")
    P = pd.read_csv("results/tables/laneR_sens_planted.csv")
    order = [s for s in SENS_LAB if s in set(S.setting)]
    fig = plt.figure(figsize=(W2, 3.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.12, left=0.19, right=0.985, bottom=0.14, top=0.9)
    y = np.arange(len(order))[::-1]
    panels = [("A5", "phantom_after", "a", "A5: phantom rate after gate", (-0.03, 0.7)),
              ("A0-A3", "retention", "b", "A0–A3: retention", (0, 1.05))]
    for k, (grp, key, letter, ttl, xl) in enumerate(panels):
        ax = fig.add_subplot(gs[0, k]); tidy(ax, grid="x")
        f = S[(S.group == grp)].set_index("setting")
        ref = f.loc["frozen", key]
        ax.axvline(ref, color=GATE, lw=1.0, ls=(0, (4, 2)))
        for yi, s in zip(y, order):
            p, lo, hi = f.loc[s, key], f.loc[s, f"{key}_lo"], f.loc[s, f"{key}_hi"]
            ax.plot([lo, hi], [yi, yi], color=GATE, lw=1.5, solid_capstyle="butt")
            ax.scatter([p], [yi], s=36, color=GATE if s == "frozen" else "white", edgecolor=GATE, lw=1.5, zorder=3)
        ax.set_yticks(y); ax.set_yticklabels([SENS_LAB[s] for s in order] if k == 0 else [])
        ax.tick_params(axis="y", length=0); ax.set_xlim(*xl); ax.set_ylim(-0.6, len(order) - 0.4)
        ax.set_xlabel("95% CI over dev seeds")
        title(ax, letter, ttl)
    ax = fig.add_subplot(gs[0, 2]); tidy(ax, grid="x")
    q = P[(P.group == "A5") & (P.category == "niche_hard")].set_index("setting")
    ref = q.loc["frozen", "pass_rate"]
    ax.axvline(ref, color=GATE, lw=1.0, ls=(0, (4, 2)))
    ax.barh(y, [q.loc[s, "pass_rate"] for s in order], 0.6, color=[GATE if s == "frozen" else MUTED for s in order])
    for yi, s in zip(y, order):
        ax.text(max(q.loc[s, "pass_rate"], ref) + 0.02, yi, f"{q.loc[s, 'pass_rate']:.2f}", va="center", fontsize=FS)
    ax.set_yticks(y); ax.set_yticklabels([]); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 0.75); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel(f"pass rate (n = {int(q.loc['frozen', 'n'])} per setting)")
    title(ax, "c", "A5 hard-niche claims passed")
    save(fig, "figA3_sensitivity")


if __name__ == "__main__":
    want = set(sys.argv[1:]) or {"fig1", "fig2", "fig3", "fig4", "figA1", "figA2", "figA3"}
    R = records()
    if "fig1" in want: fig1()
    if "fig2" in want: fig2(R)
    if "fig3" in want: fig3(R)
    if "fig4" in want: fig4(R)
    if "figA1" in want: figA1()
    if "figA2" in want: figA2()
    if "figA3" in want: figA3()
