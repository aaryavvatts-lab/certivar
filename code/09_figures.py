"""Draw every figure in the paper from the CSV files in results/. No number is typed in here.

Every figure is drawn at its printed size (6.5 in, the text width of a US Letter page with 1-inch margins) and
included at 100% scale, so all text in the figures is at least 10 pt on the page.

Palette: validated categorical slots (blue, orange, aqua, yellow), each method also has its
own marker and line style so identity never depends on colour alone.
"""
import json
import pathlib

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
import matplotlib.patches  # noqa: E402,F401

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = RES / "figures"
WIDTH = 6.5  # inches: text width of the paper
STYLE = {
    "Mondrian conformal": dict(color="#2a78d6", marker="o", ls="-"),
    "Label-shift weighted conformal": dict(color="#1baf7a", marker="s", ls="--"),
    "Pooled conformal": dict(color="#eb6834", marker="^", ls=":"),
    "Fixed 0.1/0.9 tiers": dict(color="#eda100", marker="D", ls="-."),
    "AlphaMissense classes": dict(color="#e87ba4", marker="D", ls="-"),
    "REVEL PP3/BP4 supporting": dict(color="#e34948", marker="D", ls="-"),
    "Platt (source prior)": dict(color="#eb6834", marker="^", ls=":"),
    "EM prior-adjusted": dict(color="#2a78d6", marker="o", ls="-"),
    "Oracle prior-adjusted": dict(color="#6b6b6b", marker="", ls="--"),
}
INK, MUTED = "#2b2b2b", "#8a8a8a"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 10,
                     "xtick.labelsize": 10, "ytick.labelsize": 10, "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "lines.linewidth": 1.6,
                     "lines.markersize": 5, "legend.handletextpad": 0.4, "legend.columnspacing": 1.2,
                     "figure.constrained_layout.w_pad": 0.15, "pdf.fonttype": 42})

METHODS = ["Mondrian conformal", "Pooled conformal", "Fixed 0.1/0.9 tiers", "AlphaMissense classes",
           "REVEL PP3/BP4 supporting"]
SHORT = {"Mondrian conformal": "Per-class (ours)", "Pooled conformal": "Pooled conformal",
         "Fixed 0.1/0.9 tiers": "Fixed 0.1/0.9", "AlphaMissense classes": "AM classes",
         "REVEL PP3/BP4 supporting": "REVEL PP3/BP4"}
GENES = ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]
PRED_MARK = {"AlphaMissense": "o", "REVEL": "s", "ESM-2": "^"}


def predictor_legend(fig, **kw):
    handles = [plt.Line2D([], [], marker=mk, ls="", color=MUTED, markersize=6, label=pr) for pr, mk in PRED_MARK.items()]
    fig.legend(handles=handles, frameon=False, loc="outside lower center", ncol=3, **kw)


def fig_overview():
    """Schematic of the method (the damaging shares shown in earlier drafts are given in Table 1)."""
    src_prev = json.loads((RES / "paper_numbers.json").read_text())["source_prev"]
    fig = plt.figure(figsize=(WIDTH, 1.62))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH); ax.set_ylim(0, 1.62); ax.axis("off")

    def box(x, y, w, h, text, fc):
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                                       fc=fc, ec=MUTED, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10, color=INK, linespacing=1.15)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))

    h, rows = 0.44, (1.15, 0.60, 0.04)
    box(0.03, rows[0], 1.85, h, f"ClinVar P/LP vs B/LB\n(≥2 stars, {src_prev:.0%} P/LP)", "#e9f1fb")
    box(0.03, rows[1], 1.85, h, "Assay labels, part A\n(one variant per residue)", "#e8f6f0")
    box(0.03, rows[2], 1.85, h, "Assay labels, part B\n(one variant per residue)", "#e8f6f0")
    box(2.20, rows[0], 1.95, h, "Platt scaling: common\nscore scale $q(s)$", "#f5f5f3")
    box(2.20, rows[1], 1.95, h, "Per-class thresholds\n" + r"(miss rate $\leq\alpha$)", "#f5f5f3")
    box(2.20, rows[2], 1.95, h, "Certify likelihood ratios\n(exact, simultaneous)", "#f5f5f3")
    box(4.47, 0.27, 2.00, 1.08, "Calls: LoF, functional\nor abstain, each type\nwith a certified LR\nfor the assay outcome", "#fdf3e1")
    for y in rows:
        arrow(1.90, y + h / 2, 2.18, y + h / 2)
    arrow(3.175, rows[0], 3.175, rows[1] + h + 0.01); arrow(3.175, rows[1], 3.175, rows[2] + h + 0.01)
    arrow(4.17, rows[1] + h / 2, 4.45, 0.95); arrow(4.17, rows[2] + h / 2, 4.45, 0.55)
    fig.savefig(FIG / "fig_overview.pdf")
    plt.close(fig)


def strip(ax, c, title):
    meths = [m for m in METHODS if m in set(c.method)]
    for i, m in enumerate(meths):
        s = c[c.method == m]
        for pr, mk in PRED_MARK.items():
            v = s[s.predictor == pr]
            if len(v) == 0:
                continue
            jitter = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.06 + (list(PRED_MARK).index(pr) - 1) * 0.2
            ax.scatter(v.worst_rw, i + jitter, marker=mk, s=18, color=STYLE[m]["color"], edgecolor="white",
                       lw=0.4, zorder=3)
    ax.axvline(0.90, color=INK, lw=0.8, ls="--")
    ax.set_yticks(range(len(meths)), [SHORT[m] for m in meths])
    ax.invert_yaxis()
    ax.set_xlim(0.4, 1.02)
    ax.set_xticks([0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel("Worst-class coverage")
    ax.set_title(title, loc="left")


def fig_main():
    """Per-class guarantees: method comparison, split-to-split reliability, real VUS shift."""
    conf = pd.read_csv(RES / "conformal.csv")
    c = conf[(conf.alpha == 0.10) & (conf.design == "within-gene random split")].copy()
    c = c[c.method.isin(METHODS)]
    c["worst_rw"] = c[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    rel = pd.read_csv(RES / "reliability.csv")
    rel["worst"] = rel[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    vs = pd.read_csv(RES / "vus_shift.csv")
    vs = vs[vs.part == "coverage"].groupby(["thresholds", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw"]].mean().reset_index()
    fig, ax = plt.subplots(1, 3, figsize=(WIDTH, 2.7), layout="constrained", gridspec_kw={"width_ratios": [1.15, 1, 1]})
    strip(ax[0], c, "A  By method")
    bins = np.linspace(0.70, 1.0, 31)
    for name, col, hatch, lab in (("in expectation", "#2a78d6", None, "in exp."),
                                  ("PAC, delta = 0.1", "#1baf7a", "///", "PAC")):
        v = rel[rel.thresholds == name].worst
        ax[1].hist(v, bins=bins, color=col, alpha=0.75, label=lab, hatch=hatch, edgecolor="white", lw=0.3)
    ax[1].axvline(0.90, color=INK, lw=0.8, ls="--")
    ax[1].set_xticks([0.7, 0.8, 0.9, 1.0])
    ax[1].set(xlabel="Worst-class coverage", ylabel="Splits", title="B  Per split")
    ax[1].legend(frameon=False, loc="upper left")
    for j, (name, col, mk) in enumerate((("in expectation", "#2a78d6", "o"), ("PAC, delta = 0.1", "#1baf7a", "s"))):
        v = vs[vs.thresholds == name]
        ax[2].scatter(v.coverage_func_rw, v.coverage_lof_rw, marker=mk, s=20, color=col, edgecolor="white", lw=0.4,
                      label="in exp." if j == 0 else "PAC", zorder=3)
    ax[2].axhline(0.90, color=INK, lw=0.8, ls="--"); ax[2].axvline(0.90, color=INK, lw=0.8, ls="--")
    ax[2].set(xlim=(0.80, 1.01), ylim=(0.80, 1.01), xlabel="Func. coverage", ylabel="LoF coverage",
              title="C  Test on VUS")
    ax[2].set_xticks([0.8, 0.9, 1.0]); ax[2].set_yticks([0.8, 0.9, 1.0])
    predictor_legend(fig)
    fig.savefig(FIG / "fig_main.pdf")
    plt.close(fig)


# selection rules of the oracle experiment: (key in oracle.csv, legend label, colour)
RULES = [("certified", "Certified (ours)", "#2a78d6"),
         ("Clopper-Pearson, uncorrected", "Exact, no correction", "#6250d6"),
         ("bootstrap", "Bootstrap 95%", "#eda100"),
         ("bootstrap, Bonferroni", "Bootstrap, corrected", "#1baf7a"),
         ("point estimate", "Point estimate", "#eb6834")]
# internal level keys of 13_evidence.py, printed by call type and strength (LR+ > 2.1, 4.3, 18.7; LR- < 0.48, 0.23, 0.053)
LEVELS = ["PS3-supporting", "PS3-moderate", "PS3-strong", "BS3-supporting", "BS3-moderate", "BS3-strong"]


def fig_evidence():
    """Selection rules judged against the exact population likelihood ratio (oracle design)."""
    ev = pd.read_csv(RES / "oracle.csv")
    ci = json.loads((RES / "paper_numbers.json").read_text())["oracle_validity_ci"]
    fig, ax = plt.subplots(1, 2, figsize=(WIDTH, 2.8), layout="constrained", gridspec_kw={"width_ratios": [2.6, 1]})
    x = np.arange(len(LEVELS)) + np.array([0, 0, 0, 0.5, 0.5, 0.5])
    w = 0.16
    for i, (sel, lab, col) in enumerate(RULES):
        g = ev[ev.selection == sel].groupby("level")
        claimed = g.achieved.mean().reindex(LEVELS).fillna(0).values
        held = g.apply(lambda t: (t.achieved & t.holds.fillna(False).astype(bool)).mean()).reindex(LEVELS).fillna(0).values
        off = (i - 2) * w
        ax[0].bar(x + off, held, w, color=col, edgecolor="white", lw=0, label=lab)
        ax[0].bar(x + off, claimed - held, w, bottom=held, color="white", edgecolor=col, lw=0.8, hatch="////")
    ax[0].set_xticks(x, ["supp.", "mod.", "strong"] * 2)
    for xc, lab in ((1, "LoF calls"), (4.5, "functional calls")):
        ax[0].text(xc, -0.22, lab, ha="center", va="top", transform=ax[0].get_xaxis_transform())
    ax[0].set_ylim(0, 1.05)
    ax[0].set(ylabel="Share of draws", title="A  Level claimed (bar) and truly held (solid)")
    vals = [ci[s]["ok"] / ci[s]["n"] for s, *_ in RULES]
    err = np.array([[v - ci[s]["lo"], ci[s]["hi"] - v] for v, (s, *_) in zip(vals, RULES)]).T
    ax[1].bar(range(len(RULES)), vals, color=[r[2] for r in RULES], edgecolor="white",
              yerr=np.clip(err, 0, None), error_kw=dict(lw=0.9, capsize=2, ecolor=INK))
    for i, v in enumerate(vals):
        ax[1].text(i, v / 2, f"{v:.2f}", ha="center", va="center", rotation=90, color="white", fontweight="bold")
    ax[1].axhline(0.95, color=INK, lw=0.8, ls="--")
    ax[1].set_xticks([])
    ax[1].set(ylim=(0, 1.05), yticks=[0, 0.25, 0.5, 0.75, 1.0], ylabel="True LR within bound", title="B  Validity")
    fig.legend(*ax[0].get_legend_handles_labels(), frameon=False, loc="outside lower center", ncol=3)
    fig.savefig(FIG / "fig_evidence.pdf")
    plt.close(fig)


def fig_practical():
    cs = pd.read_csv(RES / "calibration_size.csv")
    pri = pd.read_csv(RES / "prior_estimates.csv")
    conf = pd.read_csv(RES / "conformal.csv")
    fig, ax = plt.subplots(1, 3, figsize=(WIDTH, 2.85), layout="constrained", gridspec_kw={"width_ratios": [1, 1, 1.1]})
    g = cs.groupby("n_cal_residues")
    for col, lab, st in [("abstain", "abstention", dict(color="#eb6834", marker="^", ls=":")),
                         ("coverage_lof_rw", "LoF cov.", dict(color="#2a78d6", marker="o", ls="-")),
                         ("coverage_func_rw", "func. cov.", dict(color="#1baf7a", marker="s", ls="--"))]:
        mu = g[col].mean()
        ax[0].plot(mu.index, mu.values, label=lab, **st)
    ax[0].axhline(0.90, color=INK, lw=0.8, ls="--")
    ax[0].set_xscale("log"); ax[0].set_xticks([25, 50, 100, 200], ["25", "50", "100", "200"])
    ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[0].set_ylim(0.35, 1.02)
    ax[0].set(xlabel="Calibration residues", ylabel="Fraction", title="A  Assay data")
    p3 = pri[pri.predictor.isin(PRED_MARK)]
    for pr, mk in PRED_MARK.items():
        v = p3[p3.predictor == pr]
        ax[1].scatter(v.observed, v.em, marker=mk, s=22, color="#2a78d6", edgecolor="white", lw=0.4, label=pr)
    ax[1].plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls="--")
    ax[1].set(xlim=(0, 0.8), ylim=(0, 0.8), xlabel="Observed LoF share", ylabel="EM estimate",
              title="B  Estimated prior")
    ax[1].set_xticks([0, 0.4, 0.8]); ax[1].set_yticks([0, 0.4, 0.8])
    m = conf[(conf.alpha == 0.10) & (conf.method == "Mondrian conformal")].copy()
    m["worst_rw"] = m[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    des = ["within-gene random split", "within-gene region split", "leave-one-gene-out"]
    for j, dname in enumerate(des):
        s_ = m[m.design == dname]
        for pr, mk in PRED_MARK.items():
            v = s_[s_.predictor == pr]
            jitter = (np.arange(len(v)) - 2) * 0.06 + (list(PRED_MARK).index(pr) - 1) * 0.2
            ax[2].scatter(v.worst_rw, j + jitter, marker=mk, s=18, color="#2a78d6", edgecolor="white", lw=0.4, zorder=3)
    ax[2].axvline(0.90, color=INK, lw=0.8, ls="--")
    ax[2].set_yticks(range(3), ["same gene", "other\nregion", "other\ngenes"]); ax[2].invert_yaxis()
    ax[2].set_xlim(0.4, 1.02); ax[2].set_xticks([0.4, 0.6, 0.8, 1.0])
    ax[2].set(xlabel="Worst-class coverage", title="C  Calibration source")
    lines = ax[0].get_legend_handles_labels()
    marks = [plt.Line2D([], [], marker=mk, ls="", color=MUTED, markersize=6) for mk in PRED_MARK.values()]
    fig.legend(lines[0] + marks, lines[1] + list(PRED_MARK), frameon=False, loc="outside lower center", ncol=6,
               handlelength=1.6, columnspacing=0.9)
    fig.savefig(FIG / "fig_practical.pdf")
    plt.close(fig)


def fig_simulation():
    sim = pd.read_csv(RES / "simulation.csv")
    g = sim.groupby(["kind", "method", "pi_t"]).agg(["mean", "std", "count"]).reset_index()
    labels = {"Platt (source prior)": "Platt, source prior", "EM prior-adjusted": "EM-adjusted",
              "Oracle prior-adjusted": "true-prior adjusted", "Mondrian conformal": "Per-class (ours)",
              "Label-shift weighted conformal": "Weighted (EM prior)", "Pooled conformal": "Pooled"}
    fig, ax = plt.subplots(1, 3, figsize=(WIDTH, 3.1), layout="constrained")
    cal = g[g.kind == "calibration"]
    for m in ["Platt (source prior)", "EM prior-adjusted", "Oracle prior-adjusted"]:
        s = cal[cal.method == m]
        ax[0].plot(s.pi_t, s[("ece", "mean")], label=labels[m], **STYLE[m])
    ax[0].set(xlabel="Target LoF share", ylabel="Calibration error", title="A  Calibration")
    ax[0].axvline(0.45, color=MUTED, lw=0.8, ls=":")
    ax[0].set_ylim(0, 0.36)
    ax[0].legend(frameon=False, loc="upper right", handlelength=1.6, borderaxespad=0.1)
    con = g[g.kind == "conformal"]
    for j, (col, title, ylab) in enumerate([("coverage_lof", "B  LoF coverage", "Fraction covered"),
                                            ("abstain", "C  Abstention", "Fraction abstained")]):
        for m in ["Mondrian conformal", "Label-shift weighted conformal", "Pooled conformal"]:
            s = con[con.method == m]
            mean, se = s[(col, "mean")], s[(col, "std")] / np.sqrt(s[(col, "count")])
            ax[j + 1].errorbar(s.pi_t, mean, yerr=1.96 * se, label=labels[m], capsize=1.5, **STYLE[m])
        ax[j + 1].set(xlabel="Target LoF share", ylabel=ylab, title=title)
    ax[1].axhline(0.90, color=INK, lw=0.8, ls="--")
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, loc="outside lower center", ncol=3)
    fig.savefig(FIG / "fig_simulation.pdf")
    plt.close(fig)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    fig_simulation()
    fig_overview()
    fig_main()
    fig_evidence()
    fig_practical()


if __name__ == "__main__":
    main()
