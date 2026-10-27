"""Draw every figure in the paper from the CSV files in results/. No number is typed in here.

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
STYLE = {
    "Mondrian conformal": dict(color="#2a78d6", marker="o", ls="-"),
    "Label-shift weighted conformal": dict(color="#1baf7a", marker="s", ls="--"),
    "Pooled conformal": dict(color="#eb6834", marker="^", ls=":"),
    "Fixed 0.1/0.9 tiers": dict(color="#eda100", marker="D", ls="-."),
    "Platt (source prior)": dict(color="#eb6834", marker="^", ls=":"),
    "EM prior-adjusted": dict(color="#2a78d6", marker="o", ls="-"),
    "Oracle prior-adjusted": dict(color="#6b6b6b", marker="", ls="--"),
}
INK, MUTED = "#2b2b2b", "#8a8a8a"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "lines.linewidth": 1.6,
                     "lines.markersize": 4.5, "pdf.fonttype": 42})


def fig_simulation():
    sim = pd.read_csv(RES / "simulation.csv")
    g = sim.groupby(["kind", "method", "pi_t"]).agg(["mean", "std", "count"]).reset_index()
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.6), layout="constrained")
    cal = g[g.kind == "calibration"]
    for m in ["Platt (source prior)", "EM prior-adjusted", "Oracle prior-adjusted"]:
        s = cal[cal.method == m]
        ax[0].plot(s.pi_t, s[("ece", "mean")], label=m, **STYLE[m])
    ax[0].set(xlabel="LoF share among target variants", ylabel="Expected calibration error",
              title="A  Calibration under prior shift")
    ax[0].axvline(0.45, color=MUTED, lw=0.8, ls=":")
    ax[0].text(0.46, ax[0].get_ylim()[1] * 0.55, "source\nprior", fontsize=6.5, color=INK, va="center")
    ax[0].legend(frameon=False, loc="upper right")
    con = g[g.kind == "conformal"]
    for j, (col, title) in enumerate([("coverage_lof", "B  Coverage of LoF variants"),
                                      ("abstain", "C  Abstention rate")]):
        for m in ["Mondrian conformal", "Label-shift weighted conformal", "Pooled conformal"]:
            s = con[con.method == m]
            mean, se = s[(col, "mean")], s[(col, "std")] / np.sqrt(s[(col, "count")])
            ax[j + 1].errorbar(s.pi_t, mean, yerr=1.96 * se, label=m, capsize=1.5, **STYLE[m])
        ax[j + 1].set(xlabel="LoF share among target variants", title=title)
    ax[1].axhline(0.90, color=INK, lw=0.8, ls="--")
    ax[1].text(0.6, 0.902, "target 0.90", ha="right", va="bottom", fontsize=6.5, color=INK)
    ax[1].set_ylabel("Fraction of LoF variants covered")
    ax[2].set_ylabel("Fraction of variants abstained")
    h, l = ax[1].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, loc="outside lower center", ncol=3)
    fig.savefig(FIG / "fig_simulation.pdf")
    plt.close(fig)


METHODS = ["Mondrian conformal", "Pooled conformal", "Fixed 0.1/0.9 tiers", "AlphaMissense classes",
           "REVEL PP3/BP4 supporting"]
SHORT = {"Mondrian conformal": "Per-class (ours)", "Pooled conformal": "Pooled\nconformal",
         "Label-shift weighted conformal": "Weighted\n(EM prior)",
         "Label-shift weighted, true prior": "Weighted\n(true prior)", "Fixed 0.1/0.9 tiers": "Fixed\n0.1/0.9",
         "AlphaMissense classes": "AM\nclasses", "REVEL PP3/BP4 supporting": "REVEL\nPP3/BP4"}
GENES = ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]
PRED_MARK = {"AlphaMissense": "o", "REVEL": "s", "ESM-2": "^"}


def fig_overview():
    co = pd.read_csv(RES / "cohort.csv").set_index("gene").loc[GENES]
    fig = plt.figure(figsize=(7.0, 2.3), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1])
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(-0.1, 10.1); ax.set_ylim(0, 4.15); ax.axis("off")

    def box(x, y, w, h, text, fc):
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15",
                                                       fc=fc, ec=MUTED, lw=0.8))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=5.9, color=INK)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))

    src_prev = json.loads((RES / "paper_numbers.json").read_text())["source_prev"]
    box(0.0, 3.0, 2.55, 1.05, f"ClinVar labels, \u22652 stars\n{src_prev:.0%} pathogenic", "#e9f1fb")
    box(0.0, 1.55, 2.55, 1.05, "Assay labels, part A\n(one variant\nper residue)", "#e8f6f0")
    box(0.0, 0.1, 2.55, 1.05, "Assay labels, part B\n(one variant\nper residue)", "#e8f6f0")
    box(3.2, 3.0, 2.6, 1.05, "Platt score\nq(s) = P(y=1 | s)", "#f5f5f3")
    box(3.2, 1.55, 2.6, 1.05, "Per-class thresholds\n" + r"miss rate $\leq\alpha$" + "\n" + r"with prob. $\geq 1-\delta$", "#f5f5f3")
    box(3.2, 0.1, 2.6, 1.05, "Certify likelihood\nratios with exact\nbinomial bounds", "#f5f5f3")
    box(6.45, 1.0, 3.6, 2.4, "Calls: LoF, functional\nor abstain\n\neach call type with a\ncertified OddsPath level\n(vs. assay; no prior needed)", "#fdf3e1")
    arrow(2.6, 3.52, 3.15, 3.52); arrow(2.6, 2.07, 3.15, 2.07); arrow(2.6, 0.62, 3.15, 0.62)
    arrow(4.5, 3.0, 4.5, 2.65); arrow(4.5, 1.55, 4.5, 1.2)
    arrow(5.85, 2.07, 6.4, 2.2); arrow(5.85, 0.62, 6.4, 1.5)
    ax.set_title("A  Certified, prior-free variant calls", loc="left")

    ax = fig.add_subplot(gs[1])
    x = np.arange(len(GENES)); w = 0.27
    ax.bar(x - w, co.src_prev, w, color="#2a78d6", label="ClinVar labelled (pathogenic)")
    ax.bar(x, co.target_prev, w, color="#1baf7a", label="Assayed, held out (LoF)")
    ax.bar(x + w, co.vus_prev, w, color="#eda100", hatch="///", edgecolor="white", lw=0, label="Assayed VUS (LoF)")
    for i, g in enumerate(GENES):
        ax.text(i - w, co.src_prev.iloc[i] + 0.015, f"n={int(co.src_labelled.iloc[i])}", ha="center", fontsize=5.2,
                color=INK, rotation=90, va="bottom")
    ax.set_xticks(x, GENES)
    ax.set_ylabel("Share damaging"); ax.set_ylim(0, 1.0)
    ax.set_title("B  The damaging share is not stable", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    fig.savefig(FIG / "fig_overview.pdf")
    plt.close(fig)


def strip(ax, c, title, labels=True):
    meths = [m for m in METHODS if m in set(c.method)]
    for i, m in enumerate(meths):
        s = c[c.method == m]
        for pr, mk in PRED_MARK.items():
            v = s[s.predictor == pr]
            if len(v) == 0:
                continue
            jitter = (np.arange(len(v)) - (len(v) - 1) / 2) * 0.06 + (list(PRED_MARK).index(pr) - 1) * 0.2
            ax.scatter(v.worst_rw, i + jitter, marker=mk, s=13, color=STYLE[m]["color"], edgecolor="white",
                       lw=0.4, zorder=3)
    ax.axvline(0.90, color=INK, lw=0.8, ls="--")
    ax.set_yticks(range(len(meths)), [SHORT[m].replace("\n", " ") for m in meths] if labels else [])
    ax.invert_yaxis()
    ax.set_xlim(0.4, 1.01)
    ax.set_xlabel("Worst-class coverage")
    ax.set_title(title, loc="left")


def fig_main():
    """Per-class guarantees: method comparison, split-to-split reliability, real VUS shift."""
    conf = pd.read_csv(RES / "conformal.csv")
    c = conf[(conf.alpha == 0.10) & (conf.design == "within-gene random split")].copy()
    c = c[c.method.isin(["Mondrian conformal", "Pooled conformal", "Fixed 0.1/0.9 tiers",
                         "AlphaMissense classes", "REVEL PP3/BP4 supporting"])]
    c["worst_rw"] = c[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    STYLE.setdefault("AlphaMissense classes", dict(color="#e87ba4", marker="D", ls="-"))
    STYLE.setdefault("REVEL PP3/BP4 supporting", dict(color="#e34948", marker="D", ls="-"))
    rel = pd.read_csv(RES / "reliability.csv")
    rel["worst"] = rel[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    vs = pd.read_csv(RES / "vus_shift.csv")
    vs = vs[vs.part == "coverage"].groupby(["thresholds", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw"]].mean().reset_index()
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.55), layout="constrained", gridspec_kw={"width_ratios": [1.25, 0.9, 0.95]})
    strip(ax[0], c, "A  Worst-class coverage, one gene")
    bins = np.linspace(0.70, 1.0, 31)
    for name, col, hatch in (("in expectation", "#2a78d6", None), ("PAC, delta = 0.1", "#1baf7a", "///")):
        v = rel[rel.thresholds == name].worst
        ax[1].hist(v, bins=bins, color=col, alpha=0.75, label="in expectation" if name == "in expectation"
                   else r"PAC ($\delta=0.1$)", hatch=hatch, edgecolor="white", lw=0.3)
    ax[1].axvline(0.90, color=INK, lw=0.8, ls="--")
    ax[1].set(xlabel="Worst-class coverage per split", ylabel="Splits", title="B  Split-to-split reliability")
    ax[1].legend(frameon=False, fontsize=6, loc="upper left")
    for j, (name, col, mk) in enumerate((("in expectation", "#2a78d6", "o"), ("PAC, delta = 0.1", "#1baf7a", "s"))):
        v = vs[vs.thresholds == name]
        ax[2].scatter(v.coverage_func_rw, v.coverage_lof_rw, marker=mk, s=14, color=col, edgecolor="white", lw=0.4,
                      label="in expectation" if j == 0 else r"PAC ($\delta=0.1$)", zorder=3)
    ax[2].axhline(0.90, color=INK, lw=0.8, ls="--"); ax[2].axvline(0.90, color=INK, lw=0.8, ls="--")
    ax[2].set(xlim=(0.78, 1.01), ylim=(0.78, 1.01), xlabel="Functional coverage on VUS", ylabel="LoF coverage on VUS",
              title="C  Real shift: test on VUS")
    ax[2].legend(frameon=False, fontsize=6, loc="lower right", handletextpad=0.2)
    handles = [plt.Line2D([], [], marker=mk, ls="", color=MUTED, label=pr) for pr, mk in PRED_MARK.items()]
    fig.legend(handles=handles, frameon=False, loc="outside lower center", ncol=3)
    fig.savefig(FIG / "fig_main.pdf")
    plt.close(fig)


def fig_evidence():
    """Selection rules judged against the exact population likelihood ratio (oracle design)."""
    ev = pd.read_csv(RES / "oracle.csv")
    levels = ["PS3-supporting", "PS3-moderate", "PS3-strong", "BS3-supporting", "BS3-moderate", "BS3-strong"]
    sels = [("certified", "Certified (ours)", "#2a78d6", None),
            ("bootstrap", "Bootstrap 95%", "#eda100", "///"),
            ("bootstrap, Bonferroni", "Bootstrap, same error budget", "#1baf7a", "..."),
            ("point estimate", "Point estimate", "#eb6834", None)]
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.6), layout="constrained", gridspec_kw={"width_ratios": [1.2, 1.2, 0.7]})
    x = np.arange(len(levels)); w = 0.2
    for i, (sel, lab, col, hatch) in enumerate(sels):
        g = ev[ev.selection == sel]
        ach = g.groupby("level").achieved.mean().reindex(levels)
        held = g[g.achieved].groupby("level").holds.mean().reindex(levels)
        off = (i - 1.5) * w
        ax[0].bar(x + off, ach.values, w, color=col, hatch=hatch, edgecolor="white", lw=0, label=lab)
        ax[1].bar(x + off, held.fillna(0).values, w, color=col, hatch=hatch, edgecolor="white", lw=0)
    names = [l.replace("-supporting", "\nsupp.").replace("-moderate", "\nmod.").replace("-strong", "\nstrong") for l in levels]
    for a_ in ax[:2]:
        a_.set_xticks(x, names, fontsize=6); a_.set_ylim(0, 1.05)
    ax[0].set(ylabel="Share of draws", title="A  Level claimed")
    ax[1].set(ylabel="Share of claims", title="B  Level truly held")
    ax[1].axhline(0.95, color=INK, lw=0.8, ls="--")
    acc = ev[ev.achieved]
    order = [("certified", "certified"), ("Clopper-Pearson, uncorrected", "CP, no\ncorrection"),
             ("bootstrap", "bootstrap"), ("bootstrap, Bonferroni", "bootstrap,\nsame budget"), ("point estimate", "point")]
    cols = ["#2a78d6", "#6250d6", "#eda100", "#1baf7a", "#eb6834"]
    vals = [acc[acc.selection == s_].bound_ok.mean() for s_, _ in order]
    ax[2].bar(range(5), vals, color=cols, edgecolor="white")
    for i, v in enumerate(vals):
        ax[2].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=5.8, color=INK)
    ax[2].axhline(0.95, color=INK, lw=0.8, ls="--")
    ax[2].set_xticks(range(5), [l for _, l in order], fontsize=5.2, rotation=90)
    ax[2].set(ylim=(0, 1.1), ylabel="True LR respects bound", title="C  Bound validity")
    fig.legend(*ax[0].get_legend_handles_labels(), frameon=False, loc="outside lower center", ncol=4)
    fig.savefig(FIG / "fig_evidence.pdf")
    plt.close(fig)


def fig_practical():
    cs = pd.read_csv(RES / "calibration_size.csv")
    pri = pd.read_csv(RES / "prior_estimates.csv")
    conf = pd.read_csv(RES / "conformal.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3), layout="constrained")
    g = cs.groupby("n_cal_residues")
    for col, lab, st in [("abstain", "Abstention", dict(color="#eb6834", marker="^", ls=":")),
                         ("coverage_lof_rw", "LoF coverage", dict(color="#2a78d6", marker="o", ls="-")),
                         ("coverage_func_rw", "Functional coverage", dict(color="#1baf7a", marker="s", ls="--"))]:
        mu = g[col].mean()
        ax[0].plot(mu.index, mu.values, label=lab, **st)
    ax[0].axhline(0.90, color=INK, lw=0.8, ls="--")
    ax[0].set_xscale("log"); ax[0].set_xticks([25, 50, 100, 150, 200], ["25", "50", "100", "150", "200"])
    ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax[0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[0].set(xlabel="Calibration residues in the gene", ylabel="Fraction", title="A  How much assay data is needed")
    ax[0].legend(frameon=False, fontsize=6)
    p3 = pri[pri.predictor.isin(PRED_MARK)]
    for pr, mk in PRED_MARK.items():
        v = p3[p3.predictor == pr]
        ax[1].scatter(v.observed, v.em, marker=mk, s=16, color="#2a78d6", edgecolor="white", lw=0.4, label=pr)
    ax[1].plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls="--")
    ax[1].set(xlim=(0, 1), ylim=(0, 1), xlabel="Observed LoF share", ylabel="EM estimate from scores",
              title="B  Prior estimated from scores")
    ax[1].legend(frameon=False, fontsize=6, loc="upper left")
    m = conf[(conf.alpha == 0.10) & (conf.method == "Mondrian conformal")].copy()
    m["worst_rw"] = m[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    des = ["within-gene random split", "within-gene region split", "leave-one-gene-out"]
    for j, dname in enumerate(des):
        s_ = m[m.design == dname]
        for pr, mk in PRED_MARK.items():
            v = s_[s_.predictor == pr]
            jitter = (np.arange(len(v)) - 2) * 0.06 + (list(PRED_MARK).index(pr) - 1) * 0.2
            ax[2].scatter(v.worst_rw, j + jitter, marker=mk, s=13, color="#2a78d6", edgecolor="white", lw=0.4, zorder=3)
    ax[2].axvline(0.90, color=INK, lw=0.8, ls="--")
    ax[2].set_yticks(range(3), ["random\nresidues", "other\nregion", "other\ngenes"]); ax[2].invert_yaxis()
    ax[2].set_xlim(0.4, 1.01)
    ax[2].set(xlabel="Worst-class coverage", title="C  Where calibration comes from")
    fig.savefig(FIG / "fig_practical.pdf")
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
