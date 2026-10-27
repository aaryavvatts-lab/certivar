"""Write LaTeX tables from results/*.csv (no typed numbers). All tables are set at the 10 pt body size."""
import json
import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES, OUT = ROOT / "results", ROOT / "paper" / "tables"
ASSAY = {"BRCA1": ("SGE", "HAP1", "\\cite{Findlay2018}", "authors"),
         "BRCA2": ("SGE", "HAP1", "\\cite{Huang2025}", "controls"),
         "BAP1": ("SGE", "HAP1", "\\cite{Waters2024}", "authors"),
         "RAD51C": ("SGE", "HAP1", "\\cite{OlveraLen2024}", "authors"),
         "MSH2": ("DMS", "HAP1", "\\cite{Jia2021}", "authors")}
SHORT_PRED = {"AlphaMissense": "AM", "REVEL": "REVEL", "ESM-2": "ESM-2"}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    co = pd.read_csv(RES / "cohort.csv")
    rows = []
    for _, r in co.iterrows():
        a = ASSAY[r.gene]
        rows.append(f"{r.gene} & {a[0]} {a[2]} & {a[3]} & {int(r.src_labelled)} & {r.src_prev:.2f} & "
                    f"{int(r.target_n):,} & {r.target_prev:.3f} & {int(r.vus_assayed):,} & {r.vus_prev:.3f} \\\\")
    tot_n = int(co.target_n.sum())
    tex = r"""\begin{table}[t]
\centering
\caption{Datasets. \emph{ClinVar labelled}: missense single-nucleotide variants (SNVs) with $\ge$2 review stars and a P/LP or B/LB call (used
only to fit the common score scale) and their pathogenic share. All shares are variant-level (descriptive). \emph{Assayed}: missense SNVs with a binary assay label
that are not in the ClinVar labelled set, and their LoF share; \emph{assayed VUS}: the subset classified VUS in
ClinVar. Assay labels are the authors' published classes (authors) or, where none are deposited, derived from the
assay's own synonymous and nonsense controls (controls).}\label{tab:data}
\setlength{\tabcolsep}{4.5pt}
\begin{tabular}{llcrrrrrr}
\toprule
 & & & \multicolumn{2}{c}{ClinVar labelled} & \multicolumn{2}{c}{Assayed} & \multicolumn{2}{c}{Assayed VUS}\\
\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}
Gene & Assay & Labels & $n$ & P/LP & $n$ & LoF & $n$ & LoF\\
\midrule
""" + "\n".join(rows) + r"""
\midrule
Total & & & """ + f"{int(co.src_labelled.sum())}" + r""" & & """ + f"{tot_n:,}" + r""" & & """ + f"{int(co.vus_assayed.sum()):,}" + r""" & \\
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_data.tex").write_text(tex)

    # appendix: per gene x predictor Mondrian and best baseline, within-gene random split, alpha 0.1
    c = pd.read_csv(RES / "conformal.csv")
    c = c[(c.alpha == 0.10) & (c.design == "within-gene random split")]
    lines = []
    for (pr, g), s in c.groupby(["predictor", "gene"], sort=False):
        s = s.set_index("method")
        def cell(m):
            if m not in s.index:
                return "--"
            r = s.loc[m]
            return f"{r.coverage_lof_rw:.2f}/{r.coverage_func_rw:.2f}/{r.abstain:.2f}"
        lines.append(f"{pr} & {g} & {cell('Mondrian conformal')} & {cell('Pooled conformal')} & "
                     f"{cell('Fixed 0.1/0.9 tiers')} \\\\")
    tex = r"""\begin{table}[ht]
\centering
\caption{Within-gene random split, $\alpha=0.1$, mean of 50 splits. Each cell is residue-weighted LoF coverage /
functional coverage / abstention rate.}\label{tab:app_conformal}
\begin{tabular}{llccc}
\toprule
Predictor & Gene & Per-class (ours) & Pooled conformal & Fixed 0.1/0.9\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_conformal.tex").write_text(tex)

    # main text: AlphaMissense, within-gene random split, alpha 0.1
    a = c[c.predictor == "AlphaMissense"].set_index(["gene", "method"])
    meths = [("Mondrian conformal", "Mondrian (ours)"), ("Pooled conformal", "Pooled conformal"),
             ("Label-shift weighted, true prior", "Weighted, true prior"),
             ("AlphaMissense classes", "AlphaMissense classes")]
    lines = []
    for g in ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]:
        cells = []
        for m, _ in meths:
            r = a.loc[(g, m)]
            lo, fu = r.coverage_lof_rw, r.coverage_func_rw
            fmt = lambda v: f"\\textbf{{{v:.2f}}}" if v < 0.85 else f"{v:.2f}"
            cells.append(f"{fmt(lo)} & {fmt(fu)} & {r.abstain:.2f}")
        lines.append(f"{g} & " + " & ".join(cells) + " \\\\")
    head = " & ".join(f"\\multicolumn{{3}}{{c}}{{{lab}}}" for _, lab in meths)
    rules = "".join(f"\\cmidrule(lr){{{2 + 3 * i}-{4 + 3 * i}}}" for i in range(len(meths)))
    sub = " & ".join(["LoF & Fun & Abs"] * len(meths))
    tex = (r"""\begin{table}[t]
\centering
\caption{AlphaMissense calls in each gene ($\alpha=0.1$, calibration and test on disjoint halves of the gene's residues,
mean of 50 splits): residue-weighted coverage of LoF (LoF) and functional (Fun) variants, and abstention rate (Abs).
Coverage below 0.85 in bold. The published AlphaMissense classes need no calibration data.}\label{tab:main}
\footnotesize
\setlength{\tabcolsep}{3.2pt}
\begin{tabular}{l""" + "ccc" * len(meths) + r"""}
\toprule
 & """ + head + r"""\\
""" + rules + r"""
Gene & """ + sub + r"""\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
""")
    (OUT / "tab_main.tex").write_text(tex)

    pri = pd.read_csv(RES / "prior_estimates.csv")
    pri = pri[pri.predictor.isin(["AlphaMissense", "REVEL", "ESM-2"])]
    cal = pd.read_csv(RES / "calibration.csv")
    cal = cal[cal.predictor.isin(["AlphaMissense", "REVEL", "ESM-2"])].pivot_table(
        index=["predictor", "gene"], columns="method", values="ece")
    lines = []
    for _, r in pri.iterrows():
        e = cal.loc[(r.predictor, r.gene)]
        lines.append(f"{SHORT_PRED[r.predictor]} & {r.gene} & {r.observed:.3f} & {r.em:.3f} & {r.bbse:.3f} & "
                     f"{e['Platt on ClinVar']:.3f} & {e['+ EM prior']:.3f} & {e['+ oracle prior']:.3f} & "
                     f"{e['Platt on other-gene assays']:.3f} \\\\")
    tex = r"""\begin{table}[ht]
\centering
\caption{Prior estimation and calibration on assayed variants outside the ClinVar labelled set. LoF share observed and estimated from
unlabelled scores by expectation--maximization (EM) \cite{Saerens2002} and black-box shift estimation (BBSE) \cite{Lipton2018}; expected calibration error (10 equal-mass
bins) of the ClinVar-Platt posterior, after EM prior adjustment, after adjustment to the true prior, and of Platt
scaling fitted on the other four genes' assay labels. AM: AlphaMissense.}\label{tab:app_prior}
\setlength{\tabcolsep}{4.5pt}
\begin{tabular}{llrrrrrrr}
\toprule
 & & \multicolumn{3}{c}{LoF share} & \multicolumn{4}{c}{Expected calibration error}\\
\cmidrule(lr){3-5}\cmidrule(lr){6-9}
 & & & & & ClinVar & & + true & other\\
Pred. & Gene & obs. & EM & BBSE & Platt & + EM & prior & genes\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_prior.tex").write_text(tex)

    # appendix: sensitivity to the assay-label rule (14_label_sensitivity.py)
    ls = pd.read_csv(RES / "label_sensitivity.csv")
    ls["worst"] = ls[["lof_cov", "func_cov"]].min(axis=1)
    rule_name = {"main": None, "controls 0.80": "controls, 0.80", "controls 0.95": "controls, 0.95",
                 "controls 0.99": "controls, 0.99"}
    def rng(v):
        lo, hi = f"{v.min():.2f}", f"{v.max():.2f}"
        return lo if lo == hi else f"{lo}--{hi}"

    lines = []
    for (g, rule), s in ls.groupby(["gene", "rule"], sort=False):
        name = rule_name[rule] or ("controls, 0.95" if g == "BRCA2" else "authors")
        main = " (main)" if rule == "main" else ""
        lines.append(f"{g} & {name}{main} & {int(s.n.iloc[0]):,} & {s.lof_share.iloc[0]:.3f} & {s.worst.min():.3f} & "
                     f"{rng(s.abstain)} & {int(s.oracle_bound_ok.sum())}/{int(s.oracle_claims.sum())} & "
                     f"{rng(s.share_lof_supporting)} \\\\")
    tex = r"""\begin{table}[ht]
\centering
\caption{Sensitivity to the assay-label rule. Each gene is relabelled with the authors' classes or with the control
rule at the stated posterior cut-off; other genes are unchanged. Per-class: lowest mean per-class coverage over the
three predictors (50 random halvings, $\alpha=0.1$); abstention: range over predictors; bounds held: oracle selections
whose certified bound the true likelihood ratio respected (100 draws per predictor); LoF-supp.: share of draws
certifying LoF-supporting calls, range over predictors.}\label{tab:app_labels}
\setlength{\tabcolsep}{4.5pt}
\begin{tabular}{llrrrccc}
\toprule
Gene & Labels & $n$ & LoF & Per-class & Abstention & Bounds held & LoF-supp.\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_labels.tex").write_text(tex)

    # supplement: calls judged against ClinVar labels (15_clinvar_check.py)
    cc = pd.read_csv(RES / "clinvar_check.csv")

    def fmt(v):
        return "$\\infty$" if v == float("inf") else f"{v:.2f}" if v < 10 else f"{v:.1f}"

    lines = []
    for _, r in cc.iterrows():
        lines.append(f"{SHORT_PRED[r.predictor]} & {r.scope} & {r.n_path} & {r.n_benign} & "
                     f"{r.lof_calls_path}/{r.lof_calls_benign} & {fmt(r.lr_pos_lower)} & "
                     f"{r.func_calls_path}/{r.func_calls_benign} & {fmt(r.lr_neg_upper)} \\\\")
    tex = r"""\begin{table}[ht]
\centering
\caption{Calls judged against ClinVar labels (exploratory). Thresholds ($\alpha=0.1$) are set on each gene's assay
labels; one ClinVar variant with $\ge$2 stars is kept per residue. $n$: kept P/LP and B/LB variants; calls: LoF (or
functional) calls among P/LP / among B/LB variants; $L^+$, $U^-$: likelihood-ratio bounds from one-sided
Clopper--Pearson endpoints at level 0.025 each (simultaneous level 0.95 per ratio). Pooled rows sum the kept counts over
genes, each gene with its own thresholds. BAP1 and RAD51C have too few labelled residues to be informative.}
\label{tab:app_clinvar}
\setlength{\tabcolsep}{4.5pt}
\begin{tabular}{llrrcrcr}
\toprule
Pred. & Gene & $n$ P/LP & $n$ B/LB & LoF calls & $L^+$ & func.\ calls & $U^-$\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_clinvar.tex").write_text(tex)
    print("tables written")


if __name__ == "__main__":
    main()
