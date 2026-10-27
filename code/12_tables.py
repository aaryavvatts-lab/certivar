"""Write LaTeX tables from results/*.csv (no typed numbers)."""
import json
import pathlib

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES, OUT = ROOT / "results", ROOT / "paper" / "tables"
ASSAY = {"BRCA1": ("SGE", "HAP1", "\\cite{Findlay2018}", "investigator"),
         "BRCA2": ("SGE", "HAP1", "\\cite{Huang2025}", "controls"),
         "BAP1": ("SGE", "HAP1", "\\cite{Waters2024}", "investigator"),
         "RAD51C": ("SGE", "HAP1", "\\cite{OlveraLen2024}", "investigator"),
         "MSH2": ("DMS", "HAP1", "\\cite{Jia2021}", "investigator")}


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
\caption{Datasets. \emph{ClinVar labelled}: missense SNVs with $\ge$2 review stars and a P/LP or B/LB call
(the source set for Platt calibration) and their pathogenic share. \emph{Held-out assayed}: missense SNVs with a
binary assay label and not in the ClinVar labelled set, with their LoF share; \emph{VUS}: the subset classified VUS
in ClinVar. Labels come from the assay authors' published classes (investigator) or, where none are deposited,
from synonymous and nonsense controls of the same assay (controls).}\label{tab:data}
\small
\begin{tabular}{llcrrrrrr}
\toprule
 & & & \multicolumn{2}{c}{ClinVar labelled} & \multicolumn{2}{c}{Held-out assayed} & \multicolumn{2}{c}{Assayed VUS}\\
\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}
Gene & Assay & Labels & $n$ & share P & $n$ & share LoF & $n$ & share LoF\\
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
    tex = r"""\begin{table}[h]
\centering
\caption{Within-gene random split, $\alpha=0.1$, mean of 50 splits. Each cell is residue-weighted LoF coverage /
functional coverage / abstention rate.}\label{tab:app_conformal}
\scriptsize
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
        lines.append(f"{r.predictor} & {r.gene} & {r.observed:.3f} & {r.em:.3f} & {r.bbse:.3f} & "
                     f"{e['Platt on ClinVar']:.3f} & {e['+ EM prior']:.3f} & {e['+ oracle prior']:.3f} & "
                     f"{e['Platt on other-gene assays']:.3f} \\\\")
    tex = r"""\begin{table}[h]
\centering
\caption{Prior estimation and calibration on held-out assayed variants. LoF share observed and estimated from
unlabeled scores by EM \cite{Saerens2002} and BBSE \cite{Lipton2018}; expected calibration error (10 equal-mass
bins) of the ClinVar-Platt posterior, after EM prior adjustment, after adjustment to the true prior, and of Platt
scaling fitted on the other four genes' assay labels.}\label{tab:app_prior}
\scriptsize
\begin{tabular}{llrrrrrrr}
\toprule
 & & \multicolumn{3}{c}{LoF share} & \multicolumn{4}{c}{Expected calibration error}\\
\cmidrule(lr){3-5}\cmidrule(lr){6-9}
Predictor & Gene & observed & EM & BBSE & ClinVar Platt & + EM & + true prior & other-gene assays\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_prior.tex").write_text(tex)
    print("tables written")


if __name__ == "__main__":
    main()
