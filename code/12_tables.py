"""Write LaTeX tables from results/*.csv (no typed numbers). All tables are set at the 10 pt body size."""
import json
import pathlib

import numpy as np
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
    num = json.loads((RES / "paper_numbers.json").read_text())
    rows = []
    for _, r in co.iterrows():
        a = ASSAY[r.gene]
        rr = num["residues"][r.gene]
        rows.append(f"{r.gene} & {a[0]} {a[2]} & {a[3]} & {int(r.src_labelled)} & {int(round(r.src_prev * r.src_labelled))} & "
                    f"{int(r.target_n):,} & {rr['residues']:,} & {r.target_prev:.3f} & {int(r.vus_assayed):,} & "
                    f"{r.vus_prev:.3f} \\\\")
    tot_n = int(co.target_n.sum())
    tex = r"""\begin{table}[t]
\centering
\caption{Datasets. \emph{ClinVar labelled}: missense variants with $\ge$2 review stars and a P/LP or B/LB call (used
for the common score scale and an exploratory check). \emph{Assayed}: missense variants with a binary assay label outside that set;
\emph{assayed VUS}: the subset classified VUS. LoF shares are variant-level. Labels: authors' classes or a rule fitted to
the assay's controls (\cref{sec:data}).}\label{tab:data}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llcrrrrrrr}
\toprule
 & & & \multicolumn{2}{c}{ClinVar labelled} & \multicolumn{3}{c}{Assayed} & \multicolumn{2}{c}{Assayed VUS}\\
\cmidrule(lr){4-5}\cmidrule(lr){6-8}\cmidrule(lr){9-10}
Gene & Assay & Labels & $n$ & P/LP & $n$ & Residues & LoF & $n$ & LoF\\
\midrule
""" + "\n".join(rows) + r"""
\midrule
Total & & & """ + f"{int(co.src_labelled.sum())}" + r""" & """ + f"{int(round((co.src_prev * co.src_labelled).sum()))}" + r""" & """ + f"{tot_n:,}" + r""" & """ + f"{sum(v['residues'] for v in num['residues'].values()):,}" + r""" & """ + f"{(co.target_prev * co.target_n).sum() / tot_n:.3f}" + r""" & """ + f"{int(co.vus_assayed.sum()):,}" + r""" & """ + f"{(co.vus_prev * co.vus_assayed).sum() / co.vus_assayed.sum():.3f}" + r""" \\
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
            return f"{r.coverage_lof_rw:.3f}/{r.coverage_func_rw:.3f}/{r.abstain:.2f}"
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

    # main text: every gene and predictor, within-gene random split, alpha 0.1
    pub = {"AlphaMissense": "AlphaMissense classes", "REVEL": "REVEL PP3/BP4 supporting"}
    fmt = lambda v: f"\\textbf{{{v:.3f}}}" if v < 0.85 else f"{v:.3f}"
    lines = []
    for g in ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]:
        for pr in ["AlphaMissense", "REVEL", "ESM-2"]:
            s = c[(c.gene == g) & (c.predictor == pr)].set_index("method")
            r = s.loc["Mondrian conformal"]
            cells = [f"{r.coverage_lof_rw:.3f}", f"{r.coverage_func_rw:.3f}", f"{r.abstain:.2f}"]
            q = s.loc["Pooled conformal"]
            cells += [fmt(min(q.coverage_lof_rw, q.coverage_func_rw)), f"{q.abstain:.2f}"]
            if pub.get(pr) is None:
                cells += ["--", "--", "--"]
            else:
                q = s.loc[pub[pr]]
                cells += [fmt(q.coverage_lof_rw), fmt(q.coverage_func_rw), f"{q.abstain:.2f}"]
            lines.append(f"{g if pr == 'AlphaMissense' else ''} & {SHORT_PRED[pr]} & " + " & ".join(cells) + " \\\\")
        if g != "MSH2":
            lines.append("\\addlinespace[1pt]")
    m = c[c.method == "Mondrian conformal"]
    se = pd.concat([m.coverage_lof_rw_sd, m.coverage_func_rw_sd]) / m.n_rep.iloc[0] ** 0.5
    tex = (r"""\begin{table}[t]
\centering
\caption{Calls by gene and predictor ($\alpha=0.1$; calibration and test on disjoint halves of the gene's residues, mean
of 50 splits). Per-class rule: residue-weighted coverage of LoF and functional (Fun) variants and abstention rate (Abs);
Monte Carlo standard errors of the coverages are """ + f"{se.min():.3f}--{se.max():.3f}" + r""". Pooled conformal:
worst-class coverage and abstention. Published cut-offs (AlphaMissense classes for AM, REVEL PP3/BP4 thresholds for
REVEL; none exist for ESM-2): coverage of each class and abstention. Bold: below 0.85.}\label{tab:main}
\setlength{\tabcolsep}{3.6pt}
\begin{tabular}{llcccccccc}
\toprule
 & & \multicolumn{3}{c}{Per-class (ours)} & \multicolumn{2}{c}{Pooled conformal} & \multicolumn{3}{c}{Published cut-off}\\
\cmidrule(lr){3-5}\cmidrule(lr){6-7}\cmidrule(lr){8-10}
Gene & Predictor & LoF & Fun & Abs & Worst & Abs & LoF & Fun & Abs\\
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
certifying LoF-supporting calls, range over predictors. This is an independent run with its own random splits, so its
main-label rows differ from \cref{tab:app_conformal,tab:main} by Monte Carlo error.}\label{tab:app_labels}
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
    revision_tables()
    print("tables written")


def r1(v):
    """One decimal, halves rounded up (19.05 -> 19.1)."""
    return f"{np.floor(v * 10 + 0.5 + 1e-9) / 10:.1f}"


def fmt3(v):
    return "--" if v is None else f"{v:.3f}"


def pct(v):
    return "--" if v is None or v != v else f"{100 * v:.1f}"


def revision_tables():
    """Tables 2-4 of the main text and the use-case table of the appendix."""
    num = json.loads((RES / "paper_numbers.json").read_text())
    # Table 2: per-class coverage of every method, within-gene random split, alpha 0.1, and at calibration share 0.05
    c = pd.read_csv(RES / "conformal.csv")
    c = c[(c.alpha == 0.10) & (c.design == "within-gene random split")].copy()
    c["worst"] = c[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    ps = pd.read_csv(RES / "prior_shift.csv")
    ps = ps[ps.cal_share == 0.05].groupby(["method", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw"]].mean()
    ps["worst"] = ps[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    meths = [("Mondrian conformal", "Per-class (ours)", "Per-class (ours)"),
             ("Pooled conformal", "Pooled conformal", "Pooled conformal"),
             ("Label-shift weighted, true prior", "Label-shift weighted, true share", "Weighted, true prior"),
             ("Label-shift weighted conformal", "Label-shift weighted, EM share", "Weighted, EM prior"),
             ("Fixed 0.1/0.9 tiers", "Fixed tiers ($q\\ge0.9$, $q\\le0.1$)", None),
             ("AlphaMissense classes", "AlphaMissense classes", None),
             ("REVEL PP3/BP4 supporting", "REVEL PP3/BP4 thresholds", None)]
    lines = []
    for key, lab, sh in meths:
        v = c[c.method == key]
        shift = f"{ps.loc[sh].worst.min():.3f}" if sh else "--"
        lines.append(f"{lab} & {len(v)} & {v.worst.mean():.3f} & {v.worst.min():.3f} & {int((v.worst < 0.85).sum())} & "
                     f"{v.abstain.mean():.2f} & {shift} \\\\")
    tex = r"""\begin{table}[t]
\centering
\caption{Per-class coverage at $\alpha=0.1$ (residue-weighted; 50 random halvings per gene). Last column: worst pair
when the calibration half is subsampled to a LoF share of 0.05 (100 halvings); fixed cut-offs use no calibration data.}\label{tab:methods}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lcccccc}
\toprule
 & & \multicolumn{3}{c}{Worst-class coverage} & & Worst pair,\\
\cmidrule(lr){3-5}
Method & Pairs & mean & worst pair & pairs $<0.85$ & Abstention & share 0.05\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_methods.tex").write_text(tex)

    # Table 3: oracle benchmark, grid G4, every rule on the same 1,500 draws
    o = num["oracle_g4"]
    rules = [("certified", "Certified (ours), Bonferroni", "yes, all levels"),
             ("fixed sequence, simultaneous", "Certified (ours), fixed sequence", "yes, all levels"),
             ("Clopper-Pearson, uncorrected", "Clopper--Pearson, uncorrected", "no"),
             ("Katz, Bonferroni", "Katz interval, Bonferroni", "no"),
             ("bootstrap", "Percentile bootstrap", "no"),
             ("bootstrap, pseudocount", "Percentile bootstrap, $+0.5$ counts", "no"),
             ("point estimate", "Point estimate", "no")]
    lines = []
    for key, lab, guar in rules:
        r = o[key]
        a = r["achieved"]
        bounds = "--" if key == "point estimate" else pct(r["draws_all_bounds_ok"])   # a point estimate is no bound
        lines.append(f"{lab} & {guar} & {bounds} & {pct(r['draws_all_levels_held'])} & "
                     f"{pct(a['PS3-supporting'])} & {pct(a['PS3-moderate'])} & {pct(a['BS3-supporting'])} & "
                     f"{pct(a['BS3-moderate'])} \\\\")
    tex = r"""\begin{table}[t]
\centering
\caption{Selection rules against the population likelihood ratio (oracle design, $G_4$, $\delta=0.05$, the same 1,500
draws). Guarantee: finite-sample validity after selection for all levels at once. Bounds and levels held: draws in which
every claimed bound (fixed sequence certifies levels only; a point estimate is no bound) or level held (\%). Last four columns: draws certifying each
level (\%).}
\label{tab:oracle}
\setlength{\tabcolsep}{3.4pt}
\begin{tabular}{llcccccc}
\toprule
 & Finite-sample & Bounds & Levels & \multicolumn{2}{c}{LoF calls} & \multicolumn{2}{c}{Functional calls}\\
\cmidrule(lr){5-6}\cmidrule(lr){7-8}
Selection rule & guarantee & held & held & supp. & mod. & supp. & mod.\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_oracle.tex").write_text(tex)

    # Table 4: mixture-consistency check of Assumption 2
    mc = num["mixture_check"]
    sets = [("random half", "Random half of the residues", "holds by design"),
            ("assayed VUS", "Assayed VUS at other residues", "approximately"),
            ("held-out quarter", "Held-out contiguous quarter", "not assumed"),
            ("other gene", "Another gene's calibration data", "not assumed"),
            ("unassayed VUS, other residues", "Unassayed VUS, other residues", "unknown (no labels)")]
    lines = []
    for key, lab, status in sets:
        r = mc[key]
        lines.append(f"{lab} & {status} & {r['n']:,} & {pct(r['reject'])} & "
                     f"{fmt3(r.get('worst_mean_accepted'))} & {fmt3(r.get('worst_mean_rejected'))} \\\\")
    tex = r"""\begin{table}[t]
\centering
\caption{Mixture-consistency check ($p<0.05$; one target variant per residue). Coverage: mean over target sets of the
smaller per-class coverage in each set (a per-set measure, unlike \cref{tab:methods}) of in-expectation thresholds
($\alpha=0.1$) set on the same calibration data.}\label{tab:mixture}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{llcccc}
\toprule
 & & Target & Rejected & \multicolumn{2}{c}{Worst-class coverage}\\
\cmidrule(lr){5-6}
Target set & Assumption~\ref{as:shift} & sets & (\%) & accepted & rejected\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_mixture.tex").write_text(tex)

    # Appendix: ClinVar VUS by assay status, and certified calls with mixture checks on unassayed VUS at assayed residues
    uc = pd.read_csv(RES / "use_case.csv")
    rej = mc["unassayed VUS, other residues"]["by_gene"]
    genes, preds = ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"], ["AlphaMissense", "REVEL", "ESM-2"]
    lines = []
    for g in genes:
        f = uc[uc.gene == g].iloc[0]
        lines.append(f"{g} & {int(f.vus):,} & {int(f.vus_assayed):,} & {int(f.vus_unassayed_at_assayed_residues)} & "
                     f"{int(f.vus_unassayed_elsewhere):,} & {pct(rej.get(g))} \\\\")
    t = uc.groupby("gene").first()
    lines += [r"\midrule", f"Total & {int(t.vus.sum()):,} & {int(t.vus_assayed.sum()):,} & "
              f"{int(t.vus_unassayed_at_assayed_residues.sum())} & {int(t.vus_unassayed_elsewhere.sum()):,} & \\\\"]
    tex = r"""\begin{table}[ht]
\centering
\caption{ClinVar missense VUS by assay status. Assayed: VUS with an assay score, including intermediate and
splice-region variants (\cref{tab:data} counts only those with a binary label). Check: share of 20 mixture checks
rejected ($p<0.05$) for the unassayed VUS at other residues, run where they span at least 30 residues.}
\label{tab:app_usecase}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrc}
\toprule
 & & & \multicolumn{2}{c}{Unassayed} & Check rejected,\\
\cmidrule(lr){4-5}
Gene & All & Assayed & at assayed residues & elsewhere & elsewhere (\%)\\
\midrule
""" + "\n".join(lines) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_usecase.tex").write_text(tex)
    lines = []
    for g in genes:
        n = int(uc[uc.gene == g].n_target.iloc[0])
        if n == 0:
            continue
        for i, pr in enumerate(preds):
            v = uc[(uc.gene == g) & (uc.predictor == pr)]
            checked = v.check_pval.notna().any()
            chk = pct((v.check_pval < 0.05).mean()) if checked else "--"
            acc = v[v.check_pval >= 0.05] if checked else v      # calls are reported only where the check accepts
            m = acc[["PS3-supporting", "PS3-moderate", "BS3-supporting", "BS3-moderate"]].mean()
            lines.append(f"{g if i == 0 else ''} & {n if i == 0 else ''} & {pr} & {r1(m['PS3-supporting'])} & "
                         f"{r1(m['PS3-moderate'])} & {r1(m['BS3-supporting'])} & {r1(m['BS3-moderate'])} & {chk} \\\\")
        lines.append(r"\addlinespace[1pt]")
    tex = r"""\begin{table}[ht]
\centering
\caption{Certified calls on the ClinVar VUS that are unassayed but lie at assayed residues (\cref{alg:certivar}, $G_4$,
$\delta=0.05$; mean number of calls over 20 random 40/60 calibration/certification splits of all assayed residues), and
the share of mixture checks of those VUS rejected ($p<0.05$), run where they span at least 30 residues; where the check
was run, calls are averaged over the splits it accepted.}
\label{tab:app_usecase_calls}
\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrlccccc}
\toprule
 & & & \multicolumn{2}{c}{LoF calls} & \multicolumn{2}{c}{Functional calls} & Check\\
\cmidrule(lr){4-5}\cmidrule(lr){6-7}
Gene & VUS & Predictor & supp. & mod. & supp. & mod. & rejected (\%)\\
\midrule
""" + "\n".join(lines[:-1]) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    (OUT / "tab_app_usecase_calls.tex").write_text(tex)


if __name__ == "__main__":
    main()
