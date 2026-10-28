"""Main experiments on the assembled missense table (data/processed/variants.csv).

Source labels: ClinVar missense SNVs with >= 2 review stars and a P/LP or B/LB call.
Target labels: assay-derived LoF (1) / functional (0) from 06_assemble.py. Every
variant in the source set is removed from every target set (leakage-free).

E1 cohort.csv            per-gene counts and LoF prevalence
E2 discrimination.csv    AUROC/AUPRC per predictor x gene x subset, cluster-bootstrap CIs
   dual_labelled.csv     ClinVar vs assay label agreement on variants carrying both
E3 calibration.csv       ECE/Brier per predictor x gene: source Platt, EM, BBSE, oracle,
                         assay-transfer (leave-one-gene-out)
   prior_estimates.csv   estimated vs observed LoF fraction per gene and subset
E4 conformal.csv         coverage/abstention per method x gene x alpha under three designs:
                         within-gene random split, within-gene region split, leave-one-gene-out
   conformal_strata.csv  the same by difficulty tertile and for VUS (random split, alpha 0.1)
   calibration_size.csv  Mondrian coverage and abstention against number of calibration residues
E5 application.csv       decisions for unassayed VUS (no ground truth; reported as counts)
"""
import pathlib

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import certivar as sc

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC, RES = ROOT / "data" / "processed", ROOT / "results"
GENES = ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]
PRED = {"AlphaMissense": "alphamissense", "REVEL": "revel", "ESM-2": "esm_neg",
        "CADD": "cadd", "SIFT": "sift_d", "PolyPhen-2": "polyphen"}
CONF_PRED = ["AlphaMissense", "REVEL", "ESM-2"]
ALPHAS = [0.05, 0.10, 0.20]
N_BOOT, SEED = int(__import__("os").environ.get("N_BOOT", 1000)), 2027
# Published fixed cut-offs used as baselines (verified sources listed in SOURCES.md)
AM_CUTS = (0.34, 0.564, False)   # AlphaMissense: likely benign < 0.34, likely pathogenic > 0.564
REVEL_CUTS = (0.290, 0.644, True) # Pejaver et al. 2022 Table 2: BP4 <= 0.290, PP3 >= 0.644


def load():
    d = pd.read_csv(PROC / "variants.csv", low_memory=False)
    d["esm_neg"] = -d.esm_llr
    d["src"] = ((d.stars >= 2) & d.group.isin(["path", "benign"])).astype(int)
    d["y_src"] = (d.group == "path").astype(int)
    return d


def auc_ci(df, col, ycol):
    s = df[[col, ycol, "cluster"]].dropna()
    if s[ycol].nunique() < 2 or len(s) < 20:
        return dict(auroc=np.nan, lo=np.nan, hi=np.nan, auprc=np.nan, n=len(s), prev=np.nan)
    y, x, cl = s[ycol].values.astype(int), s[col].values, s.cluster.values

    def stat(i):
        return roc_auc_score(y[i], x[i]) if len(np.unique(y[i])) == 2 else np.nan

    lo, hi = sc.cluster_bootstrap(stat, cl, N_BOOT, SEED)
    return dict(auroc=roc_auc_score(y, x), lo=lo, hi=hi, auprc=average_precision_score(y, x),
                n=len(s), prev=y.mean())


def e1_cohort(d):
    rows = []
    for g in GENES:
        s = d[d.gene == g]
        t = s[(s.assay_label.isin([0, 1])) & (s.src == 0)]
        rows.append(dict(gene=g, missense=len(s), clinvar_PLP=(s.group == "path").sum(),
                         clinvar_BLB=(s.group == "benign").sum(), clinvar_VUS=(s.group == "vus").sum(),
                         src_labelled=s.src.sum(), src_prev=s[s.src == 1].y_src.mean(),
                         assayed=s.score.notna().sum(), assay_LoF=(s.assay_label == 1).sum(),
                         assay_func=(s.assay_label == 0).sum(), assay_int=(s.assay_label == -1).sum(),
                         target_n=len(t), target_prev=t.assay_label.mean(),
                         vus_assayed=len(t[t.group == "vus"]),
                         vus_prev=t[t.group == "vus"].assay_label.mean()))
    return pd.DataFrame(rows)


def e2_discrimination(d):
    rows = []
    src = d[d.src == 1]
    for name, col in PRED.items():
        rows.append(dict(predictor=name, gene="ALL", subset="ClinVar>=2*", **auc_ci(src, col, "y_src")))
        for g in GENES + ["ALL"]:
            t = d[(d.assay_label.isin([0, 1])) & (d.src == 0)]
            if g != "ALL":
                t = t[t.gene == g]
            for subset, df in [("assay", t), ("assay VUS", t[t.group == "vus"])]:
                rows.append(dict(predictor=name, gene=g, subset=subset, **auc_ci(df, col, "assay_label")))
    dual = d[(d.src == 1) & d.assay_label.isin([0, 1])]
    agree = dual.groupby("gene").apply(lambda x: pd.Series(
        dict(n=len(x), agree=(x.y_src == x.assay_label).mean()))).reset_index()
    return pd.DataFrame(rows), agree


def platt_on(df, col, ycol):
    s = df[[col, ycol]].dropna()
    return sc.Platt().fit(s[col].values, s[ycol].values.astype(int))


def e3_calibration(d):
    rows, priors = [], []
    src = d[d.src == 1]
    tgt_all = d[(d.assay_label.isin([0, 1])) & (d.src == 0)]
    for name, col in PRED.items():
        cal = platt_on(src, col, "y_src")
        q_src = cal.predict(src[col].dropna().values)
        y_src_ = src.dropna(subset=[col]).y_src.values
        for g in GENES:
            t = tgt_all[(tgt_all.gene == g)].dropna(subset=[col])
            if len(t) < 30 or t.assay_label.nunique() < 2:
                continue
            y, cl = t.assay_label.values.astype(int), t.cluster.values
            q = cal.predict(t[col].values)
            pi_em = sc.em_prior(q, cal.prior)
            pi_bbse = sc.bbse_prior(q_src, y_src_, q)
            other = tgt_all[(tgt_all.gene != g)].dropna(subset=[col])
            xfer = platt_on(other, col, "assay_label")
            q_x = xfer.predict(t[col].values)
            pi_x = sc.em_prior(q_x, xfer.prior)
            probs = {"Platt on ClinVar": q,
                     "+ EM prior": sc.adjust(q, cal.prior, pi_em),
                     "+ BBSE prior": sc.adjust(q, cal.prior, pi_bbse) if np.isfinite(pi_bbse) else None,
                     "+ oracle prior": sc.adjust(q, cal.prior, y.mean()),
                     "Platt on other-gene assays": q_x,
                     "Other-gene assays + EM prior": sc.adjust(q_x, xfer.prior, pi_x)}
            base = probs["Platt on ClinVar"]
            for m, p in probs.items():
                if p is None:
                    continue
                lo, hi = sc.cluster_bootstrap(
                    lambda i: sc.ece(base[i], y[i]) - sc.ece(p[i], y[i]), cl, N_BOOT, SEED)
                rows.append(dict(predictor=name, gene=g, method=m, n=len(y), prev=y.mean(),
                                 ece=sc.ece(p, y), brier=sc.brier(p, y),
                                 ece_gain=sc.ece(base, y) - sc.ece(p, y), gain_lo=lo, gain_hi=hi))
            v = t[t.group == "vus"]
            priors.append(dict(predictor=name, gene=g, src_prior=cal.prior, observed=y.mean(),
                               em=pi_em, bbse=pi_bbse, xfer_em=pi_x,
                               vus_observed=v.assay_label.mean() if len(v) else np.nan,
                               vus_em=sc.em_prior(cal.predict(v[col].values), cal.prior) if len(v) > 10 else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(priors)


def fixed_cut_decisions(x, cuts):
    lo, hi, inclusive = cuts
    path = x >= hi if inclusive else x > hi
    ben = x <= lo if inclusive else x < lo
    return np.where(path, 1, np.where(ben, 0, -1))


def conformal_methods(p_cal, y_cal, cl_cal, p_tst, raw_tst, name, alpha, rng, pi_t, pi_oracle=None,
                      balance=False):
    keep = sc.one_per_cluster(cl_cal, rng)
    if balance:  # shift the calibration prior to 50% LoF by subsampling the majority class
        k1, k0 = keep[y_cal[keep] == 1], keep[y_cal[keep] == 0]
        n = min(len(k1), len(k0))
        keep = np.sort(np.concatenate([rng.choice(k1, n, replace=False), rng.choice(k0, n, replace=False)]))
    pc, yc = p_cal[keep], y_cal[keep]
    out = {"Mondrian conformal": sc.decisions(p_tst, sc.mondrian_thresholds(pc, yc, alpha))}
    out["Label-shift weighted conformal"] = sc.weighted_decisions(pc, yc, p_tst, alpha, pi_t)
    if pi_oracle is not None:
        out["Label-shift weighted, true prior"] = sc.weighted_decisions(pc, yc, p_tst, alpha, pi_oracle)
    nc = np.sort(np.where(yc == 1, 1 - pc, pc))
    t = nc[min(int(np.ceil((len(nc) + 1) * (1 - alpha))), len(nc)) - 1]
    out["Pooled conformal"] = sc.decisions(p_tst, {0: t, 1: t})
    out["Fixed 0.1/0.9 tiers"] = np.where(p_tst >= 0.9, 1, np.where(p_tst <= 0.1, 0, -1))
    if name == "AlphaMissense":
        out["AlphaMissense classes"] = fixed_cut_decisions(raw_tst, AM_CUTS)
    if name == "REVEL":
        out["REVEL PP3/BP4 supporting"] = fixed_cut_decisions(raw_tst, REVEL_CUTS)
    return out


def posteriors(d, col):
    """ClinVar-Platt posterior q for every target variant, plus the per-gene EM prior estimate of q.

    The Mondrian rule uses q directly: per-class thresholds are invariant to any per-gene monotone
    re-weighting, so it needs no prior. The EM prior is used only by the label-shift weighted baseline.
    """
    src = d[d.src == 1]
    cal = platt_on(src, col, "y_src")
    q, pis = pd.Series(np.nan, index=d.index), {}
    for g in GENES:
        m = (d.gene == g) & d[col].notna() & d.assay_label.isin([0, 1]) & (d.src == 0)
        if m.sum() < 30:
            continue
        q[m] = cal.predict(d.loc[m, col].values)
        pis[g] = sc.em_prior(q[m].values, cal.prior)
    return q, pis


def region_blocks(te, n_blocks=4):
    """Assign each residue to one of n contiguous blocks of the assayed protein positions."""
    pos = te.protein_pos.values
    edges = np.quantile(np.unique(pos), np.linspace(0, 1, n_blocks + 1))
    return np.clip(np.searchsorted(edges, pos, side="right") - 1, 0, n_blocks - 1)


def report(rows, base, dec, y, cl):
    for m, dv in dec.items():
        rows.append(dict(base, method=m, **sc.coverage_report(dv, y, cl)))


def e4_conformal(d):
    rows, strata, curve = [], [], []
    rng = np.random.default_rng(SEED)
    tgt = d[(d.assay_label.isin([0, 1])) & (d.src == 0)].copy()
    tgt["tertile"] = tgt.groupby("gene").assay_margin.transform(
        lambda x: pd.qcut(x.rank(method="first"), 3, labels=["hard", "medium", "easy"]))
    for name in CONF_PRED:
        col = PRED[name]
        q, pis = posteriors(d, col)
        tgt["p"] = q.reindex(tgt.index)
        t = tgt.dropna(subset=["p"])
        for alpha in ALPHAS:
            for g in GENES:
                te, ca = t[t.gene == g], t[t.gene != g]
                if len(te) < 30:
                    continue
                y = te.assay_label.values.astype(int)
                base = dict(predictor=name, alpha=alpha, gene=g)
                # (a) leave-one-gene-out: calibrate on the other four genes
                dec = conformal_methods(ca.p.values, ca.assay_label.values.astype(int), ca.cluster.values,
                                        te.p.values, te[col].values, name, alpha, rng, pis[g], y.mean())
                report(rows, dict(base, design="leave-one-gene-out"), dec, y, te.cluster.values)
                # (b) within-gene random split of residues in half, 50 repeats
                pos_ids = te.cluster.unique()
                acc = {}
                for rep in range(50):
                    cal_pos = set(rng.choice(pos_ids, len(pos_ids) // 2, replace=False))
                    mc = te.cluster.isin(cal_pos).values
                    dec = conformal_methods(te.p.values[mc], y[mc], te.cluster.values[mc], te.p.values[~mc],
                                            te[col].values[~mc], name, alpha, rng, pis[g], y[~mc].mean())
                    for m, dv in dec.items():
                        acc.setdefault(m, []).append(sc.coverage_report(dv, y[~mc], te.cluster.values[~mc]))
                        if alpha == 0.10 and rep == 0:
                            sub = te[~mc]
                            for key, grp in (("difficulty", "tertile"), ("clinvar", "group")):
                                for lvl, idx in sub.groupby(grp, observed=True).groups.items():
                                    if key == "clinvar" and lvl != "vus":
                                        continue
                                    pos = sub.index.get_indexer(idx)
                                    strata.append(dict(base, method=m, stratum=f"{key}:{lvl}",
                                                       **sc.coverage_report(dv[pos], y[~mc][pos])))
                for m, reps in acc.items():
                    r = pd.DataFrame(reps)
                    sd = {f"{k}_sd": v for k, v in r[["coverage_lof_rw", "coverage_func_rw", "abstain"]].std().items()}
                    rows.append(dict(base, design="within-gene random split", method=m, n_rep=len(r),
                                     **r.mean(numeric_only=True).to_dict(), **sd))
                # (b2) same split, calibration prior shifted to 50% LoF (ClinVar-like), test at natural prior
                acc = {}
                for rep in range(50):
                    cal_pos = set(rng.choice(pos_ids, len(pos_ids) // 2, replace=False))
                    mc = te.cluster.isin(cal_pos).values
                    dec = conformal_methods(te.p.values[mc], y[mc], te.cluster.values[mc], te.p.values[~mc],
                                            te[col].values[~mc], name, alpha, rng, pis[g], y[~mc].mean(),
                                            balance=True)
                    for m, dv in dec.items():
                        acc.setdefault(m, []).append(sc.coverage_report(dv, y[~mc], te.cluster.values[~mc]))
                for m, reps in acc.items():
                    rows.append(dict(base, design="within-gene split, calibration at 50% LoF", method=m,
                                     **pd.DataFrame(reps).mean(numeric_only=True).to_dict()))
                # (c) within-gene region split: leave one contiguous quarter of the protein out
                blk = region_blocks(te)
                acc = {}
                for b in np.unique(blk):
                    mc = blk != b
                    dec = conformal_methods(te.p.values[mc], y[mc], te.cluster.values[mc], te.p.values[~mc],
                                            te[col].values[~mc], name, alpha, rng, pis[g], y[~mc].mean())
                    for m, dv in dec.items():
                        acc.setdefault(m, [[], [], []])
                        acc[m][0].append(dv)
                        acc[m][1].append(y[~mc])
                        acc[m][2].append(te.cluster.values[~mc])
                for m, (dvs, ys, cls) in acc.items():
                    rows.append(dict(base, design="within-gene region split", method=m,
                                     **sc.coverage_report(np.concatenate(dvs), np.concatenate(ys), np.concatenate(cls))))
                # (d) how many calibration residues are needed (Mondrian only, alpha = 0.1)
                if alpha == 0.10:
                    for n_pos in (25, 50, 100, 150, 200):  # sizes available in every gene
                        for rep in range(50):
                            cal_pos = set(rng.choice(pos_ids, n_pos, replace=False))
                            mc = te.cluster.isin(cal_pos).values
                            if len(np.unique(y[mc])) < 2:
                                continue
                            keep = sc.one_per_cluster(te.cluster.values[mc], rng)
                            thr = sc.mondrian_thresholds(te.p.values[mc][keep], y[mc][keep], alpha)
                            curve.append(dict(base, n_cal_residues=n_pos, rep=rep,
                                              **sc.coverage_report(sc.decisions(te.p.values[~mc], thr), y[~mc],
                                                                  te.cluster.values[~mc])))
    return pd.DataFrame(rows), pd.DataFrame(strata), pd.DataFrame(curve)


def e5_application(d):
    """Calibrate the Mondrian rule on a gene's assayed variants; apply it to that gene's unassayed VUS.

    These VUS lie outside the assayed region, so the guarantee does not formally apply (region shift);
    the counts are descriptive only.
    """
    rows = []
    col = PRED["AlphaMissense"]
    src = d[d.src == 1]
    cal = platt_on(src, col, "y_src")
    for g in GENES:
        lab = d[(d.gene == g) & d.assay_label.isin([0, 1]) & (d.src == 0) & d[col].notna()]
        un = d[(d.gene == g) & (d.group == "vus") & d.score.isna() & d[col].notna()]
        if len(lab) < 30 or len(un) == 0:
            continue
        conf = sc.ShiftConformal(alpha=0.10, seed=SEED).fit(
            cal.predict(lab[col].values), lab.assay_label.values.astype(int), lab.cluster.values)
        dec = conf.predict(cal.predict(un[col].values))
        rows.append(dict(gene=g, unassayed_vus=len(un), call_lof=int((dec == 1).sum()),
                         call_func=int((dec == 0).sum()), abstain=int((dec == -1).sum())))
    return pd.DataFrame(rows)


def main():
    RES.mkdir(exist_ok=True)
    d = load()
    e1_cohort(d).to_csv(RES / "cohort.csv", index=False)
    disc, agree = e2_discrimination(d)
    disc.to_csv(RES / "discrimination.csv", index=False)
    agree.to_csv(RES / "dual_labelled.csv", index=False)
    cal, pri = e3_calibration(d)
    cal.to_csv(RES / "calibration.csv", index=False)
    pri.to_csv(RES / "prior_estimates.csv", index=False)
    conf, strata, curve = e4_conformal(d)
    conf.to_csv(RES / "conformal.csv", index=False)
    strata.to_csv(RES / "conformal_strata.csv", index=False)
    curve.to_csv(RES / "calibration_size.csv", index=False)
    e5_application(d).to_csv(RES / "application.csv", index=False)
    print("done")


if __name__ == "__main__":
    main()
