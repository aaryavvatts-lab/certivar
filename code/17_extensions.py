"""Extensions to the main experiments (E11-E13).

E11 oracle_extended.csv  Oracle design of 13_evidence.py (calibration 30% and certification 50% of residues drawn with
                         replacement, 100 draws per pair) on two grids: the main grid G4 = {0.05, 0.10, 0.20, 0.30}
                         and a finer grid G20 of 20 levels from 0.02 to 0.40. Selection rules: certified with a
                         Bonferroni correction over the grid (Algorithm 1), certified by fixed-sequence testing over
                         the grid ordered by increasing alpha (Algorithm 1, fixed-sequence variant), uncorrected
                         Clopper-Pearson, one-sided 95% percentile bootstrap and point estimate. For every selection we
                         record the true residue-weighted and variant-weighted likelihood ratios on the population, the
                         population per-class coverage of the selected thresholds, and the certification counts.
E12 prior_shift.csv      Calibration share shift: in each of 100 random halvings, the calibration half (one variant per
                         residue) is subsampled to a fixed size with LoF share 0.05, 0.10, 0.20, 0.30 or 0.50, and the
                         thresholds are applied to the other half at its natural share.
E13 data_use.csv         Per-class thresholds on one variant per residue against all variants of the calibration half
                         (50 halvings), and the calls on intermediate variants of the test half.
"""
import importlib
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ev = importlib.import_module("13_evidence")
ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
SEED, DELTA, N_ORACLE, N_SHIFT, N_USE = 2027, 0.05, 100, 100, 50
GRIDS = {"G4": list(ev.GRID), "G20": [round(float(a), 2) for a in np.linspace(0.02, 0.40, 20)]}
PIS = [0.05, 0.10, 0.20, 0.30, 0.50]


def counts(dec, y):
    return dict(tp=int(((dec == 1) & (y == 1)).sum()), fp=int(((dec == 1) & (y == 0)).sum()),
                fn=int(((dec == 0) & (y == 1)).sum()), tn=int(((dec == 0) & (y == 0)).sum()),
                n1=int((y == 1).sum()), n0=int((y == 0).sum()))


def variant_rates(dec, y):
    tp, fp = np.mean(dec[y == 1] == 1), np.mean(dec[y == 0] == 1)
    fn, tn = np.mean(dec[y == 1] == 0), np.mean(dec[y == 0] == 0)
    return sc.ratio(tp, fp), sc.ratio(fn, tn)


def fixed_sequence(cert, grid, func_descending=True):
    """Fixed-sequence testing over the grid in an order fixed in advance. For LoF level tau the hypothesis at alpha is
    LR+(alpha) <= tau; it is rejected when l_TP / u_FP > tau with both Clopper-Pearson endpoints at level delta/2,
    a level-delta test. Testing stops at the first non-rejection; among the rejected alphas the one with the largest
    yield is selected. LoF levels are tested in increasing alpha. Functional levels are treated in the same way with
    u_FN / l_TN < tau', in decreasing alpha (func_descending=True) or in increasing alpha."""
    out = {}
    for name, tau in ev.LEVELS_POS.items():
        a = sc.fixed_sequence_select(grid, [cert[x]["fs_pos"] for x in grid], [cert[x]["yield_lof"] for x in grid], tau)
        out[name] = (a, cert[a]["fs_pos"]) if a is not None else None
    order = sorted(grid, reverse=True) if func_descending else list(grid)
    for name, tau in ev.LEVELS_NEG.items():
        a = sc.fixed_sequence_select(order, [cert[x]["fs_neg"] for x in order], [cert[x]["yield_func"] for x in order],
                                     tau, lof=False)
        out[name] = (a, cert[a]["fs_neg"]) if a is not None else None
    return out


def oracle_extended(d):
    rows = []
    rng = np.random.default_rng(SEED + 20)
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        for g in ev.GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            res = np.unique(cl)
            members = {r: np.flatnonzero(cl == r) for r in res}
            for draw in range(N_ORACLE):
                na, nb = int(0.3 * len(res)), int(0.5 * len(res))
                ic = np.array([rng.choice(members[r]) for r in rng.choice(res, na)])
                iz = np.array([rng.choice(members[r]) for r in rng.choice(res, nb)])
                for gname, grid in GRIDS.items():
                    level = DELTA / (2 * len(grid))
                    thrs, cert = {}, {}
                    for a in grid:
                        thrs[a] = sc.mondrian_thresholds(p[ic], y[ic], a)
                        dz = sc.decisions(p[iz], thrs[a])
                        c = sc.certify(dz, y[iz], level)
                        u = sc.certify(dz, y[iz], DELTA / 2)
                        c["lr_pos_unc"], c["lr_neg_unc"] = u["lr_pos_lower"], u["lr_neg_upper"]
                        c["fs_pos"], c["fs_neg"] = u["lr_pos_lower"], u["lr_neg_upper"]
                        c["lr_pos_boot"], c["lr_neg_boot"] = ev.boot_bounds(dz, y[iz], rng, 5.0)
                        c["counts"] = counts(dz, y[iz])
                        cert[a] = c
                    sels = {m: ev.select(cert, m) for m in
                            ("certified", "Clopper-Pearson, uncorrected", "bootstrap", "point estimate")}
                    sels["fixed sequence"] = fixed_sequence(cert, grid)
                    sels["fixed sequence, all increasing"] = fixed_sequence(cert, grid, func_descending=False)
                    for mode, chosen in sels.items():
                        for name, best in chosen.items():
                            row = dict(predictor=pr, gene=g, draw=draw, grid=gname, selection=mode, level=name,
                                       achieved=best is not None)
                            if best is not None:
                                dec = sc.decisions(p, thrs[best[0]])
                                row.update(ev.judge(name, best, ev.weighted_rates(dec, y, cl)))
                                if mode.startswith("fixed sequence"):
                                    row["bound_ok"] = np.nan   # the fixed-sequence rule certifies the level only
                                pos = name.startswith("PS3")
                                tau = (ev.LEVELS_POS if pos else ev.LEVELS_NEG)[name]
                                vp, vn = variant_rates(dec, y)
                                row["vw_lr"] = vp if pos else vn
                                row["vw_holds"] = (vp > tau) if pos else (vn < tau)
                                rep = sc.coverage_report(dec, y, cl)
                                row["cov_lof"], row["cov_func"] = rep["coverage_lof_rw"], rep["coverage_func_rw"]
                                row.update({f"cert_{k}": v for k, v in cert[best[0]]["counts"].items()})
                            rows.append(row)
            print(pr, g, flush=True)
    return pd.DataFrame(rows)


def pooled(pc, yc, p, alpha):
    nc = np.sort(np.where(yc == 1, 1 - pc, pc))
    t = nc[min(int(np.ceil((len(nc) + 1) * (1 - alpha))), len(nc)) - 1]
    return sc.decisions(p, {0: t, 1: t})


def prior_shift(d):
    rows = []
    rng = np.random.default_rng(SEED + 21)
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        for g in ev.GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            res = np.unique(cl)
            for rep in range(N_SHIFT):
                m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                keep = np.flatnonzero(m)[sc.one_per_cluster(cl[m], rng)]
                k1, k0 = keep[y[keep] == 1], keep[y[keep] == 0]
                n = int(min(len(k1) / max(PIS), len(k0) / (1 - min(PIS))))
                te = np.flatnonzero(~m)
                pi_true = ev.residue_share(y[te], cl[te])
                pi_em = sc.em_prior(p[te], cal.prior)
                for pi in PIS:
                    n1 = int(round(pi * n))
                    ix = np.concatenate([rng.choice(k1, n1, replace=False), rng.choice(k0, n - n1, replace=False)])
                    pc, yc = p[ix], y[ix]
                    dec = {"Per-class (ours)": sc.decisions(p[te], sc.mondrian_thresholds(pc, yc, 0.10)),
                           "Pooled conformal": pooled(pc, yc, p[te], 0.10),
                           "Weighted, EM prior": sc.weighted_decisions(pc, yc, p[te], 0.10, pi_em),
                           "Weighted, true prior": sc.weighted_decisions(pc, yc, p[te], 0.10, pi_true)}
                    for meth, dv in dec.items():
                        rows.append(dict(predictor=pr, gene=g, rep=rep, cal_share=pi, n_cal=n, test_share=pi_true,
                                         method=meth, **sc.coverage_report(dv, y[te], cl[te])))
            print(pr, g, flush=True)
    return pd.DataFrame(rows)


def data_use(d):
    rows = []
    rng = np.random.default_rng(SEED + 22)
    inter_all = d[(d.assay_label == -1) & (d.src == 0) & ~d.splice_region.fillna(False).astype(bool)]
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        for g in ev.GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            it = inter_all[(inter_all.gene == g) & inter_all[col].notna()]
            p_int = cal.predict(it[col].values) if len(it) else np.array([])
            cl_int = it.cluster.values
            res = np.unique(cl)
            for rep in range(N_USE):
                m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                keep = np.flatnonzero(m)[sc.one_per_cluster(cl[m], rng)]
                allc = np.flatnonzero(m)
                for name, ix in (("one per residue", keep), ("all variants", allc)):
                    thr = sc.mondrian_thresholds(p[ix], y[ix], 0.10)
                    row = dict(predictor=pr, gene=g, rep=rep, calibration=name, n_cal=len(ix),
                               **sc.coverage_report(sc.decisions(p[~m], thr), y[~m], cl[~m]))
                    if len(it):
                        mi = ~np.isin(cl_int, np.unique(cl[m]))   # intermediate variants at test-half residues
                        di = sc.decisions(p_int[mi], thr)
                        row.update(n_int=int(mi.sum()), int_lof=float(np.mean(di == 1)) if mi.any() else np.nan,
                                   int_func=float(np.mean(di == 0)) if mi.any() else np.nan)
                    rows.append(row)
    return pd.DataFrame(rows)


def main(parts=("data_use", "prior_shift", "oracle_extended")):
    d = ev.load()
    if "data_use" in parts:
        data_use(d).to_csv(RES / "data_use.csv", index=False)
    if "prior_shift" in parts:
        prior_shift(d).to_csv(RES / "prior_shift.csv", index=False)
    if "oracle_extended" in parts:
        oracle_extended(d).to_csv(RES / "oracle_extended.csv", index=False)
    print("done")


if __name__ == "__main__":
    import sys
    main(tuple(sys.argv[1:]) or ("data_use", "prior_shift", "oracle_extended"))
