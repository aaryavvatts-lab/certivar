"""Certified evidence, oracle validation, per-split reliability, residue-disjoint VUS shift and region split.

E6 evidence.csv         100 random three-way splits of each gene's residues (calibration 30%, certification 50%,
                        test 20%). For each OddsPath level, select the miscoverage level on a grid by one of five
                        rules applied to the same certification data: certified (Clopper-Pearson, Bonferroni over
                        the grid; Algorithm 1), Clopper-Pearson without correction, a one-sided 95% percentile
                        bootstrap (2,000 resamples), the same bootstrap at the Bonferroni level, and the point
                        estimate. Record the test-split likelihood ratio.
E6b oracle.csv          The same selection rules with calibration and certification residues drawn with
                        replacement from the gene's assayed residues, treated as the population, so that the true
                        residue-weighted likelihood ratio of every selected rule is known exactly.
E7 reliability.csv      Per-class coverage in 200 random halvings: in-expectation vs PAC thresholds.
E8 vus_shift.csv        Residues split into two disjoint halves; calibrate and certify on non-VUS variants of one
                        half, call the ClinVar VUS of the other half; 50 repetitions.
E9 region_evidence.csv  Calibrate and certify on three contiguous quarters of a protein, test on the fourth.
"""
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC, RES = ROOT / "data" / "processed", ROOT / "results"
GENES = ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]
PRED = {"AlphaMissense": "alphamissense", "REVEL": "revel", "ESM-2": "esm_neg"}
GRID = [0.05, 0.10, 0.20, 0.30]
# OddsPath levels (Tavtigian 2018, Brnich 2020). The keys are internal names: the paper reports them as LoF- and
# functional-supporting/moderate/strong, because the likelihood ratios are defined against the assay outcome.
LEVELS_POS = {"PS3-supporting": 2.1, "PS3-moderate": 4.3, "PS3-strong": 18.7}       # LR+ must exceed
LEVELS_NEG = {"BS3-supporting": 0.48, "BS3-moderate": 0.23, "BS3-strong": 0.053}   # LR- must be below
DELTA, SEED = 0.05, 2027
N_SPLIT, N_ORACLE, N_REL, N_SHIFT, N_BOOT = 100, 100, 200, 50, 2000
MODES = ["certified", "Clopper-Pearson, uncorrected", "bootstrap", "bootstrap, Bonferroni", "point estimate"]


def load():
    d = pd.read_csv(PROC / "variants.csv", low_memory=False)
    d["esm_neg"] = -d.esm_llr
    d["src"] = ((d.stars >= 2) & d.group.isin(["path", "benign"])).astype(int)
    d["y_src"] = (d.group == "path").astype(int)
    return d


def clinvar_posterior(d, col):
    src = d[d.src == 1].dropna(subset=[col])
    return sc.Platt().fit(src[col].values, src.y_src.values.astype(int))


def gene_data(d, g, col, cal):
    t = d[(d.gene == g) & d.assay_label.isin([0, 1]) & (d.src == 0) & d[col].notna()]
    return t, cal.predict(t[col].values), t.assay_label.values.astype(int), t.cluster.values


def weighted_rates(dec, y, cl):
    """Residue-weighted conditional call rates (draw a residue, then a variant within it)."""
    _, inv, cnt = np.unique(cl, return_inverse=True, return_counts=True)
    w = 1.0 / cnt[inv]

    def rate(mask, k):
        m = y == k
        return np.sum(w[m] * mask[m]) / np.sum(w[m])

    tp, fp, fn, tn = rate(dec == 1, 1), rate(dec == 1, 0), rate(dec == 0, 1), rate(dec == 0, 0)
    return dict(lr_pos=sc.ratio(tp, fp), lr_neg=sc.ratio(fn, tn),
                yield_lof=float(np.sum(w * (dec == 1)) / np.sum(w)),
                yield_func=float(np.sum(w * (dec == 0)) / np.sum(w)))


def residue_share(y, cl):
    """LoF share under the residue-first law: draw a residue uniformly, then one of its variants uniformly."""
    _, inv, cnt = np.unique(cl, return_inverse=True, return_counts=True)
    w = 1.0 / cnt[inv]
    return float(np.sum(w * (y == 1)) / np.sum(w))


def boot_bounds(dec, y, rng, pct, pseudo=0.0):
    """Vectorised percentile bootstrap over residues (one variant per residue) for LR+ and LR-.

    pseudo > 0 adds that count to each of the six cells of every resample (0.5 gives the Haldane-corrected bootstrap)."""
    cells = np.array([np.sum((y == 1) & (dec == 1)), np.sum((y == 1) & (dec == 0)), np.sum((y == 1) & (dec < 0)),
                      np.sum((y == 0) & (dec == 1)), np.sum((y == 0) & (dec == 0)), np.sum((y == 0) & (dec < 0))])
    c = rng.multinomial(len(y), cells / cells.sum(), size=N_BOOT).astype(float) + pseudo
    n1, n0 = c[:, :3].sum(1), c[:, 3:].sum(1)
    ok = (n1 > 0) & (n0 > 0)
    c, n1, n0 = c[ok], n1[ok], n0[ok]
    tp, fn, fp, tn = c[:, 0] / n1, c[:, 1] / n1, c[:, 3] / n0, c[:, 4] / n0
    # x / 0 = inf for x > 0; 0 / 0 is undefined and dropped from the percentile
    lp = np.where(fp > 0, tp / np.where(fp > 0, fp, 1), np.where(tp > 0, np.inf, np.nan))
    ln = np.where(tn > 0, fn / np.where(tn > 0, tn, 1), np.where(fn > 0, np.inf, np.nan))
    return float(np.nanpercentile(lp, pct)), float(np.nanpercentile(ln, 100 - pct))


def certify_all(p, y, ic, iz, rng):
    """Thresholds on calibration indices ic for each grid alpha, and all five bounds on certification indices iz."""
    level = DELTA / (2 * len(GRID))
    thrs, cert = {}, {}
    for alpha in GRID:
        thrs[alpha] = sc.mondrian_thresholds(p[ic], y[ic], alpha)
        dz = sc.decisions(p[iz], thrs[alpha])
        c = sc.certify(dz, y[iz], level)
        u = sc.certify(dz, y[iz], DELTA / 2)
        c["lr_pos_unc"], c["lr_neg_unc"] = u["lr_pos_lower"], u["lr_neg_upper"]
        c["lr_pos_boot"], c["lr_neg_boot"] = boot_bounds(dz, y[iz], rng, 5.0)
        c["lr_pos_bootb"], c["lr_neg_bootb"] = boot_bounds(dz, y[iz], rng, 100 * level)
        cert[alpha] = c
    return thrs, cert


KEYS = {"certified": ("lr_pos_lower", "lr_neg_upper"), "Clopper-Pearson, uncorrected": ("lr_pos_unc", "lr_neg_unc"),
        "bootstrap": ("lr_pos_boot", "lr_neg_boot"), "bootstrap, Bonferroni": ("lr_pos_bootb", "lr_neg_bootb"),
        "point estimate": ("lr_pos_hat", "lr_neg_hat")}


def select(cert, mode):
    """Per level: the grid alpha with the largest yield whose bound meets the level. Functional calls are chosen
    among alpha >= the LoF alpha of the same strength, which keeps the two call sets disjoint."""
    kp, kn = KEYS[mode]
    out = {}
    for name, tau in LEVELS_POS.items():
        ok = [(c["yield_lof"], a, c[kp]) for a, c in cert.items() if c[kp] > tau]
        out[name] = max(ok)[1:] if ok else None
    for name, tau in LEVELS_NEG.items():
        floor = out[name.replace("BS3", "PS3")]
        ok = [(c["yield_func"], a, c[kn]) for a, c in cert.items()
              if c[kn] < tau and (floor is None or a >= floor[0])]
        out[name] = max(ok)[1:] if ok else None
    return out


def judge(name, best, truth):
    """Compare a selected rule's bound and target with a reference likelihood ratio (test split or oracle)."""
    pos = name.startswith("PS3")
    tau = (LEVELS_POS if pos else LEVELS_NEG)[name]
    lr = truth["lr_pos"] if pos else truth["lr_neg"]
    return dict(alpha=best[0], bound=best[1], ref_lr=lr,
                holds=(lr > tau) if pos else (lr < tau),
                bound_ok=(lr >= best[1]) if pos else (lr <= best[1]),
                ref_yield=truth["yield_lof"] if pos else truth["yield_func"])


def evidence(d, oracle=False):
    rows = []
    rng = np.random.default_rng(SEED + (10 if oracle else 0))
    for pr, col in PRED.items():
        cal = clinvar_posterior(d, col)
        for g in GENES:
            t, p, y, cl = gene_data(d, g, col, cal)
            res = np.unique(cl)
            members = {r: np.flatnonzero(cl == r) for r in res}
            for split in range(N_ORACLE if oracle else N_SPLIT):
                if oracle:  # residues drawn with replacement from the finite population of assayed residues
                    na, nb = int(0.3 * len(res)), int(0.5 * len(res))
                    ic = np.array([rng.choice(members[r]) for r in rng.choice(res, na)])
                    iz = np.array([rng.choice(members[r]) for r in rng.choice(res, nb)])
                    it = np.arange(len(y))
                else:
                    perm = rng.permutation(res)
                    a, b = int(0.3 * len(res)), int(0.8 * len(res))
                    ic = np.array([rng.choice(members[r]) for r in perm[:a]])
                    iz = np.array([rng.choice(members[r]) for r in perm[a:b]])
                    it = np.flatnonzero(np.isin(cl, perm[b:]))
                thrs, cert = certify_all(p, y, ic, iz, rng)
                for mode in MODES:
                    for name, best in select(cert, mode).items():
                        row = dict(predictor=pr, gene=g, split=split, level=name, selection=mode,
                                   achieved=best is not None)
                        if best is not None:
                            truth = weighted_rates(sc.decisions(p[it], thrs[best[0]]), y[it], cl[it])
                            row.update(judge(name, best, truth))
                        rows.append(row)
    return pd.DataFrame(rows)


def reliability(d):
    rows = []
    rng = np.random.default_rng(SEED + 1)
    for pr, col in PRED.items():
        cal = clinvar_posterior(d, col)
        for g in GENES:
            t, p, y, cl = gene_data(d, g, col, cal)
            res = np.unique(cl)
            for split in range(N_REL):
                m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                keep = np.flatnonzero(m)[sc.one_per_cluster(cl[m], rng)]
                for name, delta in (("in expectation", None), ("PAC, delta = 0.1", 0.1)):
                    thr = sc.mondrian_thresholds(p[keep], y[keep], 0.10, delta=delta)
                    rep = sc.coverage_report(sc.decisions(p[~m], thr), y[~m], cl[~m])
                    rows.append(dict(predictor=pr, gene=g, split=split, thresholds=name,
                                     n_cal_lof=int((y[keep] == 1).sum()), **rep))
    return pd.DataFrame(rows)


def vus_shift(d):
    rows = []
    rng = np.random.default_rng(SEED + 2)
    for pr, col in PRED.items():
        cal = clinvar_posterior(d, col)
        for g in GENES:
            t, p, y, cl = gene_data(d, g, col, cal)
            vus = (t.group == "vus").values
            res = np.unique(cl)
            for rep in range(N_SHIFT):
                in_half = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                src_m, tgt_m = in_half & ~vus, ~in_half & vus      # residue-disjoint
                perm = rng.permutation(np.unique(cl[src_m]))
                cut = int(0.4 * len(perm))
                mc, mz = src_m & np.isin(cl, perm[:cut]), src_m & np.isin(cl, perm[cut:])
                ic = np.flatnonzero(mc)[sc.one_per_cluster(cl[mc], rng)]
                iz = np.flatnonzero(mz)[sc.one_per_cluster(cl[mz], rng)]
                iv = np.flatnonzero(tgt_m)
                if len(iv) < 30 or len(np.unique(y[iv])) < 2 or len(np.unique(y[ic])) < 2:
                    continue
                for delta, name in ((None, "in expectation"), (0.1, "PAC, delta = 0.1")):
                    thr = sc.mondrian_thresholds(p[ic], y[ic], 0.10, delta=delta)
                    rows.append(dict(predictor=pr, gene=g, rep=rep, part="coverage", thresholds=name,
                                     cal_prev=y[ic].mean(), vus_prev=y[iv].mean(), vus_prev_rw=residue_share(y[iv], cl[iv]),
                                     **sc.coverage_report(sc.decisions(p[iv], thr), y[iv], cl[iv])))
                thrs, cert = certify_all(p, y, ic, iz, rng)
                for name, best in select(cert, "certified").items():
                    row = dict(predictor=pr, gene=g, rep=rep, part="evidence", level=name, achieved=best is not None)
                    if best is not None:
                        row.update(judge(name, best, weighted_rates(sc.decisions(p[iv], thrs[best[0]]), y[iv], cl[iv])))
                    rows.append(row)
    return pd.DataFrame(rows)


def region_evidence(d):
    rows = []
    rng = np.random.default_rng(SEED + 3)
    for pr, col in PRED.items():
        cal = clinvar_posterior(d, col)
        for g in GENES:
            t, p, y, cl = gene_data(d, g, col, cal)
            pos = t.protein_pos.values
            edges = np.quantile(np.unique(pos), np.linspace(0, 1, 5))
            blk = np.clip(np.searchsorted(edges, pos, side="right") - 1, 0, 3)
            for b in range(4):
                it = np.flatnonzero(blk == b)
                if len(np.unique(y[it])) < 2:
                    continue
                tr_res = np.unique(cl[blk != b])
                for rep in range(25):
                    perm = rng.permutation(tr_res)
                    cut = int(0.4 * len(perm))
                    mc, mz = np.isin(cl, perm[:cut]), np.isin(cl, perm[cut:])
                    ic = np.flatnonzero(mc)[sc.one_per_cluster(cl[mc], rng)]
                    iz = np.flatnonzero(mz)[sc.one_per_cluster(cl[mz], rng)]
                    thrs, cert = certify_all(p, y, ic, iz, rng)
                    for name, best in select(cert, "certified").items():
                        row = dict(predictor=pr, gene=g, block=b, rep=rep, level=name, achieved=best is not None)
                        if best is not None:
                            row.update(judge(name, best, weighted_rates(sc.decisions(p[it], thrs[best[0]]), y[it], cl[it])))
                        rows.append(row)
    return pd.DataFrame(rows)


def main():
    d = load()
    evidence(d, oracle=True).to_csv(RES / "oracle.csv", index=False)
    evidence(d).to_csv(RES / "evidence.csv", index=False)
    reliability(d).to_csv(RES / "reliability.csv", index=False)
    vus_shift(d).to_csv(RES / "vus_shift.csv", index=False)
    region_evidence(d).to_csv(RES / "region_evidence.csv", index=False)
    print("done")


if __name__ == "__main__":
    main()
