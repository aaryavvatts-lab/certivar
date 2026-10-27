"""certivar: per-class (Mondrian) conformal calls for variant-effect predictors, with PAC thresholds
and certified likelihood ratios (with respect to the assay outcome) on the ACMG/AMP OddsPath scale.

Notation. y = 1 for loss of function / pathogenic, 0 otherwise. s(x) is a predictor
score (higher = more damaging). A *source* calibrator maps s to P_S(y = 1 | s) under
the source prior pi_S. Under label shift (p(s | y) fixed, p(y) changes) the target
posterior is the prior-adjusted source posterior (Saerens et al. 2002):

    P_T(y=1|s) = w1 q / (w1 q + w0 (1 - q)),  q = P_S(y=1|s),  w1 = pi_T/pi_S, w0 = (1-pi_T)/(1-pi_S).

pi_T is estimated from unlabeled target scores by the EM fixed point of the same paper.

Conformal part. Nonconformity of label k at score s is 1 - P(y = k | s). For each class
k a threshold is the finite-sample quantile of calibration nonconformities among
calibration points with y = k (Mondrian / class-conditional split conformal, Vovk et al.).
Class-conditional thresholds do not depend on the class prior, so the per-class coverage
guarantee survives any label shift. To respect dependence among substitutions at the same
residue, calibration uses one randomly chosen variant per (gene, protein position)
cluster, which makes calibration points exchangeable with a test variant drawn from a
new position. The prediction set is {k : nonconformity_k <= threshold_k}; the set {0, 1}
means abstain, and the empty set (rare) is also reported as abstain.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import beta as beta_dist
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression

EPS = 1e-6


# ---------------------------------------------------------------- calibrators
class Platt:
    """Logistic calibration of a single score (Platt 1999)."""

    def fit(self, s, y):
        self.m = LogisticRegression(C=1e6).fit(np.asarray(s).reshape(-1, 1), y)
        self.prior = float(np.mean(y))
        return self

    def predict(self, s):
        return self.m.predict_proba(np.asarray(s).reshape(-1, 1))[:, 1]


class Isotonic:
    """Isotonic calibration (Zadrozny and Elkan 2002)."""

    def fit(self, s, y):
        self.m = IsotonicRegression(y_min=EPS, y_max=1 - EPS, out_of_bounds="clip").fit(s, y)
        self.prior = float(np.mean(y))
        return self

    def predict(self, s):
        return self.m.predict(np.asarray(s))


# ---------------------------------------------------------------- label shift
def adjust(q, pi_s, pi_t):
    """Prior-adjust source posteriors q from prior pi_s to prior pi_t."""
    q = np.clip(q, EPS, 1 - EPS)
    w1, w0 = pi_t / pi_s, (1 - pi_t) / (1 - pi_s)
    return w1 * q / (w1 * q + w0 * (1 - q))


def em_prior(q, pi_s, tol=1e-8, max_iter=1000):
    """EM estimate of the target prior from unlabeled source posteriors (Saerens et al. 2002)."""
    pi = pi_s
    for _ in range(max_iter):
        new = float(np.mean(adjust(q, pi_s, pi)))
        if abs(new - pi) < tol:
            return new
        pi = new
    return pi


def bbse_prior(q_cal, y_cal, q_tgt, thr=0.5):
    """Black-box shift estimate (Lipton et al. 2018) with a thresholded binary classifier."""
    h_cal = (np.asarray(q_cal) >= thr).astype(int)
    y_cal = np.asarray(y_cal)
    tpr = h_cal[y_cal == 1].mean()
    fpr = h_cal[y_cal == 0].mean()
    rate = (np.asarray(q_tgt) >= thr).mean()
    if tpr - fpr < 1e-3:
        return np.nan
    return float(np.clip((rate - fpr) / (tpr - fpr), 0, 1))


# ---------------------------------------------------------------- conformal
def one_per_cluster(clusters, rng):
    """Indices of one randomly chosen member per cluster."""
    clusters = np.asarray(clusters)
    order = rng.permutation(len(clusters))
    _, first = np.unique(clusters[order], return_index=True)
    return np.sort(order[first])


def pac_rank(n, alpha, delta):
    """Smallest order-statistic rank j with P(coverage of the j-th smallest < 1 - alpha) <= delta.

    For n exchangeable continuous nonconformities the coverage of the j-th order statistic is
    Beta(j, n - j + 1) distributed (training-conditional validity, Vovk 2013), so the rank is the
    smallest j with BetaCDF(1 - alpha; j, n - j + 1) <= delta. Returns n + 1 (threshold = inf) if none.
    """
    j = np.arange(1, n + 1)
    ok = beta_dist.cdf(1 - alpha, j, n - j + 1) <= delta
    return int(j[ok][0]) if ok.any() else n + 1


def mondrian_thresholds(p_cal, y_cal, alpha, delta=None):
    """Per-class conformal thresholds on nonconformity 1 - P(y = k | s).

    delta=None gives the usual in-expectation guarantee; a delta in (0, 1) gives the PAC version:
    with probability at least 1 - delta over the calibration draw, per-class miscoverage <= alpha.
    `alpha` may be a number or a dict {0: alpha_0, 1: alpha_1}.
    """
    p_cal, y_cal = np.asarray(p_cal), np.asarray(y_cal)
    thr = {}
    for k in (0, 1):
        a = alpha[k] if isinstance(alpha, dict) else alpha
        sc = np.sort(1 - p_cal[y_cal == 1]) if k == 1 else np.sort(p_cal[y_cal == 0])
        n = len(sc)
        rank = int(np.ceil((n + 1) * (1 - a))) if delta is None else pac_rank(n, a, delta)
        thr[k] = np.inf if rank > n else sc[rank - 1]
    return thr


def weighted_decisions(p_cal, y_cal, p_tst, alpha, pi_t):
    """Label-shift weighted split conformal (Podkopaev and Ramdas 2021), marginal nonconformity.

    Calibration point i gets weight w(y_i) = pi_t(y_i) / pi_cal(y_i); for candidate label k the
    test point contributes mass w(k) at +infinity. Label k enters the set if its nonconformity is
    at most the weighted (1 - alpha) quantile.
    """
    p_cal, y_cal, p_tst = map(np.asarray, (p_cal, y_cal, p_tst))
    pi_c = y_cal.mean()
    w = {1: pi_t / pi_c, 0: (1 - pi_t) / (1 - pi_c)}
    nc = np.where(y_cal == 1, 1 - p_cal, p_cal)
    order = np.argsort(nc)
    nc_s, w_s = nc[order], np.where(y_cal[order] == 1, w[1], w[0])
    cum = np.cumsum(w_s)
    thr = {}
    for k in (0, 1):
        total = cum[-1] + w[k]
        idx = np.searchsorted(cum / total, 1 - alpha, side="left")
        thr[k] = np.inf if idx >= len(nc_s) else nc_s[idx]
    return decisions(p_tst, thr)


def predict_sets(p, thr):
    """Return (in1, in0) boolean arrays: is label 1 / label 0 in the prediction set."""
    p = np.asarray(p)
    return (1 - p) <= thr[1], p <= thr[0]


def decisions(p, thr):
    """1 = call LoF, 0 = call functional, -1 = abstain (set {0, 1}), -2 = empty set (also not called)."""
    in1, in0 = predict_sets(p, thr)
    out = np.full(len(p), -1)
    out[in1 & ~in0] = 1
    out[in0 & ~in1] = 0
    out[~in1 & ~in0] = -2
    return out


class ShiftConformal:
    """Mondrian conformal abstention calibrated on labelled assay data, one variant per position."""

    def __init__(self, alpha=0.1, seed=0):
        self.alpha, self.rng = alpha, np.random.default_rng(seed)

    def fit(self, p_cal, y_cal, clusters):
        keep = one_per_cluster(clusters, self.rng)
        self.n_cal = len(keep)
        self.thr = mondrian_thresholds(np.asarray(p_cal)[keep], np.asarray(y_cal)[keep], self.alpha)
        return self

    def predict(self, p):
        return decisions(p, self.thr)


# ---------------------------------------------------------------- metrics
def brier(p, y):
    return float(np.mean((np.asarray(p) - np.asarray(y)) ** 2))


def ece(p, y, bins=10):
    """Expected calibration error with equal-mass bins."""
    p, y = np.asarray(p), np.asarray(y)
    order = np.argsort(p)
    total = 0.0
    for idx in np.array_split(order, bins):
        if len(idx):
            total += len(idx) * abs(p[idx].mean() - y[idx].mean())
    return float(total / len(p))


def coverage_report(dec, y, clusters=None):
    """Coverage, per-class coverage, abstention rate and error rate among committed calls.

    With `clusters`, also reports residue-weighted per-class coverage: each variant gets weight
    1 / (number of evaluated variants at its residue), i.e. the law "draw a residue uniformly, then
    a variant uniformly within it", conditioned on the class. This is the quantity Proposition 1 bounds.
    """
    dec, y = np.asarray(dec), np.asarray(y)
    covered = (dec == -1) | (dec == y)          # the set contains the true label; empty sets (-2) do not
    committed = dec >= 0
    out = {
        "coverage": float(covered.mean()),
        "coverage_lof": float(covered[y == 1].mean()) if (y == 1).any() else np.nan,
        "coverage_func": float(covered[y == 0].mean()) if (y == 0).any() else np.nan,
        "abstain": float((~committed).mean()),
        "empty": float((dec == -2).mean()),
        "error_committed": float((dec[committed] != y[committed]).mean()) if committed.any() else np.nan,
        "n": int(len(y)),
    }
    if clusters is not None:
        clusters = np.asarray(clusters)
        _, inv, cnt = np.unique(clusters, return_inverse=True, return_counts=True)
        w_all = 1.0 / cnt[inv]
        for k, name in ((1, "coverage_lof_rw"), (0, "coverage_func_rw")):
            m = y == k
            if m.any():
                w = w_all[m]
                out[name] = float(np.sum(w * covered[m]) / np.sum(w))
            else:
                out[name] = np.nan
    return out


def cluster_bootstrap(stat, clusters, n_boot=2000, seed=0):
    """Percentile CI for stat(index_array) resampling clusters with replacement."""
    rng = np.random.default_rng(seed)
    clusters = np.asarray(clusters)
    uniq, inv = np.unique(clusters, return_inverse=True)
    members = [np.flatnonzero(inv == i) for i in range(len(uniq))]
    vals = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(uniq), len(uniq))
        idx = np.concatenate([members[i] for i in pick])
        v = stat(idx)
        if np.isfinite(v):
            vals.append(v)
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


# ---------------------------------------------------------------- certified likelihood ratios
def clopper_pearson(k, n, level):
    """Two one-sided (1 - level) Clopper-Pearson bounds for a binomial proportion k / n."""
    lo = 0.0 if k == 0 else float(beta_dist.ppf(level, k, n - k + 1))
    hi = 1.0 if k == n else float(beta_dist.ppf(1 - level, k + 1, n - k))
    return lo, hi


def certify(dec, y, level):
    """Certified bounds on the likelihood ratios of the two call types.

    dec, y: decisions and labels for one variant per residue on data independent of the thresholds.
    Returns lower bound on LR+ = P(call LoF | LoF) / P(call LoF | functional) and upper bound on
    LR- = P(call functional | LoF) / P(call functional | functional), each holding with probability
    at least 1 - 2 * level (union of two one-sided bounds), plus point estimates and yields.
    """
    dec, y = np.asarray(dec), np.asarray(y)
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    a = int(((dec == 1) & (y == 1)).sum())
    b = int(((dec == 1) & (y == 0)).sum())
    c = int(((dec == 0) & (y == 1)).sum())
    e = int(((dec == 0) & (y == 0)).sum())
    tp_lo, _ = clopper_pearson(a, n1, level)
    _, fp_hi = clopper_pearson(b, n0, level)
    _, fn_hi = clopper_pearson(c, n1, level)
    tn_lo, _ = clopper_pearson(e, n0, level)
    eps = 1e-12
    return {
        "lr_pos_lower": tp_lo / max(fp_hi, eps), "lr_pos_hat": (a / n1) / max(b / n0, eps) if n1 and n0 else np.nan,
        "lr_neg_upper": fn_hi / max(tn_lo, eps), "lr_neg_hat": (c / n1) / max(e / n0, eps) if n1 and n0 else np.nan,
        "yield_lof": float((dec == 1).mean()), "yield_func": float((dec == 0).mean()), "n1": n1, "n0": n0,
    }
