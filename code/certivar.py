"""certivar: per-class (Mondrian) conformal calls for variant-effect predictors, with PAC thresholds
and certified likelihood ratios (with respect to the assay outcome) on the ACMG/AMP OddsPath scale.

Two label types are used and never mixed:
  y_clinvar = 1 for ClinVar pathogenic/likely pathogenic, 0 for benign/likely benign. Used only by Platt (and the
              label-shift helpers below) to put scores on a common, monotone scale q(s).
  y_assay   = 1 for loss of function (LoF) in the functional assay, 0 for functional. Used by mondrian_thresholds,
              decisions, coverage_report and certify; every guarantee and certified likelihood ratio refers to y_assay.
s is a predictor score oriented so that higher = more damaging. A *source* calibrator maps s to P_S(y = 1 | s) under
the source prior pi_S. Under label shift (p(s | y) fixed, p(y) changes) the target
posterior is the prior-adjusted source posterior (Saerens et al. 2002):

    P_T(y=1|s) = w1 q / (w1 q + w0 (1 - q)),  q = P_S(y=1|s),  w1 = pi_T/pi_S, w0 = (1-pi_T)/(1-pi_S).

pi_T is estimated from unlabeled target scores by the EM fixed point of the same paper.

Conformal part. Nonconformity of label k at score s is 1 - P(y = k | s). For each class
k a threshold is the finite-sample quantile of calibration nonconformities among
calibration points with y = k (Mondrian / class-conditional split conformal, Vovk et al.).
Class-conditional thresholds do not depend on the class prior, so the per-class coverage
guarantee carries over to targets with a different class share, provided the score distribution within each class is
unchanged (class-conditional invariance, Assumption 2 of the paper). To respect dependence among substitutions at the same
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
        self.slope = float(self.m.coef_[0, 0])
        if self.slope <= 0:  # calls are invariant to the map only if it is increasing in s
            raise ValueError(f"Platt slope {self.slope:.4g} is not positive; orient the score so higher = more damaging")
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


def residue_weighted_thresholds(p_cal, y_cal, clusters, alpha):
    """Per-class thresholds from all calibration variants with residue weights.

    Each class-k variant gets weight 1 / (number of class-k calibration variants at its residue), so every residue
    carries total weight one within a class. The threshold is the smallest nonconformity whose cumulative weight
    reaches the in-expectation rank ceil((m_k + 1)(1 - alpha)), where m_k is the number of residues with a class-k
    variant (infinite if the rank exceeds m_k). This uses every label while matching the residue-first law Q.
    """
    p_cal, y_cal, clusters = map(np.asarray, (p_cal, y_cal, clusters))
    thr = {}
    for k in (0, 1):
        m = y_cal == k
        r = 1 - p_cal[m] if k == 1 else p_cal[m]
        _, inv, cnt = np.unique(clusters[m], return_inverse=True, return_counts=True)
        rank = int(np.ceil((len(cnt) + 1) * (1 - alpha)))
        if rank > len(cnt):
            thr[k] = np.inf
            continue
        order = np.argsort(r, kind="stable")
        cum = np.cumsum((1.0 / cnt[inv])[order])
        thr[k] = r[order][np.searchsorted(cum, rank - 1e-9)]
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


def ratio(num, den):
    """num / den with the conventions of the paper: x / 0 = inf for x > 0, and 0 / 0 is undefined (nan)."""
    if den > 0:
        return num / den
    return np.inf if num > 0 else np.nan


def certify(dec, y, level):
    """Certified bounds on the likelihood ratios of the two call types.

    dec, y: decisions and assay labels for one variant per residue on data independent of the thresholds.
    With one-sided Clopper-Pearson endpoints l_TP, u_FP, u_FN, l_TN at `level` each, returns
    L+ = l_TP / u_FP, a lower bound on LR+ = P(call LoF | Y=1) / P(call LoF | Y=0), and
    U- = u_FN / l_TN, an upper bound on LR- = P(call functional | Y=1) / P(call functional | Y=0)
    (U- = inf when l_TN = 0); each holds with probability at least 1 - 2 * level. Also returns the empirical
    ratios (inf for x/0 with x > 0, nan for 0/0 or an absent class) and the call yields.
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
    both = n1 > 0 and n0 > 0
    return {
        "lr_pos_lower": ratio(tp_lo, fp_hi), "lr_pos_hat": ratio(a / n1, b / n0) if both else np.nan,
        "lr_neg_upper": np.inf if tn_lo == 0 else fn_hi / tn_lo, "lr_neg_hat": ratio(c / n1, e / n0) if both else np.nan,
        "yield_lof": float((dec == 1).mean()), "yield_func": float((dec == 0).mean()), "n1": n1, "n0": n0,
    }


def katz_bounds(dec, y, level):
    """One-sided Katz log intervals (Katz et al. 1978) for LR+ = TP/FP (lower) and LR- = FN/TN (upper) at `level`.

    The usual large-sample interval for a ratio of two binomial proportions; when a count in a ratio is zero, 0.5 is
    added to every cell of that ratio (each count and its complement). Returns (L+, U-), nan if a class is absent.
    """
    from scipy.stats import norm
    dec, y = np.asarray(dec), np.asarray(y)
    n1, n0 = int((y == 1).sum()), int((y == 0).sum())
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan
    z = float(norm.ppf(1 - level))

    def one(a, b, sign):
        m1, m0 = n1, n0
        if a == 0 or b == 0:
            a, b, m1, m0 = a + 0.5, b + 0.5, n1 + 1, n0 + 1
        se = np.sqrt(max(1 / a - 1 / m1 + 1 / b - 1 / m0, 0.0))
        return (a / m1) / (b / m0) * np.exp(sign * z * se)

    a = int(((dec == 1) & (y == 1)).sum()); b = int(((dec == 1) & (y == 0)).sum())
    c = int(((dec == 0) & (y == 1)).sum()); e = int(((dec == 0) & (y == 0)).sum())
    return one(a, b, -1.0), one(c, e, 1.0)


def _mixture_distance(s1, s0, st, pis):
    """min over pi of the Kolmogorov distance between the target ECDF and pi F1 + (1 - pi) F0; returns (D, pi)."""
    x = np.unique(np.concatenate([s1, s0, st]))
    f1 = np.searchsorted(np.sort(s1), x, side="right") / len(s1)
    f0 = np.searchsorted(np.sort(s0), x, side="right") / len(s0)
    ft = np.searchsorted(np.sort(st), x, side="right") / len(st)
    dist = np.abs(ft[None, :] - pis[:, None] * f1[None, :] - (1 - pis[:, None]) * f0[None, :]).max(axis=1)
    i = int(np.argmin(dist))
    return float(dist[i]), float(pis[i])


def mixture_check(s1, s0, st, n_boot=300, rng=None, n_grid=201):
    """Necessary-condition check of class-conditional invariance (Assumption 2 of the paper).

    s1, s0: calibration scores of LoF and functional variants (one per residue); st: unlabelled target scores (one per
    residue). Under Assumption 2 the target score law is pi F1 + (1 - pi) F0 for some pi in [0, 1]. The statistic is the
    smallest Kolmogorov distance between the target empirical distribution and such a mixture. Its null distribution is
    approximated by a bootstrap that redraws both calibration classes from their empirical laws and a target sample of
    the same size from the fitted mixture. Returns (statistic, fitted pi, p-value). A small p-value shows that the target
    is not a re-weighting of the calibration classes, so the guarantees need not transfer; a large p-value does not
    prove invariance.
    """
    rng = rng if rng is not None else np.random.default_rng(0)
    s1, s0, st = map(np.asarray, (s1, s0, st))
    pis = np.linspace(0, 1, n_grid)
    d, pi = _mixture_distance(s1, s0, st, pis)
    exceed = 0
    for _ in range(n_boot):
        b1, b0 = rng.choice(s1, len(s1)), rng.choice(s0, len(s0))
        lof = rng.random(len(st)) < pi
        bt = np.where(lof, rng.choice(s1, len(st)), rng.choice(s0, len(st)))
        exceed += _mixture_distance(b1, b0, bt, pis)[0] >= d
    return d, pi, (1 + exceed) / (1 + n_boot)


def fixed_sequence_select(alphas, bounds, yields, tau, lof=True):
    """Fixed-sequence certification of one evidence level (Theorem 2 of the paper).

    alphas, bounds, yields: per threshold level, in the testing order fixed in advance; bounds are l_TP / u_FP
    (lof=True) or u_FN / l_TN (lof=False), with each Clopper-Pearson endpoint at level delta / 2. Hypotheses are tested
    in the given order and testing stops at the first level that does not meet tau (bound > tau for LoF calls,
    bound < tau for functional calls). Returns the rejected alpha with the largest yield (ties: the larger alpha), or
    None if the first level already fails. This is the rule used in code/17_extensions.py.
    """
    rejected = []
    for a, b, y in zip(alphas, bounds, yields):
        if not (b > tau if lof else b < tau):
            break
        rejected.append((y, a))
    return max(rejected)[1] if rejected else None
