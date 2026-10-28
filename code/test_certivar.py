"""Boundary-case tests for certivar (run: python code/test_certivar.py, or make test)."""
import numpy as np

import certivar as cv


def test_ratio_conventions():
    assert cv.ratio(0.2, 0.1) == 2.0
    assert cv.ratio(0.3, 0.0) == np.inf          # x / 0 with x > 0
    assert np.isnan(cv.ratio(0.0, 0.0))          # 0 / 0 is undefined


def test_clopper_pearson_endpoints():
    assert cv.clopper_pearson(0, 10, 0.025)[0] == 0.0     # x = 0: lower endpoint 0
    assert cv.clopper_pearson(10, 10, 0.025)[1] == 1.0    # x = n: upper endpoint 1
    assert cv.clopper_pearson(0, 0, 0.025) == (0.0, 1.0)  # empty class: trivial bounds
    lo, hi = cv.clopper_pearson(3, 20, 0.025)
    assert 0 < lo < 3 / 20 < hi < 1


def test_tn_zero_gives_infinite_upper_bound():
    # no functional variant is called functional: l_TN = 0, so U- must be +inf (never certified)
    y = np.array([1] * 10 + [0] * 10)
    dec = np.array([0] * 3 + [1] * 7 + [-1] * 10)
    b = cv.certify(dec, y, 0.025)
    assert b["lr_neg_upper"] == np.inf
    assert b["lr_neg_hat"] == np.inf             # 0.3 / 0


def test_fp_zero():
    # no functional variant called LoF: u_FP > 0 still, so L+ is finite; the empirical ratio is +inf
    y = np.array([1] * 20 + [0] * 20)
    dec = np.array([1] * 15 + [-1] * 5 + [0] * 20)
    b = cv.certify(dec, y, 0.025)
    assert np.isfinite(b["lr_pos_lower"]) and b["lr_pos_lower"] > 0
    assert b["lr_pos_hat"] == np.inf


def test_zero_over_zero_is_undefined():
    y = np.array([1] * 5 + [0] * 5)
    dec = np.full(10, -1)                        # no calls at all
    b = cv.certify(dec, y, 0.025)
    assert np.isnan(b["lr_pos_hat"]) and np.isnan(b["lr_neg_hat"])
    assert b["lr_pos_lower"] == 0.0 and b["lr_neg_upper"] == np.inf


def test_absent_class():
    y = np.ones(10, dtype=int)                   # n0 = 0
    dec = np.array([1] * 8 + [0] * 2)
    b = cv.certify(dec, y, 0.025)
    assert np.isnan(b["lr_pos_hat"]) and np.isnan(b["lr_neg_hat"])
    assert b["lr_neg_upper"] == np.inf           # l_TN = 0 when n0 = 0


def test_certified_bounds_bracket_truth():
    rng = np.random.default_rng(0)
    tp, fp, fn, tn = 0.6, 0.1, 0.1, 0.7
    fails = 0
    for _ in range(400):
        y = (rng.random(400) < 0.3).astype(int)
        u = rng.random(400)
        dec = np.where(y == 1, np.where(u < tp, 1, np.where(u < tp + fn, 0, -1)),
                       np.where(u < fp, 1, np.where(u < fp + tn, 0, -1)))
        b = cv.certify(dec, y, 0.025)
        fails += not (b["lr_pos_lower"] <= tp / fp and b["lr_neg_upper"] >= fn / tn)
    assert fails / 400 <= 0.10                    # guarantee: at most 2 * 2 * 0.025 = 0.10


def test_platt_requires_positive_slope():
    s = np.linspace(-2, 2, 200)
    y = (s + np.random.default_rng(1).normal(0, 0.5, 200) > 0).astype(int)
    assert cv.Platt().fit(s, y).slope > 0
    try:
        cv.Platt().fit(-s, y)
    except ValueError:
        pass
    else:
        raise AssertionError("negative slope was not rejected")


def test_disjoint_calls_across_levels():
    rng = np.random.default_rng(2)
    p, y = rng.random(300), (rng.random(300) < 0.3).astype(int)
    q = rng.random(1000)
    for a1 in (0.05, 0.1, 0.2, 0.3):
        for a2 in (0.05, 0.1, 0.2, 0.3):
            lof = cv.decisions(q, cv.mondrian_thresholds(p, y, a1)) == 1
            fun = cv.decisions(q, cv.mondrian_thresholds(p, y, a2)) == 0
            assert not np.any(lof & fun)


def test_fixed_sequence_stops_at_first_failure():
    # levels in the fixed order; the third fails, so the fourth (largest yield) must not be selected
    g = [0.05, 0.10, 0.20, 0.30]
    assert cv.fixed_sequence_select(g, [5.0, 4.5, 2.0, 9.0], [0.1, 0.2, 0.3, 0.9], tau=4.3) == 0.10
    assert cv.fixed_sequence_select(g[:2], [3.0, 9.0], [0.1, 0.2], tau=4.3) is None
    assert cv.fixed_sequence_select(g[::-1][:3], [0.2, 0.3, 0.6], [0.6, 0.5, 0.4], tau=0.48, lof=False) == 0.30


def test_fixed_sequence_level_control():
    # the level LR+ > tau is false at every threshold, so any claim is false; claims must occur at most delta of the time
    rng = np.random.default_rng(3)
    delta, tau, claims = 0.05, 4.3, 0
    tp, fp = 0.40, 0.40 / 4.3            # true LR+ equals tau exactly (worst case for the test)
    for _ in range(2000):
        alphas, bounds, yields = [0.05, 0.10, 0.20, 0.30], [], []
        for _ in range(4):
            a = rng.binomial(60, tp)
            b = rng.binomial(140, fp)
            lo, _ = cv.clopper_pearson(a, 60, delta / 2)
            _, hi = cv.clopper_pearson(b, 140, delta / 2)
            bounds.append(cv.ratio(lo, hi))
            yields.append(a + b)
        claims += cv.fixed_sequence_select(alphas, bounds, yields, tau) is not None
    assert claims / 2000 <= delta


def test_residue_weights_reduce_to_one_per_residue():
    rng = np.random.default_rng(3)
    p = rng.random(400)
    y = (rng.random(400) < p).astype(int)
    for alpha in (0.05, 0.1, 0.3):
        assert cv.residue_weighted_thresholds(p, y, np.arange(400), alpha) == cv.mondrian_thresholds(p, y, alpha)


def test_katz_bounds_bracket_point_estimate():
    dec = np.array([1] * 20 + [0] * 80 + [-1] * 50)
    y = np.array([1] * 18 + [0] * 2 + [1] * 5 + [0] * 75 + [1] * 10 + [0] * 40)
    lo, hi = cv.katz_bounds(dec, y, 0.0125)
    lr_pos = (18 / 33) / (2 / 117)
    lr_neg = (5 / 33) / (75 / 117)
    assert 0 < lo < lr_pos and hi > lr_neg
    assert all(np.isnan(cv.katz_bounds(dec, np.ones_like(y), 0.05)))


def test_mixture_check_level_and_power():
    rng = np.random.default_rng(4)
    null = power = 0
    for _ in range(40):
        s1, s0 = rng.normal(2, 1, 80), rng.normal(0, 1, 250)
        lof = rng.random(200) < 0.3
        st = np.where(lof, rng.normal(2, 1, 200), rng.normal(0, 1, 200))
        null += cv.mixture_check(s1, s0, st, n_boot=100, rng=rng)[2] < 0.05
        shifted = np.where(lof, rng.normal(1.5, 1, 200), rng.normal(0.5, 1, 200))
        power += cv.mixture_check(s1, s0, shifted, n_boot=100, rng=rng)[2] < 0.05
    assert null / 40 <= 0.15 and power / 40 >= 0.7


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
    print(f"{len(tests)} tests passed")
