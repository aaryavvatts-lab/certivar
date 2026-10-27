"""Controlled label-shift simulation.

Scores are generated per residue cluster: s = mu_y + u_c + e, with u_c ~ N(0, tau^2) a
shared position effect (so substitutions at one residue are dependent) and e ~ N(0, 1).
A calibrator is fitted on a source sample with prior pi_S = 0.45 (the ClinVar-like
setting) and applied to target samples whose prior sweeps 0.05-0.60.

Compared, over R replicates per target prior:
  calibration  - Platt as fitted (source prior), EM prior-adjusted, oracle prior-adjusted
  abstention   - pooled split conformal, label-shift weighted conformal (EM prior) and
                 Mondrian (class-conditional) conformal,
                 all calibrated on source labels with one variant per cluster.
Output: results/simulation.csv
"""
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
MU1, TAU, PER = 2.0, 0.5, 5          # class separation, position effect SD, variants per position
PI_S, N_SRC, N_TGT, ALPHA, R = 0.45, 600, 3000, 0.10, 300


def draw(n, pi, rng):
    n_pos = int(np.ceil(n / PER))
    y_pos = rng.random(n_pos) < pi                     # class is shared within a residue
    u = rng.normal(0, TAU, n_pos)
    pos = np.repeat(np.arange(n_pos), PER)[:n]
    y = y_pos[pos].astype(int)
    s = MU1 * y + u[pos] + rng.normal(0, 1, n)
    return s, y, pos


def pooled_threshold(p, y, alpha):
    """Standard (marginal) split conformal on nonconformity 1 - P(y_true | s)."""
    nc = np.sort(np.where(y == 1, 1 - p, p))
    rank = int(np.ceil((len(nc) + 1) * (1 - alpha)))
    t = nc[min(rank, len(nc)) - 1]
    return {0: t, 1: t}


def main():
    rows = []
    rng = np.random.default_rng(2027)
    for pi_t in [0.05, 0.10, 0.15, 0.20, 0.30, 0.45, 0.60]:
        for r in range(R):
            s_src, y_src, pos_src = draw(N_SRC, PI_S, rng)
            s_cal, y_cal, pos_cal = draw(N_SRC, PI_S, rng)
            s_tgt, y_tgt, _ = draw(N_TGT, pi_t, rng)
            cal = sc.Platt().fit(s_src, y_src)
            q_tgt = cal.predict(s_tgt)
            pi_hat = sc.em_prior(q_tgt, cal.prior)
            probs = {
                "Platt (source prior)": q_tgt,
                "EM prior-adjusted": sc.adjust(q_tgt, cal.prior, pi_hat),
                "Oracle prior-adjusted": sc.adjust(q_tgt, cal.prior, pi_t),
            }
            for name, p in probs.items():
                rows.append({"pi_t": pi_t, "rep": r, "kind": "calibration", "method": name,
                             "ece": sc.ece(p, y_tgt), "brier": sc.brier(p, y_tgt),
                             "pi_hat": pi_hat})
            keep = sc.one_per_cluster(pos_cal, rng)
            p_cal = cal.predict(s_cal[keep])
            decs = {
                "Pooled conformal": sc.decisions(q_tgt, pooled_threshold(p_cal, y_cal[keep], ALPHA)),
                "Label-shift weighted conformal": sc.weighted_decisions(p_cal, y_cal[keep], q_tgt, ALPHA, pi_hat),
                "Mondrian conformal": sc.decisions(q_tgt, sc.mondrian_thresholds(p_cal, y_cal[keep], ALPHA)),
            }
            for name, dec in decs.items():
                rep = sc.coverage_report(dec, y_tgt)
                rows.append({"pi_t": pi_t, "rep": r, "kind": "conformal", "method": name, **rep})
    df = pd.DataFrame(rows)
    RES.mkdir(exist_ok=True)
    df.to_csv(RES / "simulation.csv", index=False)
    summ = df.groupby(["kind", "method", "pi_t"]).mean(numeric_only=True)
    print(summ[["ece", "brier", "pi_hat", "coverage", "coverage_lof", "coverage_func", "abstain"]]
          .round(3).to_string())


if __name__ == "__main__":
    main()
