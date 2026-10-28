"""Oracle benchmark with every comparison rule on the same draws (E18, Appendix S7).

Same design as the oracle of 17_extensions.py (calibration 30% and certification 50% of residues drawn with replacement,
100 draws per pair, grids G4 and G20), with more comparison rules: Katz intervals with a Bonferroni correction over
the grid, the percentile bootstrap at the Bonferroni level and with 0.5 added to every cell, and fixed-sequence testing
simultaneous over the three evidence levels of a call type (each level at delta / 3). The draws are separate from
those of 17_extensions.py, whose output the main text reports.

Output: results/oracle_baselines.csv
"""
import importlib
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ev = importlib.import_module("13_evidence")
ex = importlib.import_module("17_extensions")
ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
SEED, DELTA, N_ORACLE = ex.SEED, ex.DELTA, ex.N_ORACLE
ev.KEYS.update({"bootstrap, Bonferroni": ("lr_pos_bootb", "lr_neg_bootb"),
                "bootstrap, pseudocount": ("lr_pos_bootpc", "lr_neg_bootpc"),
                "Katz, Bonferroni": ("lr_pos_katz", "lr_neg_katz")})
N_LEVELS = len(ev.LEVELS_POS)   # evidence levels tested per call type (supporting, moderate, strong)


def fixed_sequence(cert, grid, func_descending=True, key="fs"):
    """Fixed-sequence testing of every level (LoF levels in increasing alpha, functional levels in decreasing alpha).
    Key "fs" tests each level at delta (a per-level guarantee, Theorem 2); key "fs3" tests each level at
    delta / N_LEVELS, so that by a union bound every level of a call type holds simultaneously with probability at
    least 1 - delta."""
    out = {}
    kp, kn = f"{key}_pos", f"{key}_neg"
    for name, tau in ev.LEVELS_POS.items():
        a = sc.fixed_sequence_select(grid, [cert[x][kp] for x in grid], [cert[x]["yield_lof"] for x in grid], tau)
        out[name] = (a, cert[a][kp]) if a is not None else None
    order = sorted(grid, reverse=True) if func_descending else list(grid)
    for name, tau in ev.LEVELS_NEG.items():
        a = sc.fixed_sequence_select(order, [cert[x][kn] for x in order], [cert[x]["yield_func"] for x in order],
                                     tau, lof=False)
        out[name] = (a, cert[a][kn]) if a is not None else None
    return out


def oracle_baselines(d):
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
                for gname, grid in ex.GRIDS.items():
                    level = DELTA / (2 * len(grid))
                    thrs, cert = {}, {}
                    for a in grid:
                        thrs[a] = sc.mondrian_thresholds(p[ic], y[ic], a)
                        dz = sc.decisions(p[iz], thrs[a])
                        c = sc.certify(dz, y[iz], level)
                        u = sc.certify(dz, y[iz], DELTA / 2)
                        c["lr_pos_unc"], c["lr_neg_unc"] = u["lr_pos_lower"], u["lr_neg_upper"]
                        c["fs_pos"], c["fs_neg"] = u["lr_pos_lower"], u["lr_neg_upper"]
                        s3 = sc.certify(dz, y[iz], DELTA / (2 * N_LEVELS))
                        c["fs3_pos"], c["fs3_neg"] = s3["lr_pos_lower"], s3["lr_neg_upper"]
                        c["lr_pos_boot"], c["lr_neg_boot"] = ev.boot_bounds(dz, y[iz], rng, 5.0)
                        c["lr_pos_bootb"], c["lr_neg_bootb"] = ev.boot_bounds(dz, y[iz], rng, 100 * level)
                        c["lr_pos_bootpc"], c["lr_neg_bootpc"] = ev.boot_bounds(dz, y[iz], rng, 5.0, pseudo=0.5)
                        c["lr_pos_katz"], c["lr_neg_katz"] = sc.katz_bounds(dz, y[iz], DELTA / len(grid))
                        c["counts"] = ex.counts(dz, y[iz])
                        cert[a] = c
                    sels = {m: ev.select(cert, m) for m in
                            ("certified", "Clopper-Pearson, uncorrected", "bootstrap", "bootstrap, Bonferroni",
                             "bootstrap, pseudocount", "Katz, Bonferroni", "point estimate")}
                    sels["fixed sequence"] = fixed_sequence(cert, grid)
                    sels["fixed sequence, simultaneous"] = fixed_sequence(cert, grid, key="fs3")
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
                                vp, vn = ex.variant_rates(dec, y)
                                row["vw_lr"] = vp if pos else vn
                                row["vw_holds"] = (vp > tau) if pos else (vn < tau)
                                rep = sc.coverage_report(dec, y, cl)
                                row["cov_lof"], row["cov_func"] = rep["coverage_lof_rw"], rep["coverage_func_rw"]
                                row.update({f"cert_{k}": v for k, v in cert[best[0]]["counts"].items()})
                            rows.append(row)
            print(pr, g, flush=True)
    return pd.DataFrame(rows)


def main():
    oracle_baselines(ev.load()).to_csv(RES / "oracle_baselines.csv", index=False)
    print("done")


if __name__ == "__main__":
    main()
