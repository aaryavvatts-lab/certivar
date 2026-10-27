"""Sensitivity of the results to the assay-label rule (E10).

The BRCA2 labels come from Gaussians fitted to the assay's synonymous and nonsense controls (posterior
cut-off 0.95), the other genes' labels from the assay authors' classes. Here each gene is relabelled with
alternative rules and the two main experiments are repeated:

  BRCA2   control rule with posterior cut-off 0.80, 0.95 (main) and 0.99
  BAP1    investigator classes (main) and the control rule (cut-off 0.95)
  RAD51C  investigator classes (main) and the control rule (cut-off 0.95)

For each rule and predictor: mean per-class coverage and abstention of the in-expectation per-class rule
over 50 random halvings of the residues (alpha = 0.1), and the oracle check of the certified selection
(100 draws with replacement, as in 13_evidence.py, certified rule only).

Output: results/label_sensitivity.csv
"""
import importlib
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ev = importlib.import_module("13_evidence")
ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
N_HALF, N_ORACLE, SEED = 50, 100, 2027


def relabel(d, gene, rule):
    """Return a copy of d whose assay_label for `gene` follows `rule`; other genes are unchanged."""
    d = d.copy()
    m = d.gene == gene
    if rule != "main":
        cut = {"controls 0.80": 0.80, "controls 0.95": 0.95, "controls 0.99": 0.99}[rule]
        p = d.loc[m, "p_lof_ctrl"]
        lab = np.select([p >= cut, p <= 1 - cut], [1.0, 0.0], default=-1.0)
        lab[p.isna().values] = np.nan
        d.loc[m, "assay_label"] = lab
        d.loc[m & (d.score.isna() | d.splice_region), "assay_label"] = np.nan
    return d


def halvings(p, y, cl, rng):
    rows = []
    res = np.unique(cl)
    for _ in range(N_HALF):
        m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
        keep = np.flatnonzero(m)[sc.one_per_cluster(cl[m], rng)]
        thr = sc.mondrian_thresholds(p[keep], y[keep], 0.10)
        rows.append(sc.coverage_report(sc.decisions(p[~m], thr), y[~m], cl[~m]))
    r = pd.DataFrame(rows)
    return dict(lof_cov=r.coverage_lof_rw.mean(), func_cov=r.coverage_func_rw.mean(), abstain=r.abstain.mean())


def oracle(p, y, cl, rng):
    res = np.unique(cl)
    members = {r: np.flatnonzero(cl == r) for r in res}
    level = ev.DELTA / (2 * len(ev.GRID))
    claims = ok = draws_sup_lof = draws_sup_func = 0
    for _ in range(N_ORACLE):
        na, nb = int(0.3 * len(res)), int(0.5 * len(res))
        ic = np.array([rng.choice(members[r]) for r in rng.choice(res, na)])
        iz = np.array([rng.choice(members[r]) for r in rng.choice(res, nb)])
        thrs, cert = {}, {}
        for a in ev.GRID:
            thrs[a] = sc.mondrian_thresholds(p[ic], y[ic], a)
            cert[a] = sc.certify(sc.decisions(p[iz], thrs[a]), y[iz], level)
        for name, best in ev.select(cert, "certified").items():
            if best is None:
                continue
            truth = ev.weighted_rates(sc.decisions(p, thrs[best[0]]), y, cl)
            j = ev.judge(name, best, truth)
            claims += 1
            ok += bool(j["bound_ok"])
            draws_sup_lof += name == "PS3-supporting"
            draws_sup_func += name == "BS3-supporting"
    return dict(oracle_claims=claims, oracle_bound_ok=ok, share_lof_supporting=draws_sup_lof / N_ORACLE,
                share_func_supporting=draws_sup_func / N_ORACLE)


def main():
    d = ev.load()
    rules = [("BRCA2", r) for r in ("controls 0.80", "main", "controls 0.99")] + \
            [(g, r) for g in ("BAP1", "RAD51C") for r in ("main", "controls 0.95")]
    rows = []
    for gene, rule in rules:
        dg = relabel(d, gene, rule)
        for pr, col in ev.PRED.items():
            cal = ev.clinvar_posterior(dg, col)   # ClinVar-only fit: unaffected by the assay labels
            t, p, y, cl = ev.gene_data(dg, gene, col, cal)
            rng = np.random.default_rng(SEED)
            row = dict(gene=gene, rule=rule, predictor=pr, n=len(y), lof_share=float(y.mean()))
            row.update(halvings(p, y, cl, rng))
            row.update(oracle(p, y, cl, rng))
            rows.append(row)
            print(row, flush=True)
    pd.DataFrame(rows).to_csv(RES / "label_sensitivity.csv", index=False)


if __name__ == "__main__":
    main()
