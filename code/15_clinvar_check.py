"""Calls judged against clinical labels (E11), an exploratory check.

For each gene and predictor, per-class thresholds (alpha = 0.1, in-expectation rank) are set on all of the gene's
assay-labelled residues (one variant per residue, outside the ClinVar labelled set) and applied to the gene's ClinVar
variants with >= 2 review stars and a P/LP or B/LB call. These ClinVar variants never enter the thresholds; the Platt
map fitted on them is monotone, so it does not change which variants are called.

Estimand and interval. One ClinVar variant is kept per residue (residue-first weighting, as in the paper's law Q),
so kept variants are treated as independent draws. The four call rates against ClinVar P/LP (y = 1) and B/LB (y = 0)
get one-sided Clopper-Pearson endpoints at level 0.025 each, and the likelihood-ratio bounds are l_TP / u_FP and
u_FN / l_TN, so each ratio bound holds with probability at least 0.95 (Bonferroni over its two endpoints). Results
are given per gene and, descriptively, pooled by summing the kept counts over genes (each gene with its own
thresholds). Counts below 20 per class make a gene's bounds uninformative; they are reported, not interpreted.

Output: results/clinvar_check.csv
"""
import importlib
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ev = importlib.import_module("13_evidence")
RES = pathlib.Path(__file__).resolve().parents[1] / "results"
GENES = ev.GENES


def main():
    d = ev.load()
    rows = []
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        cells = []
        for g in GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            keep = sc.one_per_cluster(cl, np.random.default_rng(ev.SEED + 20))
            thr = sc.mondrian_thresholds(p[keep], y[keep], 0.10)
            cv = d[(d.gene == g) & (d.src == 1) & d[col].notna()]
            kc = sc.one_per_cluster(cv.cluster.values, np.random.default_rng(ev.SEED + 21))
            cv = cv.iloc[kc]
            cells.append((g, sc.decisions(cal.predict(cv[col].values), thr), cv.y_src.values.astype(int)))
        for scope in GENES + ["pooled"]:
            use = cells if scope == "pooled" else [c for c in cells if c[0] == scope]
            dec = np.concatenate([c[1] for c in use]); yc = np.concatenate([c[2] for c in use])
            b = sc.certify(dec, yc, 0.025)
            rows.append(dict(predictor=pr, scope=scope, n_path=int((yc == 1).sum()), n_benign=int((yc == 0).sum()),
                             lof_calls_path=int(((dec == 1) & (yc == 1)).sum()), lof_calls_benign=int(((dec == 1) & (yc == 0)).sum()),
                             func_calls_path=int(((dec == 0) & (yc == 1)).sum()), func_calls_benign=int(((dec == 0) & (yc == 0)).sum()),
                             lr_pos_hat=b["lr_pos_hat"], lr_pos_lower=b["lr_pos_lower"],
                             lr_neg_hat=b["lr_neg_hat"], lr_neg_upper=b["lr_neg_upper"]))
    out = pd.DataFrame(rows)
    out.to_csv(RES / "clinvar_check.csv", index=False)
    print(out.round(3).to_string())


if __name__ == "__main__":
    main()
