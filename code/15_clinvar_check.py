"""Calls judged against clinical labels (E11).

For each gene and predictor, per-class thresholds (alpha = 0.1, in-expectation rank) are set on all of the gene's
assay-labelled residues (one variant per residue, outside the ClinVar labelled set) and applied to the gene's ClinVar
variants with >= 2 review stars and a P/LP or B/LB call. These ClinVar variants never enter the thresholds; the Platt
map fitted on them is monotone, so it does not change which variants are called. Counts are pooled over genes, and
the likelihood ratios of the calls against ClinVar pathogenicity get 95% Clopper-Pearson-based bounds
(each rate at one-sided level 0.025).

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
            dec = sc.decisions(cal.predict(cv[col].values), thr)
            yc = cv.y_src.values.astype(int)
            cells.append((g, dec, yc))
        for scope in ["pooled", "BRCA1"]:
            use = cells if scope == "pooled" else [c for c in cells if c[0] == "BRCA1"]
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
