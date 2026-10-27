# CertiVar

Code and data for the RECOMB 2027 submission

> **CertiVar: Prevalence-Free Certified Error Rates and Likelihood Ratios for Missense Variant-Effect Predictor Calls**
> Aaryav Sharma (Eastview High School, Apple Valley, MN, USA)

Release `v1.0.0` is the version that matches the submitted paper. The supplementary material (full proofs,
simulation, per-gene tables, strata, label sensitivity, the ClinVar check and implementation details; items S1-S5,
Figs. S1-S2, Tables S1-S4) is in
[`supplement/CertiVar_supplement.pdf`](supplement/CertiVar_supplement.pdf).

## What CertiVar does

CertiVar turns any variant-effect predictor score into three-way calls (loss of function, functional, abstain) for a
gene with saturation genome editing (SGE) or deep mutational scanning (DMS) data:

1. **Per-class thresholds.** Class-conditional ("Mondrian") split-conformal thresholds, set from assay labels with one
   variant per residue, bound the rate at which each class is mis-called, either in expectation or with probability at
   least 1 - delta over the calibration data (PAC version).
2. **Certified likelihood ratios.** One-sided Clopper-Pearson endpoints `l_TP`, `u_FP`, `u_FN`, `l_TN` give
   `L+ = l_TP / u_FP` and `U- = u_FN / l_TN` (`U- = inf` when `l_TN = 0`), valid simultaneously over a grid of
   thresholds, for the call-averaged likelihood ratio of each call type **against the assay labels**. These are
   compared with the numerical OddsPath thresholds (2.1 / 4.3 / 18.7 and 0.48 / 0.23 / 0.053; Tavtigian et al. 2018,
   Brnich et al. 2020).
3. **No prevalence estimate.** Both certificates are class-conditional, so they hold whatever the share of damaging
   variants among the variants being called, provided that within each class those variants resemble the assayed ones.

## What it certifies, and what it does not

CertiVar provides finite-sample certificates for class-conditional call rates and call-averaged likelihood ratios with
respect to functional-assay labels, under the sampling assumptions stated in the paper. Numerical alignment with
OddsPath thresholds is not, by itself, validation of clinical PP3/BP4/PS3/BS3 use.

- The certified likelihood ratios are defined **with respect to the assay outcome** and averaged over all variants
  that receive a call. They state how strongly a call predicts the functional-assay result.
- They are **not** clinical pathogenicity evidence by themselves, and they are **not** the local, score-specific
  likelihood ratios that ClinGen uses for PP3/BP4. They become evidence about pathogenicity only to the extent that the
  assay itself is clinically validated, and they must not be combined with the same assay's own PS3/BS3 evidence for a
  variant.
- The guarantees are proved for residues exchangeable with the calibrated ones (Assumptions 1-2 of the paper). Calls
  in other protein regions or other genes are supported only by the empirical checks in the paper.
- The check against ClinVar labels in the paper is exploratory: it is reported per gene (Table S4) and is not
  independent clinical validation.
- The calls are research results, not clinical classifications.

The benchmark covers BRCA1, BRCA2, BAP1, RAD51C (SGE) and MSH2 (DMS) with AlphaMissense, REVEL and ESM-2.

## Reproduce everything

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
make reproduce
```

Tested with Python 3.13 on macOS; exact package versions are pinned in `requirements.txt`. Network access is needed for
the download and annotation steps; ESM-2 scoring downloads about 2.6 GB of weights and uses Apple MPS or CUDA when
present. All random draws use fixed seeds (2027 plus a fixed offset per experiment). Every number quoted in the paper is
read from `results/paper_numbers.json`, and `EVIDENCE.md` maps each claim to its key. SHA-256 checksums of the input
data are in `data/MANIFEST.txt`; data sources and access dates are in `SOURCES.md`.

`make test` runs the boundary-case unit tests of the method (zero denominators, empty classes, Platt slope sign,
disjoint calls, coverage of the certified bounds). To rebuild only the results, figures and tables from the processed
data already in this repository:

```bash
make test results figures
```

| Step | Script | Output |
|---|---|---|
| MaveDB score sets and class thresholds | `code/01_download.py` | `data/raw/mavedb_*.csv`, `data/processed/assays.json` |
| ClinVar and assay variant tables | `code/02_variants.py` | `data/processed/clinvar.csv`, `assay_snvs.csv` |
| Ensembl VEP annotation | `code/03_vep.py` | `data/processed/vep.csv` |
| ESM-2 log-likelihood ratios | `code/04_esm.py` | `data/processed/esm2_llr.csv` |
| gnomAD v4 allele counts | `code/05_gnomad.py` | `data/processed/gnomad.csv` |
| Join and assay labels | `code/06_assemble.py` | `data/processed/variants.csv` |
| Simulation | `code/07_simulation.py` | `results/simulation.csv` |
| Coverage, calibration and transfer experiments | `code/08_experiments.py` | `results/*.csv` |
| Certified likelihood ratios (oracle and data splits), reliability, VUS shift, region split | `code/13_evidence.py` | `results/oracle.csv`, `evidence.csv`, `reliability.csv`, `vus_shift.csv`, `region_evidence.csv` |
| Sensitivity to the assay-label rule | `code/14_label_sensitivity.py` | `results/label_sensitivity.csv` |
| Calls judged against ClinVar pathogenicity labels | `code/15_clinvar_check.py` | `results/clinvar_check.csv` |
| Numbers quoted in the paper | `code/11_numbers.py` | `results/paper_numbers.json` |
| Figures and tables | `code/09_figures.py`, `code/12_tables.py` | `results/figures/*.pdf`, LaTeX tables |
| Short-author bibliography for the paper | `code/16_short_bib.py` | `paper/refs_short.bib` |

The method itself is in `code/certivar.py` (NumPy, SciPy and scikit-learn only).

## Using the method on your own scores

Run from the `code/` directory (or add it to `PYTHONPATH`):

```python
import numpy as np
import certivar as cv

rng = np.random.default_rng(0)
cal = cv.Platt().fit(clinvar_scores, clinvar_labels)             # optional: common scale (monotone, calls unchanged)

# part A: one variant per residue -> per-class thresholds (PAC version with delta; omit delta for in-expectation)
keep_a = cv.one_per_cluster(residues_A, rng)
thr = cv.mondrian_thresholds(cal.predict(s_A[keep_a]), y_A[keep_a], alpha=0.1, delta=0.1)

calls = cv.decisions(cal.predict(s_new), thr)                     # 1 LoF, 0 functional, -1 / -2 no call

# part B (residues disjoint from A): one variant per residue -> certified bounds on LR+ and LR-
keep_b = cv.one_per_cluster(residues_B, rng)
bounds = cv.certify(cv.decisions(cal.predict(s_B[keep_b]), thr), y_B[keep_b], level=0.05 / 2)
print(bounds["lr_pos_lower"], bounds["lr_neg_upper"])
```

`level` is the one-sided level of each Clopper-Pearson bound. For a single threshold, `level = delta / 2` gives each
bound with probability at least 1 - delta; when thresholds are searched over a grid `G` (Algorithm 1 of the paper), use
`level = delta / (2 * len(G))`, as `code/13_evidence.py` does.

## Licence

MIT for the code. Derived data follow the licences and terms of their original sources (ClinVar, MaveDB score sets,
gnomAD, Ensembl VEP and the predictor score providers).
