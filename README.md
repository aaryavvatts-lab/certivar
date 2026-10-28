# CertiVar

Code, data and results for the RECOMB 2027 submission

> **CertiVar: Certified Likelihood Ratios for Variant-Effect Predictor Calls, Calibrated on Saturation Mutagenesis Data**
> Aaryav Sharma (Eastview High School, Apple Valley, MN, USA)

Release `v1.1.0` matches the submitted paper, [`manuscript/CertiVar_RECOMB2027.pdf`](manuscript/CertiVar_RECOMB2027.pdf):
ten pages of main text, the references, and an appendix (proofs, simulation, per-gene tables, label sensitivity, the
ClinVar check, the VUS counts, implementation and experimental designs). The LaTeX source is in `paper/`.

## What CertiVar does

CertiVar turns a variant-effect predictor score (oriented so that higher means more damaging) into three-way calls (loss
of function, functional, abstain) for a gene with saturation genome editing (SGE) or deep mutational scanning (DMS) data.

1. **Two cut-offs from per-class thresholds.** Class-conditional ("Mondrian") split-conformal thresholds, set from assay
   labels with one variant per residue, bound the rate at which each class is mis-called, in expectation or with
   probability at least 1 - delta_P over the calibration data (PAC version). On the score axis they give two cut-offs:
   below both, a functional call; above both, a LoF call; in between, no call.
2. **Certified likelihood ratios.** One-sided Clopper-Pearson bounds on the four call rates give a lower bound
   `L+ = lower(P(call LoF | LoF)) / upper(P(call LoF | functional))` and an upper bound
   `U- = upper(P(call functional | LoF)) / lower(P(call functional | functional))` (`U- = inf` when the denominator's
   lower bound is 0) on the call-averaged likelihood ratio of each call type **against the assay labels**. They stay
   valid after a threshold is selected from a grid, with a Bonferroni correction (Theorem 1) or with fixed-sequence
   testing, per evidence level or simultaneously over levels (Theorem 2, `certivar.fixed_sequence_select`). The levels
   are reported on the OddsPath scale (2.1 / 4.3 / 18.7 / 350 and 0.48 / 0.23 / 0.053; Tavtigian et al. 2018, Brnich
   et al. 2020).
3. **No prevalence estimate, and a check of when that holds.** Both certificates are class-conditional, so they hold
   whatever the share of damaging variants among the variants being called, provided that within each class those
   variants have the same score distribution as the assayed ones. `certivar.mixture_check` tests a necessary condition
   for this on unlabelled target scores: whether they are a mixture of the two calibration classes.

## What it certifies, and what it does not

- The certified likelihood ratios are defined **with respect to the assay outcome**, weight residues equally and are
  averaged over all variants that receive a call. They state how strongly a call predicts the functional-assay result.
- They are **not** clinical pathogenicity evidence by themselves, and they are **not** the local, score-specific
  likelihood ratios that ClinGen uses for PP3/BP4. They become evidence about pathogenicity only to the extent that the
  assay itself is clinically validated, and they must not be combined with the same assay's own PS3/BS3 evidence for a
  variant.
- The guarantees are proved for residues exchangeable with the calibrated ones (Assumptions 1-2 of the paper). Calls in
  other protein regions or other genes are supported only by the empirical checks in the paper; the mixture check
  rejects most such transfers, but passing it does not prove that the guarantees hold.
- The check against ClinVar labels is exploratory (Table S4) and is not independent clinical validation.
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
read from `results/paper_numbers.json` or from a generated table; `docs/EVIDENCE.md` maps each claim to its source.
SHA-256 checksums of the input data are in `data/MANIFEST.txt`; data sources, access dates and checked references are in
`docs/SOURCES.md`.

Pinned inputs: ClinVar `variant_summary.txt.gz` of 26 Oct 2026 (SHA-256 `657b47d0…a35703`; NCBI serves only the latest
release at that address, so `make data` stops if a fresh download has a different checksum, and the extract used in the
paper is committed as `data/processed/clinvar.csv`); MaveDB score sets urn:mavedb:00000097-0-2, 00001279-a-1,
00000662-0-1, 00000673-0-1 and 00000050-a-1 (raw files committed in `data/raw/`); Ensembl VEP REST on GRCh38
(annotations retrieved 26 Oct 2026, committed as `data/processed/vep.csv`; the service reported release 116 on
27 Oct 2026); ESM-2 checkpoint `esm2_t33_650M_UR50D` (fair-esm 2.0.0, scores committed as
`data/processed/esm2_llr.csv`).

`make test` runs the unit tests of the method (zero denominators, empty classes, Platt slope sign, disjoint calls,
coverage of the certified bounds, fixed-sequence stopping and level control, residue-weighted thresholds, Katz bounds,
and the level and power of the mixture check). To rebuild only the results, figures and tables from the processed data
in this repository (on a laptop CPU, `08_experiments.py`, `17_extensions.py` and `18_robustness.py` take about 4 to 8
minutes each):

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
| Certified likelihood ratios on data splits, reliability, VUS shift, region split | `code/13_evidence.py` | `results/oracle.csv`, `evidence.csv`, `reliability.csv`, `vus_shift.csv`, `region_evidence.csv` |
| Sensitivity to the assay-label rule | `code/14_label_sensitivity.py` | `results/label_sensitivity.csv` |
| Calls judged against ClinVar pathogenicity labels | `code/15_clinvar_check.py` | `results/clinvar_check.csv` |
| Oracle benchmark of every selection rule (Bonferroni, fixed sequence, Clopper-Pearson, Katz, bootstrap variants, point estimate) on 4- and 20-level grids; calibration-share shift; all-variant calibration | `code/17_extensions.py` | `results/oracle_extended.csv`, `prior_shift.csv`, `data_use.csv` |
| Population coverage of in-expectation and PAC thresholds; residue-weighted calibration; mixture check; ClinVar VUS use case | `code/18_robustness.py` | `results/coverage_oracle.csv`, `calibration_weighting.csv`, `mixture_check.csv`, `use_case.csv` |
| Numbers quoted in the paper | `code/11_numbers.py` | `results/paper_numbers.json` |
| Figures and tables | `code/09_figures.py`, `code/12_tables.py` | `results/figures/*.pdf`, `paper/tables/*.tex` |
| Short-author bibliography for the paper | `code/16_short_bib.py` | `paper/refs_short.bib` |

The method itself is in `code/certivar.py` (NumPy, SciPy and scikit-learn only). The paper builds with
`tectonic paper.tex` (or `latexmk -pdf paper.tex`), run from `paper/`.

## Using the method on your own scores

Run from the `code/` directory (or add it to `PYTHONPATH`):

```python
import numpy as np
import certivar as cv

rng = np.random.default_rng(0)
cal = cv.Platt().fit(clinvar_scores, clinvar_labels)             # optional: common scale (monotone, calls unchanged)

# part A: one variant per residue -> per-class thresholds (PAC version with delta; omit delta for in-expectation)
keep_a = cv.one_per_cluster(residues_A, rng)
q_a, y_a = cal.predict(s_A[keep_a]), y_A[keep_a]
thr = cv.mondrian_thresholds(q_a, y_a, alpha=0.1, delta=0.1)

# check the targets before calling them (one variant per residue); p < 0.05: do not report calls
keep_t = cv.one_per_cluster(residues_new, rng)
stat, pi_fit, p_value = cv.mixture_check(q_a[y_a == 1], q_a[y_a == 0], cal.predict(s_new[keep_t]))

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
