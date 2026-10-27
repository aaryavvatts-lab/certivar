# CertiVar: certified evidence strengths for missense variant-effect predictors

Code and data for the RECOMB 2027 submission *"CertiVar: Certified Evidence Strengths for Missense Variant-Effect
Predictors under Unknown Prevalence of Damaging Variants"* by Aaryav Sharma (Eastview High School, Apple Valley, MN).

CertiVar turns any variant-effect predictor score into three-way calls (loss of function, functional, abstain) for a
gene with saturation genome editing (SGE) or deep mutational scanning (DMS) data:

1. **Per-class thresholds** (class-conditional, "Mondrian" split conformal) set from assay labels with one variant per
   residue, in an in-expectation version and a PAC version (miss rate at most alpha with probability at least 1 - delta).
2. **Certified likelihood ratios**: exact Clopper-Pearson bounds on the likelihood ratio of the calls, valid
   simultaneously over a grid of thresholds, so the calls can be reported at an ACMG/AMP OddsPath evidence level
   (Brnich et al. 2020) with a finite-sample guarantee.
3. **No prior needed**: both guarantees condition on the true class, so they hold whatever the share of damaging
   variants among the variants being called.

The benchmark covers BRCA1, BRCA2, BAP1, RAD51C (SGE) and MSH2 (DMS), with AlphaMissense, REVEL and ESM-2.

## Reproduce everything

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
make reproduce
```

Network access is needed for the download and annotation steps; ESM-2 scoring downloads about 2.6 GB of weights and uses
Apple MPS or CUDA when present. Every number in the paper is read from `results/paper_numbers.json`.

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
| Certified evidence, reliability, VUS shift | `code/13_evidence.py` | `results/evidence.csv`, `reliability.csv`, `vus_shift.csv` |
| Numbers quoted in the paper | `code/11_numbers.py` | `results/paper_numbers.json` |
| Figures and tables | `code/09_figures.py`, `code/12_tables.py` | `results/figures/*.pdf` |

The method is in `code/certivar.py` (NumPy, SciPy and scikit-learn only).

## Using the method on your own scores

```python
import certivar as cv
cal = cv.Platt().fit(clinvar_scores, clinvar_labels)        # optional: common scale
thr = cv.mondrian_thresholds(cal.predict(s_A), y_A, alpha=0.1, delta=0.1)   # one variant per residue in part A
calls = cv.decisions(cal.predict(s_new), thr)                 # 1 LoF, 0 functional, -1 / -2 no call
bounds = cv.certify(cv.decisions(cal.predict(s_B), thr), y_B, level=0.05 / 2)  # part B, one variant per residue
```

## Intended use

The calls and certified levels are research results. Their certified likelihood ratios refer to the assay outcome;
using them as clinical PP3/BP4 evidence requires the assay itself to be clinically validated and gene-specific review.

## Licence

MIT for the code. Derived data follow the licences and terms of their original sources (ClinVar, MaveDB score sets,
gnomAD, Ensembl VEP and the predictor score providers).
