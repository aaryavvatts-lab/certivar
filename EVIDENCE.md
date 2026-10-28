# Evidence map

Every quantitative claim in the paper is read from `results/paper_numbers.json` (written by `code/11_numbers.py`) or,
where noted, computed directly from a file in `results/`. This table maps each claim to its source.

| Claim | Source (paper_numbers.json key / file) |
|---|---|
| 81% of ClinVar missense records (P/LP, B/LB or VUS) in the five genes are VUS | `clinvar_vus_share` |
| 859 ClinVar >=2-star P/LP or B/LB missense variants, 45.1% P/LP | `source_n`, `source_prev` |
| 14,477 assay-labelled variants outside the ClinVar labelled set; LoF share 6.0-25.9% by gene (variant-level) | `n_target_total`, `target_prev_min/max`, `cohort` |
| 422 dual-labelled variants, 97.4% agreement | `dual_n`, `dual_agree` |
| 729 splice-region variants excluded | `splice_excluded` |
| Control rule vs authors' classes (BRCA1, BAP1 216 of 4,023 opposite, RAD51C 77% left intermediate) | `qc/assay_labels.txt` (written by `code/06_assemble.py`) |
| AUROC 0.946 / 0.840 / 0.779 (AlphaMissense, pooled; residue-bootstrap intervals) | `auroc_AlphaMissense` |
| EM prior error 0.255 mean, 0.702 max; BBSE 0.233; ECE 0.249 / 0.259 / 0.102 | `em_abs_err_*`, `bbse_abs_err_mean`, `ece_by_method` |
| Worst-class (lower of the two mean per-class) coverage 0.896-0.908; >= 0.948 at alpha 0.05, >= 0.792 at 0.2; baselines below 0.55 | `conformal_alpha_0.1` -> `within-gene random split`, `mondrian_by_alpha`, `published_cutoffs` |
| Pooled conformal abstains on 69% (AlphaMissense, BAP1); abstention 35% vs 42% vs 53% | `pooled_abstain_AM_BAP1`, `conformal_alpha_0.1` |
| 17-93 LoF calibration residues per half; 36.4% / 43.6% vs 7.9% / 10.8% of splits below 0.90; abstention 58.4% vs 35.4% | `n_cal_lof_range`, `reliability` |
| VUS shift: residue-weighted LoF share 17.0% (calibration) and 19.8% (VUS); coverage ranges; abstention 81.7% vs 46.4% | `vus_shift_prev`, `vus_shift` |
| Oracle: certified bounds respected 2,178 / 2,181 (99.9%, interval 99.7-100%); all claimed bounds of a draw held in 1,497 / 1,500 draws | `oracle_validity_ci` -> `certified` (`ok`, `n`, `lo`, `hi`, `draws_all_ok_uncond`, `draws_total`) |
| Oracle: uncorrected Clopper-Pearson 99.2% of 3,115; bootstrap 76.5% of 5,287; corrected bootstrap 79.2% of 4,533; point estimates on the conservative side 36.4% of 6,484 | `oracle_validity`, `oracle_validity_ci` |
| Oracle: claimed level held 100% / 99.9% / 91.0% / 79.7% (certified / uncorrected / bootstrap / point) | `oracle_validity` -> `level_held` |
| Point estimates LoF-moderate claimed 98.6%, held 81.1%; bootstrap functional-strong claimed 35.7%, held 30.0%; LoF-moderate certified 13.1% vs 34.5% uncorrected | `oracle_by_level` |
| Certified yields: LoF-supporting 86.6% of draws (74-97% by gene, 25.8% called), functional-supporting 43.6% (61.1% called), LoF-moderate 13.1%, functional-moderate 2.1% | `oracle_by_level`, `evidence_by_gene` |
| Test-fifth splits: 98.7% of 2,091 (certified), 76.8% (bootstrap), 44.7% (point estimates on the conservative side) | `evidence_validity` |
| VUS evidence: LoF-supporting certified in 59.1%, met in 443 / 443; functional-supporting 123 / 125; LoF-moderate 25 / 28 | `vus_evidence` (`n` x `holds`) |
| ClinVar check (one variant per residue, CP endpoints at 0.025 each): LoF lower bounds 4.9-19.6 (BRCA1), 7.6-9.4 (BRCA2), 2.6-3.2 (MSH2); functional upper bounds 0.24-0.37 | `results/clinvar_check.csv` (also `clinvar_check`) |
| Region split: coverage 0.860-0.923; certified bounds respected 88.8% of 1,891, level held 93.0% | `conformal_alpha_0.1` -> `within-gene region split`, `region_bound_ok_all` |
| Other genes: 8 of 15 pairs below 0.85, worst 0.487 | `conformal_alpha_0.1` -> `leave-one-gene-out` |
| Calibration size: abstention 40.5% (150), 39.8% (200), 82.8% (25) | `calibration_size` |
| Strata (variant-weighted): 0.865 / 0.873 nearest the class boundary | `strata` |
| Label sensitivity: worst-class coverage 0.894-0.908 and 2,530 / 2,533 bounds held under the four alternative rules | `label_sensitivity`, `results/label_sensitivity.csv` |
| Simulation | `simulation` |
| Published cut-offs as assay classifiers: functional coverage as low as 0.540 (AlphaMissense classes, BAP1) and 0.445 (REVEL PP3/BP4, MSH2) | `conformal_alpha_0.1` -> `within-gene random split` -> `func_rw_min` |
| Every claimed bound of a draw held in 1,497 of 1,500 draws (weaker than the all-bounds event of Theorem 1) | `oracle_validity_ci` -> `certified` -> `draws_all_ok_uncond`, `draws_total` |
| Uncorrected Clopper-Pearson bounds respected in 99.2% of 3,115 selections | `oracle_validity_ci` -> `Clopper-Pearson, uncorrected` |
| ClinVar-assay agreement 97.4% on 422 variants with both a ≥2-star P/LP or B/LB call and a binary assay label | `dual_agree`; `results/dual_labelled.csv` |
| Monte Carlo standard errors of the per-class coverages 0.002-0.010 (50 splits) | `mondrian_se_range`, `mondrian_rows` |
| Abstention of the REVEL cut-offs (42%) and AlphaMissense classes (11%) | `conformal_alpha_0.1` -> `within-gene random split` -> `abstain_mean` |
| All-variant calibration: 1,446 vs 292 labels, abstention 30.4% vs 35.2%, worst-class 0.887-0.903 (13 of 15 below 0.90) vs 0.897-0.915 (4 of 15) | `data_use` |
| Calibration-share shift: per-class worst-class 0.897-0.925 at shares 0.05-0.5; abstention 37% to 83%; pooled conformal 0.113, weighted 0.53-0.60 (true share) and 0.04 (EM) in the worst pair; 36-176 calibration variants | `prior_shift`, `prior_shift_ncal` |
| Oracle draw level: every claimed bound held in 99.8% (certified), 98.8% (uncorrected), 56.5% (bootstrap), 65.9% (corrected bootstrap) of draws | `oracle_draw_level` |
| Bootstrap diagnostic (fresh 1,500 draws): 57% of failures with zero false calls; failure rate 89% with zero counts, 11% otherwise | `bootstrap_diagnostic` |
| Coverage at the selected alpha: mean 0.82 (LoF) and 0.80 (functional); 10th percentile below 0.69; mean selected alpha 0.18 / 0.24 | `coverage_at_selection` |
| Fixed-sequence vs Bonferroni (G4): LoF-moderate 29.1% vs 12.9%, functional-supporting 60.6% vs 42.5%; 99.7% of 2,705 claims held; no level falsely claimed in more than 0.4% of draws | `oracle_extended` -> `G4|fixed sequence`, `G4|certified`; `results/oracle_extended.csv` |
| G20: LoF-moderate 6.2% (Bonferroni) vs 26.7% (fixed sequence); LoF-supporting 73.8% vs 65.3%; functional levels in increasing order 25.7% (G4) and 0.2% (G20) | `oracle_extended` |
| Variant-weighted truth: certified levels held in 99.9% of claims | `oracle_extended` -> `G4|certified` -> `level_held_vw` |
| Intermediate variants: 1,274 in all, 1,262 outside the ClinVar labelled set; 19.6% called LoF, 36.6% functional | `intermediates` |
| Per-predictor LoF-supporting certification in the oracle: 96% (AlphaMissense), 86% (REVEL), 78% (ESM-2) | `results/oracle.csv` (certified, PS3-supporting, mean of `achieved` by predictor) |
| VUS LoF-moderate 25 of 28 (Clopper-Pearson 95% interval 72-98%) | `vus_evidence`; interval computed with scipy.stats.beta |
| Runtime about 5 ms per gene (20-level grid, BAP1, laptop CPU) | measured with `time.perf_counter` on `certivar.mondrian_thresholds` + `certify` |
