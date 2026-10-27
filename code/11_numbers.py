"""Collect every number quoted in the paper into results/paper_numbers.json.

The manuscript cites only values that appear in this file (checked by qc/check_numbers.py).
"""
import json
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
A = 0.10


def r3(x):
    return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 3)


def main():
    out = {}
    co = pd.read_csv(RES / "cohort.csv")
    out["cohort"] = co.round(4).to_dict(orient="records")
    out["n_target_total"] = int(co.target_n.sum())
    out["n_assayed_total"] = int(co.assayed.sum())
    out["target_prev_min"], out["target_prev_max"] = r3(co.target_prev.min()), r3(co.target_prev.max())
    out["vus_prev_min"], out["vus_prev_max"] = r3(co.vus_prev.min()), r3(co.vus_prev.max())

    disc = pd.read_csv(RES / "discrimination.csv")
    src = disc[disc.subset == "ClinVar>=2*"]
    out["source_prev"] = r3(src.prev.iloc[0])
    out["source_n"] = int(src.n.iloc[0])
    for pr in ["AlphaMissense", "REVEL", "ESM-2"]:
        s = disc[(disc.predictor == pr) & (disc.gene == "ALL")].set_index("subset")
        out[f"auroc_{pr}"] = {k: [r3(s.loc[k, "auroc"]), r3(s.loc[k, "lo"]), r3(s.loc[k, "hi"])] for k in s.index}
    dual = pd.read_csv(RES / "dual_labelled.csv")
    out["dual_n"] = int(dual.n.sum())
    out["dual_agree"] = r3((dual.n * dual.agree).sum() / dual.n.sum())

    cal = pd.read_csv(RES / "calibration.csv")
    out["ece_by_method"] = cal[cal.predictor.isin(["AlphaMissense", "REVEL", "ESM-2"])].groupby("method").ece.mean().round(3).to_dict()
    pri = pd.read_csv(RES / "prior_estimates.csv")
    main3 = pri[pri.predictor.isin(["AlphaMissense", "REVEL", "ESM-2"])]
    out["em_abs_err_mean"] = r3((main3.em - main3.observed).abs().mean())
    out["em_abs_err_max"] = r3((main3.em - main3.observed).abs().max())
    out["bbse_abs_err_mean"] = r3((main3.bbse - main3.observed).abs().mean())
    out["prior_table"] = main3[["predictor", "gene", "observed", "em", "bbse"]].round(3).to_dict(orient="records")

    conf = pd.read_csv(RES / "conformal.csv")
    c = conf[conf.alpha == A].copy()
    c["worst_rw"] = c[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    summ = {}
    for (des, m), g in c.groupby(["design", "method"]):
        summ.setdefault(des, {})[m] = {
            "worst_min": r3(g.worst_rw.min()), "worst_max": r3(g.worst_rw.max()), "worst_mean": r3(g.worst_rw.mean()),
            "lof_rw_min": r3(g.coverage_lof_rw.min()), "lof_rw_max": r3(g.coverage_lof_rw.max()),
            "func_rw_min": r3(g.coverage_func_rw.min()), "func_rw_max": r3(g.coverage_func_rw.max()),
            "abstain_min": r3(g.abstain.min()), "abstain_max": r3(g.abstain.max()),
            "abstain_mean": r3(g.abstain.mean()), "n_cells": int(len(g)),
            "n_cells_below_085": int((g.worst_rw < 0.85).sum())}
    out["conformal_alpha_0.1"] = summ
    mond = conf[(conf.method == "Mondrian conformal") & (conf.design == "within-gene random split")]
    out["mondrian_by_alpha"] = {str(a): {"lof_rw_min": r3(g.coverage_lof_rw.min()),
                                         "func_rw_min": r3(g.coverage_func_rw.min()),
                                         "abstain_min": r3(g.abstain.min()), "abstain_max": r3(g.abstain.max())}
                                for a, g in mond.groupby("alpha")}
    fixed = c[(c.design == "within-gene random split") & c.method.isin(["AlphaMissense classes", "REVEL PP3/BP4 supporting"])]
    out["published_cutoffs"] = fixed[["method", "gene", "coverage_func", "coverage_lof", "abstain"]].round(3).to_dict(orient="records")

    st = pd.read_csv(RES / "conformal_strata.csv")
    st = st[st.method == "Mondrian conformal"]
    out["strata"] = st.groupby("stratum")[["coverage_lof", "coverage_func", "abstain"]].mean().round(3).to_dict(orient="index")
    cs = pd.read_csv(RES / "calibration_size.csv")
    out["calibration_size"] = cs.groupby("n_cal_residues")[["coverage_lof_rw", "coverage_func_rw", "abstain"]].mean().round(3).to_dict(orient="index")
    out["application"] = pd.read_csv(RES / "application.csv").to_dict(orient="records")

    sim = pd.read_csv(RES / "simulation.csv")
    sg = sim.groupby(["kind", "method", "pi_t"]).mean(numeric_only=True).round(3).reset_index()
    out["simulation"] = sg.to_dict(orient="records")
    out["splice_excluded"] = int(pd.read_csv(ROOT / "data" / "processed" / "variants.csv", low_memory=False)
                                 .pipe(lambda v: (v.splice_region & v.score.notna()).sum()))
    for name in ("oracle", "evidence"):
        ev = pd.read_csv(RES / f"{name}.csv")
        acc = ev[ev.achieved]
        out[f"{name}_validity"] = {sel: {"n": int(len(g)), "bound_ok": r3(g.bound_ok.astype(float).mean()),
                                        "level_held": r3(g.holds.astype(float).mean())}
                                   for sel, g in acc.groupby("selection")}
        out[f"{name}_by_level"] = {f"{sel}|{lvl}": {"achieved": r3(g.achieved.mean()),
                                                    "holds": r3(g[g.achieved].holds.astype(float).mean()) if g.achieved.any() else None,
                                                    "n": int(g.achieved.sum()),
                                                    "yield": r3(g[g.achieved].ref_yield.mean()) if g.achieved.any() else None}
                                   for (sel, lvl), g in ev.groupby(["selection", "level"])}
    ev = pd.read_csv(RES / "oracle.csv")
    out["evidence_by_gene"] = ev[ev.selection == "certified"].groupby(["level", "gene"]).achieved.mean().round(2).unstack().to_dict(orient="index")
    reg = pd.read_csv(RES / "region_evidence.csv")
    out["region_evidence"] = {lvl: {"achieved": r3(g.achieved.mean()), "n": int(g.achieved.sum()),
                                    "holds": r3(g[g.achieved].holds.astype(float).mean()) if g.achieved.any() else None,
                                    "bound_ok": r3(g[g.achieved].bound_ok.astype(float).mean()) if g.achieved.any() else None}
                              for lvl, g in reg.groupby("level")}
    ra = reg[reg.achieved]
    out["region_bound_ok_all"] = {"n": int(len(ra)), "bound_ok": r3(ra.bound_ok.astype(float).mean()),
                                  "holds": r3(ra.holds.astype(float).mean())}
    rel = pd.read_csv(RES / "reliability.csv")
    rel["worst"] = rel[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    out["reliability"] = {t: {"frac_worst_below": r3((g.worst < 0.9).mean()), "frac_lof_below": r3((g.coverage_lof_rw < 0.9).mean()),
                              "frac_func_below": r3((g.coverage_func_rw < 0.9).mean()), "mean_worst": r3(g.worst.mean()),
                              "abstain": r3(g.abstain.mean()), "lof_sd": r3(g.groupby(["predictor", "gene"]).coverage_lof_rw.std().mean())}
                          for t, g in rel.groupby("thresholds")}
    vs = pd.read_csv(RES / "vus_shift.csv")
    vc = vs[vs.part == "coverage"]
    agg = vc.groupby(["thresholds", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw", "abstain"]].mean()
    out["vus_shift"] = {t: {"lof_min": r3(g.coverage_lof_rw.min()), "lof_max": r3(g.coverage_lof_rw.max()),
                            "func_min": r3(g.coverage_func_rw.min()), "func_max": r3(g.coverage_func_rw.max()),
                            "lof_mean": r3(g.coverage_lof_rw.mean()), "func_mean": r3(g.coverage_func_rw.mean()),
                            "abstain_mean": r3(g.abstain.mean()), "n_pairs_lof_below_09": int((g.coverage_lof_rw < 0.9).sum()),
                            "n_pairs_func_below_09": int((g.coverage_func_rw < 0.9).sum())}
                        for t, g in agg.groupby(level=0)}
    out["vus_shift_prev"] = {"cal_prev": r3(vc.cal_prev.mean()), "vus_prev": r3(vc.vus_prev.mean())}
    ve = vs[vs.part == "evidence"]
    out["vus_evidence"] = {lvl: {"achieved": r3(g.achieved.mean()), "n": int(g.achieved.sum()),
                                 "holds": r3(g[g.achieved].holds.astype(float).mean()) if g.achieved.any() else None,
                                 "bound_ok": r3(g[g.achieved].bound_ok.astype(float).mean()) if g.achieved.any() else None}
                           for lvl, g in ve.groupby("level")}
    # oracle validity with a 95% interval that resamples whole draws (selections within a draw are dependent)
    acc = pd.read_csv(RES / "oracle.csv").query("achieved")
    rng = np.random.default_rng(0)
    out["oracle_validity_ci"] = {}
    for sel, g in acc.groupby("selection"):
        per = g.groupby(["predictor", "gene", "split"]).bound_ok.agg(["sum", "count"])
        boot = []
        for _ in range(2000):
            k = per.iloc[rng.integers(0, len(per), len(per))]
            boot.append(k["sum"].sum() / k["count"].sum())
        out["oracle_validity_ci"][sel] = {"ok": int(g.bound_ok.sum()), "n": int(len(g)), "draws": int(len(per)),
                                          "draws_all_ok": r3((per["sum"] == per["count"]).mean()),
                                          "lo": r3(np.percentile(boot, 2.5)), "hi": r3(np.percentile(boot, 97.5))}
    ls = pd.read_csv(RES / "label_sensitivity.csv")
    ls["worst"] = ls[["lof_cov", "func_cov"]].min(axis=1)
    alt = ls[ls.rule != "main"]
    out["label_sensitivity"] = {
        "alt_worst_min": r3(alt.worst.min()), "alt_worst_max": r3(alt.worst.max()),
        "alt_claims": int(alt.oracle_claims.sum()), "alt_ok": int(alt.oracle_bound_ok.sum()),
        "by_rule": ls.groupby(["gene", "rule"]).agg(n=("n", "first"), lof_share=("lof_share", "first"),
                                                     worst_min=("worst", "min"), abstain_min=("abstain", "min"),
                                                     abstain_max=("abstain", "max"),
                                                     lof_sup_min=("share_lof_supporting", "min"),
                                                     func_sup_min=("share_func_supporting", "min"),
                                                     func_sup_max=("share_func_supporting", "max"),
                                                     claims=("oracle_claims", "sum"), ok=("oracle_bound_ok", "sum"))
                     .round(3).reset_index().to_dict(orient="records")}
    out["clinvar_check"] = pd.read_csv(RES / "clinvar_check.csv").round(3).to_dict(orient="records")
    (RES / "paper_numbers.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: out[k] for k in ["source_prev", "target_prev_min", "target_prev_max", "em_abs_err_mean",
                                           "em_abs_err_max", "dual_agree", "ece_by_method"]}, indent=1))
    print(json.dumps(summ["within-gene random split"], indent=1))


if __name__ == "__main__":
    main()
