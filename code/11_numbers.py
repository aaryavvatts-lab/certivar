"""Collect every number quoted in the paper into results/paper_numbers.json.

Numbers in the manuscript come from this file; a few counts derived directly from results/*.csv (e.g. the per-gene
ClinVar bounds) are listed in EVIDENCE.md with their source file.
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


def draw_level(g, col="bound_ok"):
    """Share of draws in which every claim of the draw satisfied `col` (draws without claims count as successes)."""
    ok = g.groupby(["predictor", "gene", "draw"]).apply(
        lambda t: bool((~t.achieved | t[col].fillna(False).astype(bool)).all()))
    return r3(ok.mean())


def robustness_numbers(vv, vv_src):
    """Numbers from the extended oracle (17_extensions.py) and the robustness experiments (18_robustness.py)."""
    out = {}
    # residues per gene in the evaluation set, and the residue-weighted LoF share (law Q)
    lab = vv[vv.assay_label.isin([0, 1]) & ~vv_src & vv.score.notna() & ~vv.splice_region.fillna(False).astype(bool)]
    res = {}
    for g, s in lab.groupby("gene"):
        w = 1.0 / s.groupby("cluster").cluster.transform("size")
        res[g] = {"variants": int(len(s)), "residues": int(s.cluster.nunique()),
                  "lof_share_rw": r3((w * (s.assay_label == 1)).sum() / w.sum())}
    out["residues"] = res
    # extended oracle, grid G4: every rule on the same draws
    ox = pd.read_csv(RES / "oracle_baselines.csv")
    o4 = ox[ox.grid == "G4"]
    out["oracle_g4"] = {}
    for sel, g in o4.groupby("selection"):
        a = g[g.achieved]
        out["oracle_g4"][sel] = {
            "draws": int(g.groupby(["predictor", "gene", "draw"]).ngroups), "claims": int(len(a)),
            "draws_all_bounds_ok": draw_level(g) if a.bound_ok.notna().any() else None,
            "draws_all_levels_held": draw_level(g, "holds"),
            "claims_level_held": r3(a.holds.astype(float).mean()),
            "claims_bound_ok": r3(a.bound_ok.astype(float).mean()) if a.bound_ok.notna().any() else None,
            "claims_level_held_vw": r3(a.vw_holds.astype(float).mean()),
            "achieved": {k: r3(x) for k, x in g.groupby("level").achieved.mean().items()},
            "yield": {k: r3(x) for k, x in a.groupby("level").ref_yield.mean().items()}}
    out["oracle_g20"] = {sel: {"achieved": {k: r3(x) for k, x in g.groupby("level").achieved.mean().items()},
                               "claims_level_held": r3(g[g.achieved].holds.astype(float).mean()),
                               "draws_all_levels_held": draw_level(g, "holds")}
                         for sel, g in ox[ox.grid == "G20"].groupby("selection")}
    # population (conditional) per-class coverage of in-expectation and PAC thresholds
    co = pd.read_csv(RES / "coverage_oracle.csv")
    out["coverage_oracle"] = {t: {"frac_lof_below": r3((g.coverage_lof_rw < 0.9).mean()),
                                  "frac_func_below": r3((g.coverage_func_rw < 0.9).mean()),
                                  "frac_either_below": r3(((g.coverage_lof_rw < 0.9) | (g.coverage_func_rw < 0.9)).mean()),
                                  "max_pair_frac_func_below": r3(g.groupby(["predictor", "gene"]).coverage_func_rw
                                                                 .apply(lambda x: (x < 0.9).mean()).max()),
                                  "max_pair_frac_lof_below": r3(g.groupby(["predictor", "gene"]).coverage_lof_rw
                                                                .apply(lambda x: (x < 0.9).mean()).max()),
                                  "mean_lof": r3(g.coverage_lof_rw.mean()), "mean_func": r3(g.coverage_func_rw.mean()),
                                  "abstain": r3(g.abstain.mean()), "n": int(len(g))}
                              for t, g in co.groupby("thresholds")}
    # calibration on all variants with and without residue weights
    cw = pd.read_csv(RES / "calibration_weighting.csv")
    g = cw.groupby(["calibration", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw", "abstain"]].mean()
    g["worst"] = g[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    out["calibration_weighting"] = {c: {"worst_min": r3(v.worst.min()), "worst_max": r3(v.worst.max()),
                                        "n_below_090": int((v.worst < 0.90).sum()), "abstain_mean": r3(v.abstain.mean())}
                                    for c, v in g.groupby(level=0)}
    # mixture-consistency check of Assumption 2
    mc = pd.read_csv(RES / "mixture_check.csv")
    mc["reject"] = mc.pval < 0.05
    if "coverage_lof" in mc:
        mc["worst"] = mc[["coverage_lof", "coverage_func"]].min(axis=1)
    out["mixture_check"] = {}
    for st, v in mc.groupby("setting"):
        d = {"n": int(len(v)), "reject": r3(v.reject.mean()),
             "reject_by_predictor": {k: r3(x) for k, x in v.groupby("predictor").reject.mean().items()}}
        if v.worst.notna().any():
            d["worst_mean"] = r3(v.worst.mean())
            d["worst_mean_rejected"] = r3(v[v.reject].worst.mean()) if v.reject.any() else None
            d["worst_mean_accepted"] = r3(v[~v.reject].worst.mean()) if (~v.reject).any() else None
            d["below_085_rejected"] = r3((v[v.reject].worst < 0.85).mean()) if v.reject.any() else None
            d["below_085_accepted"] = r3((v[~v.reject].worst < 0.85).mean()) if (~v.reject).any() else None
            d["pi_abs_err"] = r3((v.pi_fit - v.true_share).abs().mean())
        if "abstain" in v and v.abstain.notna().any():
            d["abstain_mean"] = r3(v.abstain.mean())
        d["by_gene"] = {k: r3(x) for k, x in v.groupby("gene").reject.mean().items()}
        out["mixture_check"][st] = d
    # use case: ClinVar VUS by assay status and certified calls on unassayed VUS at assayed residues
    uc = pd.read_csv(RES / "use_case.csv")
    lv = [c for c in uc.columns if c.startswith(("PS3", "BS3"))]
    out["use_case"] = uc.groupby(["predictor", "gene"]).agg(
        vus=("vus", "first"), vus_assayed=("vus_assayed", "first"),
        at_assayed=("vus_unassayed_at_assayed_residues", "first"), elsewhere=("vus_unassayed_elsewhere", "first"),
        n_target=("n_target", "first"), **{k: (k, "mean") for k in lv}).round(3).reset_index().to_dict(orient="records")
    return out


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
    out["published_cutoffs"] = fixed[["method", "predictor", "gene", "coverage_func_rw", "coverage_lof_rw", "coverage_func",
                                      "coverage_lof", "abstain"]].round(3).to_dict(orient="records")

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
    # residue-weighted shares for both (calibration keeps one variant per residue); variant-level VUS share kept for reference
    out["vus_shift_prev"] = {"cal_prev": r3(vc.cal_prev.mean()), "vus_prev_rw": r3(vc.vus_prev_rw.mean()),
                             "vus_prev_variant_level": r3(vc.vus_prev.mean())}
    ve = vs[vs.part == "evidence"]
    out["vus_evidence"] = {lvl: {"achieved": r3(g.achieved.mean()), "n": int(g.achieved.sum()),
                                 "holds": r3(g[g.achieved].holds.astype(float).mean()) if g.achieved.any() else None,
                                 "bound_ok": r3(g[g.achieved].bound_ok.astype(float).mean()) if g.achieved.any() else None}
                           for lvl, g in ve.groupby("level")}
    # oracle validity with a 95% interval that resamples whole draws (selections within a draw are dependent)
    orc = pd.read_csv(RES / "oracle.csv")
    acc = orc.query("achieved")
    rng = np.random.default_rng(0)
    out["oracle_validity_ci"] = {}
    for sel, g in acc.groupby("selection"):
        per = g.groupby(["predictor", "gene", "split"]).bound_ok.agg(["sum", "count"])
        boot = []
        for _ in range(2000):
            k = per.iloc[rng.integers(0, len(per), len(per))]
            boot.append(k["sum"].sum() / k["count"].sum())
        allg = orc[orc.selection == sel]
        uncond = allg.groupby(["predictor", "gene", "split"]).apply(
            lambda t: bool((~t.achieved | t.bound_ok.fillna(False).astype(bool)).all()))
        out["oracle_validity_ci"][sel] = {"ok": int(g.bound_ok.sum()), "n": int(len(g)), "draws": int(len(per)),
                                          "draws_total": int(len(uncond)), "draws_all_ok_uncond": r3(uncond.mean()),
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
    pooled = conf[(conf.alpha == A) & (conf.design == "within-gene random split") & (conf.method == "Pooled conformal")]
    out["pooled_abstain_AM_BAP1"] = r3(pooled[(pooled.predictor == "AlphaMissense") & (pooled.gene == "BAP1")].abstain.iloc[0])
    out["n_cal_lof_range"] = [int(rel.n_cal_lof.min()), int(rel.n_cal_lof.max())]
    # ---- additions for the revised paper
    tot = co[["clinvar_PLP", "clinvar_BLB", "clinvar_VUS"]].sum()
    out["clinvar_vus_share"] = r3(tot.clinvar_VUS / tot.sum())
    m = conf[(conf.alpha == A) & (conf.design == "within-gene random split") & (conf.method == "Mondrian conformal")].copy()
    if "coverage_lof_rw_sd" in m:
        se = pd.concat([m.coverage_lof_rw_sd, m.coverage_func_rw_sd]) / np.sqrt(m.n_rep.iloc[0])
        out["mondrian_se_range"] = [r3(se.min()), r3(se.max())]
        out["mondrian_rows"] = m[["predictor", "gene", "coverage_lof_rw", "coverage_lof_rw_sd", "coverage_func_rw",
                                  "coverage_func_rw_sd", "abstain", "n_rep"]].round(4).to_dict(orient="records")
    # draw-level validity: every bound claimed in a draw respected (draws without claims count as successes)
    out["oracle_draw_level"] = {}
    for sel, g in orc.groupby("selection"):
        ok = g.groupby(["predictor", "gene", "split"]).apply(
            lambda t: bool((~t.achieved | t.bound_ok.fillna(False).astype(bool)).all())).values
        lvl = g.groupby(["predictor", "gene", "split"]).apply(
            lambda t: bool((~t.achieved | t.holds.fillna(False).astype(bool)).all())).values
        boot = [ok[rng.integers(0, len(ok), len(ok))].mean() for _ in range(2000)]
        out["oracle_draw_level"][sel] = {"draws": int(len(ok)), "all_bounds_ok": r3(ok.mean()),
                                         "all_levels_held": r3(lvl.mean()),
                                         "lo": r3(np.percentile(boot, 2.5)), "hi": r3(np.percentile(boot, 97.5))}
    ps = pd.read_csv(RES / "prior_shift.csv")
    g = ps.groupby(["method", "cal_share", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw", "abstain"]].mean()
    g["worst"] = g[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    out["prior_shift"] = {f"{mt}|{sh}": {"worst_mean": r3(v.worst.mean()), "worst_min": r3(v.worst.min()),
                                         "worst_max": r3(v.worst.max()), "lof_min": r3(v.coverage_lof_rw.min()),
                                         "abstain_mean": r3(v.abstain.mean())}
                          for (mt, sh), v in g.groupby(level=[0, 1])}
    out["prior_shift_ncal"] = [int(ps.n_cal.min()), int(ps.n_cal.max())]
    du = pd.read_csv(RES / "data_use.csv")
    g = du.groupby(["calibration", "predictor", "gene"])[["coverage_lof_rw", "coverage_func_rw", "abstain", "n_cal"]].mean()
    g["worst"] = g[["coverage_lof_rw", "coverage_func_rw"]].min(axis=1)
    out["data_use"] = {c: {"worst_min": r3(v.worst.min()), "worst_max": r3(v.worst.max()), "worst_mean": r3(v.worst.mean()),
                           "n_below_090": int((v.worst < 0.90).sum()), "abstain_mean": r3(v.abstain.mean()),
                           "n_cal_mean": r3(v.n_cal.mean())} for c, v in g.groupby(level=0)}
    it = du[(du.calibration == "one per residue") & du.n_int.notna() & (du.n_int > 0)]
    vv = pd.read_csv(ROOT / "data" / "processed" / "variants.csv", low_memory=False)
    vv_src = (vv.stars >= 2) & vv.group.isin(["path", "benign"])
    out["intermediates"] = {"n_variants": int((vv.assay_label == -1).sum()),
                            "n_outside_clinvar_set": int(((vv.assay_label == -1) & ~vv_src).sum()),
                            "outside_by_gene": {g: int(n) for g, n in vv[(vv.assay_label == -1) & ~vv_src].groupby("gene").size().items()},
                            "lof": r3((it.int_lof * it.n_int).sum() / it.n_int.sum()),
                            "func": r3((it.int_func * it.n_int).sum() / it.n_int.sum())}
    ox = pd.read_csv(RES / "oracle_extended.csv")
    out["oracle_extended"] = {}
    for (gr, sel), v in ox.groupby(["grid", "selection"]):
        a = v[v.achieved]
        lvl = v.groupby(["predictor", "gene", "draw"]).apply(
            lambda t: bool((~t.achieved | t.holds.fillna(False).astype(bool)).all())).values
        out["oracle_extended"][f"{gr}|{sel}"] = {
            "claims": int(len(a)), "level_held": r3(a.holds.astype(float).mean()),
            "level_held_vw": r3(a.vw_holds.astype(float).mean()),
            "bound_ok": r3(a.bound_ok.astype(float).mean()) if a.bound_ok.notna().any() else None,
            "draws_all_levels_held": r3(lvl.mean()),
            "achieved": {k: r3(x) for k, x in v.groupby("level").achieved.mean().items()},
            "held_by_level": {k: r3(x) for k, x in a.groupby("level").holds.mean().items()},
            "n_by_level": {k: int(x) for k, x in a.groupby("level").size().items()}}
    c = ox[(ox.grid == "G4") & (ox.selection == "certified") & ox.achieved]
    lof, fun = c[c.level.str.startswith("PS3")], c[c.level.str.startswith("BS3")]
    out["coverage_at_selection"] = {"cov_lof_mean": r3(c.cov_lof.mean()), "cov_func_mean": r3(c.cov_func.mean()),
                                    "cov_lof_min": r3(c.cov_lof.min()), "cov_func_min": r3(c.cov_func.min()),
                                    "cov_lof_p10": r3(c.cov_lof.quantile(0.1)), "cov_func_p10": r3(c.cov_func.quantile(0.1)),
                                    "alpha_mean_lof_calls": r3(lof.alpha.mean()), "alpha_mean_func_calls": r3(fun.alpha.mean())}
    b = ox[(ox.grid == "G4") & (ox.selection == "bootstrap") & ox.achieved].copy()
    b["bad"] = ~b.bound_ok.astype(bool)
    b["zero"] = np.where(b.level.str.startswith("PS3"), b.cert_fp == 0, b.cert_fn == 0)
    out["bootstrap_diagnostic"] = {"n": int(len(b)), "fail": int(b.bad.sum()), "fail_rate": r3(b.bad.mean()),
                                   "share_fail_zero": r3(b[b.bad].zero.mean()), "share_zero": r3(b.zero.mean()),
                                   "fail_rate_zero": r3(b[b.zero].bad.mean()), "fail_rate_nonzero": r3(b[~b.zero].bad.mean())}
    # authors' classes against the control rule (Section 4.1)
    cmp = {}
    for g in ("BRCA1", "BAP1", "RAD51C"):
        s = vv[(vv.gene == g) & vv.score.notna()]
        both = s.assay_label.isin([0, 1]) & s.label_ctrl.isin([0, 1])
        func = s.assay_label == 0
        cmp[g] = {"both_confident": int(both.sum()),
                  "opposite": int((both & (s.assay_label != s.label_ctrl)).sum()),
                  "authors_functional_left_intermediate": r3((func & (s.label_ctrl == -1)).sum() / func.sum())}
    out["label_rule_comparison"] = cmp
    out.update(robustness_numbers(vv, vv_src))
    (RES / "paper_numbers.json").write_text(json.dumps(out, indent=1, default=str))
    print(json.dumps({k: out[k] for k in ["source_prev", "target_prev_min", "target_prev_max", "em_abs_err_mean",
                                           "em_abs_err_max", "dual_agree", "ece_by_method"]}, indent=1))
    print(json.dumps(summ["within-gene random split"], indent=1))


if __name__ == "__main__":
    main()
