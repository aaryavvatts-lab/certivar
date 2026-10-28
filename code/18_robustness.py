"""Robustness experiments (E14-E17).

E14 coverage_oracle.csv     Population (conditional) per-class coverage of in-expectation and PAC thresholds: in each of
                            200 draws per pair, half as many calibration residues as the gene has are drawn with
                            replacement (one variant per residue), thresholds are set at alpha = 0.1, and coverage is
                            computed exactly on every assayed variant of the gene under the residue-first law Q. This is
                            the quantity Propositions 1 and 2 bound; a finite test half adds sampling noise to it.
E15 calibration_weighting.csv
                            Calibration on one variant per residue, on all variants of the calibration half without
                            weights, and on all variants with residue weights (certivar.residue_weighted_thresholds);
                            50 random halvings per pair, alpha = 0.1, coverage on the other half.
E16 mixture_check.csv       Mixture-consistency check of Assumption 2 (certivar.mixture_check) in five settings, with
                            the coverage of in-expectation thresholds on the same target where labels exist:
                            random half (Assumption 2 holds by design), assayed VUS at other residues, a held-out
                            contiguous quarter of the protein, another gene's calibration data, and unassayed ClinVar VUS
                            at residues the assay did not reach.
E17 use_case.csv            ClinVar missense VUS per gene by assay status, and the certified calls (Algorithm 1, G4,
                            delta = 0.05, 20 random 40/60 residue splits of all assayed variants) on the unassayed VUS at
                            assayed residues, with the mixture check of those VUS where they span at least 30 residues.
"""
import importlib
import pathlib

import numpy as np
import pandas as pd

import certivar as sc

ev = importlib.import_module("13_evidence")
ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results"
SEED, ALPHA, N_COV, N_W, N_BOOT = 2027, 0.10, 200, 50, 300


def one_per(idx, cl, rng):
    """One random variant per residue among the indices idx."""
    idx = np.asarray(idx)
    return idx[sc.one_per_cluster(cl[idx], rng)]


def coverage_oracle(d):
    rows = []
    rng = np.random.default_rng(SEED + 30)
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        for g in ev.GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            res = np.unique(cl)
            members = {r: np.flatnonzero(cl == r) for r in res}
            for draw in range(N_COV):
                ic = np.array([rng.choice(members[r]) for r in rng.choice(res, len(res) // 2)])
                for name, delta in (("in expectation", None), ("PAC, delta = 0.1", 0.1)):
                    thr = sc.mondrian_thresholds(p[ic], y[ic], ALPHA, delta=delta)
                    rep = sc.coverage_report(sc.decisions(p, thr), y, cl)
                    rows.append(dict(predictor=pr, gene=g, draw=draw, thresholds=name, **rep))
    return pd.DataFrame(rows)


def calibration_weighting(d):
    rows = []
    rng = np.random.default_rng(SEED + 31)
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        for g in ev.GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            res = np.unique(cl)
            for rep in range(N_W):
                m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                half = np.flatnonzero(m)
                keep = one_per(half, cl, rng)
                thrs = {"one per residue": sc.mondrian_thresholds(p[keep], y[keep], ALPHA),
                        "all variants": sc.mondrian_thresholds(p[half], y[half], ALPHA),
                        "all variants, residue weights": sc.residue_weighted_thresholds(p[half], y[half], cl[half], ALPHA)}
                for name, thr in thrs.items():
                    rows.append(dict(predictor=pr, gene=g, rep=rep, calibration=name,
                                     **sc.coverage_report(sc.decisions(p[~m], thr), y[~m], cl[~m])))
    return pd.DataFrame(rows)


def check(p, y, ic, it, rng):
    """Mixture check of target indices it against calibration indices ic, and the coverage on it of in-expectation
    thresholds set on ic."""
    stat, pi, pval = sc.mixture_check(p[ic][y[ic] == 1], p[ic][y[ic] == 0], p[it], n_boot=N_BOOT, rng=rng)
    thr = sc.mondrian_thresholds(p[ic], y[ic], ALPHA)
    return dict(stat=stat, pi_fit=pi, pval=pval, n_target=len(it), n_cal=len(ic), true_share=float(y[it].mean()),
                **sc.coverage_report(sc.decisions(p[it], thr), y[it]))


def mixture_check(d):
    rows = []
    rng = np.random.default_rng(SEED + 32)
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        data = {g: ev.gene_data(d, g, col, cal) for g in ev.GENES}
        for g in ev.GENES:
            t, p, y, cl = data[g]
            res = np.unique(cl)
            vus = (t.group == "vus").values
            base = dict(predictor=pr, gene=g)
            for rep in range(100):   # (1) random half: Assumption 2 holds by design
                m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                ic, it = one_per(np.flatnonzero(m), cl, rng), one_per(np.flatnonzero(~m), cl, rng)
                rows.append(dict(base, setting="random half", rep=rep, **check(p, y, ic, it, rng)))
            for rep in range(50):    # (2) assayed VUS at other residues
                m = np.isin(cl, rng.choice(res, len(res) // 2, replace=False))
                ic, it = one_per(np.flatnonzero(m & ~vus), cl, rng), one_per(np.flatnonzero(~m & vus), cl, rng)
                if len(it) >= 30 and len(np.unique(y[it])) == 2:
                    rows.append(dict(base, setting="assayed VUS", rep=rep, **check(p, y, ic, it, rng)))
            pos = t.protein_pos.values   # (3) held-out contiguous quarter (blocks as in 13_evidence.region_evidence)
            edges = np.quantile(np.unique(pos), np.linspace(0, 1, 5))
            blk = np.clip(np.searchsorted(edges, pos, side="right") - 1, 0, 3)
            for b in range(4):
                for rep in range(25):
                    ic, it = one_per(np.flatnonzero(blk != b), cl, rng), one_per(np.flatnonzero(blk == b), cl, rng)
                    if len(np.unique(y[it])) == 2:
                        rows.append(dict(base, setting="held-out quarter", block=b, rep=rep, **check(p, y, ic, it, rng)))
            for other in ev.GENES:   # (4) calibration data from another gene
                if other == g:
                    continue
                to, po, yo, clo = data[other]
                for rep in range(5):
                    ic = sc.one_per_cluster(clo, rng)
                    it = one_per(np.arange(len(y)), cl, rng)
                    stat, pi, pval = sc.mixture_check(po[ic][yo[ic] == 1], po[ic][yo[ic] == 0], p[it], n_boot=N_BOOT,
                                                      rng=rng)
                    thr = sc.mondrian_thresholds(po[ic], yo[ic], ALPHA)
                    rows.append(dict(base, setting="other gene", source=other, rep=rep, stat=stat, pi_fit=pi, pval=pval,
                                     n_target=len(it), n_cal=len(ic), true_share=float(y[it].mean()),
                                     **sc.coverage_report(sc.decisions(p[it], thr), y[it])))
            un = d[(d.gene == g) & (d.group == "vus") & d.score.isna() & d[col].notna()]   # (5) unassayed VUS
            un = un[~un.cluster.isin(set(d[(d.gene == g) & d.score.notna()].cluster))]
            if un.cluster.nunique() >= 30:
                pu, clu = cal.predict(un[col].values), un.cluster.values
                for rep in range(20):
                    ic = sc.one_per_cluster(cl, rng)
                    iu = sc.one_per_cluster(clu, rng)
                    stat, pi, pval = sc.mixture_check(p[ic][y[ic] == 1], p[ic][y[ic] == 0], pu[iu], n_boot=N_BOOT,
                                                      rng=rng)
                    thr = sc.mondrian_thresholds(p[ic], y[ic], ALPHA)
                    du = sc.decisions(pu[iu], thr)
                    rows.append(dict(base, setting="unassayed VUS, other residues", rep=rep, stat=stat, pi_fit=pi,
                                     pval=pval, n_target=len(iu), n_cal=len(ic), abstain=float(np.mean(du < 0))))
            print(pr, g, flush=True)
    return pd.DataFrame(rows)


def use_case(d):
    rows = []
    rng = np.random.default_rng(SEED + 33)
    level = ev.DELTA / (2 * len(ev.GRID))
    for pr, col in ev.PRED.items():
        cal = ev.clinvar_posterior(d, col)
        for g in ev.GENES:
            t, p, y, cl = ev.gene_data(d, g, col, cal)
            s = d[d.gene == g]
            v = s[s.group == "vus"]
            assayed_res = set(s[s.score.notna()].cluster)
            in_res = v.cluster.isin(assayed_res)
            un_in = v[v.score.isna() & in_res & v[col].notna()]
            base = dict(predictor=pr, gene=g, vus=len(v), vus_assayed=int(v.score.notna().sum()),
                        vus_unassayed_at_assayed_residues=int((v.score.isna() & in_res).sum()),
                        vus_unassayed_elsewhere=int((v.score.isna() & ~in_res).sum()), n_target=len(un_in))
            pu = cal.predict(un_in[col].values) if len(un_in) else np.array([])
            res = np.unique(cl)
            for rep in range(20):
                perm = rng.permutation(res)
                cut = int(0.4 * len(perm))
                ic = one_per(np.flatnonzero(np.isin(cl, perm[:cut])), cl, rng)
                iz = one_per(np.flatnonzero(np.isin(cl, perm[cut:])), cl, rng)
                thrs = {a: sc.mondrian_thresholds(p[ic], y[ic], a) for a in ev.GRID}
                cert = {a: sc.certify(sc.decisions(p[iz], thrs[a]), y[iz], level) for a in ev.GRID}
                row = dict(base, rep=rep)
                clu = un_in.cluster.values
                if len(np.unique(clu)) >= 30:   # mixture check of these targets against part A (one per residue)
                    iu = sc.one_per_cluster(clu, rng)
                    row["check_pval"] = sc.mixture_check(p[ic][y[ic] == 1], p[ic][y[ic] == 0], pu[iu], n_boot=N_BOOT,
                                                         rng=rng)[2]
                for name, best in ev.select(cert, "certified").items():
                    k = 1 if name.startswith("PS3") else 0
                    row[name] = int(np.sum(sc.decisions(pu, thrs[best[0]]) == k)) if best is not None and len(pu) else 0
                rows.append(row)
    return pd.DataFrame(rows)


def main(parts=("coverage_oracle", "calibration_weighting", "mixture_check", "use_case")):
    d = ev.load()
    for part in parts:
        globals()[part](d).to_csv(RES / f"{part}.csv", index=False)
        print(part, "done", flush=True)


if __name__ == "__main__":
    import sys
    main(tuple(sys.argv[1:]) or ("coverage_oracle", "calibration_weighting", "mixture_check", "use_case"))
