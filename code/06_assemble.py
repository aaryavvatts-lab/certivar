"""Join VEP, ESM-2, gnomAD, ClinVar and assay scores into one missense table and label
assayed variants as loss of function (1), functional (0) or intermediate (-1).

Assay labels. Where the assay authors published functional classes (MaveDB
"Investigator-provided functional classes") we use them: BRCA1 score ranges,
RAD51C per-variant classes (fast or slow depleted = LoF, unchanged = functional,
enriched = excluded), BAP1 classes (BH-FDR < 0.01 with negative score = depleted,
which reproduces the deposited class counts exactly) and the MSH2 threshold. BRCA2
has no investigator classes in MaveDB, so its labels come from the assay's internal
controls: Gaussians fitted to synonymous (functional) and nonsense (LoF) SNV scores,
LoF if the equal-prior posterior is >= 0.95, functional if <= 0.05. No ClinVar label
enters any of these rules. Agreement between the control model and the investigator
classes (BRCA1, RAD51C) is written to qc/assay_labels.txt as a check.

Outputs: data/processed/variants.csv, qc/assay_labels.txt
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
QC = ROOT / "qc"
AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
       "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
       "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter"}
LOF_SIGN = {"BRCA1": -1, "BRCA2": -1, "BAP1": -1, "RAD51C": -1, "MSH2": +1}  # direction of LoF
POST = 0.95
MINUS_STRAND = {"BRCA1", "BAP1"}  # genes on the reverse strand of GRCh38 (BRCA2, RAD51C, MSH2 are forward)


def load_vep():
    rows = [json.loads(l)["r"] for l in (PROC / "vep_raw.jsonl").read_text().splitlines()]
    v = pd.DataFrame([r for r in rows if r])
    v["cdot"] = v.hgvs.str.split(":").str[1]
    v["tx"] = v.hgvs.str.split(":").str[0]
    v["wt"] = v.amino_acids.str.split("/").str[0]
    v["mt"] = v.amino_acids.str.split("/").str[1]
    # VEP reports HGVS-input alleles on the transcript strand; gnomAD uses the forward strand
    comp = str.maketrans("ACGT", "TGCA")
    minus = v.gene.isin(MINUS_STRAND)
    for col, i in (("ref", 0), ("alt", 1)):
        a = v.allele.str.split("/").str[i]
        v[col] = np.where(minus, a.str.translate(comp), a)
    return v


def control_model(scores, conseq, sign):
    """Fit Gaussians to synonymous and nonsense scores; return P(LoF | score)."""
    syn = scores[conseq == "synonymous"]
    non = scores[conseq == "nonsense"]
    m0, s0, m1, s1 = syn.mean(), syn.std(), non.mean(), non.std()

    def post(x):
        l1, l0 = norm.pdf(x, m1, s1), norm.pdf(x, m0, s0)
        p = l1 / (l1 + l0)
        # beyond the control means the ordering is monotone even where tails cross
        p = np.where(sign * (x - m1) >= 0, np.maximum(p, 0.5), p)
        p = np.where(sign * (x - m0) <= 0, np.minimum(p, 0.5), p)
        return p

    return post, {"n_syn": int(len(syn)), "n_non": int(len(non)), "mu_syn": float(m0),
                  "sd_syn": float(s0), "mu_non": float(m1), "sd_non": float(s1)}


def main():
    vep = load_vep()
    clin = pd.read_csv(PROC / "clinvar.csv", low_memory=False)
    assay = pd.read_csv(PROC / "assay_snvs.csv", low_memory=False)
    esm = pd.read_csv(PROC / "esm2_llr.csv") if (PROC / "esm2_llr.csv").exists() else None
    gn = pd.read_csv(PROC / "gnomad.csv")

    def kind(c):
        c = str(c)
        if "stop_gained" in c:
            return "nonsense"
        if "missense_variant" in c:
            return "missense"
        if c == "synonymous_variant":
            return "synonymous"
        return "other"

    vep["kind"] = vep.consequence.map(kind)
    vep = vep[vep.mane == vep.mane].copy()  # MANE Select consequences only

    # assay scores: nucleotide-level join on cDNA, MSH2 on protein change
    nt = assay[assay.cdot.notna()][["gene", "cdot", "score"]]
    d = vep.merge(nt, on=["gene", "cdot"], how="left")
    pro = assay[assay.pdot.notna()][["gene", "pdot", "score"]].rename(columns={"score": "score_pro"})
    d["pdot"] = "p." + d.wt.map(AA3) + d.protein_pos.astype("Int64").astype(str) + d.mt.map(AA3)
    syn_mask = d.kind == "synonymous"
    d.loc[syn_mask, "pdot"] = "p." + d.wt.map(AA3) + d.protein_pos.astype("Int64").astype(str) + "="
    d = d.merge(pro, on=["gene", "pdot"], how="left")
    d["score"] = d.score.fillna(d.score_pro)
    d = d.drop(columns="score_pro")

    # assay labels: investigator classes where published, internal controls otherwise
    inv = assay[assay.inv_class.notna()][["gene", "cdot", "inv_class"]]
    d = d.merge(inv, on=["gene", "cdot"], how="left")
    lines, d["p_lof_ctrl"] = [], np.nan
    for g, sign in LOF_SIGN.items():
        m = (d.gene == g) & d.score.notna()
        k = d.loc[m, "kind"].values
        if (k == "synonymous").sum() >= 30 and (k == "nonsense").sum() >= 30:
            post, info = control_model(d.loc[m, "score"].values, k, sign)
            d.loc[m, "p_lof_ctrl"] = post(d.loc[m, "score"].values)
            lines.append(f"{g} control model: {info}")
    ctrl = np.select([d.p_lof_ctrl >= POST, d.p_lof_ctrl <= 1 - POST], [1, 0], default=-1).astype(float)
    ctrl[d.p_lof_ctrl.isna().values] = np.nan
    d["label_ctrl"] = ctrl
    d["assay_label"] = np.nan
    b1 = d.gene == "BRCA1"
    d.loc[b1, "assay_label"] = np.select([d.score <= -1.328, d.score >= -0.748], [1, 0], default=-1)[b1.values]
    bp = d.gene == "BAP1"
    d.loc[bp, "assay_label"] = d.loc[bp, "inv_class"].map({"depleted": 1, "unchanged": 0, "enriched": -1})
    r51 = d.gene == "RAD51C"
    d.loc[r51, "assay_label"] = d.loc[r51, "inv_class"].map(
        {"fast depleted": 1, "slow depleted": 1, "unchanged": 0, "enriched": -1})
    ms = d.gene == "MSH2"
    d.loc[ms, "assay_label"] = np.select([d.score > 0, d.score < 0], [1, 0], default=-1)[ms.values]
    m = d.gene == "BRCA2"
    d.loc[m, "assay_label"] = d.loc[m, "label_ctrl"]
    d.loc[d.score.isna(), "assay_label"] = np.nan
    # variants that also touch a splice region are excluded: SGE measures splicing effects, DMS does not
    d["splice_region"] = d.consequence.str.contains("splice", na=False)
    d.loc[d.splice_region, "assay_label"] = np.nan
    d["label_source"] = np.where(d.gene == "BRCA2", "internal controls", "investigator classes")

    # difficulty: distance from the midpoint between class medians, scaled by their separation
    d["assay_margin"] = np.nan
    for g in LOF_SIGN:
        m = (d.gene == g) & d.assay_label.isin([0, 1])
        med1 = d.loc[m & (d.assay_label == 1), "score"].median()
        med0 = d.loc[m & (d.assay_label == 0), "score"].median()
        mg = d.gene == g
        d.loc[mg, "assay_margin"] = (d.loc[mg, "score"] - (med0 + med1) / 2).abs() / abs(med0 - med1)

    # ClinVar annotations (cDNA on the MANE transcript)
    c = clin[clin.cdot.notna()].sort_values("stars", ascending=False).drop_duplicates(["GeneSymbol", "cdot"])
    c = c.rename(columns={"GeneSymbol": "gene"})[["gene", "cdot", "group", "stars", "VariationID",
                                                    "ClinicalSignificance", "ReviewStatus"]]
    d = d.merge(c, on=["gene", "cdot"], how="left")

    # ESM-2 and gnomAD
    if esm is not None:
        d = d.merge(esm, on=["gene", "protein_pos", "wt", "mt"], how="left")
    gn["chrom"] = gn.chrom.astype(str)
    d["chrom"] = d.chrom.astype(str)
    d = d.merge(gn[["chrom", "pos", "ref", "alt", "af", "hom"]], on=["chrom", "pos", "ref", "alt"], how="left")
    d["in_gnomad"] = d.af.notna().astype(int)
    d["log_af"] = np.log10(d.af.fillna(0) + 1e-7)

    d["sift_d"] = 1 - d.sift  # higher = more damaging, like the other scores
    d["cluster"] = d.gene + ":" + d.protein_pos.astype("Int64").astype(str)

    mis = d[d.kind == "missense"].drop_duplicates(["gene", "cdot"]).copy()
    mis.to_csv(PROC / "variants.csv", index=False)

    # label agreement summary: control model vs investigator classes
    lines.append("")
    for g in ("BRCA1", "BAP1", "RAD51C"):
        sub = d[(d.gene == g) & (d.kind == "missense") & d.assay_label.isin([0, 1, -1]) & d.label_ctrl.notna()]
        if len(sub):
            tab = pd.crosstab(sub.assay_label, sub.label_ctrl)
            both = sub[sub.assay_label.isin([0, 1]) & sub.label_ctrl.isin([0, 1])]
            lines.append(f"{g}: investigator (rows) vs control model (cols), all assayed missense:\n{tab.to_string()}\n"
                         f"  agreement where both confident {(both.assay_label == both.label_ctrl).mean():.3f} on {len(both)}; "
                         f"opposite calls {int(((both.assay_label - both.label_ctrl).abs() == 1).sum())}")
    lines.append("")
    for g in LOF_SIGN:
        sub = mis[mis.gene == g]
        lab = sub.assay_label.value_counts(dropna=False).to_dict()
        lines.append(f"{g}: missense {len(sub)}, assayed {sub.score.notna().sum()}, labels {lab}, "
                     f"ClinVar groups {sub.group.value_counts().to_dict()}")
    QC.mkdir(exist_ok=True)
    (QC / "assay_labels.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
