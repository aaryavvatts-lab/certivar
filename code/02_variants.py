"""Build the variant universe: ClinVar SNVs in the five genes plus every assayed SNV.

Outputs
  data/processed/clinvar.csv   one row per ClinVar GRCh38 SNV in the genes
                               (gene, cdot, pdot, clinsig group, review stars, VCF position)
  data/processed/assay_snvs.csv one row per assayed coding SNV (gene, cdot, score)
                               MSH2 is protein-level (pdot, score) and is joined later.
  data/processed/vep_queue.txt  HGVS strings to annotate with Ensembl VEP
"""
import gzip
import pathlib
import re

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"

GENES = {
    # gene: (MANE RefSeq used for HGVS, MaveDB urn short id, assay level)
    "BRCA1": ("NM_007294.4", "00000097-0-2", "nt"),
    "BRCA2": ("NM_000059.4", "00001279-a-1", "nt"),
    "BAP1": ("NM_004656.4", "00000662-0-1", "nt"),
    "RAD51C": ("NM_058216.3", "00000673-0-1", "nt"),
    "MSH2": ("NM_000251.3", "00000050-a-1", "pro"),
}

CDOT = re.compile(r"(c\.-?\*?\d+(?:[+-]\d+)?[ACGT]>[ACGT])")
PDOT = re.compile(r"\((p\.[A-Za-z]{3}\d+[A-Za-z]{3}|p\.[A-Za-z]{3}\d+=)\)")
STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
}


def clinsig_group(s):
    s = s.lower()
    if "conflicting" in s:
        return "conflicting"
    if s.startswith("pathogenic") or s.startswith("likely pathogenic"):
        return "path"
    if s.startswith("benign") or s.startswith("likely benign"):
        return "benign"
    if s.startswith("uncertain significance"):
        return "vus"
    return "other"


def load_clinvar():
    cols = ["GeneSymbol", "Type", "Name", "ClinicalSignificance", "ReviewStatus", "Assembly",
            "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF",
            "VariationID", "LastEvaluated", "NumberSubmitters"]
    chunks = []
    with gzip.open(RAW / "variant_summary.txt.gz", "rt") as fh:
        for ch in pd.read_csv(fh, sep="\t", usecols=cols, chunksize=500_000, low_memory=False):
            ch = ch[(ch.Assembly == "GRCh38") & (ch.Type == "single nucleotide variant")
                    & ch.GeneSymbol.isin(GENES)]
            chunks.append(ch)
    cv = pd.concat(chunks)
    cv["cdot"] = cv.Name.str.extract(CDOT)[0]
    cv["pdot"] = cv.Name.str.extract(PDOT)[0]
    cv["tx"] = cv.Name.str.extract(r"^(N[MR]_\d+\.\d+)")[0]
    cv["group"] = cv.ClinicalSignificance.map(clinsig_group)
    cv["stars"] = cv.ReviewStatus.map(STARS).fillna(0).astype(int)
    return cv


def load_assay(gene, urn, level):
    df = pd.read_csv(RAW / f"mavedb_{urn}.csv", low_memory=False)
    df = df[df.score.notna()].copy()
    df["gene"] = gene
    df["inv_class"] = df["functional_classification"] if "functional_classification" in df else pd.NA
    if gene == "BAP1":  # Waters et al. classes: BH-FDR < 0.01, direction from the sign of the score
        sig = df.processed_BH_FDR_continuous < 0.01
        df["inv_class"] = np.where(sig & (df.score < 0), "depleted", np.where(sig & (df.score > 0), "enriched", "unchanged"))
    if level == "nt":
        df["cdot"] = df.hgvs_nt.str.split(":").str[1]
        df = df[df.cdot.str.fullmatch(r"c\.\d+[ACGT]>[ACGT]", na=False)]
        return df[["gene", "cdot", "score", "inv_class"]]
    df["pdot"] = df.hgvs_pro
    return df[["gene", "pdot", "score", "inv_class"]]


def main():
    cv = load_clinvar()
    cv.to_csv(PROC / "clinvar.csv", index=False)
    print("ClinVar GRCh38 SNVs per gene:\n", cv.GeneSymbol.value_counts().to_string())

    assays = []
    for gene, (tx, urn, level) in GENES.items():
        a = load_assay(gene, urn, level)
        a["tx"] = tx
        assays.append(a)
        print(gene, "assayed coding SNVs" if level == "nt" else "assayed protein changes", len(a))
    assays = pd.concat(assays, ignore_index=True)
    assays.to_csv(PROC / "assay_snvs.csv", index=False)

    # VEP queue: every exonic ClinVar SNV in the MANE transcript plus every assayed coding SNV.
    q = set()
    for gene, (tx, _, _) in GENES.items():
        sub = cv[(cv.GeneSymbol == gene) & cv.cdot.notna() & (cv.tx == tx)]
        sub = sub[sub.cdot.str.fullmatch(r"c\.\d+[ACGT]>[ACGT]")]
        q.update(f"{tx}:{c}" for c in sub.cdot)
        nt = assays[(assays.gene == gene) & assays.cdot.notna()]
        q.update(f"{tx}:{c}" for c in nt.cdot)
    (PROC / "vep_queue.txt").write_text("\n".join(sorted(q)) + "\n")
    print("VEP queue", len(q))
    print("ClinVar transcript mismatches (not MANE):",
          (~cv[cv.cdot.notna()].apply(lambda r: r.tx == GENES[r.GeneSymbol][0], axis=1)).sum())


if __name__ == "__main__":
    main()
