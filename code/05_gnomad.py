"""gnomAD v4 joint (exome + genome) allele counts for every variant in the five genes.

One GraphQL query per gene. Output: data/processed/gnomad.csv
(gene, chrom, pos, ref, alt, ac, an, af, hom). Variants absent from gnomAD are
left missing downstream and get an explicit `in_gnomad = 0` flag rather than AF = 0.
"""
import pathlib
import time

import pandas as pd
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
API = "https://gnomad.broadinstitute.org/api"
GENES = ["BRCA1", "BRCA2", "BAP1", "RAD51C", "MSH2"]
QUERY = """{ gene(gene_symbol:"%s", reference_genome:GRCh38){
  variants(dataset: gnomad_r4){ variant_id pos ref alt
    joint{ ac an homozygote_count } } } }"""


def main():
    rows = []
    for g in GENES:
        for attempt in range(5):
            r = requests.post(API, json={"query": QUERY % g}, timeout=300)
            if r.status_code == 200 and "errors" not in r.json():
                break
            time.sleep(10 * (attempt + 1))
        for v in r.json()["data"]["gene"]["variants"]:
            j = v.get("joint") or {}
            if not j.get("an"):
                continue
            chrom = v["variant_id"].split("-")[0]
            rows.append((g, chrom, v["pos"], v["ref"], v["alt"], j["ac"], j["an"],
                         j["ac"] / j["an"], j.get("homozygote_count") or 0))
        print(g, sum(1 for x in rows if x[0] == g), flush=True)
        time.sleep(3)
    pd.DataFrame(rows, columns=["gene", "chrom", "pos", "ref", "alt", "ac", "an", "af", "hom"]).to_csv(
        PROC / "gnomad.csv", index=False)


if __name__ == "__main__":
    main()
