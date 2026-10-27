"""Download MaveDB score sets and their investigator-provided functional class thresholds.

Writes data/raw/mavedb_<urn>.csv (scores) and data/raw/meta_<urn>.json (metadata),
plus data/processed/assays.json summarising each assay's class thresholds.
ClinVar variant_summary.txt.gz is fetched separately (see Makefile).
"""
import json
import pathlib
import sys

import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
API = "https://api.mavedb.org/api/v1"

ASSAYS = {
    "BRCA1": "urn:mavedb:00000097-0-2",
    "BRCA2": "urn:mavedb:00001279-a-1",
    "BAP1": "urn:mavedb:00000662-0-1",
    "RAD51C": "urn:mavedb:00000673-0-1",
    "MSH2": "urn:mavedb:00000050-a-1",
}


def classes(meta):
    """Return the investigator-provided functional classes keyed by normal/abnormal/not_specified."""
    for cal in meta.get("scoreCalibrations") or []:
        if cal.get("title", "").startswith("Investigator-provided functional classes"):
            out = {}
            for fc in cal["functionalClassifications"]:
                out.setdefault(fc["functionalClassification"], []).append({
                    "label": fc["label"],
                    "range": fc.get("range"),
                    "classValue": fc.get("class"),
                    "inclusiveLower": fc.get("inclusiveLowerBound"),
                    "inclusiveUpper": fc.get("inclusiveUpperBound"),
                })
            return out
    return None


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    summary = {}
    for gene, urn in ASSAYS.items():
        short = urn.split(":")[-1]
        meta = requests.get(f"{API}/score-sets/{urn}", timeout=120).json()
        (RAW / f"meta_{short}.json").write_text(json.dumps(meta))
        r = requests.get(f"{API}/score-sets/{urn}/scores", timeout=600)
        r.raise_for_status()
        (RAW / f"mavedb_{short}.csv").write_text(r.text)
        summary[gene] = {
            "urn": urn,
            "title": meta.get("title"),
            "numVariants": meta.get("numVariants"),
            "target": (meta.get("targetGenes") or [{}])[0].get("targetAccession"),
            "publications": [p.get("doi") for p in meta.get("primaryPublicationIdentifiers", [])],
            "classes": classes(meta),
        }
        print(gene, urn, meta.get("numVariants"), file=sys.stderr)
    (PROC / "assays.json").write_text(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
