"""Annotate every queued HGVS with Ensembl VEP REST (GRCh38).

Keeps the MANE Select transcript consequence and returns, per variant:
genomic position, consequence, protein position and amino acids, and the
predictor scores CADD (PHRED), SIFT, PolyPhen-2, REVEL and AlphaMissense.
Adapted from the earlier single-gene rebuild (fetch_vep.py): batches of 200,
retries with backoff, and a JSONL checkpoint so an interrupted run resumes.
"""
import json
import pathlib
import time
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
CKPT = PROC / "vep_raw.jsonl"
URL = "https://rest.ensembl.org/vep/human/hgvs"
HEAD = {"Content-Type": "application/json", "Accept": "application/json"}
BATCH = 200


def pick(rec):
    """Keep the consequence on the queried MANE Select transcript (never a neighbouring gene's)."""
    tx = (rec.get("input") or "").split(":")[0]
    tcs = rec.get("transcript_consequences") or []
    match = [t for t in tcs if t.get("mane_select") == tx]
    if not match:
        return None
    tc = match[0]
    return {
        "hgvs": rec.get("input"),
        "chrom": rec.get("seq_region_name"),
        "pos": rec.get("start"),
        "allele": rec.get("allele_string"),
        "transcript": tc.get("transcript_id"),
        "mane": tc.get("mane_select"),
        "gene": tc.get("gene_symbol"),
        "consequence": ",".join(tc.get("consequence_terms") or []),
        "protein_pos": tc.get("protein_start"),
        "amino_acids": tc.get("amino_acids"),
        "cadd": tc.get("cadd_phred"),
        "sift": tc.get("sift_score"),
        "polyphen": tc.get("polyphen_score"),
        "revel": tc.get("revel"),
        "alphamissense": (tc.get("alphamissense") or {}).get("am_pathogenicity")
        if isinstance(tc.get("alphamissense"), dict) else tc.get("alphamissense"),
    }


def fetch(chunk):
    body = {"hgvs_notations": chunk, "CADD": 1, "AlphaMissense": 1, "REVEL": 1, "mane": 1}
    for attempt in range(6):
        try:
            r = requests.post(URL, headers=HEAD, data=json.dumps(body), timeout=300)
            if r.status_code == 200:
                return [p for p in (pick(x) for x in r.json()) if p]
            if r.status_code in (429, 503, 504, 502):
                time.sleep(float(r.headers.get("Retry-After", 5 * (attempt + 1))))
                continue
            # 400 usually means one bad HGVS in the batch: split and retry halves
            if r.status_code == 400 and len(chunk) > 1:
                h = len(chunk) // 2
                return fetch(chunk[:h]) + fetch(chunk[h:])
            return []
        except requests.RequestException:
            time.sleep(5 * (attempt + 1))
    return []


def main():
    queue = (PROC / "vep_queue.txt").read_text().split()
    done = set()
    if CKPT.exists():
        for line in CKPT.read_text().splitlines():
            done.add(json.loads(line)["q"])
    todo = [q for q in queue if q not in done]
    batches = [todo[i:i + BATCH] for i in range(0, len(todo), BATCH)]
    print(f"queue {len(queue)}, done {len(done)}, batches left {len(batches)}", flush=True)
    with open(CKPT, "a") as out, ThreadPoolExecutor(max_workers=8) as ex:
        for i, (chunk, res) in enumerate(zip(batches, ex.map(fetch, batches))):
            got = {r["hgvs"]: r for r in res}
            for q in chunk:
                out.write(json.dumps({"q": q, "r": got.get(q)}) + "\n")
            out.flush()
            print(f"batch {i + 1}/{len(batches)} annotated {len(got)}/{len(chunk)}", flush=True)
    rows = [json.loads(l)["r"] for l in CKPT.read_text().splitlines()]
    df = pd.DataFrame([r for r in rows if r])
    df.to_csv(PROC / "vep.csv", index=False)
    print("annotated", len(df), "of", len(queue))


if __name__ == "__main__":
    main()
