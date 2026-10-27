"""Zero-shot ESM-2 (650M) log-likelihood ratios for every missense substitution.

Wild-type marginals (Meier et al. 2021): one unmasked forward pass, then
LLR(i, a->b) = log p(b | x) - log p(a | x) at position i.
Proteins longer than the 1022-residue context are tiled into overlapping
windows; each residue takes its scores from the window in which it is most central.
Output: data/processed/esm2_llr.csv (gene, protein_pos, wt, mt, esm_llr).
"""
import pathlib

import esm
import pandas as pd
import requests
import torch

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
UNIPROT = {"BRCA1": "P38398", "BRCA2": "P51587", "BAP1": "Q92560",
           "RAD51C": "O43502", "MSH2": "P43246"}
AA = "ACDEFGHIKLMNPQRSTVWY"
WIN, STRIDE = 1022, 511


def sequence(acc):
    fasta = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60).text
    return "".join(fasta.splitlines()[1:])


def windows(n):
    if n <= WIN:
        return [(0, n)]
    starts = list(range(0, n - WIN + 1, STRIDE))
    if starts[-1] + WIN < n:
        starts.append(n - WIN)
    return [(s, s + WIN) for s in starts]


@torch.no_grad()
def score(model, alphabet, device, seq):
    conv = alphabet.get_batch_converter()
    best = {}  # residue index -> (distance to window centre, log-prob row)
    for s, e in windows(len(seq)):
        _, _, toks = conv([("p", seq[s:e])])
        logp = torch.log_softmax(model(toks.to(device))["logits"][0, 1:-1].float(), dim=-1).cpu()
        centre = (s + e - 1) / 2
        for k in range(e - s):
            i = s + k
            d = abs(i - centre)
            if i not in best or d < best[i][0]:
                best[i] = (d, logp[k])
    idx = {a: alphabet.get_idx(a) for a in AA}
    rows = []
    for i, wt in enumerate(seq):
        row = best[i][1]
        for mt in AA:
            if mt != wt and wt in idx:
                rows.append((i + 1, wt, mt, float(row[idx[mt]] - row[idx[wt]])))
    return rows


def main():
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval().to(device)
    out = []
    for gene, acc in UNIPROT.items():
        seq = sequence(acc)
        rows = score(model, alphabet, device, seq)
        out += [(gene, *r) for r in rows]
        print(gene, acc, len(seq), "residues", len(rows), "substitutions", flush=True)
    pd.DataFrame(out, columns=["gene", "protein_pos", "wt", "mt", "esm_llr"]).to_csv(
        PROC / "esm2_llr.csv", index=False)


if __name__ == "__main__":
    main()
