# Sources

Venue rules: RECOMB 2027 call for papers, https://recomb.org/recomb2027/call_for_papers (read 26 Oct 2026): 10 pages US
Letter including figures, tables and captions; 10 pt minimum; 1-inch margins; title page and references not counted;
optional appendix that reviewers need not read; single-blind; public code with README linked in the paper; abstract
registration 28 Oct 2026 and full paper 4 Nov 2026 (23:59 AoE). Re-checked 27 Oct 2026: unchanged; no column format
required; LNCS optional.

Data: ClinVar variant_summary.txt.gz (NCBI FTP, Last-Modified 26 Oct 2026 13:50:58 GMT); MaveDB score sets
urn:mavedb:00000097-0-2 (BRCA1), 00001279-a-1 (BRCA2), 00000662-0-1 (BAP1), 00000673-0-1 (RAD51C), 00000050-a-1 (MSH2);
Ensembl VEP REST (GRCh38; queried 26 Oct 2026, service reported release 116 on 27 Oct 2026); gnomAD v4 GraphQL; UniProt REST; ESM-2 esm2_t33_650M_UR50D (fair-esm 2.0.0).
SHA-256 checksums are in data/MANIFEST.txt.

Thresholds and facts checked at the source (26-27 Oct 2026):
- OddsPath evidence thresholds: Brnich et al. 2020, Table 3 (PMC6938631).
- REVEL PP3 >= 0.644, BP4 <= 0.290; sliding-window local likelihood ratio with a one-sided 95% bootstrap bound
  (10,000 resamples); prior 4.41%: Pejaver et al. 2022 (PMC9748256).
- AlphaMissense classes (< 0.34, > 0.564): google-deepmind/alphamissense README and the Ensembl VEP AlphaMissense plugin.
- Geifman and El-Yaniv 2017 (NIPS 30) and Park et al. 2020 (ICLR): proceedings pages; page or article numbers of
  Cheng 2023, Livesey 2020/2023, Esposito 2019 and McLaren 2016 from Crossref (27 Oct 2026).
- Harrizi et al. 2026 conformal procedure: full text, Sections 2.9 and 4.11 (bioRxiv 10.64898/2026.08.21.746180).
- BAP1 investigator classes: deposited MaveDB class counts (5,665 / 11,912 / 531) reproduced by BH-FDR < 0.01 and the
  sign of the score.

References: every entry of paper/refs.bib resolves by DOI (Crossref) or arXiv id; checked with citation_check.py --online.
- Added references checked against Crossref on 27 Oct 2026: Richards et al. 2015 (10.1038/gim.2015.30), Chow 1970
  (10.1109/TIT.1970.1054406), Katz et al. 1978 (10.2307/2530610), Koopman 1984 (10.2307/2531405).
- Uncertainty in cited calibration workflows: Pejaver et al. 2022 and Bergquist et al. 2025 (local posterior method with
  bootstrap bounds; Bergquist abstract via Europe PMC); Badonyi and Marsh 2025 (1,000 resamples, sign-test bounds on the
  median log likelihood ratio; Bioinformatics full text, Sections 2.1-2.2).
