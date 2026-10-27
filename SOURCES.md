# Sources

Venue rules: RECOMB 2027 call for papers, https://recomb.org/recomb2027/call_for_papers (read 26 Oct 2026): 10 pages US
Letter including figures, tables and captions; 10 pt minimum; 1-inch margins; title page and references not counted;
optional appendix that reviewers need not read; single-blind; public code with README linked in the paper; abstract
registration 28 Oct 2026 and full paper 4 Nov 2026 (23:59 AoE). No AI-assistance policy stated on that page.

Data: ClinVar variant_summary.txt.gz (NCBI FTP, Last-Modified 26 Oct 2026 13:50:58 GMT); MaveDB score sets
urn:mavedb:00000097-0-2 (BRCA1), 00001279-a-1 (BRCA2), 00000662-0-1 (BAP1), 00000673-0-1 (RAD51C), 00000050-a-1 (MSH2);
Ensembl VEP REST (GRCh38); gnomAD v4 GraphQL; UniProt REST; ESM-2 esm2_t33_650M_UR50D (fair-esm 2.0.0).
SHA-256 checksums are in data/MANIFEST.txt.

Thresholds and facts checked at the source (26-27 Oct 2026):
- OddsPath evidence thresholds: Brnich et al. 2020, Table 3 (PMC6938631).
- REVEL PP3 >= 0.644, BP4 <= 0.290; sliding-window local likelihood ratio with a one-sided 95% bootstrap bound
  (10,000 resamples); prior 4.41%: Pejaver et al. 2022 (PMC9748256).
- AlphaMissense classes (< 0.34, > 0.564): google-deepmind/alphamissense README and the Ensembl VEP AlphaMissense plugin.
- Harrizi et al. 2026 conformal procedure: full text, Sections 2.9 and 4.11 (bioRxiv 10.64898/2026.08.21.746180).
- BAP1 investigator classes: deposited MaveDB class counts (5,665 / 11,912 / 531) reproduced by BH-FDR < 0.01 and the
  sign of the score.

References: every entry of paper/refs.bib resolves by DOI (Crossref) or arXiv id; checked with citation_check.py --online.
