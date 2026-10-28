PY ?= .venv/bin/python
RAW = data/raw
PROC = data/processed

.PHONY: reproduce data check-clinvar annotate results figures manifest clean-results test

reproduce: data annotate test results figures manifest

test:
	cd code && ../$(PY) test_certivar.py

$(RAW)/variant_summary.txt.gz:
	mkdir -p $(RAW)
	curl -sL -o $@ https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

# The paper uses the ClinVar release of 26 Oct 2026. The NCBI endpoint always serves the latest release, so a new
# download is rejected unless its SHA-256 matches data/MANIFEST.txt (set ALLOW_NEW_CLINVAR=1 to use a newer release;
# the numbers will then differ from the paper). `make test results figures` needs no download.
CLINVAR_SHA = 657b47d0a00ba3ad4b3ddf03e160835c8a4a3c7a2e130c92c904b70021a35703

check-clinvar: $(RAW)/variant_summary.txt.gz
	@if [ "$$(shasum -a 256 $< | cut -d' ' -f1)" != "$(CLINVAR_SHA)" ] && [ -z "$(ALLOW_NEW_CLINVAR)" ]; then \
	  echo "ClinVar file differs from the 26 Oct 2026 release used in the paper (SHA-256 mismatch)."; \
	  echo "Use the committed data/processed files (make test results figures) or set ALLOW_NEW_CLINVAR=1."; exit 1; fi

data: check-clinvar
	$(PY) code/01_download.py
	$(PY) code/02_variants.py

annotate:
	$(PY) code/03_vep.py
	$(PY) code/04_esm.py
	$(PY) code/05_gnomad.py
	$(PY) code/06_assemble.py

results:
	cd code && ../$(PY) 07_simulation.py
	cd code && ../$(PY) 08_experiments.py
	cd code && ../$(PY) 13_evidence.py
	cd code && ../$(PY) 14_label_sensitivity.py
	cd code && ../$(PY) 15_clinvar_check.py
	cd code && ../$(PY) 17_extensions.py
	cd code && ../$(PY) 18_robustness.py
	cd code && ../$(PY) 19_baselines.py
	cd code && ../$(PY) 11_numbers.py

figures:
	cd code && ../$(PY) 09_figures.py
	cd code && ../$(PY) 12_tables.py
	cd code && ../$(PY) 16_short_bib.py

clean-results:
	rm -f results/*.csv results/figures/*

manifest:
	sh code/10_manifest.sh
