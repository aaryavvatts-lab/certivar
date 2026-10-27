PY ?= .venv/bin/python
RAW = data/raw
PROC = data/processed

.PHONY: reproduce data annotate results figures manifest clean-results test

reproduce: data annotate test results figures manifest

test:
	cd code && ../$(PY) test_certivar.py

$(RAW)/variant_summary.txt.gz:
	mkdir -p $(RAW)
	curl -sL -o $@ https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

data: $(RAW)/variant_summary.txt.gz
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
	cd code && ../$(PY) 11_numbers.py

figures:
	cd code && ../$(PY) 09_figures.py
	cd code && ../$(PY) 12_tables.py
	cd code && ../$(PY) 16_short_bib.py

clean-results:
	rm -f results/*.csv results/figures/*

manifest:
	sh code/10_manifest.sh
