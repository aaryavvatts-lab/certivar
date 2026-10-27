#!/bin/sh
# Record source release dates and checksums of every input file in data/MANIFEST.txt.
cd "$(dirname "$0")/.." || exit 1
{
  echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "ClinVar variant_summary.txt.gz $(cat data/raw/variant_summary.lastmod 2>/dev/null)"
  for f in data/raw/* data/processed/*; do
    [ -f "$f" ] && echo "$(shasum -a 256 "$f" | cut -d' ' -f1)  $f"
  done
} > data/MANIFEST.txt
