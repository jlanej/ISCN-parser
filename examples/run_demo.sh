#!/usr/bin/env bash
# Demonstrates the batteries-included Docker image / installed CLI on the
# bundled example data. Produces PLINK .cnv/.fam/.cnv.map and a BED file.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-/tmp/iscn-demo}"
mkdir -p "$OUT"

echo "==> Converting microarray (arr) results to PLINK CNV"
iscn-parser --input "$HERE/array_samples.tsv" --format cnv --output "$OUT/array"

echo "==> Converting microarray (arr) results to BED"
iscn-parser --input "$HERE/array_samples.tsv" --format bed --output "$OUT/array.bed"

echo "==> Converting conventional karyotypes to PLINK CNV (GRCh37)"
iscn-parser --input "$HERE/karyotype_samples.tsv" --build GRCh37 --format cnv --output "$OUT/karyotype"

echo
echo "==> Wrote outputs to $OUT:"
ls -1 "$OUT"
echo
echo "==> Preview of $OUT/array.cnv:"
head -n 5 "$OUT/array.cnv"
