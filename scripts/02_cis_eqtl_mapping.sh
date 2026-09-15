#!/usr/bin/env bash

# 2. cis-eQTL Association Mapping
# Tests the association between HLA-DRB6 expression and cis-region
# SNPs (TSS ± 1 Mb) using PLINK linear regression (additive model).
# Significance threshold: genome-wide significance threshold P < 5e-8

set -euo pipefail

BFILE="$HOME/evlyn/eQTL/GEUVADIS/geuvadis_genotypes"
PHENO="$HOME/evlyn/eQTL/GEUVADIS/HLA-DRB6_expression.pheno"

GENE_TSS=32485154
CHR=6
WINDOW=1000000
OUT="results/HLA-DRB6_cis_eqtl"

START=$((GENE_TSS - WINDOW))
END=$((GENE_TSS + WINDOW))

mkdir -p "$(dirname "$OUT")"

plink \
  --bfile "$BFILE" \
  --chr "$CHR" \
  --from-bp "$START" \
  --to-bp "$END" \
  --pheno "$PHENO" \
  --linear \
  --ci 0.95 \
  --out "$OUT"

# Extract SNPs meeting the genome-wide significance threshold (P < 5e-8)
awk 'NR==1 || $9 < 5e-8' \
  "${OUT}.assoc.linear" \
  > "${OUT}.genome_wide_significant.tsv"

echo "cis-eQTL mapping completed -> ${OUT}.assoc.linear"
echo "Significant SNPs -> ${OUT}.genome_wide_significant.tsv"

