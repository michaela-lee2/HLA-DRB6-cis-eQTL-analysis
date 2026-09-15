#!/usr/bin/env bash

# 2. cis-eQTL Association Mapping
# Tests the association between HLA-DRB6 expression and cis-region
# SNPs (TSS ± 1 Mb) using PLINK linear regression (additive model).
# Significance threshold: genome-wide significance threshold P < 5e-8

set -euo pipefail

BFILE="/youngeun/biostat/data/geuvadis_genotypes"          # PLINK bfile prefix (.bed/.bim/.fam)
PHENO="/youngeun/biostat/data/HLA-DRB6_expression.pheno"   # PLINK phenotype file (FID IID PHENOTYPE)
GENE_TSS=32485154                        # HLA-DRB6 TSS position (GRCh38; modify if necessary)
CHR=6
WINDOW=1000000                           # ±1 Mb
OUT="/youngeun/biostat/results/HLA-DRB6_cis_eqtl"

START=$((GENE_TSS - WINDOW))
END=$((GENE_TSS + WINDOW))

mkdir -p results

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
awk 'NR==1 || $9 < 5e-8' "${OUT}.assoc.linear" > "${OUT}.genome_wide_significant.tsv"

echo "cis-eQTL mapping completed -> ${OUT}.assoc.linear"
echo "Significant SNPs -> ${OUT}.genome_wide_significant.tsv"

