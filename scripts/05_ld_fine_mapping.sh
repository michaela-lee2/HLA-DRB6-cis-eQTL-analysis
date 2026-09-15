#!/usr/bin/env bash

# 4. LD-based candidate variant analysis
#
# Calculates linkage disequilibrium (R^2) between the lead SNP
# (rs9270694) and neighboring SNPs within ±500 kb.
#
# Variants with R^2 >= 0.8 are retained as candidate variants
# for downstream functional annotation.
#
# Note:
# LD-based filtering does not establish causality.

set -euo pipefail

BFILE="$HOME/evlyn/eQTL/GEUVADIS/geuvadis_genotypes"
LEAD_SNP="rs9270694"
WINDOW_KB=500
R2_THRESHOLD=0.8
OUT="results/rs9270694_LD"

mkdir -p "$(dirname "$OUT")"

# Check required PLINK files
for ext in bed bim fam; do
    if [[ ! -f "${BFILE}.${ext}" ]]; then
        echo "Error: Missing PLINK file ${BFILE}.${ext}" >&2
        exit 1
    fi
done

plink \
  --bfile "$BFILE" \
  --r2 \
  --ld-snp "$LEAD_SNP" \
  --ld-window-kb "$WINDOW_KB" \
  --ld-window 99999 \
  --ld-window-r2 0 \
  --out "$OUT"

# Extract candidate variants with R^2 >= 0.8
awk -v thr="$R2_THRESHOLD" \
  'NR == 1 || $7 >= thr' \
  "${OUT}.ld" \
  > "${OUT}.candidates_r2_0.8.tsv"

echo "LD analysis completed -> ${OUT}.ld"
echo "Candidate variants with R^2 >= ${R2_THRESHOLD} -> ${OUT}.candidates_r2_0.8.tsv"
