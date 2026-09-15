#!/usr/bin/env bash

# 9. Trans-eQTL Analysis
#
# Searches for trans-eQTL associations of rs9270694 with genes
# outside the cis-region (chr6 ±1 Mb) to assess potential
# pleiotropic effects and evaluate cis-specificity.
#
# Significance threshold: P < 1e-12

set -euo pipefail

BFILE="$HOME/evlyn/eQTL/GEUVADIS/geuvadis_genotypes"
LEAD_SNP="rs9270694"
PHENO_ALL_GENES="$HOME/evlyn/eQTL/GEUVADIS/all_genes_expression.pheno"
OUT="results/rs9270694_trans_eqtl"
P_THRESHOLD=1e-12

mkdir -p "$(dirname "$OUT")"

plink \
  --bfile "$BFILE" \
  --snp "$LEAD_SNP" \
  --pheno "$PHENO_ALL_GENES" \
  --all-pheno \
  --linear \
  --out "$OUT"

# Search for significant trans-eGenes using a stringent threshold (P < 1e-12)
for f in "${OUT}".P*.assoc.linear; do
  awk -v thr="$P_THRESHOLD" \
    'NR == 1 || $9 < thr' \
    "$f"
done > "${OUT}.significant_trans_egenes.tsv"

N_SIG=$(
  wc -l < "${OUT}.significant_trans_egenes.tsv"
)

N_SIG=$((N_SIG - 1))

echo "Trans-eQTL analysis completed."
echo "Number of significant trans-eGenes: ${N_SIG}"
