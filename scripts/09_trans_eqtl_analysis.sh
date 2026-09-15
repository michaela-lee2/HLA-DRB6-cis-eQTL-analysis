#!/usr/bin/env bash
# 9. Trans-eQTL 분석
# 다면발현(pleiotropy)을 배제하고 시스-특이성을 확인하기 위해
# 시스-영역(chr6 ±1Mb) 밖의 모든 유전자에 대해 rs9270694의 trans-eQTL 연관성을 탐색합니다.
# 엄격한 유의성 역치: P < 1e-12

set -euo pipefail

BFILE="/youngeun/biostat/data/geuvadis_genotypes"
LEAD_SNP="rs9270694"
PHENO_ALL_GENES="/youngeun/biostat/data/all_genes_expression.pheno"  # cis 영역(HLA-DRB6, DRB1 등) 제외한 전체 유전자 표현형
OUT="/youngeun/biostat/results/rs9270694_trans_eqtl"
P_THRESHOLD=1e-12

mkdir -p results

plink \
  --bfile "$BFILE" \
  --snp "$LEAD_SNP" \
  --pheno "$PHENO_ALL_GENES" \
  --all-pheno \
  --linear \
  --out "$OUT"

# 엄격한 임계값(P < 1e-12) 초과 유의미한 trans-eGene 검색
for f in "${OUT}".P*.assoc.linear; do
  awk -v thr="$P_THRESHOLD" 'NR==1 || $9 < thr' "$f"
done > "${OUT}.significant_trans_egenes.tsv"

N_SIG=$(($(wc -l < "${OUT}.significant_trans_egenes.tsv") - 1))
echo "Trans-eQTL 분석 완료. 유의미한 trans-eGene 수: ${N_SIG}"
echo "(논문 결과와 일치하면 0개여야 함 - HLA-DRB6 시스-특이성 확인)"
