#!/usr/bin/env bash
# 4. 정밀 지도화(Fine-mapping) - LD 기반 후보 변이 추출
# 선도 SNP(rs9270694)와 ±500kb 내 인접 SNP 간 연관불균형(R^2)을 계산하고
# R^2 >= 0.8인 SNP를 인과 변이 후보로 추출합니다.

set -euo pipefail

BFILE="/youngeun/biostat/data/geuvadis_genotypes"
LEAD_SNP="rs9270694"
WINDOW_KB=500
R2_THRESHOLD=0.8
OUT="/youngeun/biostat/results/rs9270694_LD"

mkdir -p results

plink \
  --bfile "$BFILE" \
  --r2 \
  --ld-snp "$LEAD_SNP" \
  --ld-window-kb "$WINDOW_KB" \
  --ld-window 99999 \
  --ld-window-r2 0 \
  --out "$OUT"

# R^2 >= 0.8 후보 변이만 추출 (LD 블록 내 causal variant 후보)
awk -v thr="$R2_THRESHOLD" 'NR==1 || $7 >= thr' "${OUT}.ld" > "${OUT}.candidates_r2_0.8.tsv"

echo "LD 계산 완료 -> ${OUT}.ld"
echo "R^2>=${R2_THRESHOLD} 후보 변이 -> ${OUT}.candidates_r2_0.8.tsv"
