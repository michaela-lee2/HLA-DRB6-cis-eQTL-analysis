#!/usr/bin/env bash
# 2. cis-eQTL 연관성 매핑
# HLA-DRB6 발현과 시스-영역 SNP(TSS ± 1Mb) 사이의 연관성을
# PLINK 선형 회귀(가법 모델)로 테스트합니다.
# 유의성 기준: 전유전체 유의성 역치 P < 5e-8

set -euo pipefail

BFILE="/youngeun/biostat/data/geuvadis_genotypes"          # PLINK bfile prefix (.bed/.bim/.fam)
PHENO="/youngeun/biostat/data/HLA-DRB6_expression.pheno"   # PLINK 표현형 파일 (FID IID PHENOTYPE)
GENE_TSS=32485154                        # HLA-DRB6 TSS 위치 (GRCh38 기준, 필요시 수정)
CHR=6
WINDOW=1000000                           # ±1Mb
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

# 전유전체 유의성(P < 5e-8) SNP만 추출
awk 'NR==1 || $9 < 5e-8' "${OUT}.assoc.linear" > "${OUT}.genome_wide_significant.tsv"

echo "cis-eQTL 매핑 완료 -> ${OUT}.assoc.linear"
echo "유의미한 SNP -> ${OUT}.genome_wide_significant.tsv"
