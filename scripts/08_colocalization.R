# 8. Colocalization 분석
# RFX5를 통한 광범위한 HLA Class II 억제 가설을 테스트하기 위해
# HLA-DRB6 eQTL 신호와 HLA-DRB1 eQTL 신호의 공동위치 여부를 coloc 패키지로 검증합니다.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)
drb6_assoc_file <- ifelse(length(args) >= 1, args[1], "/youngeun/biostat/results/HLA-DRB6_cis_eqtl.assoc.linear")
drb1_assoc_file <- ifelse(length(args) >= 2, args[2], "/youngeun/biostat/results/HLA-DRB1_cis_eqtl.assoc.linear")
n_samples <- ifelse(length(args) >= 3, as.numeric(args[3]), 462) # GEUVADIS 표본 수 (예시)

drb6 <- fread(drb6_assoc_file)
drb1 <- fread(drb1_assoc_file)

common <- merge(drb6, drb1, by = "SNP", suffixes = c("_DRB6", "_DRB1"))

dataset1 <- list(
  beta = common$BETA_DRB6,
  varbeta = (common$SE_DRB6)^2,
  snp = common$SNP,
  type = "quant",
  N = n_samples
)
dataset2 <- list(
  beta = common$BETA_DRB1,
  varbeta = (common$SE_DRB1)^2,
  snp = common$SNP,
  type = "quant",
  N = n_samples
)

result <- coloc.abf(dataset1 = dataset1, dataset2 = dataset2)

print(result$summary)
# PP.H4 (두 신호가 동일한 인과 변이를 공유할 사후확률)가 낮으면
# rs9270694의 효과가 DRB1으로 확장되지 않는다는(=DRB6 특이적) 증거로 해석

# rs9270694 단일 SNP의 DRB1 eQTL P-value 별도 확인 (논문: P=0.7722)
drb1_target <- drb1[SNP == "rs9270694"]
cat("\nrs9270694의 HLA-DRB1 eQTL P-value:", drb1_target$P, "\n")

saveRDS(result, "results/coloc_DRB6_DRB1.rds")
