# 5. Locus Zoom Plot 생성
# rs9270694를 중심으로 한 연관성 신호 및 LD(R^2) 구조를 시각화합니다.
# (그림 1_B에 해당)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
assoc_file <- ifelse(length(args) >= 1, args[1], "/youngeun/biostat/results/HLA-DRB6_cis_eqtl.assoc.linear")
ld_file    <- ifelse(length(args) >= 2, args[2], "/youngeun/biostat/results/rs9270694_LD.ld")
lead_snp   <- "rs9270694"
out_file   <- "results/locuszoom_rs9270694.png"

assoc <- fread(assoc_file)
ld    <- fread(ld_file)

# LD 값 병합 (SNP_B가 이웃 SNP, R2 = 리드 SNP와의 LD)
merged <- merge(assoc, ld[, .(SNP_B, R2)], by.x = "SNP", by.y = "SNP_B", all.x = TRUE)
merged[SNP == lead_snp, R2 := 1]
merged[is.na(R2), R2 := 0]
merged[, negLog10P := -log10(P)]

merged[, LD_bin := cut(R2,
  breaks = c(-0.01, 0.2, 0.4, 0.6, 0.8, 1.01),
  labels = c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
)]

p <- ggplot(merged, aes(x = BP, y = negLog10P, color = LD_bin)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_point(
    data = merged[SNP == lead_snp],
    shape = 18, size = 5, color = "purple"
  ) +
  scale_color_manual(values = c(
    "0.0-0.2" = "navy", "0.2-0.4" = "skyblue",
    "0.4-0.6" = "green3", "0.6-0.8" = "orange", "0.8-1.0" = "red"
  )) +
  labs(
    title = paste("Locus Zoom Plot -", lead_snp, "(HLA-DRB6 cis-eQTL)"),
    x = "Chromosome 6 position (bp)",
    y = expression(-log[10](italic(P))),
    color = expression(R^2)
  ) +
  theme_minimal(base_size = 13)

ggsave(out_file, p, width = 8, height = 5, dpi = 300)
cat("LocusZoom 플롯 저장 완료:", out_file, "\n")
