# 13. 요약 시각화
# (a) Manhattan plot - HLA-DRB6 cis-eQTL 결과 (그림 1_A)
# (b) HLA-DRB6 / HLA-DRB1 발현 분포 히스토그램 (그림 2_B, 2_C)

suppressPackageStartupMessages({
  library(data.table)
  library(qqman)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
assoc_file <- ifelse(length(args) >= 1, args[1], "/youngeun/biostat/results/HLA-DRB6_cis_eqtl.assoc.linear")
expr_file  <- ifelse(length(args) >= 2, args[2], "/youngeun/biostat/data/log2_expression.tsv")

## (a) Manhattan plot ---------------------------------------------------
assoc <- fread(assoc_file)
assoc <- assoc[!is.na(P)]

png("results/manhattan_HLA-DRB6.png", width = 1000, height = 600, res = 120)
manhattan(
  assoc,
  chr = "CHR", bp = "BP", snp = "SNP", p = "P",
  suggestiveline = -log10(5e-8),
  genomewideline = -log10(5e-8),
  main = "HLA-DRB6 cis-eQTL Manhattan Plot",
  highlight = "rs9270694"
)
dev.off()
cat("Manhattan plot 저장 완료: results/manhattan_HLA-DRB6.png\n")

## (b) 발현 히스토그램 ----------------------------------------------------
expr <- fread(expr_file)
drb6 <- as.numeric(expr[GeneID == "HLA-DRB6", -1, with = FALSE])
drb1 <- as.numeric(expr[GeneID == "HLA-DRB1", -1, with = FALSE])

hist_df <- rbind(
  data.frame(expression = drb6, gene = "HLA-DRB6"),
  data.frame(expression = drb1, gene = "HLA-DRB1")
)

p <- ggplot(hist_df, aes(x = expression, fill = gene)) +
  geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
  facet_wrap(~gene, scales = "free") +
  labs(
    title = "HLA-DRB6 / HLA-DRB1 발현 분포",
    x = expression(log[2](TPM + 1)), y = "샘플 수"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("results/expression_histograms.png", p, width = 8, height = 4.5, dpi = 300)
cat("발현 히스토그램 저장 완료: results/expression_histograms.png\n")
