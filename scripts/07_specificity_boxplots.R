# 7. eQTL 효과 특이성 검증 - 박스플롯
# rs9270694 유전자형 그룹(0,1,2)에 따른 HLA-DRB6 및 HLA-DRB1 발현 비교
# (그림 1_C에 해당)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
expr_file <- ifelse(length(args) >= 1, args[1], "/youngeun/biostat/data/log2_expression.tsv")
geno_file <- ifelse(length(args) >= 2, args[2], "/youngeun/biostat/data/rs9270694_genotype.tsv") # columns: sample, genotype(0/1/2)
out_file  <- "results/boxplot_DRB6_vs_DRB1_by_genotype.png"

expr <- fread(expr_file)          # rows = genes, columns = samples
geno <- fread(geno_file)          # sample, genotype

drb6 <- as.numeric(expr[GeneID == "HLA-DRB6", -1, with = FALSE])
drb1 <- as.numeric(expr[GeneID == "HLA-DRB1", -1, with = FALSE])
samples <- colnames(expr)[-1]

plot_df <- rbind(
  data.frame(sample = samples, expression = drb6, gene = "HLA-DRB6"),
  data.frame(sample = samples, expression = drb1, gene = "HLA-DRB1")
)
plot_df <- merge(plot_df, geno, by = "sample")
plot_df$genotype <- factor(plot_df$genotype, levels = c(0, 1, 2))

p <- ggplot(plot_df, aes(x = genotype, y = expression, fill = genotype)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~gene, scales = "free_y") +
  labs(
    title = "rs9270694 유전자형에 따른 발현 비교",
    x = "rs9270694 genotype (효과 대립유전자 개수)",
    y = expression(log[2](TPM + 1))
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(out_file, p, width = 8, height = 5, dpi = 300)
cat("박스플롯 저장 완료:", out_file, "\n")
