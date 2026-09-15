# 10. HLA-DRB6 - HLA-DRB1 발현 피어슨 상관분석
# 대리 유전자(surrogate gene) 사용의 타당성을 확립합니다. (그림 2_D, r=0.821)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
expr_file <- ifelse(length(args) >= 1, args[1], "/youngeun/biostat/data/log2_expression.tsv")
out_plot  <- "results/correlation_DRB6_DRB1.png"

expr <- fread(expr_file)
drb6 <- as.numeric(expr[GeneID == "HLA-DRB6", -1, with = FALSE])
drb1 <- as.numeric(expr[GeneID == "HLA-DRB1", -1, with = FALSE])

cor_test <- cor.test(drb6, drb1, method = "pearson")
cat(sprintf("Pearson r = %.3f, P-value = %.3e\n", cor_test$estimate, cor_test$p.value))

df <- data.frame(DRB6 = drb6, DRB1 = drb1)

p <- ggplot(df, aes(x = DRB6, y = DRB1)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE) +
  annotate(
    "text", x = min(df$DRB6), y = max(df$DRB1),
    hjust = 0, vjust = 1,
    label = sprintf("r = %.3f", cor_test$estimate)
  ) +
  labs(
    title = "HLA-DRB6 vs HLA-DRB1 발현 상관관계",
    x = expression("HLA-DRB6" ~ log[2](TPM + 1)),
    y = expression("HLA-DRB1" ~ log[2](TPM + 1))
  ) +
  theme_minimal(base_size = 13)

ggsave(out_plot, p, width = 6, height = 5, dpi = 300)
cat("상관관계 플롯 저장 완료:", out_plot, "\n")
