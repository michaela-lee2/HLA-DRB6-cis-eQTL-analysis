# 7. eQTL Effect Specificity Validation - Boxplots
#
# Compares HLA-DRB6 and HLA-DRB1 expression across
# rs9270694 genotype groups (0, 1, and 2).
# Corresponds to Figure 1_C.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

expr_file <- ifelse(
  length(args) >= 1,
  args[1],
  file.path(
    Sys.getenv("HOME"),
    "evlyn/eQTL/GEUVADIS/log2_expression.tsv"
  )
)

geno_file <- ifelse(
  length(args) >= 2,
  args[2],
  file.path(
    Sys.getenv("HOME"),
    "evlyn/eQTL/GEUVADIS/rs9270694_genotype.tsv"
  )
)
# Columns: sample, genotype (0/1/2)

out_file <- "results/figure1c_boxplot_DRB6_vs_DRB1_by_genotype.png"

expr <- fread(expr_file)
geno <- fread(geno_file)

# Extract HLA-DRB6 and HLA-DRB1 expression values
drb6 <- as.numeric(
  expr[GeneID == "HLA-DRB6", -1, with = FALSE]
)

drb1 <- as.numeric(
  expr[GeneID == "HLA-DRB1", -1, with = FALSE]
)

samples <- colnames(expr)[-1]

plot_df <- rbind(
  data.frame(
    sample = samples,
    expression = drb6,
    gene = "HLA-DRB6"
  ),
  data.frame(
    sample = samples,
    expression = drb1,
    gene = "HLA-DRB1"
  )
)

plot_df <- merge(
  plot_df,
  geno,
  by = "sample"
)

plot_df$genotype <- factor(
  plot_df$genotype,
  levels = c(0, 1, 2)
)

p <- ggplot(
  plot_df,
  aes(
    x = genotype,
    y = expression,
    fill = genotype
  )
) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~gene, scales = "free_y") +
  labs(
    title = "Expression by rs9270694 Genotype",
    x = "rs9270694 genotype (number of effect alleles)",
    y = expression(log[2](TPM + 1))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none"
  )

ggsave(
  out_file,
  p,
  width = 8,
  height = 5,
  dpi = 300
)

cat("Boxplot saved:", out_file, "\n")
