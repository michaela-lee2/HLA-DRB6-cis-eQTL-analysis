# 3. Manhattan Plot
#
# Visualizes the HLA-DRB6 cis-eQTL association results.
# The plot is used to identify the strongest association signal
# and guide lead SNP selection.
# Corresponds to Figure 1_A.

suppressPackageStartupMessages({
  library(data.table)
  library(qqman)
})

args <- commandArgs(trailingOnly = TRUE)

assoc_file <- ifelse(
  length(args) >= 1,
  args[1],
  "results/HLA-DRB6_cis_eqtl.assoc.linear"
)

out_file <- "results/figure1a_manhattan_plot.png"

assoc <- fread(assoc_file)
assoc <- assoc[!is.na(P)]

png(
  out_file,
  width = 1000,
  height = 600,
  res = 120
)

manhattan(
  assoc,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  suggestiveline = -log10(5e-8),
  genomewideline = -log10(5e-8),
  main = "HLA-DRB6 cis-eQTL Manhattan Plot",
  highlight = "rs9270694"
)

dev.off()

cat(
  "Manhattan plot saved:",
  out_file,
  "\n"
)
