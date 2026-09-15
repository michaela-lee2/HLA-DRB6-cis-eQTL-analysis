```r
# 8. Colocalization Analysis
#
# Tests the hypothesis of broad HLA Class II suppression through RFX5
# by evaluating whether HLA-DRB6 and HLA-DRB1 eQTL signals colocalize
# using the coloc package.

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

args <- commandArgs(trailingOnly = TRUE)

drb6_assoc_file <- ifelse(
  length(args) >= 1,
  args[1],
  "results/HLA-DRB6_cis_eqtl.assoc.linear"
)

drb1_assoc_file <- ifelse(
  length(args) >= 2,
  args[2],
  "results/HLA-DRB1_cis_eqtl.assoc.linear"
)

n_samples <- ifelse(
  length(args) >= 3,
  as.numeric(args[3]),
  462
)

drb6 <- fread(drb6_assoc_file)
drb1 <- fread(drb1_assoc_file)

common <- merge(
  drb6,
  drb1,
  by = "SNP",
  suffixes = c("_DRB6", "_DRB1")
)

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

result <- coloc.abf(
  dataset1 = dataset1,
  dataset2 = dataset2
)

# Display the coloc summary
print(result$summary)

# Independently inspect the HLA-DRB1 eQTL P-value for rs9270694
drb1_target <- drb1[SNP == "rs9270694"]

cat(
  "\nHLA-DRB1 eQTL P-value for rs9270694:",
  drb1_target$P,
  "\n"
)
```

