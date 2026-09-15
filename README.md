# HLA-DRB6 cis-eQTL & Network Analysis Pipeline

This is a pipeline for identifying a cis-eQTL (rs9270694) regulating HLA-DRB6 expression using the GEUVADIS dataset, and for performing fine-mapping, specificity validation (trans-eQTL/colocalization), and network and pathway enrichment analyses.

Each script follows the same order as the **Methods** section of the paper.

## Pipeline Order

| Step | Script                        | Description                                                                            |
| ---- | ----------------------------- | -------------------------------------------------------------------------------------- |
| 1    | `01_normalize_expression.py`  | TPM/RPKM → log2(TPM+1) normalization                                                   |
| 2    | `02_cis_eqtl_mapping.sh`      | cis-eQTL mapping using PLINK linear regression (additive model; TSS ±1 Mb; P<5e-8)     |
| 3    | `03_manhattan_plot.R`         | Generate Manhattan plot and identify the strongest eQTL signal (Figure 1_A)            |
| 4    | `04_identify_lead_snp.py`     | Identify the SNP with the minimum P-value from the results (lead SNP: rs9270694)       |
| 5    | `05_ld_fine_mapping.sh`       | Calculate LD (R²) within ±500 kb of the lead SNP and extract candidates with R²≥0.8    |
| 6    | `06_locuszoom_plot.R`         | Generate a LocusZoom plot                                                              |
| 7    | `07_functional_annotation.py` | Functional annotation using the RegulomeDB / HaploReg APIs (motif changes, RFX5, etc.) |
| 8    | `08_specificity_boxplots.R`   | Boxplots of DRB6 vs DRB1 expression by genotype (0/1/2)                                |
| 9    | `09_colocalization.R`         | Colocalization (coloc) analysis of DRB6-DRB1 eQTL signals                              |
| 10   | `10_trans_eqtl_analysis.sh`   | Genome-wide trans-eQTL search (P<1e-12)                                                |
| 11   | `11_pearson_correlation.R`    | Pearson correlation analysis of DRB6-DRB1 expression (r=0.821)                         |
| 12   | `12_ppi_network_string.py`    | Construct an DRB1-centered PPI network using the STRING API                            |
| 13   | `13_enrichment_enrichr.py`    | KEGG2021/Reactome 2024 enrichment analysis using the Enrichr API (FDR<0.05)            |


## Dependencies

* PLINK 1.9/2.0 (cis-eQTL and LD analysis)
* Python 3.9+: `pandas`, `numpy`, `scipy`, `requests`
* R 4.x: `data.table`, `ggplot2`, `qqman`, `coloc`

```bash
pip install -r requirements.txt
Rscript -e 'install.packages(c("data.table","ggplot2","qqman","coloc"))'
```

## Data

Place the GEUVADIS genotype data (VCF/PLINK bfile) and RNA expression matrix (TPM) under `data/`, and modify the path variables at the top of each script as needed.

The original data are not included in the repository due to file size.

GEUVADIS: https://www.internationalgenome.org/data-portal/data-collection/geuvadis
