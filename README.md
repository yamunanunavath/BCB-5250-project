# Deep Learning–Based Gene Prediction for Alzheimer’s Disease

This repository contains a machine learning and deep learning–based pipeline for the identification of novel Alzheimer’s disease (AD) risk genes using integrated genomic and transcriptomic data. The project combines GWAS-derived SNP features with differential gene expression from human brain microarrays to train predictive models capable of identifying known and candidate AD-associated genes.

## 📊 Project Overview

Alzheimer’s disease is a complex neurodegenerative disorder. While known risk genes like APOE have been well studied, emerging evidence suggests the involvement of additional genomic factors. This project aims to uncover these using:

Genome-wide SNP significance from the NHGRI-EBI GWAS Catalog
Gene expression data from GEO datasets (GSE33000, GSE48350)
Known AD genes from DisGeNET and DISEASES databases
Our final model is a deep neural network (DNN), benchmarked against Random Forest and XGBoost, with external validation and GO enrichment analysis.


## 🔧 Technologies & Tools

Languages: R (primary)
Libraries: tidyverse, h2o, xgboost, randomForest, caret, clusterProfiler, GEOquery, SNPRelate, disgenet2r
Data:
SNP data from NHGRI-EBI GWAS Catalog
Expression data from GSE33000 (training) and GSE48350 (validation)
Known genes from DisGeNET & DISEASES

## 🧪 Pipeline Steps

1. Data Preparation
SNP Features: Minimum p-value, mean effect size, SNP count
Expression Features: log₂ fold change, p-value
Merged into a gene-level feature table
2. Modeling
Baseline Models: Random Forest and XGBoost
Primary Model: Dropout-regularized Deep Neural Network (via h2o)
Calibration: Platt scaling and isotonic regression
Validation: External test on GSE48350 (AUC = 0.706, p = 0.001)
3. Interpretation
Feature Importance via SHAP
GO Enrichment on top novel genes → cholesterol metabolism & immune signaling


## Below are some core R code snippets that replicate key stages of this project:

### 🔹 Install Required Packages

```r
install.packages(c("tidyverse", "caret", "randomForest", "xgboost", "h2o", "pROC"))
BiocManager::install(c("GEOquery", "SNPRelate", "org.Hs.eg.db", "AnnotationDbi", "disgenet2r", "clusterProfiler"))
```

## Initialize Environment

```r
library(tidyverse)
library(h2o); h2o.init(max_mem_size = "4G")
```

## SNP Feature Aggregation

```r
snps <- read.delim("ad_gwas_snps_filtered.tsv")
snps_long <- snps %>%
  separate_rows(reported_genes, sep = ",\\s*") %>%
  group_by(reported_genes) %>%
  summarise(min_p = min(p_value), mean_beta = mean(effect_size), snp_count = n(), .groups = "drop") %>%
  rename(gene = reported_genes)

```
## Process Expression Data (GSE33000)

```r
gse <- GEOquery::getGEO("GSE33000")[[1]]
expr <- exprs(gse)
# collapse to gene symbols with AnnotationDbi + org.Hs.eg.db
```
## Compute Differential Expression
```r
# Assume expr_gene has rows as gene symbols
de_stats <- tibble(
  gene = rownames(expr_gene),
  log2FC = rowMeans(expr_gene[, group1]) - rowMeans(expr_gene[, group2]),
  de_p = apply(expr_gene, 1, function(x) t.test(x[group1], x[group2])$p.value)
)
```
## Label Known Genes

```r
alz_genes <- read.delim("DISEASES_Summary_GDA_RGD_HUMAN_C0002395.tsv")$Gene
features <- snps_long %>%
  inner_join(de_stats, by = "gene") %>%
  mutate(label = if_else(gene %in% alz_genes, "Known", "Novel"))

```
## Train Deep Neural Network

```r
df_h2o <- as.h2o(features)
dnn <- h2o.deeplearning(x = c("min_p", "mean_beta", "snp_count", "log2FC", "de_p"),
                        y = "label",
                        training_frame = df_h2o,
                        activation = "RectifierWithDropout",
                        hidden = c(64, 32),
                        epochs = 50)

```
## GO Enrichment of Top Genes

```r
library(clusterProfiler)
top_genes <- features %>% arrange(desc(prob)) %>% slice_head(n = 20) %>% pull(gene)
enrich_result <- enrichGO(gene = top_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

```


## 🧠 Biological Implications

This integrative approach demonstrates that deep learning can effectively capture non-linear multi-omic interactions for disease gene prioritization. The novel candidate genes identified warrant further biological validation and may represent potential therapeutic targets.

## 📬 Contact

Yamuna Nunavath
Graduate Student, Saint Louis University
📧 yamuna.nunavath@slu.edu
📍 St. Louis, MO




