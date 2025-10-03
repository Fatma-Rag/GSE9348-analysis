# GSE9348-analysis
Microarray analysis of early stage  colorectal cancer (CRC) patient's tumor

This repository contains R scripts for preprocessing, QC, normalization, 
and differential expression analysis of the gastric cancer dataset **GSE9348**.

## Workflow
1. Download CEL files and extract metadata
2. QC before and after normalization
3. Normalize with RMA
4. Filter low-expression probes
5. Differential expression with limma

## How to Run
```R
source("scripts/01_download_data.R")
source("scripts/02_preprocessing_QC.R")
source("scripts/03_normalization.R")
source("scripts/04_filtering.R")
source("scripts/05_DE_analysis.R")

