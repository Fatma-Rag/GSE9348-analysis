if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")
# Bioconductor Packages
BiocManager::install(c("GEOquery", "affy", "limma", "arrayQualityMetrics", 
                       "AnnotationDbi", "hgu133plus2.db"))
#CRAN Packages
install.packages("dplyr")
#Load Required libraries Bioconductor
library(GEOquery)             
library(affy)                 
library(arrayQualityMetrics)  
library(dplyr)               

gse_data <- getGEO("GSE9348", GSEMatrix = TRUE)
# Extract expression data
expression_data <- exprs(gse_data$GSE9348_series_matrix.txt.gz)
#Extract feature data
feature_data <- fData(gse_data$GSE9348_series_matrix.txt.gz)
#Extract phenotype 
phenotype_data <- pData(gse_data$GSE9348_series_matrix.txt.gz)
sum(is.na(phenotype_data$source_name_ch1))
##############################################################
######Download Raw Data (Cel files)#########
#getGEOSuppFiles("GSE9348", baseDir = "Raw_Data", makeDirectory = TRUE)
# 1. Extract CEL files from the TAR archive
untar("D:\\GSE9348_RAW.tar", exdir = "D:\\GSE9348\\CEL")

# 2. Read CEL files
raw_data <- ReadAffy(celfile.path = "D:\\GSE9348\\CEL")

# 3. Quality control report
arrayQualityMetrics(expressionset = raw_data,
                    outdir = "D:\\GSE9348\\QC_before_norm",
                    force = TRUE,
                    do.logtransform = TRUE)

# 4. Normalize data using RMA
norm_data <- rma(raw_data)
# QC after data normalization 
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "D:\\GSE9348\\QC_after_norm",
                    force = TRUE)
# Extract normalized expression values into a data frame
processed_data <- as.data.frame(exprs(norm_data))

dim(processed_data)   # Dimensions: number of probes × number of samples
# Calculate median intensity per probe across samples
row_median <- rowMedians(as.matrix(processed_data))
# Visualize distribution of probe median intensities
hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")
# Set a threshold to remove low variance probes 
threshold <- 3.5 
abline(v = threshold, col = "black", lwd = 2)
# Select probes above threshold
indx <- row_median > threshold 
filtered_data <- processed_data[indx, ]
# Rename filtered expression data with sample metadata
colnames(filtered_data) <- rownames(phenotype_data)
# Overwrite processed data with filtered dataset
processed_data <- filtered_data 
# -----------------------------------
#### Phenotype Data Preparation ####
# -----------------------------------
class(phenotype_data$source_name_ch1) 
# Define experimental groups (normal vs cancer)
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("mucosa", "tumor"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)

