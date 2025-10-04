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
                    outdir = "D:\\GSE9348\\QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
# 4. Normalize data using RMA
norm_data <- rma(raw_data)
# QC after data normalization 
arrayQualityMetrics(expressionset = norm_data,
                    outdir = "D:\\GSE9348\\QC_Normalized_Data",
                    force = TRUE)
processed_data <- as.data.frame(exprs(norm_data))
dim(processed_data)
############################################################
# Filter Low Variance Probes
############################################################
row_median <- rowMedians(as.matrix(processed_data))

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

threshold <- 4
abline(v = threshold, col = "black", lwd = 2)

indx <- row_median > threshold
filtered_data <- processed_data[indx, ]
colnames(filtered_data) <- rownames(phenotype_data)
processed_data <- filtered_data

############################################################
# Define Groups (Normal vs Cancer)
############################################################
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("mucosa", "tumor"),
                 labels = c("normal", "cancer"))
# Inspect unique values in phenotype_data$source_name_ch1
unique(phenotype_data$source_name_ch1)
# Define groups based on keywords in the metadata
groups <- ifelse(grepl("tumor", phenotype_data$source_name_ch1, ignore.case = TRUE),
                 "cancer", "normal")

# Convert to factor
groups <- factor(groups, levels = c("normal", "cancer"))

# Verify
table(groups)

############################################################
# Boxplot After Normalization
############################################################
boxplot(exprs(norm_data),
        main = "Boxplot of Normalized Expression Data",
        xlab = "Samples",
        ylab = "Log2 Expression",
        col = "lightblue",
        las = 2,        # rotate x-axis labels
        outline = FALSE # ignore extreme outliers
)
############################################################
# 10. PCA Plot After Normalization
############################################################
pca <- prcomp(t(exprs(norm_data)), scale. = TRUE)

plot(pca$x[,1], pca$x[,2],
     col = as.factor(groups),
     pch = 19,
     xlab = paste0("PC1 (", round(100*summary(pca)$importance[2,1], 1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca)$importance[2,2], 1), "%)"),
     main = "PCA of Normalized Expression Data")

