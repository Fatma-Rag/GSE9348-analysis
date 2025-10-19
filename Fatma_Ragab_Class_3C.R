rm(list = ls())
gc()
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma", "AnnotationDbi", "hgu133plus2.db"), ask = FALSE)
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
# Download the dataset from GEO
gse <- getGEO("GSE9348", GSEMatrix = TRUE)
length(gse)
if (is.list(gse)) {
  gse <- gse[[1]]
}
# Now extract the data
expr_matrix <- exprs(gse)       # Expression data
phenotype_data <- pData(gse)    # Sample metadata
feature_data <- fData(gse)
cat("Expression matrix:", dim(expr_matrix), "\n")
cat("Phenotype data:", dim(phenotype_data), "\n")
# -------------------------------------------------------------
# 3. Probe-to-Gene Mapping using AnnotationDbi
# -------------------------------------------------------------
probe_ids <- rownames(expr_matrix)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

gene_map_df <- gene_symbols %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  dplyr::rename(SYMBOL = 2)

# Handle missing or duplicated genes
gene_map_df <- gene_map_df %>% filter(!is.na(SYMBOL))

# Summary of duplicate mappings
dup_summary <- gene_map_df %>%
  group_by(SYMBOL) %>%
  summarise(probes_per_gene = n()) %>%
  arrange(desc(probes_per_gene))

duplicate_genes <- dup_summary %>% filter(probes_per_gene > 1)
dup_total <- sum(duplicate_genes$probes_per_gene)
cat("Total duplicate probe mappings:", dup_total, "\n")
# -------------------------------------------------------------
# 4. Merge Mapping with Expression Data and Collapse Duplicates
# -------------------------------------------------------------
merged_df <- expr_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PROBEID") %>%
  left_join(gene_map_df, by = "PROBEID") %>%
  filter(!is.na(SYMBOL))

expr_only <- merged_df %>% select(-PROBEID, -SYMBOL)
averaged_data <- limma::avereps(expr_only, ID = merged_df$SYMBOL)
# -------------------------------------------------------------
# 5. Define Groups (Normal vs Cancer)
# -------------------------------------------------------------
groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 labels = c("normal", "cancer"))

design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
fit <- lmFit(averaged_data, design)
dim(averaged_data)
dim(design)
dim(phenotype_data)
colnames(phenotype_data)
head(phenotype_data)
table(phenotype_data$source_name_ch1)
table(phenotype_data$title)
table(phenotype_data$characteristics_ch1)
groups <- factor(
  phenotype_data$source_name_ch1,
  levels = c("gastric mucosa", "gastric adenocarcinoma"),
  labels = c("normal", "cancer")
)
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)
dim(design)
fit <- lmFit(averaged_data, design)
contrast_matrix <- makeContrasts(cancer_vs_normal = cancer - normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# View top DE genes
topTable(fit2)
deg_results <- topTable(fit2, number = Inf, adjust = "BH")
write.csv(deg_results, "Results/All_DEGs_GSE9348.csv", row.names = TRUE)
sig_genes <- subset(deg_results, adj.P.Val < 0.05 & abs(logFC) > 1)
nrow(sig_genes)  # total number of significant DEGs
# -------------------------------------------------------------
# 1️⃣ Extract all DEGs and classify them by direction
# ---------------------------------------------------------------
deg_results <- topTable(fit2, number = Inf, adjust = "BH")

# Set significance thresholds
padj_cut <- 0.05
logfc_cut <- 1

# Add significance labels
deg_results$Regulation <- "Not Sig"
deg_results$Regulation[deg_results$adj.P.Val < padj_cut & deg_results$logFC > logfc_cut] <- "Upregulated"
deg_results$Regulation[deg_results$adj.P.Val < padj_cut & deg_results$logFC < -logfc_cut] <- "Downregulated"

# ---------------------------------------------------------------
# 2️⃣ Volcano Plot
# ---------------------------------------------------------------
dir.create("Result_Plots", showWarnings = FALSE)
png("Result_Plots/Volcano_Plot_GSE9348.png", width = 2000, height = 1600, res = 300)
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Sig" = "grey70")) +
  geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "black") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: GSE9348 (CRC Tumor vs Normal)",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Regulation"
  )
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Sig" = "grey70")) +
  theme_minimal() +
  labs(title = "Volcano Plot: GSE9348 (CRC Tumor vs Normal)")
# ---------------------------------------------------------------
# 3️⃣ Heatmap of Top 25 Significant DEGs
# ---------------------------------------------------------------
# Select top 25 by adjusted p-value
top_genes <- head(rownames(deg_results[order(deg_results$adj.P.Val), ]), 25)

# Extract expression for these genes
heatmap_data <- averaged_data[top_genes, ]

# Z-score normalize each gene (row)
heatmap_scaled <- t(scale(t(heatmap_data)))

# Create annotation for sample groups
annotation_col <- data.frame(Group = groups)
rownames(annotation_col) <- colnames(heatmap_scaled)

pheatmap(
  heatmap_scaled,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top 25 DEGs (CRC Tumor vs Normal)"
)


