#Description: Using gene counts and metadata from RNA-seq experiment to perform deseq analysis. Equally cut-off gene counts less than 60, this is because gene with few reads mapping will have very large variations in fold change. though low amount of reads can either be signal of a lowly transcribed genes which maybe significant. However, there is a big chance that it's just some fluctuation of noise signal.

library(DESeq2)
library(ggplot2)
library(gridExtra) 
library(dplyr)
library(ggrepel)

# Directories
col_data <- "../results/featureCounts/metadata.csv"
counts_data <- "../results/featureCounts/counts_matrix.csv"
output_file <- "../results/featureCounts/dseq_results.txt"
output_plot <- "../results/featureCounts"

# Load the count data
countData <- read.csv(counts_data, row.names = 1)
# Convert to integers 
head(countData, 5)
#countData <- round(countData)

# Load the metadata
colData <- read.csv(col_data, row.names = 1)
colData

# Set the design formula for DESeq2
design <- "sample_name + sex + tissue + replicates"
design <- as.formula(paste("~", design))

colData$sample_name <- factor(colData$sample_name)
colData$sex <- factor(colData$sex)
colData$tissue <- factor(colData$tissue)
colData$replicates <- factor(colData$replicates)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = design)

# Remove genes with <= 60 counts in all samples
dds <- dds[rowSums(counts(dds)) > 60, ]

# Run DESeq2 Analysis
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)

# Write results to file
write.table(res, file = output_file, sep = "\t", quote = FALSE, col.names = NA)

#Dispersion plot: A dispersion plot in DESeq2 is a visual tool to check the #quality of your data before comparing gene expression levels.
#It shows how much the expression of each gene varies (dispersion) compared to #its average expression level (mean). Ideally, genes with higher expression #levels should have lower variability.
pdf(file.path(output_plot, "Dispersion_plot.pdf"))
plotDispEsts(dds, main = "Dispersion plot")
dev.off()

#Histogram of pvalues and adjusted pvalues. useful to show if there are truely #signifantly expressed genes.
pdf(file.path(output_plot, "histogram_pvalue.pdf"))
hist(res$pvalue, breaks = 20, col = "grey", main = "Histogram of p-values", xlab = "p-value")
dev.off()
pdf(file.path(output_plot, "histogram_adj_pvalue.pdf"))
hist(res$padj, breaks = 20, col = "grey", main = "Histogram of Adjusted p-values", xlab = "Adjusted p-value")
dev.off()


# Open a PDF device and output the plots
pdf(file.path(output_plot, "MA_plots_combined.pdf"), width = 14, height = 8)

# Set up a 2x3 grid layout for the plots
par(mfrow = c(2, 3), mar = c(4, 4, 4, 2)) 

# Plot 1: male vs female
res_sex <- results(dds, name = "sex_male_vs_female")
plotMA(res_sex, main = "Male vs Female", ylim = c(-5, 5))

# Plot 2: Saline vs Cocaine
res_sample <- results(dds, name = "sample_name_Saline_vs_Cocaine")
plotMA(res_sex, main = "Saline vs Cocaine", ylim = c(-5, 5))

# Plot 3: Tissue
res_tissue <- results(dds, name = "tissue_Prefrontal.Cortex_vs_nucleus.accumbens")
plotMA(res_tissue, main = "prefrontal cortex vs nucleus accumbens", ylim = c(-5, 5))

# Plot 4: Replicate 3 vs Replicate 1
res_replicate_2_vs_1 <- results(dds, name = "replicates_2_vs_1")
plotMA(res_replicate_2_vs_1, main = "replicate 2 vs Replicate 1", ylim = c(-5, 5))
dev.off()

# Perform variance stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# PCA plot using both Time and Replicates
pca_data <- plotPCA(vsd, intgroup = c("sex","tissue", "replicates"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Generate and save PCA plot
pdf(file.path(output_plot, "PCA_plot_sex_tissue_replicates.pdf"))
ggplot(pca_data, aes(x = PC1, y = PC2, color = sex, shape = replicates)) +
  geom_point(aes(size = tissue)) +
  theme_minimal() +
  labs(title = "PCA Plot: Sex, tissue and replicates",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  scale_color_manual(values = c("red", "green", "blue", "purple")) +
  scale_shape_manual(values = c(16, 17))
dev.off()

all_genes <- read.csv("../data/ref_genome/all_genes.csv") 
# Function to create a volcano plot with upregulation, downregulation, and other categories
create_volcano_plot <- function(res, title) {
  volcano_data <- as.data.frame(res)
  
  # Define significance thresholds
  log2fc_threshold <- 2 
  pval_threshold <- 0.05
  significant_threshold <- 10^(-200)

  #Categorize highly significant genes for labelling
  #volcano_data$category[volcano_data$pvalue >= label_threshold
  filtered_data <- as.data.frame(res)
  filtered_data$Geneid <- rownames(filtered_data)

  # Perform a left join to add gene_name to filtered_data
  filtered_data <- filtered_data %>%
    left_join(all_genes %>% select(gene_id, gene_name), 
              by = c("Geneid" = "gene_id")) %>%
    rename(significant_gname = gene_name)
  print(head(filtered_data))
  #filtered_data$significant_gname[is.na(filtered_data$significant_gname)] <- NA
  # Categorize genes as upregulated, downregulated, or other
  volcano_data$category <- "Not Significant"
  volcano_data$category[volcano_data$log2FoldChange >= log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Upregulated"
  volcano_data$category[volcano_data$log2FoldChange <= -log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Downregulated"
  filtered_data$category <- "Not Significant"
  filtered_data$category[volcano_data$log2FoldChange >= log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Upregulated"
  filtered_data$category[volcano_data$log2FoldChange <= -log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Downregulated"
  
  # Create the plot
  volcano_data$s_label <- ifelse(filtered_data$pvalue >= significant_threshold | filtered_data$category == "Not Significant", NA, as.character(filtered_data$significant_gname))

  volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = category)) +
    geom_point() +
    geom_text_repel(aes(label = s_label), na.rm = TRUE, box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50', color = volcano_data$label_color) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 p-value") +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-2, 2), col = "blue") +
    geom_hline(yintercept = -log10(0.05), col = "blue")
  return(volcano_plot)
}

# Generate volcano plots
volcano_sample <- create_volcano_plot(res_sample, "Saline vs Cocaine")
volcano_sex <- create_volcano_plot(res_sex, "Sex - male vs female")
volcano_tissue <- create_volcano_plot(res_tissue, "Tissue - prefrontal cortex - nuclus accumens")
volcano_replicates <- create_volcano_plot(res_replicate_2_vs_1, "Replicates 1&2")

# Combine the volcano plots into one page with a 2x3 grid layout
pdf(file.path(output_plot, "volcano_plots_combined.pdf"), width = 14, height = 8)

grid.arrange(volcano_sample, volcano_sex, volcano_tissue, volcano_replicates,
             ncol = 3)
dev.off()

# Find and label only genes identified by Wang X et al.
identify_specific_genes_expression<- function(res, title) {
  volcano_data <- as.data.frame(res)
  
  # Define significance thresholds
  log2fc_threshold <- 2 
  pval_threshold <- 0.05
  significant_threshold <- 10^(-200)

  #Categorize highly significant genes for labelling
  #volcano_data$category[volcano_data$pvalue >= label_threshold
  filtered_data <- as.data.frame(res)
  filtered_data$Geneid <- rownames(filtered_data)

  # Perform a left join to add gene_name to filtered_data
  filtered_data <- filtered_data %>%
    left_join(all_genes %>% select(gene_id, gene_name), 
              by = c("Geneid" = "gene_id")) %>%
    rename(significant_gname = gene_name)
  #filtered_data$significant_gname[is.na(filtered_data$significant_gname)] <- NA
  # Categorize genes as upregulated, downregulated, or other
  volcano_data$category <- "Not Significant"
  volcano_data$category[volcano_data$log2FoldChange >= log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Upregulated"
  volcano_data$category[volcano_data$log2FoldChange <= -log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Downregulated"
  filtered_data$category <- "Not Significant"
  filtered_data$category[volcano_data$log2FoldChange >= log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Upregulated"
  filtered_data$category[volcano_data$log2FoldChange <= -log2fc_threshold & volcano_data$pvalue < pval_threshold] <- "Downregulated"
  
  # Create the plot
  volcano_data$s_label <- ifelse(!(filtered_data$significant_gname %in% c("Fos", "Jun", "Il6", "Egr1")), NA, as.character(filtered_data$significant_gname))
  volcano_data$label_color <- ifelse(!is.na(volcano_data$s_label), "darkgreen", NA)

  volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), color = category)) +
    geom_point() +
    geom_text_repel(aes(label = s_label), na.rm = TRUE, box.padding = 0.1, point.padding = 0.1, segment.color = 'grey50', color = volcano_data$label_color) +
    theme_minimal() +
    labs(title = title, x = "Log2 Fold Change", y = "-Log10 p-value") +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    geom_vline(xintercept = c(-2, 2), col = "blue") +
    geom_hline(yintercept = -log10(0.05), col = "blue")
  return(volcano_plot)
}
volcano_tissue <- identify_specific_genes_expression()(res_tissue, "Tissue - prefrontal cortex - nuclus accumens")
# Combine the volcano plots into one page with a 2x3 grid layout
pdf(file.path(output_plot, "Identify_specific_genes.pdf"), width = 14, height = 8)
grid.arrange(volcano_tissue)
dev.off()
cat("Results and plots saved successfully to", output_plot, "\n")
