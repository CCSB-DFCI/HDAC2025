#Batch3: DESeq2 Analysis (DMSO vs TSA/E6R at 16h)
library(tximport)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

#Set paths for batch3
batch3_dir <- "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3"
output_dir <- file.path(batch3_dir, "DESeq2_results_batch3")
dir.create(output_dir, showWarnings = FALSE)

#Get sample directories
sample_dirs <- list.dirs(batch3_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[!grepl("DESeq2_results_batch3", sample_dirs)]  #Exclude results dir if it exists
sample_names <- basename(sample_dirs)

#Create file paths/metadata
files <- file.path(sample_dirs, "quant.sf")
names(files) <- sample_names

metadata <- data.frame(
  sample = sample_names,
  condition = ifelse(grepl("DMSO", sample_names), "DMSO",
                    ifelse(grepl("TSA", sample_names), "TSA",
                          ifelse(grepl("E6R", sample_names), "E6R", NA))),
  replicate = gsub(".*([0-9]+)_16hrs", "\\1", sample_names)  #Extract replicate number
)

#Remove any samples with NA (if any)
valid_samples <- !is.na(metadata$condition)
files <- files[valid_samples]
metadata <- metadata[valid_samples, ]
rownames(metadata) <- metadata$sample

print("Sample metadata for batch3:")
print(metadata)

#Import Salmon data
txi <- tximport(files, type = "salmon", txOut = TRUE)

#Create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, 
                               colData = metadata,
                               design = ~ condition)

#Pre-filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Set DMSO as reference
dds$condition <- relevel(dds$condition, ref = "DMSO")

#Run DESeq2
dds <- DESeq(dds, fitType = "local")  #Using local to avoid dispersion warnings

#Check available coefficients
print(resultsNames(dds))

#Get coefficients for shrinkage
coefficients_to_test <- c("condition_TSA_vs_DMSO", "condition_E6R_vs_DMSO")

#Generate results with apeglm shrinkage
results_list <- lapply(coefficients_to_test, function(coef) {
  res <- results(dds, name = coef, alpha = 0.05)
  res_shrunk <- lfcShrink(dds, coef = coef, type = "apeglm")
  return(list(raw = res, shrunk = res_shrunk))
})
names(results_list) <- c("TSA", "E6R")

#Convert to data frames & filter
res_TSA_df <- as.data.frame(results_list$TSA$shrunk)
res_E6R_df <- as.data.frame(results_list$E6R$shrunk)

#Filter for significant DEGs (padj < 0.05, |log2FC| >= 1)
res_TSA_sig <- res_TSA_df[!is.na(res_TSA_df$padj) & 
                          res_TSA_df$padj < 0.05 & 
                          abs(res_TSA_df$log2FoldChange) >= 1, ]

res_E6R_sig <- res_E6R_df[!is.na(res_E6R_df$padj) & 
                          res_E6R_df$padj < 0.05 & 
                          abs(res_E6R_df$log2FoldChange) >= 1, ]

#Print filtering results
cat("Batch3 filtering results:\n")
cat("TSA - Before filtering:", nrow(res_TSA_df), "After filtering:", nrow(res_TSA_sig), "\n")
cat("E6R - Before filtering:", nrow(res_E6R_df), "After filtering:", nrow(res_E6R_sig), "\n")

#Save results
write.csv(res_TSA_df, file.path(output_dir, "batch3_TSA_vs_DMSO_all.csv"))
write.csv(res_E6R_df, file.path(output_dir, "batch3_E6R_vs_DMSO_all.csv"))
write.csv(res_TSA_sig, file.path(output_dir, "batch3_TSA_vs_DMSO_filtered.csv"))
write.csv(res_E6R_sig, file.path(output_dir, "batch3_E6R_vs_DMSO_filtered.csv"))

#Generate plots
vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("Batch3: PCA by Condition") +
  theme_minimal()
ggsave(file.path(output_dir, "batch3_PCA.png"), pca_plot, width = 8, height = 6)

#Volcano plots
pdf(file.path(output_dir, "batch3_TSA_vs_DMSO_volcano.pdf"), width = 8, height = 6)
EnhancedVolcano(results_list$TSA$shrunk,
    lab = rownames(results_list$TSA$shrunk),
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Batch3: TSA vs DMSO',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 3.5,
    labSize = 0)
dev.off()

pdf(file.path(output_dir, "batch3_E6R_vs_DMSO_volcano.pdf"), width = 8, height = 6)
EnhancedVolcano(results_list$E6R$shrunk,
    lab = rownames(results_list$E6R$shrunk),
    x = 'log2FoldChange',
    y = 'padj',
    title = 'Batch3: E6R vs DMSO',
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 3.5,
    labSize = 0)
dev.off()

#MA plots
png(file.path(output_dir, "batch3_MA_plots.png"), width = 1000, height = 500)
par(mfrow = c(1,2))
plotMA(results_list$TSA$shrunk, main = "TSA vs DMSO")
plotMA(results_list$E6R$shrunk, main = "E6R vs DMSO")
dev.off()

print("Batch3 DESeq2 analysis complete!")
