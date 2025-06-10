#Batch1: DESeq2 Analysis (DMSO vs TSA/E6R at 16h)
library(tximport)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(EnhancedVolcano)

#Set paths for batch1
batch1_dir <- "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1"
output_dir <- file.path(batch1_dir, "DESeq2_results_batch1")
dir.create(output_dir, showWarnings = FALSE)

#Get sample directories
sample_dirs <- list.dirs(batch1_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[!grepl("DESeq2_results_batch1", sample_dirs)]
sample_names <- basename(sample_dirs)

#Create file paths & metadata
files <- file.path(sample_dirs, "quant.sf")
names(files) <- sample_names

metadata <- data.frame(
  sample = sample_names,
  condition = ifelse(grepl("DMSO", sample_names), "DMSO",
                    ifelse(grepl("TSA", sample_names), "TSA",
                          ifelse(grepl("E6R", sample_names), "E6R", NA))),
  replicate = gsub(".*([0-9]+)_batch1", "\\1", sample_names)
)

#Remove any samples with NA (e.g., ES samples)
valid_samples <- !is.na(metadata$condition)
files <- files[valid_samples]
metadata <- metadata[valid_samples, ]

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
dds <- DESeq(dds)

#Generate results
contrasts <- list(
  TSA = c("condition", "TSA", "DMSO"),
  E6R = c("condition", "E6R", "DMSO")
)

#Check available coefficients in model 
resultsNames(dds)

#Define coefficients (not contrasts) for apeglm 
coefficients_to_test <- c("condition_TSA_vs_DMSO", "condition_E6R_vs_DMSO")

#Generate results with apeglm shrinkage
results_list <- lapply(coefficients_to_test, function(coef) {
  res <- results(dds, name = coef, alpha = 0.05)
  res_shrunk <- lfcShrink(dds, coef = coef, type = "apeglm")
  return(list(raw = res, shrunk = res_shrunk))
})
names(results_list) <- c("TSA", "E6R")

#Save results
for(name in names(results_list)){
  write.csv(as.data.frame(results_list[[name]]$shrunk),
            file.path(output_dir, paste0("batch1_", name, "_vs_DMSO.csv")))
}

#results_list is a named list with two elements: "TSA" & "E6R".
#To retrieve the shrunken results for TSA 
res_TSA_shrunk <- results_list$TSA$shrunk
res_TSA_df <- as.data.frame(res_TSA_shrunk)

#Filter for significant DEGs (padj < 0.05, |log2FC| >= 1)
res_TSA_sig <- res_TSA_df[!is.na(res_TSA_df$padj) & 
                           res_TSA_df$padj < 0.05 & 
                           abs(res_TSA_df$log2FoldChange) >= 1, ]


#To retrieve the shrunken results for E6R
res_E6R_shrunk <- results_list$E6R$shrunk
res_E6R_df <- as.data.frame(res_E6R_shrunk)

#Filter for significant DEGs (padj < 0.05, |log2FC| >= 1)
res_E6R_sig <- res_E6R_df[!is.na(res_E6R_df$padj) & 
                         res_E6R_df$padj < 0.05 & 
                         abs(res_E6R_df$log2FoldChange) >= 1, ]


#Save filtered TSA vs DMSO results
write.csv(res_TSA_sig, file = file.path(output_dir, "batch1_TSA_vs_DMSO_filtered.csv"), row.names = TRUE)

#Save filtered E6R vs DMSO results
write.csv(res_E6R_sig, file = file.path(output_dir, "batch1_E6R_vs_DMSO_filtered.csv"), row.names = TRUE)


#Diagnostic plots
vsd <- vst(dds, blind = FALSE)

#PCA plot colored by condition
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("Batch1: PCA by Condition") +
  theme_minimal()

ggsave(file.path(output_dir, "batch1_PCA.png"), pca_plot, width = 8, height = 6)

#Volcano plots
plot_volcano <- function(res, title) {
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange',
                  y = 'padj',
                  title = title,
                  pCutoff = 0.05,
                  FCcutoff = 1)
}

walk2(results_list, names(results_list), ~{
  ggsave(file.path(output_dir, paste0("batch1_volcano_", .y, ".png")),
         plot_volcano(.x$shrunk, paste("Batch1:", .y, "vs DMSO")),
         width = 8, height = 6)
})

#MA plots
png(file.path(output_dir, "batch1_MA_plots.png"), width = 1000, height = 500)
par(mfrow = c(1,2))
plotMA(results_list$TSA$shrunk, main = "TSA vs DMSO")
plotMA(results_list$E6R$shrunk, main = "E6R vs DMSO")
dev.off()

#Session info
sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()
