library(tximport)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(EnhancedVolcano)
library(purrr)

#Set paths for batch2
batch2_dir <- "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2"
output_dir <- file.path(batch2_dir, "DESeq2_results_batch2")
dir.create(output_dir, showWarnings = FALSE)

#Get sample directories
sample_dirs <- list.dirs(batch2_dir, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[!grepl("DESeq2_results_batch2", sample_dirs)]
sample_names <- basename(sample_dirs)

#Create file paths & metadata
files <- file.path(sample_dirs, "quant.sf")
names(files) <- sample_names

metadata <- data.frame(
  sample = sample_names,
  condition = ifelse(grepl("DMSO", sample_names), "DMSO",
                    ifelse(grepl("TSA001", sample_names), "TSA001",
                      ifelse(grepl("TSA1", sample_names), "TSA1",
                        ifelse(grepl("E6R", sample_names), "E6R", NA)))),
  replicate = gsub(".*([0-9]+)$", "\\1", sample_names)
)

#Remove any samples with NA (e.g., ES samples)
valid_samples <- !is.na(metadata$condition)
files <- files[valid_samples]
metadata <- metadata[valid_samples, ]

#Import Salmon data (transcript-level)
txi <- tximport(files, type = "salmon", txOut = TRUE)

#Create DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, 
                               colData = metadata,
                               design = ~ condition)

#Pre-filter low-count transcripts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Set DMSO as reference
dds$condition <- relevel(dds$condition, ref = "DMSO")

#Run DESeq2
dds <- DESeq(dds)

#Generate results for both TSA groups vs DMSO
contrasts <- list(
  TSA001 = c("condition", "TSA001", "DMSO"),
  TSA1 = c("condition", "TSA1", "DMSO"),
  E6R = c("condition", "E6R", "DMSO")
)

#Check available coefficients in model
resultsNames(dds)

#Define coefficients (not contrasts) for apeglm
coefficients_to_test <- c("condition_TSA001_vs_DMSO", "condition_TSA1_vs_DMSO", "condition_E6R_vs_DMSO")

#Generate results with apeglm shrinkage
results_list <- lapply(coefficients_to_test, function(coef) {
  res <- results(dds, name = coef, alpha = 0.05)
  res_shrunk <- lfcShrink(dds, coef = coef, type = "apeglm")
  return(list(raw = res, shrunk = res_shrunk))
})
names(results_list) <- c("TSA001", "TSA1", "E6R")

#Save results
for(name in names(results_list)){
  write.csv(as.data.frame(results_list[[name]]$shrunk),
            file.path(output_dir, paste0("batch2_", name, "_vs_DMSO.csv")))
}

#Example: Filter significant DEGs for TSA001 vs DMSO
res_TSA001_shrunk <- results_list$TSA001$shrunk
res_TSA001_df <- as.data.frame(res_TSA001_shrunk)
res_TSA001_sig <- res_TSA001_df[!is.na(res_TSA001_df$padj) &
                                res_TSA001_df$padj < 0.05 &
                                abs(res_TSA001_df$log2FoldChange) >= 1, ]
write.csv(res_TSA001_sig, file = file.path(output_dir, "batch2_TSA001_vs_DMSO_filtered.csv"), row.names = TRUE)

# TSA1 vs DMSO
res_TSA1_shrunk <- results_list$TSA1$shrunk
res_TSA1_df <- as.data.frame(res_TSA1_shrunk)
res_TSA1_sig <- res_TSA1_df[!is.na(res_TSA1_df$padj) &
                            res_TSA1_df$padj < 0.05 &
                            abs(res_TSA1_df$log2FoldChange) >= 1, ]
write.csv(res_TSA1_sig, file = file.path(output_dir, "batch2_TSA1_vs_DMSO_filtered.csv"), row.names = TRUE)

#E6R vs DMSO
res_E6R_shrunk <- results_list$E6R$shrunk
res_E6R_df <- as.data.frame(res_E6R_shrunk)
res_E6R_sig <- res_E6R_df[!is.na(res_E6R_df$padj) &
                          res_E6R_df$padj < 0.05 &
                          abs(res_E6R_df$log2FoldChange) >= 1, ]
write.csv(res_E6R_sig, file = file.path(output_dir, "batch2_E6R_vs_DMSO_filtered.csv"), row.names = TRUE)

#Define file names for each comparison
comparisons <- c("TSA001", "TSA1", "E6R")
unfiltered_files <- file.path(output_dir, paste0("batch2_", comparisons, "_vs_DMSO.csv"))
filtered_files   <- file.path(output_dir, paste0("batch2_", comparisons, "_vs_DMSO_filtered.csv"))

#Count total transcripts in each file:
#Total transcripts (unfiltered): 94,861 (all tested)
for(i in seq_along(comparisons)) {
  unfiltered_df <- read.csv(unfiltered_files[i], row.names=1)
  filtered_df   <- read.csv(filtered_files[i], row.names=1)
  
  cat(sprintf("%s vs DMSO:\n", comparisons[i]))
  cat(sprintf("  Total transcripts (unfiltered): %d\n", nrow(unfiltered_df)))
  cat(sprintf("  Significant transcripts (filtered): %d\n", nrow(filtered_df)))
  cat(sprintf("  Percentage significant: %.2f%%\n\n", (nrow(filtered_df)/nrow(unfiltered_df))*100))
}

#TSA1 is an internal control (same as TSA from batch1)- hence, we see an extremely high overlap & very similar #DEGs. 

#Diagnostic plots for batch 2
vsd <- vst(dds, blind = FALSE)

#PCA plot colored by condition
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  ggtitle("Batch2: PCA by Condition") +
  theme_minimal()

ggsave(file.path(output_dir, "batch2_PCA.jpeg"), pca_plot, width = 8, height = 6)


#Define a volcano plot function to depict the DEGs per condition (TSA001, TSA1 & E6R) 
plot_volcano <- function(res, title) {
  EnhancedVolcano(
    res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    title = title,
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0
  )
}

#Loop & save as PDF for each comparison
walk2(results_list, names(results_list), ~{
  pdf(file.path(output_dir, paste0("batch2_volcano_", .y, ".pdf")), width = 8, height = 6)
  print(plot_volcano(.x$shrunk, paste("Batch2:", .y, "vs DMSO")))
  dev.off()
})

#MA plots for batch 2
png(file.path(output_dir, "batch2_MA_plots.png"), width = 1000, height = 500)
par(mfrow = c(1,3))
plotMA(results_list$TSA001$shrunk, main = "TSA001 vs DMSO")
plotMA(results_list$TSA1$shrunk, main = "TSA1 vs DMSO")
plotMA(results_list$E6R$shrunk, main = "E6R vs DMSO")
dev.off()

#Session info
sink(file.path(output_dir, "session_info_batch2.txt"))
sessionInfo()
sink()



