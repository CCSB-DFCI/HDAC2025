#Libraries
library(tximport)
library(DESeq2)
library(readr)

#Batch 2 samples only
batch2_dir <- "/n/scratch/users/o/old656/RNA_seq/yeast/salmon_output/"
output_dir <- "/n/scratch/users/o/old656/RNA_seq/yeast/salmon_output/DESeq2_results/pairwise_normalization/"

#Load sample files
samples_batch2 <- list.files(path = batch2_dir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)

#Filter for ume6 & M208
samples_batch2 <- samples_batch2[grepl("ume6|M208", samples_batch2)]

#Sample names
sample_names <- gsub("_quant/quant.sf", "", basename(dirname(samples_batch2)))
condition <- ifelse(grepl("ume6", sample_names), "ume6", "M208")

#Metadata
coldata <- data.frame(sample = sample_names, condition = factor(condition))
rownames(coldata) <- sample_names

#tx2gene path
tx2gene <- read_csv("/n/scratch/users/o/old656/RNA_seq/yeast/reference_genome/tx2gene.csv")

#Import
txi <- tximport(samples_batch2, type = "salmon", tx2gene = tx2gene)

#Run DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

#Get results
res <- results(dds, contrast = c("condition", "ume6", "M208"))
res_df <- as.data.frame(res)

#Save raw
write.csv(res_df, file = paste0(output_dir, "DESeq2_ume6_vs_M208_raw_batch2.csv"), row.names = TRUE)

#Filter
padj_threshold <- 0.05
lfc_threshold <- 1
filtered_df <- res_df[!is.na(res_df$padj) & res_df$padj < padj_threshold & abs(res_df$log2FoldChange) >= lfc_threshold, ]
write.csv(filtered_df, file = paste0(output_dir, "DESeq2_ume6_vs_M208_filtered_batch2.csv"), row.names = TRUE)

#Print summary
cat(paste0("ume6_vs_M208 (Batch 2) — Total: ", nrow(filtered_df),
           " | Up: ", sum(filtered_df$log2FoldChange > 0),
           " | Down: ", sum(filtered_df$log2FoldChange < 0), "\n"))
