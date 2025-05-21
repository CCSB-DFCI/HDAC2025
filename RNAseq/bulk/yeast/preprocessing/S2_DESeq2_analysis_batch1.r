#Load required libraries
library(tximport)
library(DESeq2)
library(readr)

#Define directories
batch1_dir <- "/n/scratch/users/o/old656/yeast_31032025_Tina/salmon_output/"
tx2gene <- read_csv("/n/scratch/users/o/old656/RNA_seq/yeast/reference_genome/tx2gene.csv")
output_dir <- "/n/scratch/users/o/old656/RNA_seq/yeast/salmon_output/DESeq2_results/pairwise_normalization/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#Get all quant.sf paths from batch1
samples_batch1 <- list.files(path = batch1_dir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)

#Extract sample names and conditions
sample_names <- gsub("_quant/quant.sf", "", basename(dirname(samples_batch1)))
condition <- gsub("_\\d+_quant$", "", sample_names)

#Create coldata for batch1
coldata <- data.frame(sample = sample_names, condition = condition)
rownames(coldata) <- coldata$sample
coldata$condition <- factor(coldata$condition)

#Import quantifications
txi <- tximport(samples_batch1, type = "salmon", tx2gene = tx2gene)

#Run DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

#Define comparisons
comparisons <- list(
  E6R_vs_DMSO = c("condition", "E6R", "DMSO"),
  TSA_vs_DMSO = c("condition", "TSA", "DMSO"),
  rpd3_vs_DMSO = c("condition", "rpd3D", "DMSO"),
  sin3_vs_DMSO = c("condition", "sin3D", "DMSO"),
  rpd3_vs_sin3 = c("condition", "rpd3D", "sin3D")
)

#Set thresholds
padj_threshold <- 0.05 #only genes with adjusted p-value < 0.05 are considered significant
lfc_threshold <- 1 #only genes with |log2FoldChange| ≥ 1 are considered differentially expressed 

#Run your DEG saving block
summary_list <- list()

for (name in names(comparisons)) {
  contrast <- comparisons[[name]]
  res <- results(dds, contrast = contrast)
  res_df <- as.data.frame(res)

  #Save raw DEGs
  raw_file <- paste0(output_dir, "DESeq2_", name, "_raw_batch1.csv")
  write.csv(res_df, file = raw_file, row.names = TRUE)

  #Filter
  filtered_df <- res_df[!is.na(res_df$padj) & res_df$padj < padj_threshold & abs(res_df$log2FoldChange) >= lfc_threshold, ]
  filtered_file <- paste0(output_dir, "DESeq2_", name, "_filtered_batch1.csv")
  write.csv(filtered_df, file = filtered_file, row.names = TRUE)

  #DEG summary
  up <- sum(filtered_df$log2FoldChange > 0)
  down <- sum(filtered_df$log2FoldChange < 0)
  total <- nrow(filtered_df)

  summary_list[[name]] <- data.frame(
    Comparison = name,
    Total_DEGs = total,
    Upregulated = up,
    Downregulated = down
  )

  cat(paste0("", name, " - Total: ", total, " | Up: ", up, " | Down: ", down, "\n"))
}
