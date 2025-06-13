#------------------------------------------------------------------
#Script: DESeq2_gene_level_batchwise.R
#
#Description:
#This script runs gene-level differential expression analysis
#using DESeq2 for each RNA-seq batch based on gene-aggregated
#Salmon quantifications. It uses:
#   - DESeq2_counts_matrix_by_gene.tsv
#   - sample_table.csv (with biological condition labels: DMSO, TSA, E6R)
#
#The script performs:
#1. DESeq2 modeling
#2. log2FC shrinkage using apeglm
#3. Result export (raw, shrunk & statistically filtered)
#
#Author: Olivia D.
#Date: 2025-06-13
#-------------------------------------------------------------------

library(DESeq2)
library(tidyverse)
library(apeglm)

#Define batch directories
batch_dirs <- c(
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4"
)

#Subdirectory to store DESeq2 results
out_subdir <- "DESeq2_analysis_by_gene_13062025"

#Navigate through the batch directories (1-4)
walk(batch_dirs, function(batch_dir) {
  message("\n🔬 Running DESeq2 for batch: ", basename(batch_dir))

  counts_path <- file.path(batch_dir, "DESeq2_counts_matrix_by_gene.tsv")
  sample_path <- file.path(batch_dir, "sample_table.csv")
  out_dir <- file.path(batch_dir, out_subdir)
  dir.create(out_dir, showWarnings = FALSE)

  if (!file.exists(counts_path) | !file.exists(sample_path)) {
    message("Skipping ", basename(batch_dir), ": missing counts or sample table.")
    return(NULL)
  }

  counts <- read_tsv(counts_path, show_col_types = FALSE) %>%
    column_to_rownames("Name") %>%
    as.data.frame()

  sample_table <- read_csv(sample_path, show_col_types = FALSE) %>%
    column_to_rownames("sample")

  stopifnot(all(colnames(counts) == rownames(sample_table)))

  sample_table$condition <- factor(sample_table$condition)
  dds <- DESeqDataSetFromMatrix(countData = round(counts),
                                colData = sample_table,
                                design = ~ condition)

  dds$condition <- relevel(dds$condition, ref = "DMSO")
  dds <- DESeq(dds)

  results_list <- list()
  for (treat in setdiff(levels(dds$condition), "DMSO")) {
    contrast <- c("condition", treat, "DMSO")
    coef_name <- paste0("condition_", treat, "_vs_DMSO")

    if (!(coef_name %in% resultsNames(dds))) {
      message("Skipping contrast not found in model: ", coef_name)
      next
    }

    res <- results(dds, contrast = contrast)
    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

    res_df <- as.data.frame(res)
    res_shrunk_df <- as.data.frame(res_shrunk)

    res_sig <- res_shrunk_df %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1) #Filtered subset of the shrunk results → only genes with padj < 0.05 & LFC cut-off 

    #Save all results within individual batch dir
    write.csv(res_df, file.path(out_dir, paste0("DESeq2_", treat, "_vs_DMSO_raw.csv")))
    write.csv(res_shrunk_df, file.path(out_dir, paste0("DESeq2_", treat, "_vs_DMSO_shrunk.csv")))
    write.csv(res_sig, file.path(out_dir, paste0("DESeq2_", treat, "_vs_DMSO_DEGs_FDR05_LFC1.csv")))

    results_list[[treat]] <- list(raw = res_df, shrunk = res_shrunk_df, filtered = res_sig)

    message("Saved: ", treat, " vs DMSO")
  }

  message("Completed batch: ", basename(batch_dir))
})


#Create DEG summary: Compare filtered vs shrunk (as shrunk.csv is already denoised)
#Batch-level DESeq2 result folders
batch_dirs <- c(
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1/DESeq2_analysis_by_gene_13062025",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2/DESeq2_analysis_by_gene_13062025",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3/DESeq2_analysis_by_gene_13062025",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4/DESeq2_analysis_by_gene_13062025"
)

deg_summary <- map_dfr(batch_dirs, function(batch_dir) {
  shrunk_paths <- list.files(batch_dir, pattern = "_shrunk\\.csv$", full.names = TRUE)
  filtered_paths <- list.files(batch_dir, pattern = "_DEGs_FDR05_LFC1\\.csv$", full.names = TRUE)
  
  #Extract comparison labels like TSA, E6R etc.
  shrunk_tbl <- tibble(
    comparison = shrunk_paths %>%
      basename() %>%
      str_extract("(?<=DESeq2_).+(?=_vs_DMSO_shrunk\\.csv)"),
    shrunk_path = shrunk_paths
  )
  
  filtered_tbl <- tibble(
    comparison = filtered_paths %>%
      basename() %>%
      str_extract("(?<=DESeq2_).+(?=_vs_DMSO_DEGs_FDR05_LFC1\\.csv)"),
    filtered_path = filtered_paths
  )

  full_join(shrunk_tbl, filtered_tbl, by = "comparison") %>%
    mutate(
      shrunk_genes = map_int(shrunk_path, ~ nrow(read_csv(.x, show_col_types = FALSE))),
      filtered_DEGs = map_int(filtered_path, ~ if (!is.na(.x) && file.exists(.x)) nrow(read_csv(.x, show_col_types = FALSE)) else 0),
      batch = basename(dirname(dirname(batch_dir)))
    ) %>%
    select(batch, comparison, shrunk_genes, filtered_DEGs)
})

#Save CSV
write_csv(deg_summary, "DESeq2_DEG_summary_by_gene_13062025.csv")
message("DEG summary saved to DESeq2_DEG_summary_by_gene_13062025.csv")


