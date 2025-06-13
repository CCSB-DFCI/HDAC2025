# ---------------------------------------------------------------
# Script: DESeq2_counts_matrix_by_gene_per_batch.R
#
# Description:
# This script processes transcript-level Salmon quantification 
# files (quant.sf) from multiple RNA-seq experiment batches, 
# aggregates transcript-level TPM & counts to gene level using 
# tximport and a cleaned tx2gene map, and outputs:
#
# 1. quant_aggr_by_gene.tsv per sample (with counts, TPM, length)
# 2. DESeq2_counts_matrix_by_gene.tsv per batch (gene × sample)
# 3. sample_table.csv per batch (for use in DESeq2 design)
#
# Each batch is handled independently to allow for batch-specific
# DESeq2 comparisons (e.g., DMSO vs. treated within a batch).
#
# Inputs:
# - quant.sf files (from Salmon)
# - tx2gene_from_index_clean.csv (transcript → gene mapping)
#
# Outputs (per batch):
# - quant_aggr_by_gene.tsv in each sample directory
# - DESeq2_counts_matrix_by_gene.tsv
# - sample_table.csv
#
# Author: [Olivia D.]
# Date: [2025-06-13]
# ---------------------------------------------------------------

library(tidyverse)

#Our batch directories consisting of Salmon quant_aggr_by_gene.tsv files: 
batch_dirs <- c(
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4"
)

walk(batch_dirs, function(batch_dir) {
  message("\n Merging counts in: ", basename(batch_dir))
  
  files <- list.files(batch_dir, pattern = "quant_aggr_by_gene.tsv$", recursive = TRUE, full.names = TRUE)
  
  if (length(files) == 0) {
    message("No quant_aggr_by_gene.tsv found in ", batch_dir)
    return(NULL)
  }
  
  count_list <- map(files, ~ read_tsv(.x, show_col_types = FALSE) %>%
                      select(Name, NumReads) %>%
                      rename(!!basename(dirname(.x)) := NumReads))
  
  counts_matrix <- reduce(count_list, full_join, by = "Name") %>%
    arrange(Name)
  
  counts_matrix[is.na(counts_matrix)] <- 0
  
  output_path <- file.path(batch_dir, "DESeq2_counts_matrix_by_gene.tsv")
  write_tsv(counts_matrix, output_path)
  
  message("Saved: ", output_path)
})

#In bash, double-checked each file:
#less  /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1/DESeq2_counts_matrix_by_gene.tsv
#less /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2/DESeq2_counts_matrix_by_gene.tsv
#less /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3/DESeq2_counts_matrix_by_gene.tsv
#less /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4/DESeq2_counts_matrix_by_gene.tsv
