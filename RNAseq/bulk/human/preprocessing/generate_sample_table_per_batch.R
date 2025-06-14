# -----------------------------------------------------------------
# Script: generate_sample_table_per_batch.R
#
# Description:
# This script generates `sample_table.csv` files for each RNA-seq 
# batch based on the previously merged gene-level count matrices 
# (`DESeq2_counts_matrix_by_gene.tsv`). These sample tables are 
# essential input for DESeq2, providing metadata on each sample's 
# biological condition (e.g., DMSO, TSA, E6R).
#
# The script:
# 1. Loads the merged gene-by-sample count matrix per batch
# 2. Extracts sample names from column headers
# 3. Assigns experimental condition based on name pattern:
#    - "DMSO" = DMSO
#    - "TSA"  = TSA
#    - "E6R"  = E6R
#    - Other  = NA (excluded)
# 4. Writes `sample_table.csv` into the corresponding batch folder
#
# Inputs:
# - DESeq2_counts_matrix_by_gene.tsv (one per batch folder)
#
# Outputs:
# - sample_table.csv (sample, condition) per batch folder
#
# Author: Olivia D.
# Date: 2025-06-13
# -----------------------------------------------------------------

library(tidyverse)

batch_dirs <- c(
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4"
)

walk(batch_dirs, function(batch_dir) {
  counts_path <- file.path(batch_dir, "DESeq2_counts_matrix_by_gene.tsv")
  
  if (!file.exists(counts_path)) {
    message("Counts matrix not found in ", batch_dir)
    return(NULL)
  }
  
  counts <- read_tsv(counts_path, show_col_types = FALSE)
  sample_names <- colnames(counts)[-1]  # Drop 'Name' column
  
  sample_table <- tibble(
    sample = sample_names,
    condition = case_when(
      str_detect(sample, regex("DMSO", ignore_case = TRUE)) ~ "DMSO",
      str_detect(sample, regex("TSA", ignore_case = TRUE)) ~ "TSA",
      str_detect(sample, regex("E6R", ignore_case = TRUE)) ~ "E6R",
      TRUE ~ NA_character_
    )
  ) %>% filter(!is.na(condition))

  output_path <- file.path(batch_dir, "sample_table.csv")
  write_csv(sample_table, output_path)
  message("Sample table saved to: ", output_path)
})

#Output sample tables written to individual batch directories: 
#Sample table saved to: /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1/sample_table.csv
#Sample table saved to: /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2/sample_table.csv
#Sample table saved to: /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3/sample_table.csv
#Sample table saved to: /n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4/sample_table.csv
