#Unified DESeq2 analysis for Batches 1, 3, & 4 (inlude batch2 later)
library(tximport)
library(DESeq2)
library(EnhancedVolcano)

#Batch configuration -----------------------------------------------------
batch_config <- list(
  batch1 = list(
    dir = "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1",
    concentrations = c(TSA = "1uM", E6R = "25uM"),
    timepoint = "16h",
    replicate_pattern = "_batch1"
  ),
  batch3 = list(
    dir = "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3",
    concentrations = c(TSA = "0.1uM", E6R = "1uM"), 
    timepoint = "16h",
    replicate_pattern = "_16hrs"
  ),
  batch4 = list(
    dir = "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4",
    concentrations = c(TSA = "0.01uM", E6R = "0.25uM"),
    timepoint = "120h",
    replicate_pattern = "_5days"
  )
)

#Analysis pipeline -------------------------------------------------------
process_batch <- function(batch_name, config) {
  #Set paths
  output_dir <- file.path(config$dir, paste0("DESeq2_results_", batch_name))
  dir.create(output_dir, showWarnings = FALSE)
  
  #Get sample metadata
  sample_dirs <- list.dirs(config$dir, recursive = FALSE, full.names = TRUE)
  sample_dirs <- sample_dirs[!grepl(paste0("DESeq2_results_", batch_name), sample_dirs)]
  sample_names <- basename(sample_dirs)
  
  #Create file paths & metadata
  files <- file.path(sample_dirs, "quant.sf")
  names(files) <- sample_names
  
  metadata <- data.frame(
    sample = sample_names,
    condition = ifelse(grepl("DMSO", sample_names), "DMSO",
                      ifelse(grepl("TSA", sample_names), "TSA",
                            ifelse(grepl("E6R", sample_names), "E6R", NA))),
    replicate = gsub(paste0(".*([0-9]+)", config$replicate_pattern), "\\1", sample_names),
    timepoint = config$timepoint,
    concentration = ifelse(grepl("TSA", sample_names), config$concentrations["TSA"],
                          ifelse(grepl("E6R", sample_names), config$concentrations["E6R"], "DMSO"))
  )
  
  #Remove any invalid samples
  valid_samples <- !is.na(metadata$condition)
  files <- files[valid_samples]
  metadata <- metadata[valid_samples, ]
  rownames(metadata) <- metadata$sample
  
  #DESeq2 pipeline
  txi <- tximport(files, type = "salmon", txOut = TRUE)
  dds <- DESeqDataSetFromTximport(txi, metadata, ~condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds$condition <- relevel(dds$condition, ref = "DMSO")
  dds <- DESeq(dds, fitType = "local")
  
  #Generate results
  coefficients <- paste0("condition_", c("TSA", "E6R"), "_vs_DMSO")
  
  results_list <- lapply(coefficients, function(coef) {
    res <- results(dds, name = coef, alpha = 0.05)
    res_shrunk <- lfcShrink(dds, coef = coef, type = "apeglm")
    list(raw = res, shrunk = res_shrunk)
  })
  names(results_list) <- c("TSA", "E6R")
  
  #Save results
  save_results <- function(treatment) {
    df <- as.data.frame(results_list[[treatment]]$shrunk)
    df_sig <- df[!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) >= 1, ]
    
    write.csv(df, file.path(output_dir, paste0(batch_name, "_", treatment, "_vs_DMSO_all.csv")))
    write.csv(df_sig, file.path(output_dir, paste0(batch_name, "_", treatment, "_vs_DMSO_filtered.csv")))
    
    #Visualization- Volcano Plots 
    pdf(file.path(output_dir, paste0(batch_name, "_", treatment, "_volcano.pdf")))
    print(EnhancedVolcano(df_sig,
                          lab = rownames(df_sig),
                          x = "log2FoldChange",
                          y = "padj",
                          title = paste(batch_name, treatment, "vs DMSO"),
                          pCutoff = 0.05,
                          FCcutoff = 1))
    dev.off()
  }
  
  lapply(c("TSA", "E6R"), save_results)
  
  #Save session info
  sink(file.path(output_dir, paste0("session_info_", batch_name, ".txt")))
  print(sessionInfo())
  sink()
}

#Execute analysis --------------------------------------------------------
purrr::walk(names(batch_config), ~process_batch(.x, batch_config[[.x]]))
