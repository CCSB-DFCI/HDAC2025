library(tximport)
library(readr)
library(tidyverse)

# ------------------------------
# 0. Generate tx2gene from FASTA
# ------------------------------
#Run this once to generate the mapping if not already done
#Replace path if different
fasta_path <- "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_reference_genome/Homo_sapiens.GRCh38.cdna.all.fa.gz"

tx2gene <- read_lines(fasta_path) %>%
  str_subset("^>") %>%
  tibble(header = .) %>%
  mutate(
    TXNAME = str_extract(header, "^>\\S+") %>% str_remove("^>") %>% str_remove("\\.[0-9]+$"),
    GENEID = str_extract(header, "gene:\\S+") %>% str_remove("gene:") %>% str_remove("\\.[0-9]+$")
  ) %>%
  filter(!is.na(TXNAME) & !is.na(GENEID)) %>%
  distinct(TXNAME, GENEID)

write_csv(tx2gene, "tx2gene_from_index_clean.csv")  #here we stripped off the suffixes 

# ------------------------------
# 1. Load cleaned tx2gene
# ------------------------------
tx2gene <- read_csv("tx2gene_from_index_clean.csv")

# ------------------------------
# 2. Define batch directories
# ------------------------------
batch_dirs <- c(
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch1/salmon_output_batch1",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/neuroblastoma_batch2/salmon_output_batch2",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_16hrs_batch3_27022025/salmon_output_batch3",
  "/n/data2/dfci/genetics/vidal/HDAC_project/RNA_seq/human_SKNBE2C_5days_batch4/salmon_output_batch4"
)

# ---------------------------------------
# 3. Helper: safely coerce vector/matrix
# --------------------------------------
to_matrix <- function(x, sample_name) {
  if (is.null(x)) return(NULL)
  if (is.vector(x)) {
    matrix(x, ncol = 1, dimnames = list(names(x), sample_name))
  } else {
    matrix(x, ncol = 1, dimnames = list(rownames(x), sample_name))
  }
}

# ------------------------------
# 4. Loop over samples
# ------------------------------
walk(batch_dirs, function(batch_dir) {
  message("\n Processing batch: ", basename(batch_dir))

  #Salmon quantification files per sample per batch dirs (1-4)
  quant_files <- list.files(batch_dir, pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)
  
  walk(quant_files, function(qf) {
    sample_dir <- dirname(qf)
    sample_name <- basename(sample_dir)
    output_file <- file.path(sample_dir, "quant_aggr_by_gene.tsv")
    
    tryCatch({
      txi <- tximport(
        files = setNames(qf, sample_name),
        type = "salmon",
        tx2gene = tx2gene,
        ignoreTxVersion = TRUE,
        dropInfReps = TRUE,
        countsFromAbundance = "no"  # You can set to "lengthScaledTPM" if needed
      )
      
      #Coerce outputs (vector → matrix if needed)
      txi$counts <- to_matrix(txi$counts, sample_name)
      txi$abundance <- to_matrix(txi$abundance, sample_name)
      txi$length <- to_matrix(txi$length, sample_name)
      txi$effLen <- to_matrix(txi$effLen, sample_name)
      
      # ------------------------------
      # Why do we check for txi$effLen?!
      # ------------------------------
      #In some cases (e.g., if the EffectiveLength column is missing/NA in quant.sf, 
      #or all values are zero), tximport will return effLen = NULL.
      #If we attempt to access txi$effLen[, 1] when it's NULL, R will throw:
      #❌ Error in matrix(x, ...) : 'data' must be of a vector type, was 'NULL'
      #To avoid this, we check with `is.null()` & populate NA if absent.
      #This allows us to still export a valid table and run DESeq2.
      has_efflen <- !is.null(txi$effLen)
      
      gene_quant <- data.frame(
        Name = rownames(txi$counts),
        Length = round(txi$length[, 1], 2),
        TPM = round(txi$abundance[, 1], 3),
        NumReads = round(txi$counts[, 1], 2),
        EffectiveLength = if (has_efflen) round(txi$effLen[, 1], 2) else NA
      )
      
      write_tsv(gene_quant, output_file)
      message("✅ Saved: ", output_file)
    }, error = function(e) {
      message("❌ Failed: ", sample_name, " — ", e$message)
    })
  })
})


#Scan our batch directories, check for quant_aggr_by_gene.tsv files & print file size/number of genes per file:
#Get info on each aggregated file: 

file_info <- map_dfr(batch_dirs, function(batch_dir) {
  tsv_files <- list.files(batch_dir, pattern = "quant_aggr_by_gene.tsv$", recursive = TRUE, full.names = TRUE)
  
  map_dfr(tsv_files, function(f) {
    tryCatch({
      df <- read_tsv(f, show_col_types = FALSE)
      tibble(
        sample = basename(dirname(f)),
        file = f,
        size_kb = round(file.info(f)$size / 1024, 1),
        n_genes = nrow(df)
      )
    }, error = function(e) {
      tibble(sample = basename(dirname(f)), file = f, size_kb = NA, n_genes = NA)
    })
  })
})

#Print result
print(file_info, n = Inf) #37951 genes 

#Quick sanity check: Every sample directory with a quant.sf also has a quant_aggr_by_gene.tsv
#Gather quant.sf & quant_aggr_by_gene.tsv paths
check_df <- map_dfr(batch_dirs, function(batch_dir) {
  quant_paths <- list.files(batch_dir, pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)
  
  tibble(
    sample = basename(dirname(quant_paths)),
    batch = basename(batch_dir),
    quant_sf = quant_paths,
    quant_aggr = file.path(dirname(quant_paths), "quant_aggr_by_gene.tsv"),
    aggr_exists = file.exists(file.path(dirname(quant_paths), "quant_aggr_by_gene.tsv"))
  )
})

#Show missing ones (if any)
missing_aggr <- filter(check_df, !aggr_exists)

#Report
if (nrow(missing_aggr) == 0) {
  message("All quant.sf files have corresponding quant_aggr_by_gene.tsv files.")
} else {
  message("Missing quant_aggr_by_gene.tsv for the following samples:")
  print(missing_aggr)
}

write_csv(check_df, "quant_aggr_sanity_check.csv")
